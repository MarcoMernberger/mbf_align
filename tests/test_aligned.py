import pytest
import shutil
from pathlib import Path
import pypipegraph as ppg
from .shared import data_path
import mbf_align
import pysam
from mbf_qualitycontrol.testing import assert_image_equal


@pytest.mark.usefixtures("new_pipegraph")
class TestAligned:
    def test_from_existing_bam(self):
        bam_path = Path(data_path) / "ex2.bam"
        bam_job = ppg.FileInvariant(bam_path)
        genome = object()
        lane = mbf_align.AlignedSample("test_lane", bam_job, genome, False, "AA123")
        assert lane.name == "test_lane"
        assert lane.load()[0] is bam_job
        assert isinstance(lane.load()[1], ppg.FileInvariant)
        assert lane.genome is genome
        assert not lane.is_paired
        assert lane.vid == "AA123"

        with pytest.raises(ValueError):
            mbf_align.AlignedSample("test_lane", bam_job, genome, False, "AA123")
        lane2 = mbf_align.AlignedSample("test_lane2", bam_job, genome, True, "AA123")
        assert lane2.is_paired

        b = lane.get_bam()
        assert isinstance(b, pysam.Samfile)
        b = lane.get_unique_aligned_bam()
        assert isinstance(b, pysam.Samfile)

        assert lane.mapped_reads() == 8
        assert lane.unmapped_reads() == 0

    def test_lane_invariants_on_non_accepted_value(self):
        genome = object()
        with pytest.raises(ValueError):
            mbf_align.AlignedSample("test_lane", 123, genome, False, "AA123")

    def test_lane_raises_on_multifilegeneratingJobWithTwoBams(self):
        mfg = ppg.MultiFileGeneratingJob(["a.bam", "b.bam"], lambda: 5)
        genome = object()
        with pytest.raises(ValueError):
            mbf_align.AlignedSample("test_lane", mfg, genome, False, "AA123")

    def test_lane_raises_on_multifilegeneratingJobWithTwoBais(self):
        mfg = ppg.MultiFileGeneratingJob(["a.bam", "a.bam.bai", "b.bam.bai"], lambda: 5)
        genome = object()
        with pytest.raises(ValueError):
            mbf_align.AlignedSample("test_lane", mfg, genome, False, "AA123")

    def test_lane_raises_on_multifilegeneratingJobWithNoBAM(self):
        mfg = ppg.MultiFileGeneratingJob(["a.sam"], lambda: 5)
        genome = object()
        with pytest.raises(ValueError):
            mbf_align.AlignedSample("test_lane", mfg, genome, False, "AA123")

    def test_lane_invariants_on_string(self):
        bam_path = Path(data_path) / "ex2.bam"
        genome = object()
        lane = mbf_align.AlignedSample("test_lane", bam_path, genome, False, "AA123")
        assert isinstance(lane.load()[0], ppg.FileInvariant)

    def test_missing_index_file(self):
        bam_path = Path(data_path) / "ex2.bam"
        no_index = "noindex.bam"
        shutil.copy(bam_path, no_index)
        genome = object()
        lane = mbf_align.AlignedSample("test_lane", no_index, genome, False, "AA123")
        assert isinstance(lane.load()[0], ppg.FileInvariant)
        assert isinstance(lane.load()[1], ppg.FileGeneratingJob)
        assert lane.load()[1].job_id != "noindex.bam.bai"
        assert lane.load()[0] in lane.load()[1].prerequisites
        with pytest.raises(FileNotFoundError):
            lane.mapped_reads()
        ppg.run_pipegraph()
        assert lane.mapped_reads() == 8

    def test_creating_index_for_fg_job(self):
        def gen():
            shutil.copy(Path(data_path) / "ex2.bam", "sample.bam")

        ppg.util.global_pipegraph.quiet = False

        job = ppg.FileGeneratingJob("sample.bam", gen)
        genome = object()
        lane = mbf_align.AlignedSample("test_lane", job, genome, False, "AA123")
        assert isinstance(lane.load()[1], ppg.FileGeneratingJob)
        assert lane.load()[0] in lane.load()[1].prerequisites
        ppg.run_pipegraph()
        assert Path("sample.bam").exists()
        assert Path("sample.bam.bai").exists()

    def test_creating_index_for_mfg(self):
        def gen():
            shutil.copy(Path(data_path) / "ex2.bam", "sample.bam")

        ppg.util.global_pipegraph.quiet = False

        job = ppg.MultiFileGeneratingJob(["sample.bam"], gen)
        genome = object()
        lane = mbf_align.AlignedSample("test_lane", job, genome, False, "AA123")
        assert isinstance(lane.load()[1], ppg.FileGeneratingJob)
        assert lane.load()[0] in lane.load()[1].prerequisites
        ppg.run_pipegraph()
        assert Path("sample.bam").exists()
        assert Path("sample.bam.bai").exists()


@pytest.mark.usefixtures("new_pipegraph")
class TestQualityControl:
    def test_complexity(self):
        from mbf_qualitycontrol import do_qc

        genome = object()
        lane = mbf_align.AlignedSample(
            "test_lane", Path(data_path) / "rnaseq_spliced.bam", genome, False, "AA123"
        )
        do_qc()
        ppg.run_pipegraph()
        p = lane.result_dir / "complexity.png"
        assert_image_equal(p)
