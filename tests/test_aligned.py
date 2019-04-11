import pytest
import shutil
from pathlib import Path
import pypipegraph as ppg
import mbf_align
import pysam
from mbf_qualitycontrol.testing import assert_image_equal
from mbf_sampledata import get_sample_data, get_sample_path


@pytest.mark.usefixtures("new_pipegraph")
class TestAligned:
    def test_from_existing_bam(self):
        bam_path = get_sample_data(Path("mbf_align/ex2.bam"))
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
        bam_path = get_sample_data(Path("mbf_align/ex2.bam"))
        genome = object()
        lane = mbf_align.AlignedSample("test_lane", bam_path, genome, False, "AA123")
        assert isinstance(lane.load()[0], ppg.FileInvariant)

    def test_missing_index_file(self):
        bam_path = get_sample_data(Path("mbf_align/ex2.bam"))
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
            shutil.copy(get_sample_data(Path("mbf_align/ex2.bam")), "sample.bam")

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
            shutil.copy(get_sample_data(Path("mbf_align/ex2.bam")), "sample.bam")

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
    def test_qc_plots(self, new_pipegraph):
        new_pipegraph.new_pipegraph(quiet=False)
        from mbf_qualitycontrol import do_qc
        from mbf_sampledata import get_human_22_fake_genome

        # straight from chr22 of the human genome
        genome = get_human_22_fake_genome()

        lane = mbf_align.AlignedSample(
            "test_lane",
            get_sample_data(Path("mbf_align/rnaseq_spliced_chr22.bam")),
            genome,
            False,
            "AA123",
        )  # index creation is automatic
        assert "bam_indices" in lane.index_job.job_id
        do_qc()
        ppg.run_pipegraph()
        p = lane.result_dir / "complexity.png"
        assert_image_equal(p, suffix="_complexity")
        p = lane.result_dir / "strandedness.png"
        assert_image_equal(p, suffix="_strandedness")
        p = lane.result_dir / "reads_per_biotype.png"
        assert_image_equal(p, suffix="_biotypes")
        p = lane.result_dir / ".." / "alignment_statistics.png"
        assert_image_equal(p, suffix="_alignment_statistics")
        p = lane.result_dir / f"subchromosomal_distribution.png"
        assert_image_equal(p, suffix="_subchromosomal_distribution")

    def test_alignment_stats(self):
        genome = object()
        lane = mbf_align.AlignedSample(
            "test_lane",
            get_sample_data(Path("mbf_align/rnaseq_spliced_chr22.bam")),
            genome,
            False,
            "AA123",
        )  # index creation is automatic
        counts = {"get_bam": 0}

        def get_bam():
            counts["get_bam"] += 1

            class DummySam:
                mapped = 5
                unmapped = 10

                def __enter__(self):
                    return self

                def __exit__(self, *args):
                    pass

            return DummySam()

        lane.get_bam = get_bam
        assert lane.get_alignment_stats() == {"Mapped": 5, "Unmapped": 10}
        assert counts["get_bam"] == 1

        class DummyAlignerWithout:
            pass

        lane = mbf_align.AlignedSample(
            "test_lane2",
            get_sample_data(Path("mbf_align/rnaseq_spliced_chr22.bam")),
            genome,
            False,
            "AA123",
            aligner=DummyAlignerWithout(),
        )  # index creation is automatic
        lane.get_bam = get_bam
        assert counts["get_bam"] == 1
        assert lane.get_alignment_stats() == {"Mapped": 5, "Unmapped": 10}
        assert counts["get_bam"] == 2

        class DummyAlignerWith:
            def get_alignment_stats(self, bam_filename):
                assert (
                    Path(bam_filename).resolve()
                    == get_sample_path("mbf_align/rnaseq_spliced_chr22.bam").resolve()
                )
                return {"Hello": 23}

        lane = mbf_align.AlignedSample(
            "test_lane3",
            get_sample_data("mbf_align/rnaseq_spliced_chr22.bam"),
            genome,
            False,
            "AA123",
            aligner=DummyAlignerWith(),
        )  # index creation is automatic
        lane.get_bam = get_bam
        assert counts["get_bam"] == 2
        assert lane.get_alignment_stats() == {"Hello": 23}
        assert counts["get_bam"] == 2
