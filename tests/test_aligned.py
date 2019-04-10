import pytest
import shutil
from pathlib import Path
import pypipegraph as ppg
import mbf_align
import pysam
import pandas as pd
from mbf_qualitycontrol.testing import assert_image_equal
from mbf_sampledata import get_sample_data


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


def get_human_22_fake_genome():
    from mbf_genomics.testing import MockGenome
    import gzip

    genes = pd.read_msgpack(
        gzip.GzipFile(get_sample_data(Path("mbf_align/hs_22_genes.msgpack.gz")))
    ).reset_index()
    tr = pd.read_msgpack(
        gzip.GzipFile(get_sample_data(Path("mbf_align/hs_22_transcripts.msgpack.gz")))
    ).reset_index()
    return MockGenome(df_genes=genes, df_transcripts=tr, chr_lengths={"22": 50818468})


@pytest.mark.usefixtures("new_pipegraph")
class TestQualityControl:
    def test_qc_plots(self, new_pipegraph):
        new_pipegraph.new_pipegraph(quiet=False)
        from mbf_qualitycontrol import do_qc

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
