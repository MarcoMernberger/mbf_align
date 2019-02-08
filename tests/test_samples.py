import pytest
import gzip
from pathlib import Path
from pypipegraph import FileGeneratingJob, MultiFileGeneratingJob
import requests_mock
import pypipegraph as ppg

from mbf_align import (
    FASTQsFromFile,
    FASTQsFromFolder,
    FASTQsFromJob,
    FASTQsFromURLs,
    FASTQsFromAccession,
    build_fastq_strategy,
)
from mbf_align import Sample
from mbf_align import PairingError
from mbf_align import fastq2

from .shared import data_path


def test_FASTQsFromFile():
    fn = data_path / "sample_a" / ".." / "sample_a" / "a.fastq"
    o = FASTQsFromFile(fn)
    assert o() == [(fn.resolve(),)]


def test_FASTQsFromFileRaisesOnMissing():
    fn = data_path / "sample_a" / "a.fastq.nosuchfile"
    with pytest.raises(IOError):
        FASTQsFromFile(fn)


def test_FASTQsFromFilePaired():
    fn = data_path / "sample_b" / "a_R1_.fastq.gz"
    fn2 = data_path / "sample_b" / ".." / "sample_b" / "a_R2_.fastq.gz"
    o = FASTQsFromFile(fn, fn2)
    assert o() == [(fn.resolve(), fn2.resolve())]


def test_FASTQsFromFilePairedMissingR2():
    fn = data_path / "sample_b" / "a_R1_.fastq.gz"
    fn2 = data_path / "sample_b" / "a_R2_.fastq.gz.nosuchfile"
    with pytest.raises(IOError):
        FASTQsFromFile(fn, fn2)


def test_FASTQsFromFolder():
    folder = data_path / "sample_a"
    o = FASTQsFromFolder(folder)
    import pprint

    pprint.pprint(o())
    assert o() == [
        ((folder / "a.fastq").resolve(),),
        ((folder / "b.fastq.gz").resolve(),),
    ]


def test_FASTQsFromFolder_raises_on_non_existing():
    with pytest.raises(IOError):
        FASTQsFromFolder("shu")


def test_FASTQsFromFolder_raises_on_no_fastqs():
    with pytest.raises(ValueError):
        FASTQsFromFolder(data_path / "sample_f")


def test_FASTQsFromFolder_raises_on_not_a_folder():
    with pytest.raises(ValueError):
        FASTQsFromFolder(data_path / "sample_a" / "a.fastq")


def test_FASTQsFromFolderPaired():
    folder = data_path / "sample_b"
    o = FASTQsFromFolder(folder)
    assert o() == [
        ((folder / "a_R1_.fastq.gz").resolve(), (folder / "a_R2_.fastq.gz").resolve())
    ]


def test_FASTQsFromFolderR2_but_missing_any_r1():
    folder = data_path / "sample_c"
    o = FASTQsFromFolder(folder)
    with pytest.raises(ValueError):
        o()


def test_FASTQsFromFolder_pairing_files_fails():
    folder = data_path / "sample_d"
    o = FASTQsFromFolder(folder)
    with pytest.raises(ValueError):
        o()


def test_FASTQsFromFolder_pairing_files_fails2():
    folder = data_path / "sample_e"
    o = FASTQsFromFolder(folder)
    with pytest.raises(ValueError):
        o()


@pytest.mark.usefixtures("new_pipegraph")
class TestSamples:
    def test_FASTQsFromJob(self):
        job = FileGeneratingJob("test.fastq.gz", lambda of: None)
        o = FASTQsFromJob(job)
        assert o() == [(Path("test.fastq.gz").resolve(),)]

    def test_FASTQsFromJob_R1_ok(self):
        job = FileGeneratingJob("test_R1_.fastq.gz", lambda of: None)
        o = FASTQsFromJob(job)
        assert o() == [(Path("test_R1_.fastq.gz").resolve(),)]

    def test_FASTQsFromJob_Multiple_Unpaired(self):
        job = MultiFileGeneratingJob(
            ["test.fastq.gz", "test2.fastq.gz"], lambda of: None
        )
        o = FASTQsFromJob(job)
        assert o() == [
            (Path("test.fastq.gz").resolve(),),
            (Path("test2.fastq.gz").resolve(),),
        ]

    def test_FASTQsFromJob_Multiple_Unpaired_R1(self):
        job = MultiFileGeneratingJob(
            ["test_R1_.fastq.gz", "test2_R1_.fastq.gz"], lambda of: None
        )
        o = FASTQsFromJob(job)
        assert o() == [
            # 2 sorts before _
            (Path("test2_R1_.fastq.gz").resolve(),),
            (Path("test_R1_.fastq.gz").resolve(),),
        ]

    def test_FASTQsFromJob_Multiple_Paired(self):
        job = MultiFileGeneratingJob(
            [
                "test_R1_.fastq.gz",
                "test2_R1_.fastq.gz",
                "test_R2_.fastq.gz",
                "test2_R2_.fastq.gz",
            ],
            lambda of: None,
        )
        o = FASTQsFromJob(job)
        assert set(o()) == set(
            [
                # 2 sorts before _
                (
                    Path("test2_R1_.fastq.gz").resolve(),
                    Path("test2_R2_.fastq.gz").resolve(),
                ),
                (
                    Path("test_R1_.fastq.gz").resolve(),
                    Path("test_R2_.fastq.gz").resolve(),
                ),
            ]
        )

    def test_build_fastq_strategy(self):
        # single filename
        assert build_fastq_strategy(
            data_path / "sample_b" / ".." / "sample_b" / "a_R1_.fastq.gz"
        )() == [((data_path / "sample_b" / "a_R1_.fastq.gz").resolve(),)]
        assert build_fastq_strategy(
            str(data_path / "sample_b" / ".." / "sample_b" / "a_R1_.fastq.gz")
        )() == [((data_path / "sample_b" / "a_R1_.fastq.gz").resolve(),)]
        # multiple files
        assert build_fastq_strategy(
            [
                data_path / "sample_b" / ".." / "sample_b" / "a_R1_.fastq.gz",
                data_path / "sample_b" / ".." / "sample_b" / "a_R2_.fastq.gz",
            ]
        )() == [
            ((data_path / "sample_b" / "a_R1_.fastq.gz").resolve(),),
            ((data_path / "sample_b" / "a_R2_.fastq.gz").resolve(),),
        ]
        # folder
        assert build_fastq_strategy(data_path / "sample_b" / ".." / "sample_b")() == [
            (
                (data_path / "sample_b" / "a_R1_.fastq.gz").resolve(),
                (data_path / "sample_b" / "a_R2_.fastq.gz").resolve(),
            )
        ]
        # job
        assert isinstance(
            build_fastq_strategy(FileGeneratingJob("test.fastq", lambda of: None)),
            FASTQsFromJob,
        )

        # pass through
        fn = data_path / "sample_a" / ".." / "sample_a" / "a.fastq"
        o = FASTQsFromFile(fn)
        assert build_fastq_strategy(o) is o

        with pytest.raises(ValueError):
            build_fastq_strategy(55)

    def test_lane(self):

        lane = Sample("Sample_a", data_path / "sample_a", False, vid="VA000")
        assert lane.vid == "VA000"
        temp_job = lane.prepare_input()
        real_job = lane.save_input()
        ppg.run_pipegraph()
        assert not Path(temp_job.filenames[0]).exists()
        assert Path(real_job.filenames[0]).exists()
        with gzip.GzipFile(real_job.filenames[0], "r") as op:
            lines = op.readlines()
            assert len(lines) == 20 + 20

    def test_paired_modes(self):
        with pytest.raises(PairingError):
            lane = Sample("Sample_a", data_path / "sample_b", False, vid="VA000")
            lane.prepare_input()

    def test_lane_paired_straight(self):

        lane = Sample(
            "Sample_a", data_path / "sample_b", False, vid="VA000", pairing="paired"
        )
        assert lane.vid == "VA000"
        temp_job = lane.prepare_input()
        real_job = lane.save_input()
        ppg.run_pipegraph()
        assert not Path(temp_job.filenames[0]).exists()
        assert not Path(temp_job.filenames[1]).exists()
        assert Path(real_job.filenames[0]).exists()
        assert Path(real_job.filenames[1]).exists()
        assert "_R1_" in real_job.filenames[0]
        assert "_R2_" in real_job.filenames[1]
        assert ".fastq.gz" in real_job.filenames[0]
        assert ".fastq.gz" in real_job.filenames[1]

        for input_fn, output_fn in zip(
            [
                (data_path / "sample_b" / "a_R1_.fastq.gz"),
                (data_path / "sample_b" / "a_R2_.fastq.gz"),
            ],
            real_job.filenames,
        ):
            with gzip.GzipFile(output_fn, "r") as op:
                actual = op.read()
            with gzip.GzipFile(input_fn, "r") as op:
                should = op.read()
            assert actual == should

    def test_lane_paired_filtered(self):

        lane = Sample(
            "Sample_a",
            data_path / "sample_b",
            False,
            vid="VA000",
            pairing="paired",
            fastq_processor=fastq2.Paired_Filtered(lambda *args: True),
        )
        assert lane.vid == "VA000"
        temp_job = lane.prepare_input()
        real_job = lane.save_input()
        ppg.run_pipegraph()
        assert not Path(temp_job.filenames[0]).exists()
        assert not Path(temp_job.filenames[1]).exists()
        assert Path(real_job.filenames[0]).exists()
        assert Path(real_job.filenames[1]).exists()
        assert "_R1_" in real_job.filenames[0]
        assert "_R2_" in real_job.filenames[1]
        assert ".fastq.gz" in real_job.filenames[0]
        assert ".fastq.gz" in real_job.filenames[1]

        for input_fn, output_fn in zip(
            [
                (data_path / "sample_b" / "a_R1_.fastq.gz"),
                (data_path / "sample_b" / "a_R2_.fastq.gz"),
            ],
            real_job.filenames,
        ):
            with gzip.GzipFile(output_fn, "r") as op:
                actual = op.read()
            with gzip.GzipFile(input_fn, "r") as op:
                should = op.read()
            assert actual == should

    def test_lane_paired_as_single(self):

        lane = Sample(
            "Sample_a",
            data_path / "sample_b",
            False,
            vid="VA000",
            pairing="paired_as_single",
        )
        assert lane.vid == "VA000"
        temp_job = lane.prepare_input()
        real_job = lane.save_input()
        ppg.run_pipegraph()
        assert not Path(temp_job.filenames[0]).exists()
        assert len(temp_job.filenames) == 1
        assert Path(real_job.filenames[0]).exists()
        assert len(real_job.filenames) == 1
        assert not "_R1_" in real_job.filenames[0]
        assert ".fastq.gz" in real_job.filenames[0]

        should = b""
        for input_fn in [
            (data_path / "sample_b" / "a_R1_.fastq.gz"),
            (data_path / "sample_b" / "a_R2_.fastq.gz"),
        ]:
            with gzip.GzipFile(input_fn, "r") as op:
                should += op.read()
        with gzip.GzipFile(real_job.filenames[0], "r") as op:
            actual = op.read()
        assert actual == should

    def test_lane_paired_missing_R2(self):

        lane = Sample(
            "Sample_a", data_path / "sample_a", False, vid="VA000", pairing="paired"
        )
        with pytest.raises(PairingError):
            lane.prepare_input()

    def test_lane_paired_only_first(self):

        lane = Sample(
            "Sample_a", data_path / "sample_b", False, vid="VA000", pairing="only_first"
        )
        assert lane.vid == "VA000"
        temp_job = lane.prepare_input()
        real_job = lane.save_input()
        ppg.run_pipegraph()
        assert not Path(temp_job.filenames[0]).exists()
        assert len(temp_job.filenames) == 1
        assert Path(real_job.filenames[0]).exists()
        assert len(real_job.filenames) == 1
        assert not "_R1_" in real_job.filenames[0]
        assert ".fastq.gz" in real_job.filenames[0]

        should = b""
        for input_fn in [
            (data_path / "sample_b" / "a_R1_.fastq.gz"),
            # (data_path / "sample_b" / "a_R2_.fastq.gz"),
        ]:
            with gzip.GzipFile(input_fn, "r") as op:
                should += op.read()
        with gzip.GzipFile(real_job.filenames[0], "r") as op:
            actual = op.read()
        assert actual == should

    def test_lane_paired_only_second(self):

        lane = Sample(
            "Sample_a",
            data_path / "sample_b",
            False,
            vid="VA000",
            pairing="only_second",
        )
        assert lane.vid == "VA000"
        temp_job = lane.prepare_input()
        real_job = lane.save_input()
        ppg.run_pipegraph()
        assert not Path(temp_job.filenames[0]).exists()
        assert len(temp_job.filenames) == 1
        assert Path(real_job.filenames[0]).exists()
        assert len(real_job.filenames) == 1
        assert not "_R1_" in real_job.filenames[0]
        assert ".fastq.gz" in real_job.filenames[0]

        should = b""
        for input_fn in [
            # (data_path / "sample_b" / "a_R1_.fastq.gz"),
            (data_path / "sample_b" / "a_R2_.fastq.gz")
        ]:
            with gzip.GzipFile(input_fn, "r") as op:
                should += op.read()
        with gzip.GzipFile(real_job.filenames[0], "r") as op:
            actual = op.read()
        assert actual == should

    def test_pairing_invalid_value(self):
        with pytest.raises(ValueError):
            Sample(
                "Sample_a", data_path / "sample_a", False, pairing="do_what_you_want"
            )
        with pytest.raises(ValueError):
            Sample("Sample_a", data_path / "sample_a", False, pairing=False)
        with pytest.raises(ValueError):
            Sample("Sample_a", data_path / "sample_a", False, pairing=None)
        with pytest.raises(ValueError):
            Sample("Sample_a", data_path / "sample_a", False, pairing=[5])

    def test_lane_raises_on_pe_as_se(self):
        lane = Sample("Sample_a", data_path / "sample_b", False)
        with pytest.raises(PairingError):
            lane.prepare_input()

    def test_lane_with_job_generating_fastq(self):
        def gen_fastq(fn):
            with open(fn, "wb") as op:
                op.write(b"@shu\nAGTC\n+\nasdf")

        job = FileGeneratingJob("input.fastq", gen_fastq)

        lane = Sample("Sample_a", job, False, vid="VA000")
        assert lane.vid == "VA000"
        temp_job = lane.prepare_input()
        assert job in temp_job.prerequisites
        real_job = lane.save_input()
        ppg.run_pipegraph()
        assert not Path(temp_job.filenames[0]).exists()
        assert Path(real_job.filenames[0]).exists()
        with gzip.GzipFile(real_job.filenames[0], "r") as op:
            lines = op.readlines()
        assert len(lines) == 4

    def test_align(self, local_store):
        import json
        import gzip

        class FakeGenome:
            name = "FakeGenome"

            def build_index(self, aligner, fasta_to_use=None, gtf_to_use=None):
                return ppg.FileGeneratingJob(
                    "fake_index", lambda: Path("fake_index").write_text("hello")
                )

        class FakeAligner:
            name = "FakeAligner"

            def align_job(
                self,
                input_fastq,
                paired_end_filename,
                index_basename,
                output_bam_filename,
                parameters,
            ):
                def align():
                    with open(output_bam_filename, "w") as op:
                        json.dump(
                            [
                                open(input_fastq).read(200),
                                open(paired_end_filename).read(200)
                                if paired_end_filename
                                else "",
                                index_basename,
                                str(parameters),
                            ],
                            op,
                        )
                    with open(str(output_bam_filename) + ".bai", "w") as op:
                        op.write("Done")

                return ppg.MultiFileGeneratingJob(
                    [output_bam_filename, str(output_bam_filename) + ".bai"], align
                )

        aligner = FakeAligner()
        lane = Sample(
            "Sample_a", data_path / "sample_b", False, vid="VA000", pairing="paired"
        )
        genome = FakeGenome()
        params = {"shu": 123}
        aligned_lane = lane.align(aligner, genome, params)
        ppg.run_pipegraph()
        assert Path("fake_index").exists()
        assert Path("fake_index").read_text() == "hello"
        assert aligned_lane.load()[0].filenames[0].endswith(lane.name + ".bam")
        assert aligned_lane.load()[0].filenames[1].endswith(lane.name + ".bam.bai")
        assert Path(aligned_lane.load()[0].filenames[0]).exists()
        with open(aligned_lane.load()[0].filenames[0]) as op:
            actual = json.load(op)
        with gzip.GzipFile(data_path / "sample_b" / "a_R1_.fastq.gz") as op:
            should_0 = op.read(200).decode("utf-8")
        with gzip.GzipFile(data_path / "sample_b" / "a_R2_.fastq.gz") as op:
            should_1 = op.read(200).decode("utf-8")

        assert actual[0] == should_0
        assert actual[1] == should_1
        assert actual[2] == "fake_index"
        assert actual[3] == str(params)

    def test_from_url(self):
        import requests_mock

        with requests_mock.mock() as m:
            url = "https://www.imt.uni-marburg.de/sample.fastq.gz"
            m.get(url, text="hello_world")
            o = FASTQsFromURLs(url)
            ppg.run_pipegraph()
            assert len(o.target_files) == 1
            assert len(o.dependencies) == 2
            assert Path(o.dependencies[0].filenames[0]).read_text() == "hello_world"
            assert Path(o.dependencies[0].filenames[1]).read_text() == url
            assert o() == [(Path(o.dependencies[0].filenames[0]).absolute(),)]

    def test_from_url_paired(self):
        import requests_mock

        with requests_mock.mock() as m:
            url1 = "https://www.imt.uni-marburg.de/sample_R1_.fastq.gz"
            url2 = "https://www.imt.uni-marburg.de/sample_R2_.fastq.gz"
            m.get(url1, text="hello_world1")
            m.get(url2, text="hello_world2")
            o = FASTQsFromURLs([url1, url2])
            ppg.run_pipegraph()
            assert len(o.target_files) == 2
            assert len(o.dependencies) == 3
            assert Path(o.dependencies[0].filenames[0]).read_text() == "hello_world1"
            assert "_R1_.fastq.gz" in o.dependencies[0].filenames[0]
            assert Path(o.dependencies[0].filenames[1]).read_text() == url1
            assert Path(o.dependencies[1].filenames[0]).read_text() == "hello_world2"
            assert Path(o.dependencies[1].filenames[1]).read_text() == url2
            assert o() == [
                (
                    Path(o.dependencies[0].filenames[0]).absolute(),
                    Path(o.dependencies[1].filenames[0]).absolute(),
                )
            ]

    def test_from_url_detects_404(self):

        with requests_mock.mock() as m:
            url = "https://www.imt.uni-marburg.de/sample.fastq.gz"
            m.get(url, text="hello_world", status_code=404)
            o = FASTQsFromURLs(url)
            with pytest.raises(ppg.RuntimeError):
                ppg.run_pipegraph()
            assert "Error return" in str(o.dependencies[0].exception)
            assert url in str(o.dependencies[0].exception)

    def test_fastqs_from_err(self):
        with requests_mock.mock() as m:
            m.get(
                "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=ERR2223563&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes",
                text="""run_accession	fastq_ftp	fastq_md5	fastq_bytes
ERR2223563	ftp.sra.ebi.ac.uk/vol1/fastq/ERR222/003/ERR2223563/ERR2223563_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR222/003/ERR2223563/ERR2223563_2.fastq.gz	0e29c053bcd31072c8bed9eddece1cec;5d848b65379c195fe158a5d7324b4a18	1170312089;1246298835""",
            )
            o = FASTQsFromAccession("ERR2223563")
            print(o.urls)
            assert (
                o.urls[0]
                == "http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR222/003/ERR2223563/ERR2223563_1.fastq.gz"
            )
            assert (
                o.urls[1]
                == "http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR222/003/ERR2223563/ERR2223563_2.fastq.gz"
            )
