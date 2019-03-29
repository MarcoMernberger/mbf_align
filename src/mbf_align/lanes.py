# forwards for compatibility with old chipseq code

from .raw import Sample
import pypipegraph as ppg
import pysam
from pathlib import Path
import pandas as pd
from dppd import dppd
import dppd_plotnine  # noqa:F401 -

dppd, X = dppd()


class AlignedSample:
    def __init__(self, name, alignment_job, genome, is_paired, vid, result_dir=None):
        """
        Create an aligned sample from a BAM producing job.
        See Sample.align()

        Parameters:
            alignment_job FileGeneratingJob, FileInvariant, str, pathlib.Path
            Where does the BAM come from?
            str and Path get's converted into a FileInvariant
        """

        self.alignment_job, self.index_job, bam_name, index_fn = self._parse_alignment_job_input(
            alignment_job
        )
        self.name = name
        ppg.util.assert_uniqueness_of_object(self)
        self.result_dir = (
            Path(result_dir)
            if result_dir
            else (Path("results") / "aligned" / self.name)
        )
        self.result_dir.mkdir(exist_ok=True, parents=True)
        self.genome = genome
        self.is_paired = is_paired
        self.vid = vid
        self.bam_filename = bam_name
        self.index_filename = index_fn
        self.register_qc()

    def _parse_alignment_job_input(self, alignment_job):
        if isinstance(alignment_job, (str, Path)):
            alignment_job = ppg.FileInvariant(alignment_job)
        if not isinstance(alignment_job, (ppg.FileInvariant, ppg.FileGeneratingJob)):
            raise ValueError(
                "alignment_job must be a ppg.FileGeneratingJob or FileChecksumInvariant"
                "was %s" % (type(alignment_job))
            )
        bam_name = None
        bai_name = None
        for fn in alignment_job.filenames:
            if str(fn).endswith(".bam"):
                if bam_name is None:
                    bam_name = str(fn)
                else:
                    raise ValueError(
                        "Job passed to AlignedSample had multiple .bam filenames"
                    )
            elif str(fn).endswith(".bai"):
                if bai_name is None:
                    bai_name = str(fn)
                else:
                    raise ValueError(
                        "Job passed to AlignedSample had multiple .bai filenames"
                    )

        if bam_name is None:
            raise ValueError("Job passed to AlignedSample had no .bam filenames")

        if isinstance(alignment_job, ppg.MultiFileGeneratingJob):
            if bai_name is None:
                index_fn = bam_name + ".bai"
                index_job = ppg.FileGeneratingJob(
                    index_fn, self._index(bam_name, index_fn)
                )
                index_job.depends_on(alignment_job)

            else:
                index_fn = bai_name
                index_job = alignment_job

        elif isinstance(alignment_job, ppg.FileGeneratingJob):
            index_fn = bam_name + ".bai"
            index_job = ppg.FileGeneratingJob(index_fn, self._index(bam_name, index_fn))
            index_job.depends_on(alignment_job)
        elif isinstance(alignment_job, ppg.FileInvariant):
            index_fn = bam_name + ".bai"
            if Path(index_fn).exists():
                index_job = ppg.FileInvariant(index_fn)
            else:
                cache_dir = Path(ppg.util.global_pipegraph.cache_folder) / "bam_indices"
                cache_dir.mkdir(exist_ok=True)
                index_fn = cache_dir / (Path(bam_name).name + ".bai")
                index_job = ppg.FileGeneratingJob(
                    index_fn, self._index(bam_name, index_fn)
                )
                index_job.depends_on(alignment_job)
        else:
            raise NotImplementedError("Should not happe / covered by earlier if")
        return alignment_job, index_job, bam_name, index_fn

    def load(self):
        return self.alignment_job, self.index_job

    def _index(self, input_fn, output_fn):
        def do_index():
            print(str(input_fn), Path(input_fn).exists())
            print(str(output_fn), Path(output_fn).exists())
            pysam.index(str(input_fn), str(output_fn))

        return do_index

    def get_bam(self):
        import multiprocessing

        return pysam.Samfile(
            self.bam_filename,
            index_filename=self.index_filename,
            threads=multiprocessing.cpu_count(),
        )

    def get_unique_aligned_bam(self):
        """Deprecated compability with older pipeline"""
        return self.get_bam()

    def _parse_idxstat(self):
        by_chr = self.get_bam().get_index_statistics()
        mapped = 0
        unmapped = 0
        for record in by_chr:
            mapped += record.mapped
            unmapped += record.unmapped
        return mapped, unmapped

    def mapped_reads(self):
        """How many mapped entrys are in the bam?"""
        return self._parse_idxstat()[0]

    def unmapped_reads(self):
        """How many unmapped entrys are in the bam?"""
        return self._parse_idxstat()[1]

    def register_qc(self):
        self.register_qc_complexity()

    def register_qc_complexity(self):
        from mbf_qualitycontrol import register_qc, QCCallback

        output_filename = self.result_dir / "complexity.png"

        def build():
            def calc():
                import mbf_bam

                counts = mbf_bam.calculate_duplicate_distribution(
                    self.bam_filename
                )  # , self.index_filename
                return pd.DataFrame(
                    {
                        "source": self.name,
                        "Repetition count": list(counts.keys()),
                        "Count": list(counts.values()),
                    }
                )

            def plot(df):
                unique_count = df["Count"].sum()
                total_count = (df["Count"] * df["Repetition count"]).sum()
                severity = "severe"
                pcb = float(unique_count) / total_count
                if pcb >= 0.9:
                    severity = "none"
                elif pcb >= 0.8:
                    severity = "mild"
                elif pcb >= 0.5:
                    severity = "moderate"
                title = (
                    "Genomic positions with repetition count reads\nTotal read count: %i\nPCR Bottleneck coefficient: %.2f (%s)"
                    % (total_count, pcb, severity)
                )
                return (
                    dppd(df)
                    .p9()
                    .add_point("Repetition count", "Count")
                    .add_line("Repetition count", "Count")
                    .scale_y_continuous(trans="log2")
                    .title(title)
                    .pd
                )

            return ppg.PlotJob(output_filename, calc, plot).depends_on(self.load)

        register_qc(output_filename, QCCallback(build))
        return output_filename


__all__ = [Sample, AlignedSample]
