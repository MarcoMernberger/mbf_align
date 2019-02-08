# forwards for compatibility with old chipseq code

from .raw import Sample
import pypipegraph as ppg
import pysam
from pathlib import Path


class AlignedSample:
    def __init__(self, name, alignment_job, genome, is_paired, vid):
        """
        Create an aligned sample from a BAM producing job.
        See Sample.align()

        Parameters:
            alignment_job FileGeneratingJob, FileInvariant, str, pathlib.Path
            Where does the BAM come from?
            str and Path get's converted into a FileInvariant
        """
        if isinstance(alignment_job, (str, Path)):
            alignment_job = ppg.FileInvariant(alignment_job)
        if not isinstance(alignment_job, (ppg.FileInvariant, ppg.FileGeneratingJob)):
            raise ValueError(
                "alignment_job must be a ppg.FileGeneratingJob or FileChecksumInvariant"
                "was %s" % (type(alignment_job))
            )
        bam_name = alignment_job.filenames[0]
        if isinstance(alignment_job, ppg.MultiFileGeneratingJob):
            index_fn = alignment_job.filenames[1]
            if index_fn != bam_name + ".bai": # pragma: no cover
                raise ValueError(
                    "MultiFileGeneratingJob as alignment_job "
                    "had not .bam and .bam.bai as first two filenames"
                )
            index_job = alignment_job
        elif isinstance(alignment_job, ppg.FileGeneratingJob):
            index_fn = alignment_job.job_id + ".bai"
            index_job = ppg.FileGeneratingJob(
                index_fn, self._index(alignment_job.job_id, index_fn)
            )
            index_job.depends_on(alignment_job)
        elif isinstance(alignment_job, ppg.FileInvariant):
            index_fn = alignment_job.job_id + ".bai"
            if Path(index_fn).exists():
                index_job = ppg.FileInvariant(index_fn)
            else:
                cache_dir = Path(ppg.util.global_pipegraph.cache_folder) / "bam_indices"
                cache_dir.mkdir(exist_ok=True)
                index_fn = cache_dir / (Path(alignment_job.job_id).name + ".bai")
                index_job = ppg.FileGeneratingJob(
                    index_fn, self._index(alignment_job.job_id, index_fn)
                )
                index_job.depends_on(alignment_job)

        self.alignment_job = alignment_job
        self.index_job = index_job
        self.name = name
        ppg.util.assert_uniqueness_of_object(self)
        self.genome = genome
        self.is_paired = is_paired
        self.vid = vid
        self.bam_filename = bam_name
        self.index_filename = index_fn

    def load(self):
        return self.alignment_job, self.index_job

    def _index(self, input_fn, output_fn):
        def do_index():
            pysam.index(str(input_fn), str(output_fn))

        return do_index

    def get_bam(self):
        import multiprocessing
        return pysam.Samfile(self.bam_filename, 
                             index_filename=self.index_filename,
                             threads=multiprocessing.cpu_count())

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


__all__ = [Sample, AlignedSample]
