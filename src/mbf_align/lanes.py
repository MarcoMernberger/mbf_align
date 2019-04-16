# forwards for compatibility with old chipseq code

from .raw import Sample
import pypipegraph as ppg
import pysam
from pathlib import Path
import pandas as pd
import collections
from dppd import dppd
import dppd_plotnine  # noqa:F401 -
from mbf_qualitycontrol import register_qc, QCCallback, get_qc

dp, X = dppd()


class AlignedSample:
    def __init__(
        self, name, alignment_job, genome, is_paired, vid, result_dir=None, aligner=None
    ):
        """
        Create an aligned sample from a BAM producing job.
        See Sample.align()

        Parameters:
            alignment_job FileGeneratingJob, FileInvariant, str, pathlib.Path
            Where does the BAM come from?
            str and Path get's converted into a FileInvariant
        """

        self.name = name
        ppg.util.assert_uniqueness_of_object(self)
        self.alignment_job, self.index_job, bam_name, index_fn = self._parse_alignment_job_input(
            alignment_job
        )
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
        self.aligner = aligner
        self.register_qc()

    def __hash__(self):
        return hash(self.__class__.__name__ + self.name)

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
                    index_fn = str(fn)
                    bai_name = index_fn
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
                index_fn = cache_dir / (self.name + "_" + Path(bam_name).name + ".bai")
                index_job = ppg.FileGeneratingJob(
                    index_fn, self._index(bam_name, index_fn)
                )
                index_job.depends_on(alignment_job)
        else:
            raise NotImplementedError("Should not happe / covered by earlier if")
        return alignment_job, index_job, Path(bam_name), Path(index_fn)

    def load(self):
        return self.alignment_job, self.index_job

    def _index(self, input_fn, output_fn):
        def do_index():
            pysam.index(str(Path(input_fn).absolute()), str(Path(output_fn).absolute()))

        return do_index

    def get_bam(self):
        import multiprocessing

        return pysam.Samfile(
            self.bam_filename,
            index_filename=str(self.index_filename),
            threads=multiprocessing.cpu_count(),
        )

    def get_bam_names(self):
        """Retrieve the bam filename and index name as strings"""
        return (str(self.bam_filename), str(self.index_filename))

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

    def get_alignment_stats(self):
        if self.aligner is not None and hasattr(self.aligner, "get_alignment_stats"):
            return self.aligner.get_alignment_stats(Path(self.bam_filename))
        else:
            with self.get_bam() as f:
                return {"Mapped": f.mapped, "Unmapped": f.unmapped}

    def register_qc(self):
        self.register_qc_complexity()
        self.register_qc_gene_strandedness()
        self.register_qc_biotypes()
        self.register_qc_alignment_stats()
        self.register_qc_subchromosomal()
        self.register_qc_splicing()

    def register_qc_complexity(self):

        output_filename = self.result_dir / "complexity.png"

        def build():
            def calc():
                import mbf_bam

                counts = mbf_bam.calculate_duplicate_distribution(
                    str(self.bam_filename), str(self.index_filename)
                )
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
                pcb = float(unique_count) / total_count
                if pcb >= 0.9:  # pragma: no cover
                    severity = "none"
                elif pcb >= 0.8:  # pragma: no cover
                    severity = "mild"
                elif pcb >= 0.5:  # pragma: no cover
                    severity = "moderate"
                else:
                    severity = "severe"
                title = (
                    "Genomic positions with repetition count reads\nTotal read count: %i\nPCR Bottleneck coefficient: %.2f (%s)"
                    % (total_count, pcb, severity)
                )
                return (
                    dp(df)
                    .p9()
                    .theme_bw()
                    .add_point("Repetition count", "Count")
                    .add_line("Repetition count", "Count")
                    .scale_y_continuous(trans="log2")
                    .title(title)
                    .pd
                )

            return (
                ppg.PlotJob(output_filename, calc, plot)
                .depends_on(self.load)
                .use_cores(-1)
            )

        register_qc(output_filename, QCCallback(build))
        return output_filename

    def register_qc_gene_strandedness(self):
        from mbf_genomics.genes.anno_tag_counts import _IntervalStrategy

        class IntervalStrategyExonIntronClassification(_IntervalStrategy):
            """For QC purposes, defines all intron/exon intervals tagged
            with nothing but intron/exon

            See mbf_align.lanes.AlignedLane.register_qc_gene_strandedness
            
            """

            def _get_interval_tuples_by_chr(self, genome):
                from mbf_nested_intervals import IntervalSet

                coll = {chr: [] for chr in genome.get_chromosome_lengths()}
                ii = 0
                for g in genome.genes.values():
                    exons = g.exons_overlapping
                    if len(exons[0]) == 0:
                        exons = g.exons_merged
                    for start, stop in zip(*exons):
                        coll[g.chr].append(
                            (start, stop, 0b0101 if g.strand == 1 else 0b0110)
                        )
                    for start, stop in zip(*g.introns):
                        coll[g.chr].append(
                            (start, stop, 0b1001 if g.strand == 1 else 0b1010)
                        )
                result = {}
                for chr, tups in coll.items():
                    iset = IntervalSet.from_tuples_with_id(tups)
                    # iset = iset.merge_split()
                    iset = iset.merge_hull()
                    if iset.any_overlapping():
                        raise NotImplementedError("Should not be reached")
                    result[chr] = []
                    for start, stop, ids in iset.to_tuples_with_id():
                        ids = set(ids)
                        if len(ids) == 1:
                            id = list(ids)[0]
                            if id == 0b0101:
                                tag = "exon"
                                strand = +1
                            elif id == 0b0110:
                                tag = "exon"
                                strand = -1
                            elif id == 0b1001:
                                tag = "intron"
                                strand = +1
                            elif id == 0b1010:
                                tag = "intron"
                                strand = -1
                            else:  # pragma: no cover
                                raise NotImplementedError("Should not be reached")
                        else:
                            down = 0
                            for i in ids:
                                down |= i
                            if down & 0b1100 == 0b1100:
                                tag = "both"
                            elif down & 0b0100 == 0b0100:
                                tag = "exon"
                            else:
                                tag = "intron"
                            if down & 0b11 == 0b11:
                                tag += "_undecidable"
                                strand = (
                                    1
                                )  # doesn't matter, but must be one or the other
                            elif down & 0b01:
                                strand = 1
                            else:
                                strand -= 1

                        result[chr].append((tag, strand, [start], [stop]))
                return result

        output_filename = self.result_dir / "strandedness.png"

        def build():
            def calc():
                from mbf_genomics.genes.anno_tag_counts import IntervalStrategyGene
                from mbf_bam import count_reads_stranded

                interval_strategy = IntervalStrategyExonIntronClassification()
                intervals = interval_strategy._get_interval_tuples_by_chr(self.genome)

                bam_filename, bam_index_name = self.get_bam_names()
                forward, reverse = count_reads_stranded(
                    bam_filename,
                    bam_index_name,
                    intervals,
                    IntervalStrategyGene()._get_interval_tuples_by_chr(self.genome),
                    each_read_counts_once=True,
                )
                result = {"what": [], "count": [], "sample": self.name}
                for k in forward.keys() | reverse.keys():
                    if k.endswith("_undecidable"):
                        result["what"].append(k)
                        result["count"].append(forward.get(k, 0) + reverse.get(k, 0))
                    elif not k.startswith("_"):
                        result["what"].append(k + "_correct")
                        result["count"].append(forward.get(k, 0))
                        result["what"].append(k + "_reversed")
                        result["count"].append(reverse.get(k, 0))
                    elif k == "_outside":
                        result["what"].append("outside")
                        result["count"].append(forward.get(k, 0))

                return pd.DataFrame(result)

            def plot(df):
                return (
                    dp(df)
                    .mutate(
                        what=pd.Categorical(
                            df["what"],
                            [
                                "exon_correct",
                                "exon_reversed",
                                "exon_undecidable",
                                "intron_correct",
                                "intron_reversed",
                                "intron_undecidable",
                                "both_correct",
                                "both_reversed",
                                "both_undecidable",
                                "outside",
                            ],
                        )
                    )
                    .p9()
                    .add_bar("sample", "count", fill="what", position="dodge")
                    .turn_x_axis_labels()
                    .pd
                )

            return (
                ppg.PlotJob(output_filename, calc, plot)
                .depends_on(self.load)
                .use_cores(-1)
            )

        register_qc(output_filename, QCCallback(build))
        return output_filename

    def register_qc_biotypes(self):
        output_filename = self.result_dir / "reads_per_biotype.png"

        def build():
            from mbf_genomics.genes import Genes
            from mbf_genomics.genes.anno_tag_counts import GeneUnstranded

            genes = Genes(self.genome)
            anno = GeneUnstranded(self)
            anno.register_qc = lambda self: None  # can't register qc when doing qc...

            def plot(output_filename):
                return (
                    dp(genes.df)
                    .groupby("biotype")
                    .summarize((anno.columns[0], lambda x: x.sum(), "read count"))
                    .mutate(sample=self.name)
                    .p9()
                    .theme_bw()
                    .annotation_stripes()
                    .add_bar("biotype", "read count", stat="identity")
                    # .turn_x_axis_labels()
                    .coord_flip()
                    .title(self.name)
                    .render(
                        output_filename,
                        width=6,
                        height=2 + len(genes.df.biotype.unique()) * 0.25,
                    )
                )

            return (
                ppg.FileGeneratingJob(output_filename, plot)
                .depends_on(genes.add_annotator(anno))
                .use_cores(-1)
            )

        register_qc(output_filename, QCCallback(build))

    def register_qc_alignment_stats(self):
        output_filename = self.result_dir / ".." / "alignment_statistics.png"
        try:
            q = get_qc(output_filename)
        except KeyError:

            class AlignmentStatQC:
                def __init__(self):
                    self.lanes = set()

                def get_qc_job(self):
                    def calc():
                        parts = []
                        for l in self.lanes:
                            p = l.get_alignment_stats()
                            parts.append(
                                pd.DataFrame(
                                    {
                                        "what": list(p.keys()),
                                        "count": list(p.values()),
                                        "sample": l.name,
                                    }
                                )
                            )
                        return pd.concat(parts)

                    def plot(df):
                        return (
                            dp(df)
                            .p9()
                            .theme_bw()
                            .annotation_stripes()
                            .add_bar(
                                "sample",
                                "count",
                                fill="what",
                                position="stack",
                                stat="identity",
                            )
                        )

                    return ppg.PlotJob(output_filename, calc, plot).depends_on(
                        [x.load() for x in self.lanes]
                    )

            q = AlignmentStatQC()
            register_qc(output_filename, q)
        q.lanes.add(self)

    def register_qc_subchromosomal(self):
        """Subchromosom distribution plot - good to detect amplified regions
        or ancient virus awakening"""
        import mbf_genomics

        output_filename = self.result_dir / "subchromosomal_distribution.png"

        class IntervalStrategyWindows(
            mbf_genomics.genes.anno_tag_counts._IntervalStrategy
        ):
            """For QC purposes, spawn all chromosomes with 
            windows of the definied size

            See mbf_align.lanes.AlignedLane.register_qc_subchromosomal
            
            """

            def __init__(self, window_size):
                self.window_size = window_size

            def _get_interval_tuples_by_chr(self, genome):
                result = {}
                for chr, length in genome.get_chromosome_lengths().items():
                    result[chr] = []
                    for ii in range(0, length, self.window_size):
                        result[chr].append(
                            ("%s_%i" % (chr, ii), 0, [ii], [ii + self.window_size])
                        )
                return result

        def build():
            def calc():
                from mbf_bam import count_reads_unstranded

                interval_strategy = IntervalStrategyWindows(250_000)
                intervals = interval_strategy._get_interval_tuples_by_chr(self.genome)

                bam_filename, bam_index_name = self.get_bam_names()
                counts = count_reads_unstranded(
                    bam_filename,
                    bam_index_name,
                    intervals,
                    intervals,
                    each_read_counts_once=True,
                )
                true_chromosomes = set(self.genome.get_true_chromosomes())
                result = {"chr": [], "window": [], "count": []}
                for key, count in counts.items():
                    if not key.startswith("_"):
                        chr, window = key.split("_", 2)
                    if chr in true_chromosomes:
                        window = int(window)
                        result["chr"].append(chr)
                        result["window"].append(window)
                        result["count"].append(count)
                return pd.DataFrame(result)

            def plot(df):
                return (
                    dp(df)
                    .p9()
                    .theme_bw()
                    .add_line("window", "count")
                    .scale_y_log10()
                    .facet_wrap("chr", scales="free", ncol=1)
                    .title(self.name)
                    .render(
                        output_filename,
                        width=6,
                        height=2 + len(df["chr"].unique()) * 0.75,
                    )
                )

            return (
                ppg.PlotJob(output_filename, calc, plot)
                .depends_on(self.load())
                .use_cores(-1)
            )

        register_qc(output_filename, QCCallback(build))

    def register_qc_splicing(self):
        """How many reads were spliced? How many of those splices were known splice sites,
        how many were novel"""
        output_filename = self.result_dir / "splice_sites.png"

        def build():
            def calc():
                from mbf_bam import count_introns

                bam_filename, bam_index_name = self.get_bam_names()
                counts_per_chromosome = count_introns(bam_filename, bam_index_name)
                known_splice_sites_by_chr = {
                    chr: set() for chr in self.genome.get_chromosome_lengths()
                }
                for gene in self.genome.genes.values():
                    for start, stop in zip(*gene.introns):
                        known_splice_sites_by_chr[gene.chr].add((start, stop))
                total_counts = collections.Counter()
                known_count = 0
                unknown_count = 0
                for chr, counts in counts_per_chromosome.items():
                    for k, v in counts.items():
                        if k[0] == 0xFFFFFFFF:
                            intron_counts = 0xFFFFFFFF - k[1]
                            total_counts[intron_counts] += v
                        else:
                            if k in known_splice_sites_by_chr[chr]:
                                known_count += v
                            else:
                                unknown_count += v
                result = {"side": [], "x": [], "count": []}
                result["side"].append("splice sites")
                result["x"].append("unknown")
                result["count"].append(unknown_count)
                result["side"].append("splice sites")
                result["x"].append("known")
                result["count"].append(known_count)

                for x, count in total_counts.items():
                    result["side"].append("reads with x splices")
                    result["x"].append(x)
                    result["count"].append(count)

                return pd.DataFrame(result)

            def plot(df):
                return (
                    dp(df)
                    .p9()
                    .theme_bw()
                    .add_bar("x", "count", stat="identity")
                    .facet_wrap("side", scales="free", ncol=1)
                    .title(self.name)
                    .theme(panel_spacing_y=0.2)
                    .render(output_filename)
                )

            return (
                ppg.PlotJob(output_filename, calc, plot)
                .depends_on(self.load())
                .use_cores(-1)
            )

        register_qc(output_filename, QCCallback(build))


__all__ = [Sample, AlignedSample]
