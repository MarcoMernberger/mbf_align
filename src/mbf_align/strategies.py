from pathlib import Path
import pypipegraph as ppg


def build_fastq_strategy(input_strategy):
    """smartly find fastq input files

    Parameters
    ----------
        input_strategy - varied
            FASTQsFrom* object - return the object
            str/Path - file: treat as single fastq(.gz) file
            str/Path - folder: treat as folder
            pypipegraph job: extract filenames
            list - recurvively apply build_fastq_strategy
        """

    if isinstance(input_strategy, _FASTQsBase):
        return input_strategy
    elif isinstance(input_strategy, list):
        return _FASTQsJoin([build_fastq_strategy(x) for x in input_strategy])
    elif isinstance(input_strategy, str) or isinstance(input_strategy, Path):
        p = Path(input_strategy)
        if p.is_dir():
            input_strategy = FASTQsFromFolder(p)
        else:
            input_strategy = FASTQsFromFile(p)
    elif isinstance(input_strategy, ppg.FileGeneratingJob):
        input_strategy = FASTQsFromJob(input_strategy)
    else:
        raise ValueError("Could not parse input_strategy")
    return input_strategy


class _FASTQsBase:

    """
    All FASTQs* strategies
    return a list of tuples of (read1, read2, read3, ...)
    filenames upon being called.

    """

    def __call__(self):
        raise NotImplementedError()  # pragma: no cover

    def _combine_r1_r2(self, forward, reverse):
        results = []
        if not forward:
            raise ValueError(f"No _R1_*.fastq*  files found in {self}")
        if not reverse:
            results.extend(zip(forward))
        if reverse:
            expected_reverse = [Path(str(f).replace("_R1_", "_R2_")) for f in forward]
            if expected_reverse != reverse:
                raise ValueError(
                    f"Error pairing forward/reverse files.\nF:{forward}\nR:{reverse}\nE:{expected_reverse}"
                )
            results.extend(zip(forward, reverse))
        return results


class _FASTQsJoin(_FASTQsBase):
    """join files from multiple strategies"""

    def __init__(self, strategies):
        self.strategies = strategies

    def __call__(self):
        res = []
        for s in self.strategies:
            res.extend(s())
        return res


class FASTQsFromFile(_FASTQsBase):
    """Use as single file (or two for paired end) as input"""

    def __init__(self, r1_filename, r2_filename=None):
        self.r1_filename = Path(r1_filename)
        self.r2_filename = Path(r2_filename) if r2_filename else None
        if not self.r1_filename.exists():
            raise IOError(f"file {self.r1_filename} not found")
        if self.r2_filename and not self.r2_filename.exists():
            raise IOError(f"file {self.r2_filename} not found")

    def __call__(self):
        if self.r2_filename:
            return [(self.r1_filename.resolve(), self.r2_filename.resolve())]
        else:
            return [(self.r1_filename.resolve(),)]


class FASTQsFromFolder(_FASTQsBase):
    """Discover fastqs(.gz) in a single folder,
    with automatic pairing based on R1/R2
    """

    def __init__(self, folder):
        self.folder = Path(folder)
        if not self.folder.exists():
            raise IOError(f"folder {self.folder} not found")
        if not any(self.folder.glob("*.fastq.gz")) and not any(
            self.folder.glob("*.fastq")
        ):
            raise ValueError(f"No *.fastq or *.fastq.gz in {self.folder} found")

    def __call__(self):
        fastqs = [x.resolve() for x in self.folder.glob("*.fastq*")]
        forward = sorted([x for x in fastqs if "_R1_" in x.name])
        reverse = sorted([x for x in fastqs if "_R2_" in x.name])
        if not forward and not reverse and fastqs:  # no R1 or R2, but fastqs present
            return sorted(zip(fastqs))
        else:
            return self._combine_r1_r2(forward, reverse)

    def __str__(self):
        return f"FASTQsFromFolder({self.folder})"


class FASTQsFromJob(_FASTQsBase):
    def __init__(self, job):
        self.dependencies = job
        self.job = job

    def __call__(self):
        forward = [Path(x).resolve() for x in self.job.filenames if "_R1_" in x]
        reverse = [Path(x).resolve() for x in self.job.filenames if "_R2_" in x]
        if forward or reverse:
            return self._combine_r1_r2(forward, reverse)
        else:  # no R1, R2, assume single end
            return [(Path(p).resolve(),) for p in self.job.filenames]

    def __str__(self):
        return f"FASTQsFromJob({self.job})"  # pragma: no cover
