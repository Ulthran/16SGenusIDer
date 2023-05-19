import logging
import os
import shutil
from pathlib import Path
from .parsers import parse_fasta


class OutputFile:
    """
    Base class to define a file in the output directory and associated operations
    """

    def __init__(self, fp: Path) -> None:
        self.fp = fp
        self.fn = fp.name

    def get(self) -> Path:
        return self.fp


class QueryFile(OutputFile):
    """
    Holds and creates a the query sequence in fasta format
    :param str seq: either a filepath to the query in fasta format or the sequence itself as a string
    :param root_fp Path: the Path to the output root
    """

    def __init__(self, seq: str, root_fp: Path) -> None:
        try:
            if Path(seq).exists():
                self.fp = Path(seq)
                with open(self.fp) as f:
                    counter = 0
                    for header_str, seq_str in parse_fasta(f):
                        self.seq = seq_str
                        counter += 1
                        if counter > 1:
                            logging.warn("Multiple sequences in query file, only the last entry will be used")
            else:
                # Assuming seq is the literal sequence, not a file
                self.fp = root_fp / "query.fasta"
                self.seq = seq
                self.write_query()
        except OSError as e:
            logging.info(f"--seq arg is too long to be a path, interpreting as literal sequence")
            self.fp = root_fp / "query.fasta"
            self.seq = seq
            self.write_query()
        
        super().__init__(self.fp)
        self.desc = "UNKNOWN"
        
    def write_query(self):
        with open(self.fp, "w") as f:
            f.write(f">{self.desc}\n")
            f.write(f"{self.seq}\n")


class NearestSeqsFile(OutputFile):
    """
    Holds and operates on the nearest seqs file from vsearch and a reduced version if necessary
    """

    def __init__(self, fp: Path, query: QueryFile) -> None:
        super().__init__(fp)
        self.query = query
        self.reduced = False
    
    def get(self) -> Path:
        if self.fp.exists():
            self.reduce_subtree()
        return super().get()
    
    def reduce_subtree(self):
        if not self.reduced:
            ids = []
            temp_fp = self.fp.parents[0] / "reduced_alignment.fasta"
            with open(self.fp) as f_in, open(temp_fp, "w") as f_out:
                for id, seq in parse_fasta(f_in):
                    if len(id.strip()) > 2:
                        ids.append(id)
                        f_out.write(f">{id}\n")
                        f_out.write(f"{seq}\n")
                    if len(ids) >= 50:
                        break

                logging.debug(f"Reduced accessions list: {str(ids)}")

            with open(temp_fp, "a+") as f:
                f.write(f">{self.query.desc}\n")
                f.write(f"{self.query.seq}\n")

            os.remove(self.fp)
            os.rename(temp_fp, self.fp)


class ProbabilityFile(OutputFile):
    """
    Holds and writes the probability output file
    """

    def __init__(self, fp: Path) -> None:
        super().__init__(fp)
    
    def write_probs(self, probs: dict, header: str = "", write_mode: str = "a+"):
        with open(self.fp, write_mode) as f:
            f.write(f"\n{header}\n\n")
            for s, p in probs.items():
                if p > 0.0001:
                    f.write(f"{s.split(' ')[0]}\t{round(p * 100, 5)}\n")


class OutputDir:
    """
    Controller for all of GenusFinder's output files
    Interface directly with internal objects\n
    i.e. call QueryFile's get() method directly
    """

    def __init__(self, fp: Path, seq: str, overwrite: bool) -> None:
        self.root_fp = Path(fp)

        self.overwriteQ = overwrite
        if overwrite and self.root_fp.exists():
            logging.warning(f"Found existing output dir {self.root_fp}, overwriting...")
            shutil.rmtree(self.root_fp)

        os.makedirs(self.root_fp, exist_ok=True)

        self.query = QueryFile(seq, self.root_fp)
        self.nearest_seqs = NearestSeqsFile(self.root_fp / "nearest_seqs.fasta", self.query)
        self.nearest_seqs_aligned = OutputFile(self.root_fp / "nearest_seqs_aligned.fasta")
        self.bootstraps = OutputFile(self.root_fp / "RAxML_bootstrap.subtree1")
        self.base_tree = OutputFile(self.root_fp / "RAxML_bestTree.subtree2")
        self.bootstrapped_tree = OutputFile(self.root_fp / "RAxML_bipartitions.final")
        self.combined_alignment = OutputFile(self.root_fp / "combined_alignment.fasta")
        self.combined_tree = OutputFile(self.root_fp / "RAxML_bestTree.combined")
        self.probs = ProbabilityFile(self.root_fp / "probabilities.tsv")
