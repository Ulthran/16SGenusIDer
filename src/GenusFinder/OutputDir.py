import logging
import os
import shutil
from pathlib import Path


class OutputDir:
    """
    Controller for all of GenusFinder's output files
    """

    def __init__(self, fp: Path, seq: str, overwrite: bool) -> None:
        self.root_fp = Path(fp)

        self.overwriteQ = overwrite
        if overwrite and self.root_fp.exists():
            logging.warning(f"Found existing output dir {self.root_fp}, overwriting...")
            shutil.rmtree(self.root_fp)

        os.makedirs(self.root_fp, exist_ok=True)

        if Path(seq).exists():
            self.query_fp = Path(seq)
        else:
            # Assuming seq is the literal sequence, not a file
            self.query_fp = self.root_fp / "query"
            self.write_query(seq)

        self.bootstraps_fp = self.root_fp / "RAxML_bestTree.genus1"
        self.base_tree_fp = self.root_fp / "RAxML_bestTree.genus2"
        self.bootstrapped_tree_fp = self.root_fp / "RAxML_bipartitions.final"

        self.combined_alignment_fp = self.root_fp / "combined_alignment.fasta"
        self.combined_tree_fp = self.root_fp / "RAxML_bestTree.combined"

        self.nearest_seqs_fp = self.root_fp / "nearest_seqs.fasta"
        self.nearest_seqs_aligned_fp = self.root_fp / "nearest_seqs_aligned.fasta"

        self.probs_fp = self.root_fp / "probabilities.tsv"

    def get_bootstraps(self) -> Path:
        return self.bootstraps_fp

    def get_base_tree(self) -> Path:
        return self.base_tree_fp

    def get_bootstrapped_tree(self) -> Path:
        return self.bootstrapped_tree_fp

    def get_combined_alignment(self) -> Path:
        return self.combined_alignment_fp

    def get_combined_tree(self) -> Path:
        return self.combined_tree_fp
    
    def get_nearest_seqs(self) -> Path:
        return self.nearest_seqs_fp

    def get_nearest_seqs_aligned(self) -> Path:
        return self.nearest_seqs_aligned_fp

    def get_query(self) -> Path:
        return self.query_fp

    def write_probs(self, probs: dict, header: str = ""):
        with open(self.probs_fp, "a+") as f:
            f.write(f"{header}\n")
            for s, p in probs.items():
                f.write(f"{s.split(' ')[0]}\t{round(p * 100, 5)}\n")

    def write_query(self, query: str):
        with open(self.query_fp, "w") as f:
            f.write(query)
