import logging
import os
from pathlib import Path


class OutputDir:
    """
    Controller for all of GenusFinder's output files
    """

    def __init__(self, fp: Path, overwrite: bool, seq: str) -> None:
        self.root_fp = Path(fp)
        os.makedirs(self.root_fp, exist_ok=True)

        self.overwriteQ = overwrite

        if Path(seq).exists():
            self.query_fp = Path(seq)
        else:
            # Assuming seq is the literal sequence, not a file
            self.query_fp = self.root_fp / "query"
            self.write_query(seq)

    def get_combined_alignment(self) -> Path:
        return self.root_fp / "combined_alignment.fasta"

    def get_combined_tree(self) -> Path:
        return self.root_fp / "RAxML_bestTree.combined"

    def get_query(self) -> Path:
        return self.query_fp

    def write_probs(self, probs: dict):
        with open(self.root_fp / "probabilities.tsv", "w") as f:
            for s, p in probs.items():
                f.write(f"{s.split(' ')[0]}\t{p}\n")

    def write_query(self, query: str):
        with open(self.query_fp, "w") as f:
            f.write(query)
