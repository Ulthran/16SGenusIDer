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
        
        self.probs_fp = self.root_fp / "probabilities.tsv"

    def get_combined_alignment(self) -> Path:
        return self.root_fp / "combined_alignment.fasta"

    def get_combined_tree(self) -> Path:
        return self.root_fp / "RAxML_bestTree.combined"

    def get_query(self) -> Path:
        return self.query_fp

    def write_probs(self, probs: dict):
        with open(self.probs_fp, "w") as f:
            for s, p in probs.items():
                f.write(f"{s.split(' ')[0]}\t{p}\n")

    def write_query(self, query: str):
        with open(self.query_fp, "w") as f:
            f.write(query)
