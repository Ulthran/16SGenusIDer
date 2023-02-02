import logging
import os
from pathlib import Path

class OutputDir():
    """
    Controller for all of GenusFinder's output files
    """
    def __init__(self, fp: Path, overwrite: bool) -> None:
        self.root_fp = fp
        os.makedirs(self.root_fp, exist_ok=True)
        self.overwriteQ = overwrite
    

    def get_combined_alignment(self) -> Path:
        return self.root_fp / "combined_alignment.fasta"
    

    def get_combined_tree(self) -> Path:
        return self.root_fp / "RAxML_bestTree.combined"

    
    def get_query(self) -> Path:
        return self.root_fp / "query"

    
    def write_probs(self, probs: dict):
        with open(self.root_fp / "probabilities.tsv", "w") as f:
            for s, p in probs.items():
                f.write(f"{s.split(' ')[0]}\t{p}\n")


    def write_query(self, query: str):
        with open(self.root_fp / "query", "w") as f:
            f.write(query)