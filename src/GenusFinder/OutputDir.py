import logging
import os
import shutil
from pathlib import Path
from . import parse_fasta


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

        try:
            if Path(seq).exists():
                self.query_fp = Path(seq)
            else:
                # Assuming seq is the literal sequence, not a file
                self.query_fp = self.root_fp / "query.fasta"
                self.write_query(seq)
        except OSError as e:
            logging.info(f"--seq arg is too long to be a path, interpreting as literal")
            self.query_fp = self.root_fp / "query.fasta"
            self.write_query(seq)

        self.bootstraps_fp = self.root_fp / "RAxML_bootstrap.subtree1"
        self.base_tree_fp = self.root_fp / "RAxML_bestTree.subtree2"
        self.bootstrapped_tree_fp = self.root_fp / "RAxML_bipartitions.final"

        self.combined_alignment_fp = self.root_fp / "combined_alignment.fasta"
        self.combined_tree_fp = self.root_fp / "RAxML_bestTree.combined"

        self.nearest_seqs_fp = self.root_fp / "nearest_seqs.fasta"
        self.temp_nearest_seqs_fp = self.root_fp / "temp_nearest_seqs.fasta"
        self.nearest_seqs_reduced_fp = self.root_fp / "nearest_seqs_reduced.fasta"
        self.nearest_seqs_aligned_fp = self.root_fp / "nearest_seqs_aligned.fasta"

        self.probs_fp = self.root_fp / "probabilities.tsv"

    ### Getters

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
    
    def get_nearest_reduced_seqs(self) -> Path:
        if not self.nearest_seqs_reduced_fp.exists():
            self.reduce_subtree()
        return self.nearest_seqs_reduced_fp

    def get_nearest_seqs_aligned(self) -> Path:
        return self.nearest_seqs_aligned_fp

    def get_query(self) -> Path:
        return self.query_fp
    
    ### Utilities
    
    def reduce_subtree(self):
        ids = []
        with open(self.nearest_seqs_fp) as f_in, open(self.nearest_seqs_reduced_fp, "w") as f_out:
            for id, seq in parse_fasta(f_in):
                if len(id.strip()) > 2:
                    ids.append(id)
                    f_out.write(f"> {id}\n")
                    f_out.write(f"{seq}\n")
                if len(ids) >= 50:
                    break

            logging.debug(f"Reduced accessions list: {str(ids)}")

        with open(self.nearest_seqs_reduced_fp, "a+") as f, open(self.query_fp) as query_f:
            f.write(f"{query_f.readline().strip()}\n")
            f.write(f"{query_f.readline().strip()}\n")

    ### Writers

    def write_probs(self, probs: dict, header: str = "", write_mode: str = "a+"):
        with open(self.probs_fp, write_mode) as f:
            f.write(f"\n{header}\n\n")
            for s, p in probs.items():
                if p > 0.0001:
                    f.write(f"{s.split(' ')[0]}\t{round(p * 100, 5)}\n")

    def write_query(self, query: str):
        with open(self.query_fp, "w") as f:
            f.write("> UNKNOWN\n")
            f.write(f"{query}\n")
