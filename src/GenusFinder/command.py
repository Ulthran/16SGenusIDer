import argparse
import logging
import os

from src.GenusFinder.CLI import MuscleAligner, RAxMLTreeBuilder
from src.GenusFinder.DB import DB
from src.GenusFinder.OutputDir import OutputDir
from src.GenusFinder.Trainer import Trainer


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument(
        "--seq", help="the 16S sequence to be identified or a file containing it"
    )
    p.add_argument("--id", help="the identity value to use with vsearch", default="0.9")
    p.add_argument(
        "--ncbi_api_key",
        help="your NCBI API key for making more esearch requests per second",
        default="",
    )
    p.add_argument(
        "--output", help="the directory in which to put all output", default="output/"
    )
    p.add_argument(
        "--db", help="the directory in which to put all database files", default="db/"
    )
    p.add_argument(
        "--overwrite",
        help="overwrites any existing files of the same name",
        action="store_true",
    )
    p.add_argument(
        "--log_level",
        type=int,
        help="Sets the log level, default is info, 10 for debug (Default: 20)",
        default=20,
    )

    args = p.parse_args(argv)
    logging.basicConfig()
    logging.getLogger().setLevel(args.log_level)

    out = OutputDir(args.output, args.seq, args.overwrite)

    db = DB(args.db, args.ncbi_api_key)

    aligner = MuscleAligner()
    aligner.call(
        True, db.get_LTP_aligned(), out.get_query(), out.get_combined_alignment()
    )

    tree_builder = RAxMLTreeBuilder()
    tree_builder.call(
        "GTRCAT",
        "combined",
        10000,
        "y",
        out.get_combined_alignment(),
        db.get_LTP_tree(),
        out.root_fp,
    )

    trainer = Trainer(out.get_combined_tree())
    out.write_probs(trainer.train(db.get_LTP_tree()))

    # identify_genus(
    #    build_tree(find_similar(args.seq, args.id), args.seq, args.keep_output),
    #    args.seq,
    # )
