import argparse
import logging

from src.GenusFinder.CLI import MuscleAligner, RAxMLTreeBuilder, VsearchSearcher
from src.GenusFinder.DBDir import DBDir
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
        "--subtree_only",
        help="only use the subtree method, not more computationally intensive full tree alignment",
        action="store_false",
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

    db = DBDir(args.db, args.ncbi_api_key)

    ### Subtree alignment method ###

    searcher = VsearchSearcher()
    searcher.call(
        db.get_type_species(),
        out.get_query(),
        0.9,
        out.get_nearest_seqs()
    )

    aligner = MuscleAligner()
    aligner.call_simple(out.get_nearest_seqs(), out.get_nearest_seqs_aligned())

    tree_builder = RAxMLTreeBuilder()
    # Create 100 bootstrap trees
    tree_builder.call(
        392781,
        None,
        100,
        "GTRCAT",
        "genus1",
        10000,
        out.get_nearest_seqs_aligned(),
        None,
        out.root_fp,
        None,
    )
    # Create the base tree to use the bootstrapping trees with
    tree_builder.call(
        None,
        None,
        None,
        "GTRCAT",
        "genus2",
        10000,
        out.get_nearest_seqs_aligned(),
        None,
        out.root_fp,
        None,
    )
    # Create bootstrapped tree
    tree_builder.call(
        None,
        "b",
        None,
        "PROTGAMMAILG",
        "final",
        None,
        None,
        out.get_base_tree(),
        out.root_fp,
        out.get_bootstraps(),
    )

    

    ### Full tree alignment method ###

    if not args.subtree_only:
        aligner.call_profile(
            True, db.get_LTP_aligned(), out.get_query(), out.get_combined_alignment()
        )

        tree_builder.call(
            None, "y", None,
            "GTRCAT",
            "combined",
            10000,
            out.get_combined_alignment(),
            db.get_LTP_tree(),
            out.root_fp, None,
        )

        trainer = Trainer(out.get_combined_tree())
        out.write_probs(trainer.train(db.get_LTP_tree()))

    # identify_genus(
    #    build_tree(find_similar(args.seq, args.id), args.seq, args.keep_output),
    #    args.seq,
    # )
