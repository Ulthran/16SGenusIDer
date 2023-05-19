import argparse
import logging
import sys

from .CLI import MuscleAligner, RAxMLTreeBuilder, VsearchSearcher
from .DBDir import DBDir
from .OutputDir import OutputDir
from .Algorithms import Algorithms


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
        help="overwrites any existing files of the same name, otherwise assumes existing files should be used as they are",
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
    if not args.seq:
        p.print_help(sys.stderr)
        sys.exit(1)
    logging.basicConfig()
    logging.getLogger().setLevel(args.log_level)

    out = OutputDir(args.output, args.seq, args.overwrite)

    db = DBDir(args.db, args.ncbi_api_key)

    ### Subtree alignment method ###

    searcher = VsearchSearcher()
    if not (args.overwrite or out.get_nearest_seqs().exists()):
        searcher.call(db.get("type_species"), out.get_query(), 0.9, out.get_nearest_seqs())

    aligner = MuscleAligner()
    if not (args.overwrite or out.get_nearest_seqs_aligned().exists()):
        aligner.call_simple(out.get_nearest_reduced_seqs(), out.get_nearest_seqs_aligned())

    tree_builder = RAxMLTreeBuilder()
    # Create 100 bootstrap trees
    tree_builder.call(
        b=392781,
        N=100,
        m="GTRCAT",
        n="subtree1",
        p=10000,
        s=out.get_nearest_seqs_aligned(),
        w=out.root_fp,
    )
    # Create the base tree to use the bootstrapping trees with
    tree_builder.call(
        m="GTRCAT",
        n="subtree2",
        p=10000,
        s=out.get_nearest_seqs_aligned(),
        w=out.root_fp,
    )
    # Create bootstrapped tree
    tree_builder.call(
        f="b",
        m="PROTGAMMAILG",
        n="final",
        t=out.get_base_tree(),
        w=out.root_fp,
        z=out.get_bootstraps(),
    )

    algorithms = Algorithms(
        out.get_bootstrapped_tree(), db.get("type_species"), out.get_query()
    )
    # Set write_mode to "w" to clear any existing output
    out.write_probs(algorithms.distance_probs(), "Distance-based subtree probabilities", "w")
    out.write_probs(
        algorithms.bootstrap_probs(), "Bootstrap-based subtree probabilities"
    )

    logging.info(f"Subtree method finished! Check {out.probs_fp} for results.")

    ### Full tree alignment method ###

    if args.subtree_only:
        if not (args.overwrite or out.get_combined_alignment().exists()):
            aligner.call_profile(
                True, db.get("alignment"), out.get_query(), out.get_combined_alignment()
            )

        # Add the unknown onto the LTP tree using the combined alignment
        tree_builder.call(
            None,
            "y",
            None,
            "GTRCAT",
            "combined",
            10000,
            out.get_combined_alignment(),
            db.get("tree"),
            out.root_fp,
            None,
        )

        algorithms = Algorithms(db.get("tree"), db.get("type_species"), out.get_query())
        out.write_probs(
            algorithms.train(), "Full tree alignment probabilities"
        )

        logging.info(f"Full tree method finished! Check {out.probs_fp} for results.")