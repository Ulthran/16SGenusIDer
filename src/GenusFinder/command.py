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
    if not (args.overwrite or out.nearest_seqs.get().exists()):
        searcher.call(db.type_species.get(), out.query.get(), 0.9, out.nearest_seqs.get())

    aligner = MuscleAligner()
    if not (args.overwrite or out.nearest_seqs_aligned.get().exists()):
        aligner.call_simple(out.nearest_seqs.get(), out.nearest_seqs_aligned.get())

    tree_builder = RAxMLTreeBuilder()
    # Create 100 bootstrap trees
    tree_builder.call(
        b=392781,
        N=100,
        m="GTRCAT",
        n="subtree1",
        p=10000,
        s=out.nearest_seqs_aligned.get(),
        w=out.root_fp,
    )
    # Create the base tree to use the bootstrapping trees with
    tree_builder.call(
        m="GTRCAT",
        n="subtree2",
        p=10000,
        s=out.nearest_seqs_aligned.get(),
        w=out.root_fp,
    )
    # Create bootstrapped tree
    tree_builder.call(
        f="b",
        m="PROTGAMMAILG",
        n="final",
        t=out.base_tree.get(),
        w=out.root_fp,
        z=out.bootstraps.get(),
    )

    algorithms = Algorithms(
        out.bootstrapped_tree.get(), db.type_species.get(), out.query.get()
    )
    # Set write_mode to "w" to clear any existing output
    out.probs.get().write_probs(algorithms.distance_probs(), "Distance-based subtree probabilities", "w")
    out.probs.get().write_probs(
        algorithms.bootstrap_probs(), "Bootstrap-based subtree probabilities"
    )

    logging.info(f"Subtree method finished! Check {out.probs_fp} for results.")

    ### Full tree alignment method ###

    if args.subtree_only:
        if not (args.overwrite or out.combined_alignment.get().exists()):
            aligner.call_profile(
                True, db.LTP_aligned.get(), out.query.get(), out.combined_alignment.get()
            )

        # Add the unknown onto the LTP tree using the combined alignment
        tree_builder.call(
            None,
            "y",
            None,
            "GTRCAT",
            "combined",
            10000,
            out.combined_alignment.get(),
            db.LTP_tree.get(),
            out.root_fp,
            None,
        )

        algorithms = Algorithms(db.LTP_tree.get(), db.type_species.get(), out.query.get())
        out.probs.get().write_probs(
            algorithms.train(), "Full tree alignment probabilities"
        )

        logging.info(f"Full tree method finished! Check {out.probs.get()} for results.")