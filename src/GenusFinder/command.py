import argparse
import logging
import os

from .db import create_db, find_similar
from .tree import build_tree, identify_genus


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--seq", help="the 16S sequence to be identified")
    p.add_argument("--id", help="the identity value to use with vsearch", default="0.9")
    p.add_argument(
        "--work_dir", help="the directory in which to put all output", default="output/"
    )
    p.add_argument(
        "--keep_output", help="switch to not remove output files", action="store_true"
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

    work_dir = args.work_dir
    if not os.path.exists(work_dir):
        logging.info(f"{work_dir} not found, creating...")
        os.makedirs(work_dir)

    # identify_genus(
    #    build_tree(find_similar(args.seq, args.id), args.seq, args.keep_output),
    #    args.seq,
    # )
