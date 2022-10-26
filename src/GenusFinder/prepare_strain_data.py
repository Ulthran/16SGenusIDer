###
# Taken from unassigner
# https://github.com/kylebittinger/unassigner
###

import argparse
import sys
import os

from GenusFinder.download import (
    LTP_TREE_URL, get_url, clean,
    LTP_METADATA_URL, LTP_SEQS_URL, LTP_ALIGN_URL,
    #process_ltp_seqs,
    )


def use_or_download(optional_fp, url, db_dir):
    if optional_fp is None:
        return get_url(url, os.path.join(db_dir, os.path.basename(url)))
    else:
        return optional_fp


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--ltp_metadata_fp", help=(
        "Filepath for LTP metadata (.csv file) "
        "[default: download from LTP website]"))
    p.add_argument("--ltp_seqs_fp", help=(
        "Filepath for unaligned 16S sequences from LTP (.fasta file) "
        "[default: download from LTP website]"))
    p.add_argument("--ltp_align_fp", help=(
        "Filepath for aligned 16S sequences from LTP (.fasta file) "
        "[default: download from LTP website]"))
    p.add_argument("--clean", action="store_true", help=(
        "Remove all downloaded and processed files."))
    p.add_argument("--db-dir", help=(
        "Filepath to download the files to "
        "[default: db/]"))
    args = p.parse_args(argv)

    if args.db_dir:
        db_dir = args.db_dir
        if not os.path.exists(args.db_dir) and not args.clean:
            os.mkdir(args.db_dir)
    else:
        db_dir = os.path.join(os.getcwd(), "db/")

    if args.clean is True:
        clean(db_dir)
        sys.exit(0)

    ltp_metadata_fp = use_or_download(
        args.ltp_metadata_fp, LTP_METADATA_URL, db_dir)
    ltp_seqs_fp = use_or_download(
        args.ltp_seqs_fp, LTP_SEQS_URL, db_dir)
    #process_ltp_seqs(ltp_seqs_fp, db_dir)
    ltp_align_fp = use_or_download(
        args.ltp_align_fp, LTP_ALIGN_URL, db_dir)
    ltp_tree_fp = use_or_download(
        None, LTP_TREE_URL, db_dir)

