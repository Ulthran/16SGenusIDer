import logging
import re
from io import TextIOWrapper
from itertools import groupby

def parse_fasta(f: TextIOWrapper):
    """
    Lazily return description/sequence pairs from a fasta file
    """
    faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        header_str = header.__next__()[1:].strip()

        # join all sequence lines to one
        seq_str = "".join(s.strip() for s in faiter.__next__())

        yield (header_str, seq_str)

def parse_desc(desc: str) -> tuple:
    """
    Return the accession and species name from a fasta sequence description\n
    Otherwise return None, None
    """
    try:
        accession = re.findall(r"\[accession=(.*?)\]", desc)[0]
        species_name = re.findall(r"\[organism=(.*?)\]", desc)[0]
    except IndexError as e:
        logging.error(f"Couldn't find accession and/or organism identifier in {desc}")
        logging.error(f"Skipping this sequence...")
        return None, None
    return accession, species_name