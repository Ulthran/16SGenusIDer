import collections
import eutils
import logging
import os
import re
import requests
import shutil
import sys
import tempfile
from .CLI import MuscleAligner
from io import StringIO, TextIOWrapper
from pathlib import Path
from tqdm import tqdm
from urllib.request import urlopen
from xml.etree import ElementTree as ET


class DBDir:
    """
    Controller for all of GenusFinder's database files\n
    Maintains a 16S db made from an NCBI eutils query and mulitiple LTP files
    """

    def __init__(self, fp: Path, esearch_api_key: str) -> None:
        self.root_fp = Path(fp)
        os.makedirs(self.root_fp, exist_ok=True)

        self.key = esearch_api_key

        self.LTP_VERSION = "06_2022"
        self.LTP_URL = f"https://imedea.uib-csic.es/mmg/ltp/wp-content/uploads/ltp/"

        self._16S_db = self.root_fp / "16S.db"
        self.LTP_aligned_fp = self.root_fp / f"LTP_{self.LTP_VERSION}_aligned.fasta"
        self.LTP_blastdb_fp = self.root_fp / f"LTP_{self.LTP_VERSION}_blastdb.fasta"
        self.LTP_tree_fp = self.root_fp / f"LTP_all_{self.LTP_VERSION}.ntree"
        self.LTP_csv_fp = self.root_fp / f"LTP_{self.LTP_VERSION}.csv"
        self.type_species_fp = self.root_fp / "type_species.fasta"

    def get_16S_db(self) -> Path:
        if not self._16S_db.exists():
            logging.info(f"Creating {self._16S_db}...")
            self._create_16S_db()
        else:
            logging.info(f"Found {self._16S_db}, skipping download...")

        return self._16S_db

    def get_LTP_aligned(self) -> Path:
        ret = self._get_LTP(self.LTP_aligned_fp, self.LTP_aligned_fp.name)
        if not self.verify_alignment():
            sys.exit()
        return ret

    def get_LTP_blastdb(self) -> Path:
        return self._get_LTP(self.LTP_blastdb_fp, self.LTP_blastdb_fp.name)

    def get_LTP_tree(self) -> Path:
        return self._get_LTP(self.LTP_tree_fp, self.LTP_tree_fp.name)

    def get_LTP_csv(self) -> Path:
        return self._get_LTP(self.LTP_csv_fp, self.LTP_csv_fp.name)

    def get_type_species(self) -> Path:
        if not self.type_species_fp.exists():
            logging.info(f"Creating {self.type_species_fp}...")
            self._generate_type_species()
        else:
            logging.info(f"Found {self.type_species_fp}, skipping creation...")

        return self.type_species_fp
    
    def _generate_type_species(self):
        accession_cts = collections.defaultdict(int)
        with open(self.get_LTP_blastdb()) as f_in:
            with open(self.type_species_fp, "w") as f_out:
                for desc, seq in self._parse_fasta(f_in):
                    accession, species_name = self._parse_desc(desc)
                    if not accession or not species_name:
                        continue
                    # Some accessions refer to genomes with more than one 16S gene
                    # So accessions can be legitiamtely repeated with distinct gene sequences
                    accession_times_previously_seen = accession_cts[accession]
                    accession_cts[accession] += 1
                    if accession_times_previously_seen > 0:
                        accession = "{0}_repeat{1}".format(
                            accession, accession_times_previously_seen
                        )
                    f_out.write(">{0}\t{1}\n{2}\n".format(accession, species_name, seq))

    def _get_LTP(self, fp: Path, name: str) -> Path:
        if not fp.exists():
            url = self.url_for(name)
            logging.info(f"Fetching {url}...")
            with urlopen(url) as resp, open(fp, "wb") as f:
                shutil.copyfileobj(resp, f)
            
            if "aligned" in name:
                self.clean_alignment()
        else:
            logging.info(f"Found {fp}, skipping download...")

        return fp

    def url_for(self, name: str) -> str:
        return f"{self.LTP_URL}{name}"

    def _create_16S_db(self):
        def chunker(seq, size):
            return (seq[pos : pos + size] for pos in range(0, len(seq), size))

        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=33175%5BBioProject%5D%20OR%2033317%5BBioProject%5D&retmax=25000"
        search_response = requests.get(search_url)
        search_tree = ET.fromstring(search_response.content)

        ids = list()
        for id in search_tree[3]:
            ids.append(id.text)

        ec = eutils.Client(api_key=self.key)

        with open(self._16S_db, "w") as db, tqdm(total=round(len(ids) / 250)) as pbar:
            for group in chunker(ids, 250):
                pbar.update(1)
                egs = ec.efetch(db="nuccore", id=",".join(group))
                for seq in egs:
                    db.write(f">{str(seq)[6:-1]} {seq.organism}\n")
                    db.write(f"{seq.sequence}\n")

    def clean_alignment(self):
        replacements_map = {
            " ": "", # LTP's weird syntax
            ".": "----",
            "U": "T---",
            "R": "AG--", # Expand these abbreviations
            "Y": "TC--",
            "M": "CA--",
            "K": "TG--",
            "S": "CG--",
            "W": "TA--",
            "H": "TCA-",
            "B": "TCG-",
            "V": "CAG-",
            "D": "TAG-",
            "N": "TCAG",
            "A": "A---", # Keep what's already good
            "C": "C---",
            "G": "G---",
            "T": "T---",
            "-": "----",
            "\n": "\n",
        }

        logging.info("Cleaning LTP alignment...")
        temp_fp = self.root_fp / "temp_alignment.fasta"
        with open(temp_fp, "w") as f_temp, open(self.LTP_aligned_fp) as f_align:
            with tqdm(total=len(f_align.readlines())) as pbar:
                for line in f_align.readlines():
                    pbar.update(1)
                    if line[0] == ">":
                        f_temp.write(f"{line}")
                    else:
                        f_temp.write("".join([replacements_map[c] for c in line]))
                
        os.remove(self.LTP_aligned_fp)
        os.rename(temp_fp, self.LTP_aligned_fp)
    
    def verify_alignment(self) -> bool:
        with open(self.LTP_aligned_fp) as f:
            last = False # False: seq last, True: annotation last
            seq_len = 0 # Get seq len on first pass
            for line in f.readlines():
                if last:
                    if line[0] != ">":
                        logging.error("Annotation line doesn't start with '>'")
                        return False
                    if len(line.strip()) < 2:
                        logging.error("Annotation line empty")
                        return False
                else:
                    if seq_len == 0:
                        seq_len = len(line)
                    elif seq_len != len(line):
                        logging.error("Ragged alignment")
                        return False
                    if set(list(line)).issubset(set(["A", "C", "G", "T", "-", "\n"])):
                        logging.error(f"{set(list(line))} is not subset of {['A', 'C', 'G', 'T', '-', '\n']}")
                        return False
                last = not last
            return True
    
    @staticmethod
    def _parse_fasta(f: TextIOWrapper, trim_desc = False):
        f = iter(f)
        try:
            desc = next(f).strip()[1:]
            if trim_desc:
                desc = desc.split()[0]
        except StopIteration:
            return
        seq = StringIO()
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                yield desc, seq.getvalue()
                desc = line[1:]
                if trim_desc:
                    desc = desc.split()[0]
                seq = StringIO()
            else:
                seq.write(line.replace(" ", "").replace("U", "T"))
        yield desc, seq.getvalue()
    
    @staticmethod
    def _parse_desc(desc: str) -> tuple:
        try:
            accession = re.findall(r"\[accession=(.*?)\]", desc)[0]
            species_name = re.findall(r"\[organism=(.*?)\]", desc)[0]
        except IndexError as e:
            logging.error(f"Couldn't find accession and/or organism identifier in {desc}")
            logging.error(f"Skipping this sequence...")
            return None, None
        return accession, species_name


# Finds similar sequences to the one given with vsearch
# @param seq is the query sequence
# @param id is the identity value for vsearch
# @return is a list of the closest sequences to the query
# def find_similar(seq: str, id: str) -> list:
#    # return ['AB681979', 'HG972968', 'KP319034', 'EF643377', 'MK942857', 'KM089834', 'AF029227', 'DQ453129', 'EU014688', 'AF530475', 'EU014685', 'EU014684', 'EU014680', 'AF008580', 'AE006468', 'KJ787692', 'X96963', 'FR870445', 'MK734184', 'HF558388', 'X80725', 'AJ508775', 'LR745848', 'FR870441', 'X96966', 'AF025371', 'HQ992945', 'AB682276', 'AJ508303', 'AY373829', 'MN603664', 'MK040622', 'AF130982', 'MK040621', 'HG933296', 'AF025363', 'Y17657', 'AF025364', 'AJ417484', 'HG933295', 'AF009171', 'X87276', 'HQ651841', 'HQ888848', 'JF795013', 'EU672801', 'EF488759', 'AJ251469', 'CP010523']
#    with open("output/query.fasta", "w") as f:
#        f.write(f">UNKNOWN\n")
#        f.write(f"{seq}\n")

#    sp.run(
#        [
#            "vsearch",
#            "--usearch_global",
#            "db/type_species.fasta",
#            "--db",
#            "output/query.fasta",
#            "--id",
#            id,
#            "--fastapairs",
#            "output/nearest_res.fasta",
#        ]
#    )

#    ids = list()
#    with open("output/nearest_res.fasta") as f:
#        for l in f.readlines():
#            if l[0] == ">" and l.strip() != ">UNKNOWN":
#                ids.append(l[1:].strip())

# if len(ids) > 50:
#    ids = ids[:49]

#    print(ids)
#    return ids
