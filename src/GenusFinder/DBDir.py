import collections
import eutils
import logging
import os
import re
import requests
import shutil
import sys
from ete3 import Tree
from ete3.parser.newick import NewickError
from io import TextIOWrapper
from itertools import groupby
from pathlib import Path
from tqdm import tqdm
from urllib.request import urlopen
from xml.etree import ElementTree as ET
from .parsers import parse_desc, parse_fasta


class DBFile:
    """
    Base class to define a file in the database and associated operations
    """

    def __init__(self, fp: Path) -> None:
        self.fp = fp
        self.fn = fp.name

    def get(self) -> Path:
        return self.fp


class LTPFile(DBFile):
    """
    Base class to define an LTP file in the database and associated operations
    """

    def __init__(self, fp: Path) -> None:
        super().__init__(fp)

        self.LTP_URL = "https://imedea.uib-csic.es/mmg/ltp/wp-content/uploads/ltp/"

    def get(self) -> Path:
        if not self.fp.exists():
            url = self.url_for(self.fn)
            logging.info(f"Fetching {url}...")
            with urlopen(url) as resp, open(self.fp, "wb") as f:
                shutil.copyfileobj(resp, f)
        else:
            logging.info(f"Found {self.fp}, skipping download...")

        return self.fp

    def url_for(self, name: str) -> str:
        return f"{self.LTP_URL}{name}"


class LTPAlignment(LTPFile):
    """
    Holds and operates on the LTP alignment
    """

    def __init__(self, fp: Path) -> None:
        super().__init__(fp)
    
    def get(self) -> Path:
        if self.fp.exists():
            return super().get()
        else:
            ret = super().get()
            self.clean()
            return ret
    
    def clean(self):
        replacements_map = {
            " ": "", # LTP's fun syntax
            "\n": "",
            ".": "-",
            "U": "T",
            "R": "A", # Choose a letter per abbreviation and stick with it
            "Y": "T",
            "M": "C",
            "K": "G",
            "S": "C",
            "W": "T",
            "H": "A",
            "B": "G",
            "V": "C",
            "D": "T",
            "N": "A",
            "A": "A", # Keep what's already good
            "C": "C",
            "G": "G",
            "T": "T",
            "-": "-",
        }

        logging.info("Cleaning LTP alignment...")
        temp_fp = self.fp.parents[0] / "temp_alignment.fasta"
        with open(self.fp, "rbU") as f:
            num_lines = sum(1 for _ in f)
        with open(temp_fp, "w") as f_temp, open(self.fp) as f_align:
            with tqdm(total=num_lines) as pbar:
                for header_str, seq_str in parse_fasta(f_align):
                    pbar.update(1)
                    f_temp.write(f">{header_str}")
                    f_temp.write("".join([replacements_map[c] for c in seq_str]))
                
        os.remove(self.fp)
        os.rename(temp_fp, self.fp)


class LTPTree(LTPFile):
    """
    Holds and operates on the LTP tree
    """
    
    def __init__(self, fp: Path) -> None:
        super().__init__(fp)
    
    def get(self) -> Path:
        if self.fp.exists():
            return super().get()
        else:
            ret = super().get()
            self.clean()
            return ret
    
    def clean(self):
        logging.info("Cleaning LTP tree...")
        temp_fp = self.fp.parents[0] / "temp_tree.ntree"
        with open(temp_fp, "w") as f_temp, open(self.fp) as f_tree:
            header = False
            for line in f_tree.readlines():
                if line[0] == "[":
                    header = True
                if header:
                    if line[0] == "]":
                        header = False
                elif "'" in line and "," in line:
                    new_id = line[line.index("'")+1:line.index(",", line.index("'"))]
                    f_temp.write(line[:line.index("'")] + new_id + line[line.index("'", line.index("'")+1)+1:])
                elif "'" in line:
                    f_temp.write(line.replace("'", "").replace(" ", ""))
                else:
                    f_temp.write(line)

        os.remove(self.fp)
        os.rename(temp_fp, self.fp)


class LTPBlast(LTPFile):
    """
    Holds and operates on the LTP blast db
    """

    def __init__(self, fp: Path) -> None:
        super().__init__(fp)
    
    def get(self) -> Path:
        return super().get()


class TypeSpecies(DBFile):
    """
    Holds and creates the type_species fasta file
    """

    def __init__(self, fp: Path, blastdb: LTPBlast) -> None:
        super().__init__(fp)
        self.blastdb = blastdb
    
    def get(self) -> Path:
        assert self.blastdb.get().exists()

        if not self.fp.exists():
            logging.info(f"Creating {self.fp}...")
            self._generate_type_species()
        else:
            logging.info(f"Found {self.fp}, skipping creation...")

        return super().get()
    
    def _generate_type_species(self):
        accession_cts = collections.defaultdict(int)
        with open(self.blastdb.get()) as f_in:
            with open(self.fp, "w") as f_out:
                for desc, seq in parse_fasta(f_in):
                    accession, species_name = parse_desc(desc)
                    if not accession or not species_name:
                        continue
                    # Some accessions refer to genomes with more than one 16S gene
                    # So accessions can be legitimately repeated with distinct gene sequences
                    accession_times_previously_seen = accession_cts[accession]
                    accession_cts[accession] += 1
                    if accession_times_previously_seen > 0:
                        accession = "{0}_repeat{1}".format(
                            accession, accession_times_previously_seen
                        )
                    f_out.write(">{0}\t{1}\n{2}\n".format(accession, species_name, seq.replace(" ", "").replace("U", "T")))


class DBDir:
    """
    Controller for all of GenusFinder's database files\n
    Maintains a 16S db made from an NCBI eutils query and mulitiple LTP files\n
    Interface directly with internal objects\n
    i.e. call LTPAlignment's get() method directly
    """

    def __init__(self, fp: Path, esearch_api_key: str) -> None:
        self.root_fp = Path(fp)
        os.makedirs(self.root_fp, exist_ok=True)

        self.key = esearch_api_key
        self.LTP_VERSION = "06_2022"

        #self._16S_db = self.root_fp / "16S.db"
        self.LTP_aligned = LTPAlignment(self.root_fp / f"LTP_{self.LTP_VERSION}_aligned.fasta")
        self.LTP_blastdb = LTPBlast(self.root_fp / f"LTP_{self.LTP_VERSION}_blastdb.fasta")
        self.LTP_tree = LTPTree(self.root_fp / f"LTP_all_{self.LTP_VERSION}.ntree")
        #self.LTP_csv_fp = self.root_fp / f"LTP_{self.LTP_VERSION}.csv"
        self.type_species = TypeSpecies(self.root_fp / "type_species.fasta", self.LTP_blastdb)

    
    
    
    
    
    
    
    
    
#    def get_16S_db(self) -> Path:
#        if not self._16S_db.exists():
#            logging.info(f"Creating {self._16S_db}...")
#            self._create_16S_db()
#        else:
#            logging.info(f"Found {self._16S_db}, skipping download...")
#
#        return self._16S_db
#
#    def _create_16S_db(self):
#        def chunker(seq, size):
#            return (seq[pos : pos + size] for pos in range(0, len(seq), size))
#
#        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=33175%5BBioProject%5D%20OR%2033317%5BBioProject%5D&retmax=25000"
#        search_response = requests.get(search_url)
#        search_tree = ET.fromstring(search_response.content)
#
#        ids = list()
#        for id in search_tree[3]:
#            ids.append(id.text)
#
#        ec = eutils.Client(api_key=self.key)
#
#        with open(self._16S_db, "w") as db, tqdm(total=round(len(ids) / 250)) as pbar:
#            for group in chunker(ids, 250):
#                pbar.update(1)
#                egs = ec.efetch(db="nuccore", id=",".join(group))
#                for seq in egs:
#                    db.write(f">{str(seq)[6:-1]} {seq.organism}\n")
#                    db.write(f"{seq.sequence}\n")
    


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
