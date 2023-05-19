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


class LTPAlignment(DBFile):
    """
    Holds and operates on the LTP alignment
    """

    def __init__(self, fp: Path) -> None:
        super().__init__(fp)


class LTPTree(DBFile):
    """
    Holds and operates on the LTP tree
    """
    
    def __init__(self, fp: Path) -> None:
        super().__init__(fp)


class LTPBlast(DBFile):
    """
    Holds and operates on the LTP blast db
    """

    def __init__(self, fp: Path) -> None:
        super().__init__(fp)


class TypeSpecies(DBFile):
    """
    Holds and creates the type_species fasta file
    """

    def __init__(self, fp: Path, blastdb_fp: LTPBlast) -> None:
        super().__init__(fp)
        self.blastdb_fp = blastdb_fp
    
    def get(self) -> Path:
        assert self.blastdb_fp.get().exists()

        if not self.fp.exists():
            logging.info(f"Creating {self.fp}...")
            self._generate_type_species()
        else:
            logging.info(f"Found {self.fp}, skipping creation...")

        return super().get()
    
    def _generate_type_species(self):
        accession_cts = collections.defaultdict(int)
        with open(self.blastdb_fp.get()) as f_in:
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
        self.rooted_LTP_tree_fp = self.root_fp / "RAxML_bestTree.rooted"
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
        #ret = self._get_LTP(self.LTP_aligned_fp, self.LTP_aligned_fp.name)
        shutil.copyfile(self.root_fp / "LTP", self.LTP_aligned_fp)
        self.clean_alignment()
        if not self.verify_alignment():
            sys.exit()
        ret = self.LTP_aligned_fp
        return ret

    def get_LTP_blastdb(self) -> Path:
        return self._get_LTP(self.LTP_blastdb_fp, self.LTP_blastdb_fp.name)

    def get_LTP_tree(self) -> Path:
        ret = self._get_LTP(self.LTP_tree_fp, self.LTP_tree_fp.name)
        #if not self.verify_tree():
        #    sys.exit()
        return ret
    
    def get_rooted_LTP_tree(self) -> Path:
        return 

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
                    f_out.write(">{0}\t{1}\n{2}\n".format(accession, species_name, seq.replace(" ", "").replace("U", "T")))

    def _get_LTP(self, fp: Path, name: str) -> Path:
        if not fp.exists():
            url = self.url_for(name)
            logging.info(f"Fetching {url}...")
            with urlopen(url) as resp, open(fp, "wb") as f:
                shutil.copyfileobj(resp, f)
            
            if "aligned" in name:
                self.clean_alignment()
            elif "tree" in name:
                self.clean_tree()
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
                f_align.seek(0)
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
            seq_len = 0 # Get seq len on first pass
            for desc, seq in self._parse_fasta(f):
                if len(desc.strip()) < 2:
                    logging.error("Annotation line empty")
                    return False
                acceptable_chars = set(["A", "C", "G", "T", "-"])
                if seq_len == 0:
                    seq_len = len(seq)
                elif seq_len != len(seq):
                    logging.error(f"Ragged alignment of length {len(seq)} instead of {seq_len}: {desc}")
                    return False
                if not set(list(seq)).issubset(acceptable_chars):
                    logging.error(f"{set(list(seq))} is not a subset of {acceptable_chars}")
                    return False
            return True
    
    def clean_tree(self):
        logging.info("Cleaning LTP tree...")
        temp_fp = self.root_fp / "temp_tree.ntree"
        with open(temp_fp, "w") as f_temp, open(self.LTP_tree_fp) as f_tree:
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

        os.remove(self.LTP_tree_fp)
        os.rename(temp_fp, self.LTP_tree_fp)
    
    def verify_tree(self) -> bool:
        try:
            with open(self.LTP_tree_fp) as f:
                Tree("".join([l.strip() for l in f.readlines()]))
        except NewickError as e:
            logging.error("Error reading tree into ete3")
            return False
        
        return True
    
    


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
