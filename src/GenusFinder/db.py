import eutils
import logging
import os
import requests
import shutil
import subprocess as sp
from pathlib import Path
from tqdm import tqdm
from urllib.request import urlopen
from xml.etree import ElementTree as ET


class DB:
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
        self.LTP_tree_fp = self.root_fp / f"Tree_LTP_all_{self.LTP_VERSION}.ntree"
        self.LTP_csv_fp = self.root_fp / f"LTP_{self.LTP_VERSION}.csv"

    def get_16S_db(self) -> Path:
        if not self._16S_db.exists():
            logging.info(f"Creating {self._16S_db}...")
            self._create_16S_db()
        else:
            logging.info(f"Found {self._16S_db}, skipping download...")

        return self._16S_db

    def get_LTP_aligned(self) -> Path:
        return self._get_LTP(self.LTP_aligned_fp, self.LTP_aligned_fp.name)

    def get_LTP_blastdb(self) -> Path:
        return self._get_LTP(self.LTP_blastdb_fp, self.LTP_blastdb_fp.name)

    def get_LTP_tree(self) -> Path:
        return self._get_LTP(self.LTP_tree_fp, self.LTP_tree_fp.name)

    def get_LTP_csv(self) -> Path:
        return self._get_LTP(self.LTP_csv_fp, self.LTP_csv_fp.name)

    def _get_LTP(self, fp: Path, name: str) -> Path:
        if not fp.exists():
            url = self.url_for(name)
            logging.info(f"Fetching {url}...")
            with urlopen(url) as resp, open(fp, "wb") as f:
                shutil.copyfileobj(resp, f)
        else:
            logging.info(f"Found {fp}, skipping download...")

        return fp

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

    def url_for(self, name: str) -> str:
        return f"{self.LTP_URL}{name}"


# Finds similar sequences to the one given with vsearch
# @param seq is the query sequence
# @param id is the identity value for vsearch
# @return is a list of the closest sequences to the query
def find_similar(seq: str, id: str) -> list:
    # return ['AB681979', 'HG972968', 'KP319034', 'EF643377', 'MK942857', 'KM089834', 'AF029227', 'DQ453129', 'EU014688', 'AF530475', 'EU014685', 'EU014684', 'EU014680', 'AF008580', 'AE006468', 'KJ787692', 'X96963', 'FR870445', 'MK734184', 'HF558388', 'X80725', 'AJ508775', 'LR745848', 'FR870441', 'X96966', 'AF025371', 'HQ992945', 'AB682276', 'AJ508303', 'AY373829', 'MN603664', 'MK040622', 'AF130982', 'MK040621', 'HG933296', 'AF025363', 'Y17657', 'AF025364', 'AJ417484', 'HG933295', 'AF009171', 'X87276', 'HQ651841', 'HQ888848', 'JF795013', 'EU672801', 'EF488759', 'AJ251469', 'CP010523']
    with open("output/query.fasta", "w") as f:
        f.write(f">UNKNOWN\n")
        f.write(f"{seq}\n")

    sp.run(
        [
            "vsearch",
            "--usearch_global",
            "db/type_species.fasta",
            "--db",
            "output/query.fasta",
            "--id",
            id,
            "--fastapairs",
            "output/nearest_res.fasta",
        ]
    )

    ids = list()
    with open("output/nearest_res.fasta") as f:
        for l in f.readlines():
            if l[0] == ">" and l.strip() != ">UNKNOWN":
                ids.append(l[1:].strip())

    # if len(ids) > 50:
    #    ids = ids[:49]

    print(ids)
    return ids
