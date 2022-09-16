import eutils
import os
import requests
import subprocess
from tqdm import tqdm
from xml.etree import ElementTree as ET

# https://stackoverflow.com/questions/434287/how-to-iterate-over-a-list-in-chunks
def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

# Creates db with eutils
def create_db(api_key: str):
    if os.path.isfile("16S.db"):
        print("Existing 16S.db file found, skipping download...")
        return None

    search_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=33175%5BBioProject%5D%20OR%2033317%5BBioProject%5D&retmax=25000'
    search_response = requests.get(search_url)
    search_tree = ET.fromstring(search_response.content)
    
    ids = list()
    for id in search_tree[3]:
        ids.append(id.text)

    ec = eutils.Client(api_key=api_key)

    with open("16S.db", "w") as db, tqdm(total=round(len(ids)/250)) as pbar:
        for group in chunker(ids, 250):
            pbar.update(1)
            egs = ec.efetch(db='nuccore', id=",".join(group))
            for seq in egs:
                db.write(f">{str(seq)[6:-1]} {seq.organism}\n")
                db.write(f"{seq.sequence}\n")

    try:
        os.remove("esearch.fcgi")
    except OSError:
        pass

# Finds similar sequences to the one given with vsearch
# @param seq is the query sequence
# @return is a list of the closest sequences to the query
def find_similar(seq: str) -> list:
    #return ['AB681979', 'HG972968', 'KP319034', 'EF643377', 'MK942857', 'KM089834', 'AF029227', 'DQ453129', 'EU014688', 'AF530475', 'EU014685', 'EU014684', 'EU014680', 'AF008580', 'AE006468', 'KJ787692', 'X96963', 'FR870445', 'MK734184', 'HF558388', 'X80725', 'AJ508775', 'LR745848', 'FR870441', 'X96966', 'AF025371', 'HQ992945', 'AB682276', 'AJ508303', 'AY373829', 'MN603664', 'MK040622', 'AF130982', 'MK040621', 'HG933296', 'AF025363', 'Y17657', 'AF025364', 'AJ417484', 'HG933295', 'AF009171', 'X87276', 'HQ651841', 'HQ888848', 'JF795013', 'EU672801', 'EF488759', 'AJ251469', 'CP010523']
    with open("output/query.fasta", "w") as f:
        f.write(f">UNKNOWN\n")
        f.write(f"{seq}\n")
    
    subprocess.run(["vsearch",
    "--usearch_global", "db/type_species.fasta",
    "--db", "output/query.fasta",
    "--id", "0.9",
    "--fastapairs", "output/nearest_res.fasta"])

    ids = list()
    with open("output/nearest_res.fasta") as f:
        for l in f.readlines():
            if l[0] == ">" and l.strip() != ">UNKNOWN":
                ids.append(l[1:].strip())
    
    #if len(ids) > 50:
    #    ids = ids[:49]

    print(ids)
    return ids