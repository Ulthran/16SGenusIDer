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
    with open("query.fasta", "w") as f:
        f.write(f"> UNKNOWN\n")
        f.write(f"{seq}\n")
    
    subprocess.run(["vsearch",
    "--usearch_global", "16S.db",
    "--db", "query.fasta",
    "--id", "0.9",
    #"--uc_allhits",
    #"--maxaccepts", "50",
    "--fastapairs", "nearest.fasta"])

    ids = list()
    with open("nearest.fasta") as f:
        for l in f.readlines():
            if l[0] == ">" and len(l) > 1:
                ids.append(l[1:].strip())

    try:
        os.remove("query.fasta")
        os.remove("nearest.fasta")
    except OSError:
        pass

    return ids