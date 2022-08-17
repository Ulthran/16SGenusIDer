import eutils
import os
import requests
from tqdm import tqdm
from xml.etree import ElementTree as ET

# https://stackoverflow.com/questions/434287/how-to-iterate-over-a-list-in-chunks
def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

# Creates db with eutils
def create_db():
    search_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=33175%5BBioProject%5D%20OR%2033317%5BBioProject%5D&retmax=25000'
    search_response = requests.get(search_url)
    search_tree = ET.fromstring(search_response.content)
    
    ids = list()
    for id in search_tree[3]:
        ids.append(id.text)

    ec = eutils.Client()

    with open("16S.db", "w") as db, tqdm(total=round(len(ids)/250)) as pbar:
        for group in chunker(ids, 250):
            pbar.update(1)
            egs = ec.efetch(db='nuccore', id=",".join(group))
            for seq in egs:
                db.write(f"> {seq.organism} {str(seq)[6:-1]}\n")
                db.write(f"{seq.sequence}\n")

    os.remove("esearch.fcgi")

# Finds similar sequences to the one given with vsearch
# @param seq is the query sequence
# @return is a list of the closest sequences to the query
def find_similar(seq: str) -> list:
    return list()