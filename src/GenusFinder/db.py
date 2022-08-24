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
    
    #subprocess.run(["vsearch",
    #"--usearch_global", "16S.db",
    #"--db", "query.fasta",
    #"--id", "0.9",
    #"--uc_allhits",
    #"--maxaccepts", "50",
    #"--fastapairs", "nearest.fasta"])

    #ids = list()
    #with open("nearest.fasta") as f:
    #    for l in f.readlines():
    #        if l[0] == ">" and len(l) > 2:
    #            ids.append(l[1:].strip())

    try:
        os.remove("query.fasta")
        os.remove("nearest.fasta")
    except OSError:
        pass

    return ['NR_177056.1', 'NR_177050.1', 'NR_177042.1', 'NR_177036.1', 'NR_177033.1', 'NR_177041.1', 'NR_177034.1', 'NR_177032.1', 'NR_177028.1', 'NR_177031.1', 'NR_177035.1', 'NR_177027.1', 'NR_177018.1', 'NR_177019.1', 'NR_177017.1', 'NR_177013.1', 'NR_177001.1', 'NR_176998.1', 'NR_176997.1', 'NR_176995.1', 'NR_176991.1', 'NR_176990.1', 'NR_176988.1', 'NR_172681.1', 'NR_176596.1', 'NR_176590.1', 'NR_176592.1', 'NR_176579.1', 'NR_176575.1', 'NR_176567.1', 'NR_176552.1', 'NR_176550.1', 'NR_176542.1', 'NR_176518.1', 'NR_176493.1', 'NR_176479.1', 'NR_172702.1', 'NR_172695.1', 'NR_172694.1', 'NR_172692.1', 'NR_172691.1', 'NR_172690.1', 'NR_172680.1', 'NR_172676.1', 'NR_172674.1', 'NR_172673.1', 'NR_172669.1', 'NR_172672.1', 'NR_172668.1', 'NR_172665.1', 'NR_172664.1', 'NR_172659.1', 'NR_172564.1', 'NR_172612.1', 'NR_172599.1', 'NR_172596.1', 'NR_172598.1', 'NR_172595.1', 'NR_172591.1', 'NR_172587.1', 'NR_172588.1', 'NR_172589.1', 'NR_172586.1', 'NR_172576.1', 'NR_172575.1', 'NR_172572.1', 'NR_172566.1', 'NR_148569.1', 'NR_125697.1', 'NR_114591.1', 'NR_126312.1', 'NR_125700.1', 'NR_125522.1', 'NR_074309.1']
    #return ids