import re
import requests
import time
import xml.etree.ElementTree as ET
import urllib.request
import gzip
import shutil

counter = 0 
cids = set()
prots = set()
with open("../data/raw/Deloitte_DrugDiscovery_dataset.csv","r") as f:
    next(f)
    for line in f:
        counter+=1
        prot, cid, _,_ = line.split(',')
        if cid:
            cids.add(int(float(cid)))
        if prot:
            prots.add(prot)

total_cids = len(cids)
total_prots = len(prots)

print("Total #unique cids: {}".format(total_cids))
print("Total #unqiue prots: {}".format(total_prots))


import polars as pl

dat = pl.read_csv("../data/raw/Deloitte_DrugDiscovery_dataset.csv")

#dat = pl.read_csv("../data/raw/smoke_dataset.csv")
dat = dat.with_columns(pl.col("pubchem_cid").cast(pl.Int64, strict = True))


print("data read ok")

print("Total len: {}".format(dat.shape[0]))

smiles = set(filter(None,dat["pubchem_cid"].to_list()))


if smiles == total_cids:
    print("sets are equal")
else:
    print(len(smiles))


def pullSmiles(sset):
    """
    Convert list of PubChem CIDs to SMILES using their format conversion service
    in the power user gateway
    """


    # Create the submission file
    with open("../data/external/pubchem_id_exchange_template.xml", "r") as fin, open("../data/interim/pubchem_conversion_sub.xml", "w") as fout:
        for line in fin:
            if line.strip() == "*"*10:
                for cid in sset:
                    ws = re.match(r"\s*",line).group()
                    fout.write(ws+"<PCT-ID-List_uids_E>"+str(cid)+"</PCT-ID-List_uids_E>\n")
                continue
            fout.write(line)

    # Submit the file
    with open("../data/interim/pubchem_conversion_sub.xml", "r") as f:
        payload = f.read()
    
    status = requests.post("https://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi", payload, headers = {"ContentType":"application/xml"})
    
    if str(status) != "<Response [200]>":
        print("Connection error. Try again later. Shutting down")
        quit()

    print("Submitted conversion job to cluster")
        

    status_str = ET.fromstring(status.text)

    try:
        req_id = status_str[0][0][1][0][0][0].text
    except:
        print(status.text)
        print("Terminating")
        quit()

    # Grab poll template
    with open("../data/external/pubchem_poll_template.xml","r") as f:
        poll = f.read()

    poll = poll.replace("*"*10, req_id)
        
    time.sleep(5)
    while status_str[0][0][1][0][0].tag == 'PCT-Waiting':
        status = requests.post("https://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi", poll, headers = {"ContentType":"application/xml"})
        status_str = ET.fromstring(status.text)
        print("Waiting...")
        time.sleep(5)


    print("Got response!")

    link = status_str[0][0][1][0][0][0].text.strip().replace("ftp://", "https://")


    # Download file

    with urllib.request.urlopen(link) as compressed_in, open("../data/interim/smiles_mapping.txt", "wb") as smiles_out:
        with gzip.GzipFile(fileobj = compressed_in) as smiles_in:
            shutil.copyfileobj(smiles_in, smiles_out)

    print("saved output to file")
    print("done")
    
                

pullSmiles(smiles)

