import re
import requests
import time
from requests.adapters import HTTPAdapter, Retry
import xml.etree.ElementTree as ET
import polars as pl
import ast
import urllib.request
import gzip
import shutil
import argparse


# The block of code below is provided by uniprot to help with accessing their api and was originally founds here: https://www.uniprot.org/help/api_queries
# It is provided without further comment. 
### START ### 
re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)
### END ###


# Curl requests to the uniprot endpoint should be of the form: 
#curl --request POST 'https://rest.uniprot.org/idmapping/run' --form 'ids="P21802,P12345"' --form 'from="UniProtKB_AC-ID"' --form 'to="UniRef90"'

def convertUniProt(uniprot_ids):
    """
    Convert list of UniProt IDs to sequences using their format conversion service.
    """

    data = {
    "ids" : ','.join(uniprot_ids),
    "from" : "UniProtKB_AC-ID",
    "to" : "UniProtKB"
    }

    print("Sending request to UniProt...")
    initial_request = requests.post("https://rest.uniprot.org/idmapping/run", data = data)

    job_id = ast.literal_eval(initial_request.text)["jobId"]

    print("Done. Job ID is {}".format(job_id))

    # Poll until response is complete

    status = ''
    while status != "FINISHED":
        response = requests.get('https://rest.uniprot.org/idmapping/status/' + job_id, allow_redirects = False)
        try:
            status = response.json()["jobStatus"]
        except:
            print("Failed to get response")
            print(response.text)
            quit()
        print("Waiting...")
        time.sleep(5)

    print("Got response!")

    print("Attempting to download and write to file")

    url = "https://rest.uniprot.org/idmapping/uniprotkb/results/" + job_id + "?format=json&size=500&fields=sequence"

    progress = 0

    testout = "../data/interim/TEST_prots.csv"
    with open(testout, "w") as f:
        for batch, total in get_batch(url):
            if not progress:
                f.write("UniProt_ID,Sequence\n")
            data = batch.json()
            for element in data["results"]:
                f.write(element["from"] + "," + element["to"]["sequence"]["value"] + "\n")
                progress += 1
                print('{}/{} Downloaded'.format(progress,len(uniprot_ids)))

    print("Done. Proteins written to file at {} ".format(testout))
    


def convertPubChem(pubchem_ids):
    """
    Convert list of PubChem CIDs to SMILES using their format conversion service
    in the power user gateway
    """

    # Create the submission file
    with open("../data/external/pubchem_id_exchange_template.xml", "r") as fin, open("../data/interim/TEST_pubchem_conversion_sub.xml", "w") as fout:
        for line in fin:
            if line.strip() == "*"*10:
                for cid in pubchem_ids:
                    ws = re.match(r"\s*",line).group()
                    fout.write(ws+"<PCT-ID-List_uids_E>"+str(cid)+"</PCT-ID-List_uids_E>\n")
                continue
            fout.write(line)

    # Submit the file
    with open("../data/interim/TEST_pubchem_conversion_sub.xml", "r") as f:
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
        print("Unhandled exception. Try again later. Shutting down")
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

    # FTP is increasingly unsupported
    link = status_str[0][0][1][0][0][0].text.strip().replace("ftp://", "https://")


    # Download file

    test_smiles_location = "../data/interim/TEST_smiles_mapping.txt"
    with urllib.request.urlopen(link) as compressed_in, open(test_smiles_location, "wb") as smiles_out:
        with gzip.GzipFile(fileobj = compressed_in) as smiles_in:
            shutil.copyfileobj(smiles_in, smiles_out)

    print("saved smiles mapping to {}".format(test_smiles_location))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Lookup SMILES and protein sequences from a csv containing columns for UniProt_ID, pubchem_cid, kiba_score, and kiba_score estimated. Columns need not be in the same order. Processed files will be written to predetermined locations for ease of testing")
    parser.add_argument("filename", help = "Absolute path to test file")

    args = parser.parse_args()

    dat = pl.read_csv(args.filename)

    dat = dat.with_columns(pl.col("pubchem_cid").cast(pl.Int64, strict = True))

    pub_id = set(dat["pubchem_cid"].to_list())

    convertPubChem(pub_id)

    uni_id = set(dat["UniProt_ID"].to_list())

    convertUniProt(uni_id)

    print("SMILES and Sequences saved to file")

    quit()




