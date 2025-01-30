import re
import requests
import time
from requests.adapters import HTTPAdapter, Retry

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


# Curl is of the form: %
#curl --request POST 'https://rest.uniprot.org/idmapping/run' --form 'ids="P21802,P12345"' --form 'from="UniProtKB_AC-ID"' --form 'to="UniRef90"'

import polars as pl
import ast

#dat = pl.read_csv("../data/raw/smoke_dataset2.csv")
dat = pl.read_csv("../data/raw/Deloitte_DrugDiscovery_dataset.csv")


prot_ids = set(filter(None, dat["UniProt_ID"].to_list()))

data = {
    "ids" : ','.join(prot_ids),
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
#        status = ast.literal_eval(response.headers)["jobStatus"]
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

with open("../data/processed/all_prots.csv", "w") as f:
    for batch, total in get_batch(url):
        if not progress:
            f.write("UniProt_ID, Sequence\n")
        data = batch.json()
        for element in data["results"]:
            f.write(element["from"] + ", " + element["to"]["sequence"]["value"] + "\n")
            progress += 1
            print('{}/{} Downloaded'.format(progress,len(prot_ids)))

print("done. Proteins written to file")
        



#{"jobStatus":"FINISHED"}

# url = 'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Ccc_interaction&format=tsv&query=Insulin%20AND%20%28reviewed%3Atrue%29&size=500'
# progress = 0
# with open('insulin-interactions.tsv', 'w') as f:
#     for batch, total in get_batch(url):
#         lines = batch.text.splitlines()
#         if not progress:
#             print(lines[0], file=f)
#         for line in lines[1:]:
#             print(line, file=f)
#         progress += len(lines[1:])
#         print(f'{progress} / {total}')


