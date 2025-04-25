import requests 
from concurrent.futures import ThreadPoolExecutor, as_completed
import tqdm
import os 
import json

def susdat(mols, config):
    """
    This function searches the SUSDAT database for each NS ID in mols and returns a list with the SUSDAT string instead. 
    Input is a list of NS ID strings. 
    Output is a list of SUSDAT strings.
    IF CACHE the function will cache downloaded json files in SUSPATH and search SUSPATH for previously downloaded files. 
    """
    SUSDAT = config['SUSDAT']
    THREADS = config['THREADS']
    CACHE = config['CACHE']
    FILEPATH = config['SUSPATH']

    if CACHE:
        os.makedirs(FILEPATH, exist_ok=True)

    def parsedb(m):
        if f"{m}.json" in savedfiles and CACHE:
            with open(f"{FILEPATH}{m}.json", 'r') as f:
                jsonfile = json.load(f)
        else: 
                # convert mols to IUPAC using susdat look up
                url = f"https://www.norman-network.com/nds/api/susdat/nsid/{m}/JSON"
                r = requests.get(url)
                jsonfile = r.json()
                #save the json file for future caching 
                with open(f"{FILEPATH}{m}.json", 'w') as f:
                    json.dump(jsonfile, f)
                savedfiles.append(f"{m}.json")
        try:
            sus.append(jsonfile[SUSDAT])
        except KeyError:
            sus.append(m)

    sus = []
    savedfiles = os.listdir(FILEPATH)
    pbar = tqdm.tqdm(total=len(mols), desc='susdat')  # Init pbar
    with ThreadPoolExecutor(max_workers=THREADS) as executor:
        futures = [executor.submit(parsedb, m) for m in mols]
        for _ in as_completed(futures):
            pbar.update(n=1)
    return sus
