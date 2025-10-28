#libraries
import requests
import pandas as pd
from io import StringIO
import yaml
import argparse
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
import tqdm

def main(config):
    FILEPATH = config['METAPATH']
    THREADS = config['THREADS']
    METACACHE = config['METACACHE']

    os.makedirs(FILEPATH, exist_ok=True)
    cachedFiles = os.listdir(FILEPATH)
    #1 get collections
    url="https://dsfp.norman-data.eu/api/1/metastore/schemas/dataset/all"
    r = requests.get(url)
    json = r.json()
    dataframe = pd.DataFrame(json)
    print(f"There are {len(dataframe['internal_id'].unique())} collections.")
    print(f"There are {len(dataframe['participating_labs'].unique())} participating labs.")
    dataframe

    # loop
    # 2 get collection files
    urls = []
    iids = []
    mtypes = []
    fails = []
    mlists = {}

    for iid in dataframe['internal_id'].unique():
        matrixTypes=[]
        # check if we've already downloaded the file
        if METACACHE:
            if f"{iid}_metadata.csv" in cachedFiles:
                continue
        # figure out the matrix type 
        try: 
            r = requests.get(f'https://dsfp.norman-data.eu/data/{iid}/samples.json')
            samplesjson = pd.json_normalize(r.json())
            matrixTypes = list(samplesjson['matrix'].unique())
            mlists[iid] = matrixTypes
        except KeyError:
            fails.append(iid)
            mlists[iid] = []
            continue
        for m in matrixTypes:
            mf = m.replace(' ','-').lower()
            urls.append(f'https://dsfp.norman-data.eu/data/{iid}/samples-{mf}.csv')
            iids.append(iid)
            mtypes.append(mf)

    pbar = tqdm.tqdm(total=len(urls), desc='download')  # Init pbar
    with ThreadPoolExecutor(max_workers=THREADS) as executor:
        futures = [executor.submit(download_file, url, iid, mf, FILEPATH) for url, iid, mf in zip(urls,iids,mtypes)]
        for _ in as_completed(futures):
            pbar.update(n=1)

    print('consolidating files ...')
    # now go through and consolidate each collections files into a single .csv file 
    for iid in mlists.keys():
        mtypes = mlists[iid]
        baseDf = pd.DataFrame()
        for mt in mtypes:
            mt = mt.replace(' ','-').lower()
            try: 
                mtdf = pd.read_csv(f'{FILEPATH}/{iid}_{mt}_metadata.csv')
                mtdf['matrix'] = mt 
                baseDf = pd.concat([baseDf,mtdf])
            except FileNotFoundError:
                pass
        if len(baseDf) > 0:
            baseDf.to_csv(f'{FILEPATH}/{iid}_metadata.csv', index=False)
    
    print(f'Failed to download the following collections: {fails}.')
  
#FUNCTIONS 
def download_file(url, iid, mf, FILEPATH):
    r = requests.get(url)
    csv_file = StringIO(r.content.decode('utf-8'))
    samples = pd.read_csv(csv_file)
    samples.to_csv(f'{FILEPATH}/{iid}_{mf}_metadata.csv', index=False)

#COMMAND LINE
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download the metadata associated with DSFP collection.')
    parser.add_argument('config', type=str, help='Path to the YAML configuration file')
    
    args = parser.parse_args()
    
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    
    main(config)


