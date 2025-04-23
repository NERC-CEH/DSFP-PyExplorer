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

    os.makedirs(FILEPATH, exist_ok=True)

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

    for iid in dataframe['internal_id'].unique():
        urls.append(f'https://dsfp.norman-data.eu/data/{iid}/samples-biota.csv')
        iids.append(iid)

    pbar = tqdm.tqdm(total=len(urls), desc='download')  # Init pbar
    with ThreadPoolExecutor(max_workers=THREADS) as executor:
        futures = [executor.submit(download_file, url, iid, FILEPATH) for url, iid in zip(urls,iids)]
        for _ in as_completed(futures):
            pbar.update(n=1)
  
#FUNCTIONS 
def download_file(url, iid, FILEPATH):
        r = requests.get(url)
        csv_file = StringIO(r.content.decode('utf-8'))
        samples = pd.read_csv(csv_file)
        samples.to_csv(f'{FILEPATH}/{iid}_metadata.csv')

#COMMAND LINE
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download all the DSFP detections in a collection and create an ordination matrix.')
    parser.add_argument('config', type=str, help='Path to the YAML configuration file')
    
    args = parser.parse_args()
    
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    
    main(config)


