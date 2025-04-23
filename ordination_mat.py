# LIBRARIES
import os
import requests
import json
import pandas as pd
import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import yaml
import argparse

# MAIN
def main(config):
    DOWNLOAD = config['DOWNLOAD']
    METHOD = config['METHOD']
    THREADS = config['THREADS']
    FILEPATH = config['ORDPATH']
    DOWNLOAD_DIR = f"{config['DOWNLOAD_DIR']}{config['COLLECTION_ID']}"
    COLLECTION_ID = config['COLLECTION_ID']
    NS_MIN = config['NS_MIN']
    NS_MAX = config['NS_MAX']

    compounds = [f"NS{str(i).zfill(8)}" for i in range(NS_MIN, NS_MAX)]# full 

    # Download files
    url = "https://files.dsfp.norman-data.eu/detections/"
    urls = [f"{url}{COLLECTION_ID}/{COLLECTION_ID}_{compound}.json" for compound in compounds]
    os.makedirs(DOWNLOAD_DIR, exist_ok=True)
    os.makedirs(FILEPATH, exist_ok=True)
    if DOWNLOAD:
        pbar = tqdm.tqdm(total=len(urls), desc='download')  # Init pbar
        with ThreadPoolExecutor(max_workers=THREADS) as executor:
            futures = [executor.submit(download_file, url, DOWNLOAD_DIR) for url in urls]
            for _ in as_completed(futures):
                pbar.update(n=1)

    # Read files and save in a list object
    local_urls = [os.path.join(DOWNLOAD_DIR, file) for file in os.listdir(DOWNLOAD_DIR)]
    print(local_urls)
    data_as_list = []
    pbar = tqdm.tqdm(total=len(local_urls), desc='load')  # Init pbar
    with ThreadPoolExecutor(max_workers=THREADS) as executor:
        futures = [executor.submit(load_file, local_url, data_as_list) for local_url in local_urls]
        for _ in as_completed(futures):
            pbar.update(n=1)

    #create the ordination matrix
    samples = list(set([item['sample_id'] for sublist in data_as_list for item in sublist]))
    compounds = list(set([item['substance_id'] for sublist in data_as_list for item in sublist]))
    screening_results = pd.DataFrame(columns=compounds,index=samples)
    pbar = tqdm.tqdm(total=len(data_as_list), desc='ordinate')  # Init pbar
    with ThreadPoolExecutor(max_workers=THREADS) as executor:
        futures = [executor.submit(ordinate, detections,METHOD, screening_results) for detections in data_as_list]
        for _ in as_completed(futures):
            pbar.update(n=1)
    screening_results.to_csv(f"{FILEPATH}/{COLLECTION_ID}_ordination.csv")

#FUNCTIONS 
def download_file(url, DOWNLOAD_DIR):
    file_name = os.path.join(DOWNLOAD_DIR, os.path.basename(url))
    response = requests.get(url)
    data = response.json()
    # only write the data if there is content
    if len(data) > 0:
        with open(file_name, 'wb') as f:
            f.write(response.content)
        return file_name
    else:
         return 0
    
def load_file(local_url, data_as_list):
    with open(local_url, 'r') as f:
        data = json.load(f)
    if len(data) > 0:
        data_as_list.append(data)

def ordinate(detections,METHOD,screening_results):
    detections = pd.DataFrame(detections)
    for _,sample in detections.iterrows(): 
                #check there is at least one match 
                if len(sample["matches"]) > 0:
                    if METHOD == 'sum':
                        screening_results.at[sample['sample_id'],sample['substance_id']] = 0
                        for m in range(0,len(sample["matches"])):
                            screening_results.at[sample['sample_id'],sample['substance_id']] = screening_results.at[sample['sample_id'],sample['substance_id']] + sample["matches"][m]["peak_area"]  
                    elif METHOD == 'max':
                        areas = []
                        for m in range(0,len(sample["matches"])):
                            areas.append(sample["matches"][m]["peak_area"])
                        screening_results.at[sample['sample_id'],sample['substance_id']] = max(areas)

#COMMAND LINE
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download all the DSFP detections in a collection and create an ordination matrix.')
    parser.add_argument('config', type=str, help='Path to the YAML configuration file')
    
    args = parser.parse_args()
    
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    
    main(config)
