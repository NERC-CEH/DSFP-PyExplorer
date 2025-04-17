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
    FILEPATH = config['FILEPATH']
    DOWNLOAD_DIR = config['DOWNLOAD_DIR']
    COLLECTION_ID = config['COLLECTION_ID']
    NS_MIN = config['NS_MIN']
    NS_MAX = config['NS_MAX']

    compounds = [f"NS{str(i).zfill(8)}" for i in range(NS_MIN, NS_MAX)]# full 

    # Download files
    url = "https://files.dsfp.norman-data.eu/detections/"
    urls = [f"{url}{COLLECTION_ID}/{COLLECTION_ID}_{compound}.json" for compound in compounds]
    os.makedirs(DOWNLOAD_DIR, exist_ok=True)

    if DOWNLOAD:
        pbar = tqdm.tqdm(total=len(urls), desc='download')  # Init pbar
        with ThreadPoolExecutor(max_workers=THREADS) as executor:
            futures = [executor.submit(download_file, url, DOWNLOAD_DIR) for url in urls]
            for _ in as_completed(futures):
                pbar.update(n=1)

    # Read files and save in a list object
    local_urls = [os.path.join(DOWNLOAD_DIR, file) for file in os.listdir(DOWNLOAD_DIR)]
    data_as_list = []

    for idx, local_url in enumerate(local_urls):
            with open(local_url, 'r') as f:
                data = json.load(f)
            if len(data) > 0:
                data_as_list.append(data)

    samples = list(set([item['short_name'] for sublist in data_as_list for item in sublist]))
    compounds = list(set([item['substance_id'] for sublist in data_as_list for item in sublist]))

    screening_results = pd.DataFrame(columns=compounds,index=samples)

    for detections in data_as_list:
                detections = pd.DataFrame(detections)
                for _,sample in detections.iterrows(): 
                            #check there is at least one match 
                            if len(sample["matches"]) > 0:
                                if METHOD == 'sum':
                                    screening_results.at[sample['short_name'],sample['substance_id']] = 0
                                    for m in range(0,len(sample["matches"])):
                                        screening_results.at[sample['short_name'],sample['substance_id']] = screening_results.at[sample['short_name'],sample['substance_id']] + sample["matches"][m]["peak_area"]  
                                elif METHOD == 'max':
                                    areas = []
                                    for m in range(0,len(sample["matches"])):
                                        areas.append(sample["matches"][m]["peak_area"])
                                    screening_results.at[sample['short_name'],sample['substance_id']] = max(areas)

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

#COMMAND LINE
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download all the DSFP detections in a collection and create an ordination matrix.')
    parser.add_argument('config', type=str, help='Path to the YAML configuration file')
    
    args = parser.parse_args()
    
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    
    main(config)
