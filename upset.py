import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import argparse
import yaml
import os
from susdat import susdat
import upsetplot


def main(config):

    FILEPATH = f"{config['UPSET_PATH']}{config['COLLECTION_ID']}"
    COLLECTION_ID = config['COLLECTION_ID']
    ORDPATH = config["ORDPATH"]
    METAPATH = config["METAPATH"]
    SUBSET = config["SUBSET"]
    PLOT = config['PLOT']
    SUSDAT = config["SUSDAT"]
    GROUPING = config["GROUPING"]
    MAX_SUBSET_RANK = config["MAX_SUBSET_RANK"]
    MIN_SUBSET_SIZE = config["MIN_SUBSET_SIZE"]
    if MAX_SUBSET_RANK == False:
        MAX_SUBSET_RANK = None
    if MIN_SUBSET_SIZE == False:
        MIN_SUBSET_SIZE = None

    def upset_ROUTINE(g=''):
        if (SUSDAT != False): 
            sus = susdat(ordinationData.columns,config)
            ordinationData.columns = sus
        else:
            pass

        #merge the data
        if GROUPING != "ID":
            mergeData = pd.merge(ordinationData,metaData.filter(['ID',GROUPING]),how='left',left_index=True,right_on='ID')
            mergeData.drop(columns=["ID"], inplace=True)

        else: 
            mergeData = pd.merge(ordinationData,metaData.filter([GROUPING]),how='left',left_index=True,right_on='ID')

        #make a dictionary with a set of the unique compounds for each species group
        groupDict = {}
        for group in mergeData[GROUPING].unique():
            if type(group) == float:
                if np.isnan(group):
                    groupData = mergeData.loc[metaData[GROUPING].isnull(),:]
                else:
                    groupData = mergeData[mergeData[GROUPING] == group]
            else:
                groupData = mergeData[mergeData[GROUPING] == group]
            groupData = groupData.drop(columns=[GROUPING])
            groupDict[group] = set(groupData.columns[groupData.sum(axis=0) > 0])
        #do the upset
        fromContents = upsetplot.from_contents(groupDict)
        fromContents.to_csv(f'{FILEPATH}/{g}fromContents.csv')
        upset = upsetplot.UpSet(fromContents, show_counts=True, show_percentages = True, max_subset_rank=MAX_SUBSET_RANK, min_subset_size=MIN_SUBSET_SIZE)
        if PLOT: 
            upset.plot()
            plt.savefig(f'{FILEPATH}/{g}upset_plot.svg', bbox_inches='tight')

    os.makedirs(FILEPATH, exist_ok=True)
    #subset the data 
    if SUBSET != 'none':
        #read the metadata
        metaData = pd.read_csv(f"{METAPATH}/{COLLECTION_ID}_metadata.csv",index_col=0)
        #determine the groupings
        groups = set(metaData[SUBSET])
        for g in groups:
            #make dirs
            os.makedirs(f"{FILEPATH}/{g}/", exist_ok=True)
            #read the ordination data 
            ordinationData = pd.read_csv(f"{ORDPATH}/{COLLECTION_ID}_ordination.csv",index_col=0)
            #re-read the metadata
            metaData = pd.read_csv(f"{METAPATH}/{COLLECTION_ID}_metadata.csv",index_col=0)
            #subset the metadata
            metaData = metaData[metaData[SUBSET] == g]
            #subset the ordination data 
            ordinationData = ordinationData[ordinationData.index.isin(metaData["ID"])]
            #do the kmeans
            upset_ROUTINE(g=f"{g}/")
    
    else: 
        #read the metadata
        metaData = pd.read_csv(f"{METAPATH}/{COLLECTION_ID}_metadata.csv",index_col=0)
        #read the ordination data 
        ordinationData = pd.read_csv(f"{ORDPATH}/{COLLECTION_ID}_ordination.csv",index_col=0)
        #do the kmeans
        upset_ROUTINE()

#COMMAND LINE
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform and visualise upsetplot analysis')
    parser.add_argument('config', type=str, help='Path to the YAML configuration file')
    
    args = parser.parse_args()
    
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    
    main(config)

