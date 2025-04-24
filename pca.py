import pandas as pd
import pykrev as pk
import numpy as np
import seaborn as sns
from sklearn import decomposition
from matplotlib import pyplot as plt
import argparse
import yaml
import os

def main(config):

    FILEPATH = f"{config['PCAPATH']}{config['COLLECTION_ID']}"
    COLLECTION_ID = config['COLLECTION_ID']
    ORDPATH = config["ORDPATH"]
    METAPATH = config["METAPATH"]
    NORM_METHOD = config["NORM_METHOD"]
    NORM_TRANSFORM = config["NORM_TRANSFORM"]
    COMPONENTS = config['COMPONENTS']
    PLOT = config['PLOT']

    os.makedirs(FILEPATH, exist_ok=True)

    #read the ordination data 
    ordinationData = pd.read_csv(f"{ORDPATH}/{COLLECTION_ID}_ordination.csv",index_col=0)

    #zero-fill the ordination data 
    zeroFill = ordinationData.fillna(0)

    #normalise the ordination data
    normalise = pk.normalise_intensity(zeroFill,norm_method=NORM_METHOD,norm_transform=NORM_TRANSFORM)

    #read the metadata
    metaData = pd.read_csv(f"{METAPATH}/{COLLECTION_ID}_metadata.csv",index_col=0)

    #create an empty dataframe to make the pca plots in
    plotData = pd.DataFrame()
    plotData['sample_id'] = ordinationData.index
    plotData = pd.merge(plotData,metaData, how='left', left_on='sample_id', right_on='ID')

    #do the PCA
    PCA = decomposition.PCA()
    pca_model = PCA.fit(normalise)
    pca_result = pca_model.transform(normalise)

    #plot data update 
    for i in range (0,COMPONENTS):
        plotData[f'PC{i+1}'] = pca_result[:,i]
    plotData.to_csv(f"{FILEPATH}/pca_data.csv")

    #save svg biplots
    if PLOT: 
        for i in range (1,COMPONENTS):
            for j in range(i+1,COMPONENTS+1):
                ax = sns.scatterplot(x=f'PC{i}',y=f'PC{j}',data=plotData, hue='Instrument setup used', style='Species group')
                ax.set_xlabel(f'PC{i} ({np.round(pca_model.explained_variance_ratio_[i-1],3)}%)')
                ax.set_ylabel(f'PC{j} ({np.round(pca_model.explained_variance_ratio_[j-1],3)}%)')
                ax.set_title(f'{COLLECTION_ID} PCA plot: PC{i}/PC{j}')
                plt.savefig(f'{FILEPATH}/PC{i}_PC{j}.svg', bbox_inches='tight')
                plt.close()

#COMMAND LINE
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download all the DSFP detections in a collection and create an ordination matrix.')
    parser.add_argument('config', type=str, help='Path to the YAML configuration file')
    
    args = parser.parse_args()
    
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    
    main(config)

