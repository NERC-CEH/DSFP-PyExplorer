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

    FILEPATH = config['PCAPATH']
    COLLECTION_ID = config['COLLECTION_ID']
    ORDPATH = config["ORDPATH"]
    METAPATH = config["METAPATH"]
    NORM_METHOD = config["NORM_METHOD"]
    NORM_TRANSFORM = config["NORM_TRANSFORM"]
    COMPONENTS = config['COMPONENTS']

    os.makedirs(FILEPATH, exist_ok=True)

    #read the ordination data 
    ordinationData = pd.read_csv(f"{ORDPATH/{COLLECTION_ID}}_ordination.csv",index_col=0)

    #zero-fill the ordination data 
    zeroFill = ordinationData.fillna(0)

    #normalise the ordination data
    normalise = pk.normalise_intensity(zeroFill,NORM_METHOD=NORM_METHOD,norm_transform=NORM_TRANSFORM)

    #read the metadata
    metaData = pd.read_csv(f"{METAPATH/{COLLECTION_ID}}_metadata.csv",index_col=0)

    #create an empty dataframe to make the pca plots in
    plotData = pd.DataFrame()
    plotData['sample_id'] = ordinationData.index
    plotData = pd.merge(plotData,metaData, how='left', left_on='sample_id', right_on='ID')

    #do the PCA
    PCA = decomposition.PCA()
    pca_model = PCA.fit(normalise)
    pca_result = pca_model.transform(normalise)

    #plot data update 
    for i in range (0,len(COMPONENTS)):
        plotData[f'PC{i+1}'] = pca_result[:,i]
    plotData.to_csv(f"{FILEPATH}/{COLLECTION_ID}_pca_data.csv")

    ax = sns.scatterplot(x='PC1',y='PC2',data=plotData, hue='Instrument setup used', style='Species group')
    ax.set_xlabel(f'PC1 ({np.round(pca_model.explained_variance_ratio_[0],3)}%)')
    ax.set_ylabel(f'PC2 ({np.round(pca_model.explained_variance_ratio_[1],3)}%)')
    plt.savefig(f'{FILEPATH}{COLLECTION_ID}_pca_plot.svg')


#COMMAND LINE
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download all the DSFP detections in a collection and create an ordination matrix.')
    parser.add_argument('config', type=str, help='Path to the YAML configuration file')
    
    args = parser.parse_args()
    
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    
    main(config)

