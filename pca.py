import pandas as pd
import pykrev as pk
import numpy as np
import seaborn as sns
from sklearn import decomposition
from matplotlib import pyplot as plt
import argparse
import yaml
import os
from susdat import susdat

def main(config):

    FILEPATH = f"{config['PCAPATH']}{config['COLLECTION_ID']}"
    COLLECTION_ID = config['COLLECTION_ID']
    ORDPATH = config["ORDPATH"]
    METAPATH = config["METAPATH"]
    SUBSET = config["SUBSET"]
    NORM_ORDER = config["NORM_ORDER"]
    ROW_METHOD = config["ROW_METHOD"]
    ROW_TRANSFORM = config["ROW_TRANSFORM"]
    COL_METHOD = config["COL_METHOD"]
    COL_TRANSFORM = config["COL_TRANSFORM"]
    COMPONENTS = config['COMPONENTS']
    LOADINGS = config['LOADINGS']
    SUSDAT = config['SUSDAT']
    PLOT = config['PLOT']
    HUE = config['HUE']
    STYLE = config['STYLE']

    def PCA_ROUTINE(g=''):
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
        plotData.to_csv(f"{FILEPATH}/{g}pca_data.csv")

        #save svg biplots
        if PLOT: 
            for i in range (1,COMPONENTS):
                for j in range(i+1,COMPONENTS+1):
                    ax = sns.scatterplot(x=f'PC{i}',y=f'PC{j}',data=plotData, hue=HUE, style=STYLE)
                    ax.set_xlabel(f'PC{i} ({np.round(pca_model.explained_variance_ratio_[i-1]*100,3)}%)')
                    ax.set_ylabel(f'PC{j} ({np.round(pca_model.explained_variance_ratio_[j-1]*100,3)}%)')
                    ax.set_title(f'{COLLECTION_ID} PCA plot: PC{i}/PC{j}')
                    plt.savefig(f'{FILEPATH}/{g}PC{i}_PC{j}.svg', bbox_inches='tight')
                    plt.close()

        #save loadings plots
        for i in range(0,COMPONENTS):
            top_N_idx = np.argsort(abs(pca_model.components_[i,:]))[-LOADINGS:]
            scores = pca_model.components_[i,:][top_N_idx]
            mols = ordinationData.columns[top_N_idx]
            argidx = np.argsort(scores)
            scores = scores[argidx]
            mols = mols[argidx]
            loadings= pd.DataFrame()
            loadings['scores'] = scores

            if (SUSDAT != False): 
                sus = susdat(mols,config)
                loadings['mols'] = sus
            else:
                loadings['mols'] = mols

            loadings.to_csv(f"{FILEPATH}/{g}PC{i+1}_loadings.csv")
            if PLOT:
                if i % 2 == 0:
                    ax = sns.stripplot(x=f'PC{i+1}',data=plotData, alpha=.4, marker='X', legend=False,c='grey')
                    sns.stripplot(x="scores",ax=ax, hue="mols", palette='coolwarm',legend=True,data=loadings)
                else:
                    ax = sns.stripplot(y=f'PC{i+1}',data=plotData, alpha=.4, marker='X', legend=False,c='grey')
                    sns.stripplot(y="scores",ax=ax, hue="mols", palette='coolwarm',legend=True, data=loadings)
                # Adjusting the legend
                handles, labels = ax.get_legend_handles_labels()
                ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, title='loadings')
                ax.set_title(f'PC{i+1} loading plot')
                plt.savefig(f'{FILEPATH}/{g}PC{i+1}_loadings.svg', bbox_inches='tight')
                plt.close()

    
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
            #zero-fill the ordination data 
            zeroFill = ordinationData.fillna(0)
            #normalise the ordination data
            if NORM_ORDER == 'None':
                normalise = zeroFill
            elif NORM_ORDER == 'row':
                normalise = pk.normalise_intensity(zeroFill,norm_method=ROW_METHOD,norm_transform=ROW_TRANSFORM,norm_direction='rows')
            elif NORM_ORDER == 'col':
                normalise = pk.normalise_intensity(zeroFill,norm_method=COL_METHOD,norm_transform=COL_TRANSFORM,norm_direction='columns')
            elif NORM_ORDER == 'rowcol':
                normalise = pk.normalise_intensity(zeroFill,norm_method=ROW_METHOD,norm_transform=ROW_TRANSFORM,norm_direction='rows')
                normalise = pk.normalise_intensity(zeroFill,norm_method=COL_METHOD,norm_transform=COL_TRANSFORM,norm_direction='columns')
            elif NORM_ORDER == 'colrow':
                normalise = pk.normalise_intensity(zeroFill,norm_method=COL_METHOD,norm_transform=COL_TRANSFORM,norm_direction='columns')
                normalise = pk.normalise_intensity(zeroFill,norm_method=ROW_METHOD,norm_transform=ROW_TRANSFORM,norm_direction='rows')
            #do the PCA
            PCA_ROUTINE(g=f"{g}/")

    else: 
        #read the metadata
        metaData = pd.read_csv(f"{METAPATH}/{COLLECTION_ID}_metadata.csv",index_col=0)
        #read the ordination data 
        ordinationData = pd.read_csv(f"{ORDPATH}/{COLLECTION_ID}_ordination.csv",index_col=0)
        #zero-fill the ordination data 
        zeroFill = ordinationData.fillna(0)
        #normalise the ordination data
        if NORM_ORDER == 'None':
            normalise = zeroFill
        elif NORM_ORDER == 'row':
            normalise = pk.normalise_intensity(zeroFill,norm_method=ROW_METHOD,norm_transform=ROW_TRANSFORM,norm_direction='rows')
        elif NORM_ORDER == 'col':
            normalise = pk.normalise_intensity(zeroFill,norm_method=COL_METHOD,norm_transform=COL_TRANSFORM,norm_direction='columns')
        elif NORM_ORDER == 'rowcol':
            normalise = pk.normalise_intensity(zeroFill,norm_method=ROW_METHOD,norm_transform=ROW_TRANSFORM,norm_direction='rows')
            normalise = pk.normalise_intensity(zeroFill,norm_method=COL_METHOD,norm_transform=COL_TRANSFORM,norm_direction='columns')
        elif NORM_ORDER == 'colrow':
            normalise = pk.normalise_intensity(zeroFill,norm_method=COL_METHOD,norm_transform=COL_TRANSFORM,norm_direction='columns')
            normalise = pk.normalise_intensity(zeroFill,norm_method=ROW_METHOD,norm_transform=ROW_TRANSFORM,norm_direction='rows')
        #do the PCA
        PCA_ROUTINE()

#COMMAND LINE
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download all the DSFP detections in a collection and create an ordination matrix.')
    parser.add_argument('config', type=str, help='Path to the YAML configuration file')
    
    args = parser.parse_args()
    
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    
    main(config)

