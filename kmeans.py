import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import argparse
import yaml
import os
from susdat import susdat
import math
from sklearn.cluster import KMeans
import matplotlib.patches as mpatches


def main(config):

    FILEPATH = f"{config['KMEANS_PATH']}{config['COLLECTION_ID']}"
    COLLECTION_ID = config['COLLECTION_ID']
    ORDPATH = config["ORDPATH"]
    METAPATH = config["METAPATH"]
    SUBSET = config["SUBSET"]
    PLOT = config['PLOT']
    HUE = config['HUE']
    SUSDAT = config["SUSDAT"]
    KMEANS_P = config["KMEANS_P"]
    KMEANS_CLUSTERS = config["KMEANS_CLUSTERS"]
    LEGEND_BBOX_X = config["LEGEND_BBOX_X"]
    LEGEND_BBOX_Y = config["LEGEND_BBOX_Y"]

    def kmeans_ROUTINE(g=''):
        #filter out empty or completely full columns
        N = ordinationData.shape[0]
        P = KMEANS_P # the required proportion of data
        M = int(np.floor(N * P)) #the minimum number of recorded values for a sample (round down)
        U =  int(np.ceil(N * (1-P))) #the maximum number of recorded values for a sample (round up)
        boolA = (ordinationData.notnull().sum(axis=0) >= M) & (ordinationData.notnull().sum(axis=0) <= U)
        filterData = ordinationData.loc[:,boolA]
        C = filterData.shape[1]
        columns = filterData.columns


        plotData = pd.DataFrame()
        plotData['sample_id'] = ordinationData.index
        plotData = pd.merge(plotData,metaData, how='left', left_on='sample_id', right_on='ID')

        #figure out how many pairwise combinations are possible 
        C_combinations = math.factorial(C)/(2*math.factorial(C-2))
        shuffles = C_combinations/C
        #figure out how many times you need to shuffle the isnull matrix, round up 
        shuffles = int(np.round(shuffles))

        #now we are going to filter out any compounds that don't have a complementary match
        notnull = np.array(filterData.notnull(),dtype=int)
        notnull_shift = np.array(filterData.notnull(),dtype=int)
        exclusives = set()
        for i in range (0,shuffles):
            notnull_shift = np.concatenate((notnull_shift[:, -1:], notnull_shift[:, :-1]), axis=1)
            boolB = np.max(notnull + notnull_shift,axis=0) == 1
            exclusives.update(columns[boolB])
        exclusiveData = ordinationData.loc[:,list(exclusives)]

        #now we are going to cluster this filtered data using kmeans
        kmeans = KMeans(n_clusters=KMEANS_CLUSTERS)
        clusters = kmeans.fit_predict(exclusiveData.notnull().transpose())

        #create a CSV file with the clusters in it
        clusterCSV = pd.DataFrame({"compounds":exclusiveData.columns,"cluster":clusters})
        if (SUSDAT != False): 
            sus = susdat(clusterCSV['compounds'],config)
            clusterCSV['compound_name'] = sus
        else:
            pass
        clusterCSV.to_csv(f'{FILEPATH}/{g}kmeans_clusters.csv', index=False)

        #save svg biplots
        if PLOT: 
            # Sort columns based on cluster membership 
            sorted_indices = np.argsort(clusters)
            sorted_matrix = exclusiveData.iloc[:, sorted_indices]
            #if hue not specified then do a black and white plot
            if HUE == 'none':
                plt.imshow(sorted_matrix.isnull(),cmap='gray', interpolation='None')
                plt.xlabel('Compound')
                plt.ylabel('Sample')
                black_patch = mpatches.Patch(color='black', label='Present')
                white_patch = mpatches.Patch(color='white', label='Missing')
                plt.legend(handles=[black_patch, white_patch], bbox_to_anchor=(LEGEND_BBOX_X, LEGEND_BBOX_Y))

                # Define the number of columns for each label and annotate the plots 
                for i in range(0,KMEANS_CLUSTERS):
                        X = sum(clusters < i)
                        N = sum(clusters == i)
                        plt.text(X + N/2 - 0.5, -20, f'Cluster {i}', ha='center', va='center', fontsize=10)

                
                # Save the plot
                plt.savefig(f'{FILEPATH}/{g}kmeans_clusters.svg', bbox_inches='tight')
                plt.close()

            #otherwise colour the rows based on hue
            else: 
                 # Create a color map for the HUEs
                hue_groups = plotData[HUE].unique()
                colors = plt.cm.tab10(np.linspace(0, 1, len(hue_groups)))

                # Create a color dictionary for HUEs
                color_dict = {hue_groups[i]: colors[i] for i in range(len(hue_groups))}

                # Create a matrix to hold the colors
                colored_matrix = np.empty(shape=(sorted_matrix.shape[0],sorted_matrix.shape[1],3))

                # Apply colors to the matrix based on HUE
                for i in range(sorted_matrix.shape[0]):
                    for j in range(sorted_matrix.shape[1]):
                        if pd.notnull(sorted_matrix.iloc[i, j]):
                            colored_matrix[i, j, 0] = color_dict[plotData[HUE][i]][0]
                            colored_matrix[i, j, 1] = color_dict[plotData[HUE][i]][1]
                            colored_matrix[i, j, 2] = color_dict[plotData[HUE][i]][2]
                        else:
                            colored_matrix[i, j, 0] = 1
                            colored_matrix[i, j, 1] = 1
                            colored_matrix[i, j, 2] = 1

                # Plot the matrix with the categorical colors
                plt.imshow(colored_matrix, interpolation='None')
                plt.xlabel('Compound')
                plt.ylabel('Sample')

                # Add a legend slightly off the figure to one side
                legend_patches = [mpatches.Patch(color='white', label='Missing')] + [mpatches.Patch(color=color_dict[group], label=group) for group in hue_groups]
                plt.legend(handles=legend_patches, bbox_to_anchor=(LEGEND_BBOX_X, LEGEND_BBOX_Y))

                # Define the number of columns for each label and annotate the plots 
                for i in range(0,KMEANS_CLUSTERS):
                        X = sum(clusters < i)
                        N = sum(clusters == i)
                        plt.text(X + N/2 - 0.5, -20, f'Cluster {i}', ha='center', va='center', fontsize=10)

                # Save the plot
                plt.savefig(f'{FILEPATH}/{g}kmeans_clusters.svg', bbox_inches='tight')
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
            #do the kmeans
            kmeans_ROUTINE(g=f"{g}/")
    
    else: 
        #read the metadata
        metaData = pd.read_csv(f"{METAPATH}/{COLLECTION_ID}_metadata.csv",index_col=0)
        #read the ordination data 
        ordinationData = pd.read_csv(f"{ORDPATH}/{COLLECTION_ID}_ordination.csv",index_col=0)
        #do the kmeans
        kmeans_ROUTINE()

#COMMAND LINE
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform and visualise KMEANS clustered missing data')
    parser.add_argument('config', type=str, help='Path to the YAML configuration file')
    
    args = parser.parse_args()
    
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    
    main(config)

