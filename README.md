# EEx-NTA
A repository for storing mass-spectrometry non-target analysis (NTA) code linked to the Environmental Exposure Hub. 

## ```python ordination_mat.py config.yaml```

This script is adapted from ```instant_search_processing.R``` by Nikiforos Alygizakis. It downloads all of the detection files in a collections and collates them into a N by M matrix. Where N is the number of samples in the collection and M is the number of detected compounds. Each sample compound pair has its peak area as a value in the matrix. This is ultimately written to a .csv file. 

The ```config.yaml``` file contains the following relevant parameters: 
* ```DOWNLOAD: True``` - Whether to download detection json files
* ```METHOD: 'max'``` - How to deal with compounds that match to multiple peaks in a sample. ```max``` takes the maximum value, ```sum```, sums the values. 
* ```THREADS: 10``` - How many parallel threads to use when downloading/loading/ordinating data
* ```ORDPATH: "ordination/"``` - Filepath to where output data should be written. Output data is written as ```{COLLECTION_ID}_ordination.csv```. 
* ```DOWNLOAD_DIR: "downloads"``` - Where to store/read downloaded json files
* ```COLLECTION_ID: 2380``` - Which collection ID to search
* ```NS_MIN: 1``` - Compounds are searched from MIN:MAX in the Normal SusDat datbase
* ```NS_MAX: 120000```- Compounds are searched from MIN:MAX in the Normal SusDat datbase

## ```python metadata.py config.yaml```

This script downloads metadata for all collections. 

The ```config.yaml``` file contains the following relevant parameters: 
* ```METAPATH: "metadata/"``` - Filepath to where output data should be written. Output data is written as ```{COLLECTION_ID}_metadata.csv```. 

## ```python pca.py config.yaml``` 

This script performs PCA on the ordination matrix produced by ```ordination_mat.py``` and collates the results with the metadata downloaded by ```metadata.py``` to allow PCA biplots to be created. 

The ```config.yaml``` file contains the following relevant parameters: 
* ```COLLECTION_ID: 2380``` - Which collection ID to plot
* ```NORM_METHOD: "sum"``` - How to row-normalise the ordination data before doing PCA. See PyKrev.normalise_intensity() for more info. 
* ```NORM_TRANSFORM: "power3"```- How to row-transform the ordination data before doing PCA. See PyKrev.normalise_intensity() for more info.
* ```COMPONENTS: 4``` - Number of principal components to keep.
* ```LOADINGS: 10``` - Number of loadings to keep/plot. 
* ```PLOT: True``` - Save biplot svgs?
* ```HUE: 'Instrument setup used'``` - Column of the metadata file to colour the points by. 
* ```STYLE: 'Species group'``` - Column of the metadata file to style the points by. 
* ```PCAPATH: "pca/"``` - File path to save the pca data/plots. 

## To do
* cache download of norman susdat data
* parallelise susdat download
* missingno analysis
* random forest predictions