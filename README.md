# DSFP-PyExplorer
A Python package for doing exploratory data analysis of collections on the [NORMAN Digital Sample Freezing Platform](https://dsfp.norman-data.eu/). This work is linked to the UKCEH Environmental Exposure Hub Project. 

![logo](dsfp-pylogo.png)

## ``config.yaml``

The parameters for the scripts are all contained in the config.yaml file. A description of some generic parameters used by most scripts is given here:
* ```COLLECTION_ID: 2380``` - Which collection ID to search
* ```PLOT: True``` - Create svg plots?
* ```HUE: 'Instrument setup used'``` - Column of the metadata file to colour the data by in plots. 
* ```STYLE: 'Species group'``` - Column of the metadata file to colour the data by in plots. 
* ```SUBSET: "Instrument type"``` - Whether to subset the data before performing an analysis. Set to ```none``` to do analysis on all the data. 

## ```python ordination_mat.py config.yaml```

This script is adapted from ```instant_search_processing.R``` by Nikiforos Alygizakis. It downloads all of the detection files in a collections and collates them into a N by M matrix. Where N is the number of samples in the collection and M is the number of detected compounds. Each sample compound pair has its peak area as a value in the matrix. This is ultimately written to a .csv file. 

The ```config.yaml``` file contains the following relevant parameters: 
* ```DOWNLOAD: True``` - Whether to download detection json files
* ```METHOD: 'max'``` - How to deal with compounds that match to multiple peaks in a sample. ```max``` takes the maximum value, ```sum```, sums the values. 
* ```THREADS: 10``` - How many parallel threads to use when downloading/loading/ordinating data
* ```ORDPATH: "ordination/"``` - Filepath to where output data should be written. Output data is written as ```{COLLECTION_ID}_ordination.csv```. 
* ```DOWNLOAD_DIR: "downloads"``` - Where to store/read downloaded json files
* ```NS_MIN: 1``` - Compounds are searched from MIN:MAX in the Normal SusDat datbase
* ```NS_MAX: 120000```- Compounds are searched from MIN:MAX in the Normal SusDat datbase

## ```python metadata.py config.yaml```

This script downloads metadata for all collections. 

The ```config.yaml``` file contains the following relevant parameters: 
* ```METAPATH: "metadata/"``` - Filepath to where output data should be written. Output data is written as ```{COLLECTION_ID}_metadata.csv```. 

## ```python pca.py config.yaml``` 

This script performs principal component analysis (PCA) on the ordination matrix produced by ```ordination_mat.py``` and collates the results with the metadata downloaded by ```metadata.py``` to allow PCA biplots to be created. 

The ```config.yaml``` file contains the following relevant parameters: 
* ```NORM_ORDER: "rowcol"``` - What order to normalise the data in. One of ```row```, ```col```, ```rowcol```, ```colrow``` or ```none```.
* ```ROW_METHOD: "sum"``` - How to row-normalise the ordination data before doing PCA. See PyKrev.normalise_intensity() for more info. 
* ```ROW_TRANSFORM: "power3"```- How to row-transform the ordination data before doing PCA. See PyKrev.normalise_intensity() for more info.
* ```COL_METHOD: "sum"``` - How to col-normalise the ordination data before doing PCA. See PyKrev.normalise_intensity() for more info. 
* ```COL_TRANSFORM: "power3"```- How to col-transform the ordination data before doing PCA. See PyKrev.normalise_intensity() for more info.
* ```COMPONENTS: 4``` - Number of principal components to keep.
* ```LOADINGS: 10``` - Number of loadings to keep/plot. 
* ```PLOT: True``` - Save biplot svgs?
* ```HUE: 'Instrument setup used'``` - Column of the metadata file to colour the points by. 
* ```STYLE: 'Species group'``` - Column of the metadata file to style the points by. 
* ```PCAPATH: "pca/"``` - File path to save the pca data/plots. 
* ```SUBSET: "Instrument type"``` - Whether to subset the data before performing PCA. Set to ```none``` to do PCA on all the data. 

## ```python lda.py config.yaml``` 

This script performs linear discriminant analysis (LDA) on the ordination matrix and collates the results with the metadata downloaded by ```metadata.py``` to allow LDA biplots to be created. 

In addition to the parameters shown above for PCA, the ```config.yaml``` file contains the following relevant parameters: 
* ```LDAPATH: "lda/"``` - File path to save the lda data/plots
* ```LDACLASS: "Species group"``` - Categorical variable to use as the basis of the linear discriminant analysis. 

## ```python kmeans.py config.yaml```

This script clusters compounds based on their presence or absence across samples within a collection. The algorithm first filters compounds based on a data completeness threshold, and then checks that complements of these filtered compounds exist in the data (i.e. there is at least one pair (A,B) such that when compound A is present in samples compound B is not and vice versa). This vastly reduces the size of the dataset, allowing computationally intensive kmeans clustering to be performed on this filtered data. The script returns a csv and graphic of these filtered compounds and their cluster groups. 

The ```config.yaml``` file contains the following relevant parameters: 
* ```KMEANS_PATH: "kmeans/"``` - File path to save the csv and plots to.
* ```KMEANS_P: 0.4``` - The completeness threshold to use. The threshold X is set as (P < X < 1-P). Reducing the value of P will increase the amount of compound in the analysis. Max value is 0.5. 
* ```KMEANS_CLUSTERS: 2``` - The number of kmeans clusters to fit
* ```HUE: 'Instrument setup used'``` - Column of the metadata file to colour the rows in the plot by. 
* ```SUBSET: "Instrument type"``` - Whether to subset the data before performing kmeans. Set to ```none``` to do kmeans on all the data. 

## ```python upset.py config.yaml``` 

This script creates an upset plot showing the intersections of compounds between groups.

The ```config.yaml``` file contains the following relevant parameters: 
* ```UPSET_PATH: "upset/"``` - Where to save the outputs
* ```GROUPING: "Tissue"```- What grouping variable to use from the metadata



## ```susdat.py```

This Python module can convert NS Ids into compound names or any other ID returned by the susdat database. It is used by ```pca.py``` and ```lda.py``` to rename the loadings.
The ```config.yaml``` - file contains the following relevant parameters: 
* ```SUSDAT: "Compound name"``` - which field in the SUSDAT json file to replace the NSID with. Set to "None" if you don't want to change the IDs. 
* ```CACHE: True``` - whether or not to use the downloaded susdat json files cached locally
* ```SUSPATH: "susdata/"``` - where to save/read the json files


## To do
* random forest predictions
* better documentation
* requirements file