# EEx-NTA
A repository for storing mass-spectrometry non-target analysis (NTA) code linked to the Environmental Exposure Hub. 

## ordination_mat.py 

This script is adapted from ```instant_search_processing.R``` by Nikiforos Alygizakis. It downloads all of the detection files in a collections and collates them into a N by M matrix. Where N is the number of samples in the collection and M is the number of detected compounds. Each sample compound pair has its peak area as a value in the matrix. This is ultimately written to a .csv file. 

The script should be run from the command line as follows ```python ordination_mat.py config.yaml```. The ```config.yaml``` file contains the following parameters: 
* ```DOWNLOAD: True``` - Whether to download detection json files
* ```METHOD: 'max'``` - How to deal with compounds that match to multiple peaks in a sample. ```max``` takes the maximum value, ```sum```, sums the values. 
* ```THREADS: 100``` - How many parallel threads to use when downloading data
* ```FILEPATH: "."``` - Filepath to where output data should be written. Output data is written as ```{COLLECTION_ID}_ordination.csv```. 
* ```DOWNLOAD_DIR: "downloads"``` - Where to store/read downloaded json files
* ```COLLECTION_ID: 2380``` - Which collection ID to search
* ```NS_MIN: 1``` - Compounds are searched from MIN:MAX in the Normal SusDat datbase
* ```NS_MAX: 120000```- Compounds are searched from MIN:MAX in the Normal SusDat datbase
