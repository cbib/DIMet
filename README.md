DIMet: Differential Isotope-labeled targeted Metabolomics
===

# Introduction

DIMet is a bioinformatics pipeline for differential analysis of targeted isotope-labelled data.

DIMet supports the analysis of full metabolite abundances and isotopologue contributions, and allows to perform it either in the differential comparison mode or as a time-series analysis. As input, the DIMet accepts three types of measures: a) isotopologues’ contributions, b) fractional contributions (also known as mean enrichment), c) full metabolites’ abundances. Specific functions process each of the three types of measures separately.

Note: DIMet is intended for downstream analysis of tracer metabolomics data that has been corrected for the presence of natural isotopologues. Make sure you that the metabolomics platform provides you the output of the correction procedure before using this pipeline.

# Provided datasets

Datasets, configuration and bash commands corresponding to the results presented in the manuscript "DIMet: An open-source tool for Differential analysis of targeted Isotope-labeled Metabolomics data
" are available on [Zenodo](https://sandbox.zenodo.org/record/1224020).
Download and uncompress the file `datasets_manuscript_DIMet`.

You can use the files in `data/Cycloserine_data/raw`, `data/LDHAB-Control_data/raw` and
`data/LDHAB-Control_data/integration_files` and perform the analyses through the Galaxy version of DIMet, which contains the general instructions to do so.

Otherwise, after downloading and uncompressing the datasets folder, place it inside any directory NOT being the `DIMet` directory. Perform the analyses in a command line environment, by placing yourself inside `datasets_manuscript_DIMet` and running the `.sh` executable files to reproduce the analyses. 
The rest of the present document is therefore addressed to users 
desiring to work in command line environment. 


-----------------------------------------------------------------------------------------------------


# Installing DIMet 

## System Requirements
DIMet installation requires an Unix environment with [python 3.9](http://www.python.org/). 
It was tested under Linux and MacOS environments.


## Installation 

The full installation process should take less than 15 minutes on a standard computer.

Via pip command:
`pip install dimet`

Or if you are a developer working in a local cloned version, you can install:
`pip install -e .`


## Code organization

* `src/processing` directory contains the implemented high-level analysis scripts that produced the tables for the DIMet paper
* `src/visualization` : directory contains the implemented high-level scripts that produced the figures in the DIMet paper
* `src/data` directory contains the python classes for data initialization  
* `src/method` directory contains the python classes for configuration handling
* `src/tests` directory contains unit tests
* `tools` directory contains venv setup scripts. Alternatively, the yml file is also provided in this directory, to perform the installation via conda.

## Running unit tests 

* with pytest, by running `pytest` from `DIMet`
* Place yourself in `DIMet/tests` and execute `python -m unittest` 


# Using DIMet 

DIMet has a Galaxy version. It also runs in a command line environment. The runtime is dependent on the hardware, the number of metabolites in the dataset, and the selected options.


## Running analyses on the manuscript data

Make sure you have activated your virtual environment, or DIMet conda environment.
After downloading the data from Zenodo, you have a folder named `datasets_manuscript_DIMet`.
You have inside, two .sh files. Each .sh file will run sequentially, different analyses for its respective dataset.
Make them executable
```
chmod a+x *.sh
```
and finally:
```
./run_LDHAB-Control.sh
./run_Cycloserine_timeseries.sh
```


## Running analyses on your data

To run any analysis on your own data, first your files have to be preprocessed.
Preprocessing includes, in this order: 

-  *(a)* The correction by the abundance of naturally occurring isotopologues by external tools such IsoCor, ElMaven, AccuCor, etc). 

-  *(b)* The normalization by internal standard and/or the normalization by amount of material and/or the formatting of the quantification tables.

For the *(b)* type of preprocessing we offer our accompanying tool [Tracegroomer](https://github.com/johaGL/Tracegroomer). 
Only use DIMet when your data has already underwent the complete preprocessing, as it will not be carried out by DIMet.

The following sections will help you to understand the general expected organization of your 
input data and configuration files. Further minimal examples can also be downloaded from [Zenodo](https://sandbox.zenodo.org/record/1224019)


## General organization of the input data and configuration (when in command line environment)

Input data and configuration files that were used to generate the figures in the DIMet paper are provided as examples of how the DIMet pipeline can be used. 


## Input format

DIMet, in its command line environment,  takes a folder of hierarchically organized folders containing the data and the configuration files.
Below there is a minimal structure example, `MYPROJECT` must be replaced by your project name, for example "Brain" or your field of study. The word `DATANAME1` represents an experiment, and must be replaced by a meaningful noun.
More than one subfolder in the `data` folder can be added, but make sure to elaborate the respective -and unambiguous- configuration files.

```
MYPROJECT
├── config
│   ├── analysis
│   │   ├── dataset
│   │   │   └── DATANAME1_data.yaml
│   │   ├── differential_analysis_DATANAME1.yaml
│   │   ├── ... # ---> other 'analysis configuration' yml files
│   │   └── pca_analysis_DATANAME1.yaml
│   ├── general_config_differential_analysis_DATANAME1.yaml
│   ├──  ... # ---> other 'general configuration' yml files
│   └── general_config_pca_analysis_DATANAME1.yaml
└── data
    └── DATANAME1_data
        └── raw
            ├── AbundanceCorrected.csv
            ├── Isotopologues_proportions.csv
            ├── MeanEnrichment.csv
            └── DATANAME1_metadata.csv
```

### Data

Using the aliases shown in the example structure above, the `data` folder contains:

* `DATANAME1_data` that is a folder named accordingly with the experiment  (for example `Cycloserine_data` or `LDHAB-Control_data` or `example1_data`). In turn this `DATANAME1_data`, contains:                           
    -  One `raw` subfolder where your quantifications and metadata are located.
  Note that for practical reasons we used the name "raw", but as explained above, this data has been already preprocessed before using DIMet.
    -  After running any analysis, one `processed` subfolder is generated, which contains the tables split by cellular compartment, cleaned from rows containing only NaN, and also cleaned from rows containing only 0. 


### Configuration files

There are three types of configuration files, all of them in json format: 

1. `dataset configuration` files

2. `analysis configuration` files

3. `general configuration` files

 
A given `dataset configuration` file refers to the corresponding quantifications and metadata files that are present in the `data` folder, and is located in the `config/analysis/dataset/` folder.

The `analysis configuration` file  should be located in the `config/analysis/` folder
Thus, an `analysis configuration` file indicates the name of the `dataset configuration` file, which in turn points to the input data in the `data` folder. 

The `general configuration` file points to the `analysis configuration` file, and is located in `config`.

Further examples of these configuration files, with their respective minimal datasets are provided on [Zenodo](https://sandbox.zenodo.org/record/1224019)


# Getting help

For any information or help running DIMet, you can get in touch with: 

* [Macha Nikolski](mailto:macha.nikolski[AT]u-bordeaux.fr)
* [Benjamin Dartigues](mailto:benjamin.dartigues[AT]u-bordeaux.fr)
* [Johanna Galvis](mailto:deisy-johanna.galvis-rodriguez[AT]u-bordeaux.fr)

# LICENSE MIT

    Copyright (c) 2023 
    
    Johanna Galvis (1,3)    (deisy-johanna.galvis-rodriguez@u-bordeaux.fr)
    Joris Guyon (2)         (joris.guyon@u-bordeaux.fr)
    Benjamin Dartigues (1)	(benjamin.dartigues@u-bordeaux.fr)	 
    Florian Specque (3)	  	(florian.specque@u-bordeaux.fr)
    Slim Karkar (1,3)       (slim.karkar@u-bordeaux.fr)
    Helge Hecht             (helge.hecht@recetox.muni.cz)
    Bjorn Gruening          (bjoern.gruening@gmail.com)
    Thomas Daubon (3)       (thomas.daubon@u-bordeaux.fr)  
    Macha Nikolski (1,3)    (macha.nikolski@u-bordeaux.fr)

    (1) CBiB - University of Bordeaux,
    146, rue Leo Saignat, 33076 Bordeaux, France
    
    (2) University of Bordeaux, INSERM, BPH U1219, Bordeaux, France

    (3) CNRS, IBGC - University of Bordeaux,
    1, rue Camille Saint-Saens, 33077 Bordeaux, France

    (4) Galaxy Europe
    

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
