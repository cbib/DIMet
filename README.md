DIMet: Differential analysis Isotope-labeled targeted Metabolomics data
===

# Introduction

DIMet is a bioinformatics pipeline for differential and time-course analysis of targeted isotope-labelled data.

DIMet supports the analysis of full metabolite abundances and isotopologue contributions, and allows to perform it either in the differential comparison mode or as a time-series analysis. As input, the DIMet accepts three types of measures: a) isotopologues’ contributions, b) fractional contributions (also known as mean enrichment), c) full metabolites’ abundances. Specific functions process each of the three types of measures separately.

Note: DIMet is intended for downstream analysis of tracer metabolomics data that has been corrected for the presence of natural isotopologues. 

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

-----------------------------------------------------------------------------------------------------
# Organising your data for the analysis

## Input data folder structure

DIMet relies on a specific folder organisation for input data and the accompanying configuration files.
Data should be structured following the example `MYPROJECT` below. 

```
MYPROJECT
├── config
│   ├── analysis
│   │   ├── dataset
│   │   │   └── # --->'dataset configuration' yml files
│   │   ├── # --->'analysis configuration' yml files
│   ├── # ---> 'general configuration' yml files
└── data
    └── DATANAME1_data
        └── raw
            ├── # ---> raw .csv files

```

`DATANAME1` represents an experiment, you should replace it by a meaningful name.
More than one subfolder in the `data` folder can be present.

### Configuration files

As it is shown in the folder structure above, there are three types of configuration files, all of them in yaml format: 

1. `dataset configuration` files:  provide the quantifications and metadata file names that are themselves present in the `data` directory.

2. `analysis configuration` files: indicates which analysis (e.g. "differential analysis") is to be run for which data file and with which parameters.

4. `general configuration` files: configures the analysis to be run by both pointing to the `analysis configuration` and providing general information such as e.g. output folder names etc.

Examples of these configuration files, with their respective minimal datasets are provided on [Zenodo](https://sandbox.zenodo.org/record/1224019)

## Provided datasets

Datasets, configuration and bash commands corresponding to the results presented in the manuscript "DIMet: An open-source tool for Differential analysis of Isotope-labeled targeted Metabolomics data" by J. Galvis et al. are available at [Zenodo](https://sandbox.zenodo.org/record/1224340). 

Download and uncompress the file `datasets_manuscript_DIMet.zip`.
<details>
<summary>Following the folder structure above, the `data` folder contains</summary>
 ```
├── config
│   ├── analysis
│   │   ├── abundance_plot_Cycloserine.yaml
│   │   ├── abundance_plot_LDHAB-Control.yaml
│   │   ├── dataset
│   │   │   ├── Cycloserine_data.yaml
│   │   │   ├── LDHAB-Control_data_integrate.yaml
│   │   │   └── LDHAB-Control_data.yaml
│   │   ├── differential_analysis_pairwise_LDHAB-Control.yaml
│   │   ├── enrichment_lineplot_Cycloserine.yaml
│   │   ├── isotopologues_plot_Cycloserine.yaml
│   │   ├── isotopologues_plot_LDHAB-Control.yaml
│   │   ├── metabologram_abundance_LDHAB-Control.yaml
│   │   ├── metabologram_enrichment_LDHAB-Control.yaml
│   │   ├── pca_plot_LDHAB-Control.yaml
│   │   ├── pca_tables_Cycloserine.yaml
│   │   └── timecourse_analysis_Cycloserine.yaml
│   ├── general_config_abundance_plot_Cycloserine.yaml
│   ├── general_config_abundance_plot_LDHAB-Control.yaml
│   ├── general_config_differential_analysis_LDHAB-Control.yaml
│   ├── general_config_enrichment_lineplot_Cycloserine.yaml
│   ├── general_config_isotopologues_plot_Cycloserine.yaml
│   ├── general_config_isotopologues_plot_LDHAB-Control.yaml
│   ├── general_config_metabologram_abundance_LDHAB-Control.yaml
│   ├── general_config_metabologram_enrichment_LDHAB-Control.yaml
│   ├── general_config_pca_plot_LDHAB-Control.yaml
│   ├── general_config_pca_tables_Cycloserine.yaml
│   └── general_config_timecourse_analysis_Cycloserine.yaml
├── data
│   ├── Cycloserine_data
│   │   └── raw
│   │       ├── CorrectedIsotopologues.csv
│   │       ├── FracContribution_C.csv
│   │       ├── metadata_cycloser.csv
│   │       └── rawAbundances.csv
│   └── LDHAB-Control_data
│       ├── integration_files
│       │   ├── DEG_Control_LDHAB.csv
│       │   ├── pathways_kegg_metabolites.csv
│       │   ├── pathways_kegg_transcripts.csv
│       │   └── readme.txt
│       └── raw
│           ├── AbundanceCorrected.csv
│           ├── IsotopologuesAbs.csv
│           ├── IsotopologuesProp.csv
│           ├── MeanEnrichment13C.csv
│           └── metadata_endo_ldh.csv
├── run_Cycloserine_timeseries.sh
└── run_LDHAB-Control.sh
 ```
</details>

-----------------------------------------------------------------------------------------------------

# Using DIMet 

DIMet runs in the command line environment. 

## Running analyses on the provided datasets

Make sure you have activated your virtual environment. In the `datasets_manuscript_DIMet` folder you have two `.sh` files: `run_LDHAB-Control.sh` and `run_Cycloserine_timeseries.sh`; Make them executable:
```
chmod a+x *.sh
```
and finally run:
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


# Getting help

For any information or help running DIMet, you can get in touch with: 

* [Johanna Galvis](mailto:deisy-johanna.galvis-rodriguez[AT]u-bordeaux.fr)
* [Macha Nikolski](mailto:macha.nikolski[AT]u-bordeaux.fr)
* [Benjamin Dartigues](mailto:benjamin.dartigues[AT]u-bordeaux.fr)

# LICENSE MIT

    Copyright (c) 2023 
    
    Johanna Galvis (1,2)    deisy-johanna.galvis-rodriguez@u-bordeaux.fr
    Benjamin Dartigues (2)	 benjamin.dartigues@u-bordeaux.fr
    Florian Specque (1,2)	  florian.specque@u-bordeaux.fr
    Helge Hecht (3,5)       helge.hecht@recetox.muni.cz
    Bjorn Gruening (4,5)    bjoern.gruening@gmail.com
    Hayssam Soueidan (2)    massyah@gmail.com
    Macha Nikolski (1,2)    macha.nikolski@u-bordeaux.fr
    
    (2) CNRS, IBGC - University of Bordeaux,
    1, rue Camille Saint-Saens, Bordeaux, France

    (2) CBiB - University of Bordeaux,
    146, rue Leo Saignat, Bordeaux, France

    (3) Spectrometric Data Processing and Analysis,
    Masaryk University, Brno, Czech Republic
    
    (4) University of Freiburg, 
    Freiburg, Germany
    
    (5) Galaxy Europe
    
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
