DIMet: Differential analysis of Isotope-labeled targeted Metabolomics data
===

# Introduction

DIMet is a bioinformatics pipeline for **differential and time-course analysis of targeted isotope-labeled data**.

DIMet supports the analysis of full metabolite abundances and isotopologue contributions, and allows to perform it either in the differential comparison mode or as a time-series analysis. As input, the DIMet accepts three types of measures: a) isotopologues’ contributions, b) fractional contributions (also known as mean enrichment), c) full metabolites’ abundances. Specific functions process each of the three types of measures separately.

_Note_: DIMet is intended for downstream analysis of tracer metabolomics data that has been corrected for the presence of natural isotopologues. 

_Formatting and normalisation helper_: scripts for formatting and normalization are provided in [Tracegroomer](https://github.com/johaGL/Tracegroomer).

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
* `tools` directory contains venv setup scripts.


## Running unit tests 

* with pytest, by running `pytest` from `DIMet`
* Place yourself in `DIMet/tests` and execute `python -m unittest` 

-----------------------------------------------------------------------------------------------

# Using DIMet 

DIMet runs in the command line environment. 

### Using DIMet with provided datasets

To **test the use of DIMet**, we provide datasets, configuration and bash scripts corresponding to the results presented in the manuscript "DIMet: An open-source tool for Differential analysis of Isotope-labeled targeted Metabolomics data" by J. Galvis *et al*. 
are available at [Zenodo (manuscript_data)](https://sandbox.zenodo.org/record/1234735):

 - Download and uncompress the file `datasets_manuscript_DIMet.zip`.    
    
    <details>
     
    <summary>Whole structure of the downloaded folder <sub><sup>(click to show/hide)</sup></sub>
    </summary>
    
    ```
    datasets_manuscript_DIMet
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
    
* Make sure you have activated your virtual environment. In the `datasets_manuscript_DIMet` folder you have two `.sh` files: `run_LDHAB-Control.sh` and `run_Cycloserine_timeseries.sh`; Make them executable:
  ```
  chmod a+x *.sh
  ```
* and finally run:
  ```
  ./run_LDHAB-Control.sh
  ./run_Cycloserine_timeseries.sh
  ```


## Available analyses

- _pca_analysis_ computes the PCA and outputs tables with principal components and explained variances
- _pca_plot_ generates classical PCA plots
- _abundance_plot_ plots with bars of total metabolite abundances
- _mean_enrichment_line_plot_ generates lineplots of mean enrichment
- _isotopologue_proportions_plot_ generates stacked bars of isotopologue proportions
- _differential_analysis_ runs differential analysis and computes the corresponding statistics
- _multi_group_comparison_ same as differential analysis before, but for > 2 groups
- _time_course_analysis_ runs differential analysis for time-course experiments in pairwise fashion for consecutive time points
- _metabologram_integration_ pathway-based integration between *labeled targeted metabolomic* and *trascriptomic* data, resulting in metabologram plots

All the analyses must be in a specific folder structure. Once the structure is ready, the **generic command** for running each analysis is:

```commandline
python -m dimet -cd config -cn GENERAL_CONFIGURATION_FILENAME
```

## Documentation

The documentation of DIMet is found in its [Wiki](https://github.com/cbib/DIMet/wiki) page, which contains all the information that the user needs to build her/his own folder structure of data and configuration files, and to successfully run DIMet. The section **organising your data for the analysis** and its respective subsections in the [Wiki](https://github.com/cbib/DIMet/wiki) provide explanations and examples to the user.


-----------------------------------------------



  
# Getting help

For any information or help running DIMet, you can get in touch with: 

* [Johanna Galvis](mailto:deisy-johanna.galvis-rodriguez[AT]u-bordeaux.fr)
* [Macha Nikolski](mailto:macha.nikolski[AT]u-bordeaux.fr)
* [Benjamin Dartigues](mailto:benjamin.dartigues[AT]u-bordeaux.fr)

# LICENSE MIT

    Copyright (c) 2023 
    
    Johanna Galvis (1,2)    deisy-johanna.galvis-rodriguez@u-bordeaux.fr
    Benjamin Dartigues (2)	benjamin.dartigues@u-bordeaux.fr
    Florian Specque (1,2)   florian.specque@u-bordeaux.fr
    Slim Karkar (1,2)       slim.karkar@u-bordeaux.fr
    Helge Hecht (3,5)       helge.hecht@recetox.muni.cz
    Bjorn Gruening (4,5)    bjoern.gruening@gmail.com
    Hayssam Soueidan (2)    massyah@gmail.com
    Macha Nikolski (1,2)    macha.nikolski@u-bordeaux.fr
    
    (1) CNRS, IBGC - University of Bordeaux,
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
