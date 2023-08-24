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

All the analyses must be in a structure that is explained in the section [Organising your data for the analysis](#organising-your-data-for-the-analysis).
Once the structure is ready, the generic command for running each analysis is:

```commandline
python -m dimet -cd config -cn GENERAL_CONFIGURATION_FILENAME
```


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

### Data files

The user has to provide his data files within the data folder. Original data files themselves
have to be placed in the `raw` subfolder along with a metadata file that contains the experimental setup corresponding to the data. 

Both quantification and metadata files must be provided as in .csv 'tab' delimited format.

_Note_:  if you provide at least one type of measure, you can still run some of the analyses, by making sure that the data you provide is suitable for the analysis that you choose.

#### a. Metadata file
The structure of the **metadata file** has to contain 6 columns named 
<code>name_to_plot</code>, <code>timepoint</code>, 
<code>timenum</code>, <code>condition</code>, 
<code>compartment</code>, <code>original_name</code>. 

<details><!--Metadata format-->

 <summary>Explanation of the metadata format
 <sub><sup>(click to show/hide)</sup></sub>
 </summary>
  
Here is the semantics of the columns:

- <code>name_to_plot</code> is the string that will appear on the figures produced by DIMet
- <code>condition</code> is the experimental condition
- <code>timepoint</code> is the sampling time as it is defined in your experimental setup
(it is an arbitary string that can contain non numerical characters)
- <code>timenum</code> is the numerical encoding of the <code>timepoint</code>
- <code>compartment</code> is the name of the cellular compartment for which the measuring
has been done (e.g. "endo", "endocellular", "cyto", etc)
- <code>original_name</code> contains the column names that are provided in the quantification files

_Example_:

| name_to_plot | condition | timepoint | timenum | compartment | original_name |
|--------------|-----------|-----------|---------|-------------|---------------|
| Cond1 T0     | cond1     | T0        | 0       | comp_name   | T0_cond_1     |
| Cond1 T24    | cond1     | T24       | 24      | comp_name   | T24_cond_1    |
| Cond2 T0     | cond2     | T0        | 0       | comp_name   | T0_cond_2     |
| Cond3 T24    | cond2     | T24       | 24      | comp_name   | T24_cond_2    |

</details>

#### b. Quantification files

Each quantification file is expected to correspond to one type of measure. Supported measure types are:

1. Isotopologue absolute values
2. Total metabolite abundances
3. Mean enrichment (also called Fractional contribution) 
4. Isotopologue proportions

<details><!--Quantifications format-->
 <summary>Expected format of the quantification files with examples
 <sub><sup>(click to show/hide)</sup></sub>
 </summary>

Each row in the quantification files contains measurements for a given metabolite. Expected columns are the following:
- <code>ID</code> contains the molecule identifiers
- All the other columns contain measures in numeric format (no letters or symbols, only numbers).

_Note 1_: quantification columns' names have to match with the column <code>original_name</code> in the **metadata** file.
_Note 2_: For the isotopologues, the </code>ID</code> must follow the convention: `metaboliteID_m+X` (for example: `AMP_m+4`, `cit_m+0`, `cit_m+1`)

The total metabolites' Abundances file:
    
| ID       | T0_cond_1  | T24_cond_1  | T0_cond_2  | T24_cond_2 |
|----------|------------|-------------|------------|------------|
| PEP      | 3364610.46 | 10250098.25 | 1124772.29 | 1035932.25 |
| citrate  | 5783654.51 | 5934305.65  | 3546334.99 | 3460334.88 |
| fumarate | 354387.74  | 360087.74   | 334287.74  | 350387.74  |
| OA       | 9435186.33 | 9435186.33  | 9435186.33 | 9435186.33 |
    
    
The Mean enrichment (also called Fractional contribution) file:


| ID       | T0_cond_1 | T24_cond_1 | T0_cond_2 | T24_cond_2 |
|----------|-----------|------------|-----------|------------|
| PEP      | 0.5603    | 0.6391     | 0.9591    | 0.9553     |
| citrate  | 0.8057    | 0.8870     | 0.7809    | 0.6918     |
| fumarate | 0.001     | 0          | 0.1508    | 0.1511     |
| OA       | 0.7030    | 0.7006     | 0.001     | 0          |
    
  
The Isotopologue absolute values file:

| ID      | T0_cond_1  | T24_cond_1 | T0_cond_2 | T24_cond_2 |
|---------|------------|------------|-----------|------------|
| PEP_m+0 | 357354.66  | 387054.66  | 0         | 0          |
| PEP_m+1 | 965435.68  | 975030.68  | 668.91    | 568.87     |
| PEP_m+2 | 1435050.95 | 7987654.66 | 136749.05 | 137709.05  |
| PEP_m+3 | 606769.17  | 900358.25  | 987354.33 | 897654.33  |


The Isotopologue proportions file :


| ID      | T0_cond_1 | T24_cond_1 | T0_cond_2 | T24_cond_2 |
|---------|-----------|------------|-----------|------------|
| PEP_m+0 | 0.106     | 0.038      | 0.000     | 0.000      |
| PEP_m+1 | 0.287     | 0.095      | 0.001     | 0.001      |
| PEP_m+2 | 0.427     | 0.779      | 0.122     | 0.133      |
| PEP_m+3 | 0.180     | 0.088      | 0.878     | 0.867      |


</details>


#### c. Data files for the omics integration (optional)

DIMet offers the possitibilty of pathway-based integration of the metabolome and the transcriptome though **metabolograms**. 

<details>
  <summary>Data files required for omics integration<sup><sub> (click to show/hide)</sub></sup>
  </summary>
  Two data types are required:

   1. Metabolite quantification files in the <code>raw</code> directory.
   2. Results, provided by the user, of the differential analysis of the transcriptome data placed in the folder <code>integration_files</code>

 <code>integration_files</code> folder contains files with differentially expressed genes provided by the user as well as the pathways files.
 
 Thus the expected project data structure becomes:
  
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
            ├── integration_files  # ---> NEW FOLDER
            │   ├── # ---> .csv files required for omics integration
            └── raw
                ├── # ---> raw .csv files

  ```
   
  * Files for differentially expressed genes (DEGs)

Files for differentially expressed genes (DEGs) must be provided in the tab delimited .csv format. For each file:
1. The rows represent the genes (except the first one, which is the header having the names of the columns)
2. The columns provide the information to be integrated, two columns are compulsory:
   1. the gene names, given as strings
   2. the Fold-Changes (or the log2 Fold-Changes) in numeric format (no letters or symbols, only numbers)

Formatting example of differentially expressed genes files:
    
| log2FoldChange    | gene_symbol  |
|-------------------|--------------|
| -16.1660338229612 | GPI          |
| 3.32192809488736  | HK1          |
| 2.32192809488736  | RPIA         |
| 0.807354922057604 | PFKL         |

   
  * The *metabolites per pathway* and *genes or transcripts per pathway* files

 These files contain the user-provided metabolites and genes for each pathway. It is allowed for a metabolite or gene to appear in several pathways.
 Identifiers must match with those appearing in the quantification files in the <code>raw</code> folder. Gene names must match with those appearing in the DEGs file 

_Example_ for metabolites per pathway:

| GLYCOLYSIS | PENTOSE_PHOSPHATE | ... |
|------------|-------------------|-----|
| Glucose_6P | Ribose_5P         | ... |
| Pyruvate   | Xylulose_5P       | ... |
| PEP        | Glucose_6P        | ... |
| ...        | ...               | ... |


_Example_ for genes per pathway:

| GLYCOLYSIS | PENTOSE_PHOSPHATE | ... |
|------------|-------------------|-----|
| GPI        | RPIA              | ... |
| HK1        | PGD               | ... |
| PKFL       | RBKS              | ... |
| ...        | ...               | ... |

All these files must be provided in the tab delimited .csv format.

</details>


### Configuration files

As it is shown in the folder structure above, there are three types of configuration files, 
all of them in `.yaml` format: 

1. `dataset configuration` files:  provide the quantifications and metadata file names that are themselves present in the `data` directory.

2. `analysis configuration` files: indicates which analysis (e.g. "differential analysis") is to be run for which data file and with which parameters.

3. `general configuration` files: configures the analysis to be run by both pointing to the `analysis configuration` and providing general information such as e.g. output folder names etc.


<details><!-- section organize the configuration files -->

<summary><b> How to organize the configuration files</b> <sub><sup>(click to show/hide)</sup></sub></summary>

  The user must place the configuration files in the `config` folder, as follows: 
* the *dataset configuration* files, inside the `config/analysis/dataset/` subfolder
* the *analysis configuration* files,  inside the `config/analysis/` subfolder
* the *general configuration* files, directly inside the `config` folder
*   <details>
    <summary>Example of a structure of the  <code>config</code>  folder  
     <sub><sup>(click to show/hide)</sup></sub>
    </summary>
    
      ```
      config
      ├── analysis
      │   ├── dataset
      │   │   ├── dataname1_data.yaml
      │   │   └── dataname1_data_integrate.yaml  # <-- if suitable
      │   ├── differential_analysis_pairwise_dataname1.yaml
      │   ├── enrichment_lineplot_dataname1.yaml
      │   ├── isotopologues_plot_dataname1.yaml
      │   ├── metabologram_abundance_dataname1.yaml
      │   ├── pca_plot_dataname1.yaml
      │   └── timecourse_analysis_dataname1.yaml
      ├── general_config_abundance_plot_dataname1.yaml
      ├── general_config_differential_analysis_dataname1.yaml
      ├── general_config_enrichment_lineplot_dataname1.yaml
      ├── general_config_isotopologues_plot_dataname1.yaml
      ├── general_config_metabologram_abundance_dataname1.yaml
      ├── general_config_pca_plot_dataname1.yaml
      └── general_config_timecourse_analysis_dataname1.yaml
    
      ```
     </details> 
   


The next sections explain each type of these files. Furthermore, when analyzing your own data, you can
use these configuration examples as templates and fill them coherently.

The `# <-` comment inside the examples indicates the parts that the user must fill. 

 <details><!--section dataset config -->
 <summary><b>The dataset configuration </b></summary>  
  
   For each given dataset, there must exist a corresponding *dataset configuration* file, which is located inside the `config/analysis/dataset/` folder.
   This file must specify:
  - the name of the metadata file
  - the names of the quantification files
  - the "order" of the conditions, where the *control* is the first to be listed 
(note that the name assigned to the control condition can be different than the one
we used in our data).

  - Note that the parameter `_target_` must **not** be changed.
  - If you do not have all the four types of quantification that are shown below, just delete
    the entire line making reference to the type of quantification that is absent.
  
  ```
  _target_: dimet.data.DatasetConfig

  label:  # <- name of your dataset, fill after the colon
  name:  # <- short description of your dataset, fill after the colon
  subfolder:  DATANAME1_data  # <- subfolder name in the data folder, change after the colon
    
  # ALWAYS WITHIN THE RAW SUBFOLDER
  metadata:    # <- file name, fill after the colon
  abundances:    # <- file name, fill after the colon
  mean_enrichment:    # <- file name, fill after the colon
  isotopologue_proportions:    # <- file name, fill after the colon
  isotopologues:    # <- file name, fill after the colon
    
  conditions :
   - cond1 # <- first must be control, fill
   - cond2 # <- the rest of the conditions vertically listed
  ```
  This dataset configuration will be used by all the types of analysis we offer except the omics integration, see the section *The integration dataset configuration* for that purpose.



   <details><!-- section integration dataset config-->
   <summary>Special case of dataset configuration: the integration dataset configuration (for <b>omics integration</b>)</summary>
  
   When performing the **omics integration** , 
   there must also exist a *dataset configuration* file specific for the integration, 
   that will point to `integration_files` and `raw` subfolders.
       
   Note that the parameter `_target_` must **not** be changed.
       
    ```
    _target_: dimet.data.DataIntegrationConfig  

    label: integrate_DATANAME1   # <- change after the colon
    name:     # <- short description of your dataset
    subfolder: DATANAME1_data  # <- subfolder name in the data folder, change after the colon
    
    conditions :
     - cond1 # <- first must be control
     - cond2 # <- the rest of the conditions vertically listed
    
    # ALWAYS WITHIN THE RAW SUBFOLDER
    metadata:    # <- file name, fill after the colon
    abundances:    # <- file name, fill after the colon
    mean_enrichment:    # <- file name, fill after the colon
    
    # WITHIN THE INTEGRATION_FILES SUBFOLDER
    transcripts:  
      - myDEG_1   # <- file name
    
    pathways: 
      metabolites:   # <- file name, change fill the colon
      transcripts:    # <- file name, change fill the colon

    ```

   </details><!--closes  section integration dataset config-->

 
  </details><!--closes section dataset config-->



  <details><!--section analysis config-->
  <summary><b>The analysis configuration </b></summary>   
   
   The *analysis configuration* is located inside the `config/analysis/` folder. 
   For each analysis to be performed, one analysis configuration file must exist. 
   It indicates which is the type of analysis we want to run, 
   on which dataset this analysis will be applied, and the parameters that are 
   specific to that analysis.

   
     ```
     label: differential-analysis-DATANAME1   # <- change after the colon

     defaults:
       - dataset:      # <- name of the dataset cofiguration file
       - method: differential_analysis  # <- see 'Available analyses'
     
     comparisons :
       - [[cond2, T24], [cond1, T24]]   # <-      
     
     statistical_test:
       abundances:    # <- see statistic test options, fill after the colon
       mean_enrichment: prm-scipy    # <- see options, fill after the colon
       isotopologues: disfit    # <- see options, fill after the colon
       isotopologue_proportions: ranksum    # <- see options, fill after the colon

     ```
   
   Currently, the **statistic test options** are:

   `MW` (Mann-Whitney test), `KW` (Kruskall-Wallis test), 
   `ranksum` (Wilcoxon's rank sum test), `Wcox` (Wilcoxon's signed rank test),
    `BrMu` (Brunner-Munzel test) , `prm-scipy` (permutation testing),
    `Tt` (t-test) and
   `disfit` (Fitting of the distribution of the z-scores).
   
   Note that the parameters in the *analysis configuration* must be coherent with our dataset.
 
        
   <details><!-- section analysis config for integration-->
   <summary> The analysis configuration for the <b>omics integration</b> </summary> 

   The integration performed by DIMet displays the 'differential omics' in the form of Metabologram(s). 
   The analysis configuration file:
    
       ```
       label: metabologram-using-abundance-DATANAME1  # <- change after the colon
   
       defaults:
         - dataset:   # <- integration dataset config
         - method: metabologram_integration
       
       comparisons :
         - [[cond2, T24], [cond1, T24]]   # <-  see note (**) 
       
       # running for total abundances
       statistical_test:
         abundances:    # <-  can be abundances OR mean_enrichment
       
       columns_metabolites:
         ID : metabolite
         values :    # <- log2FC or FC, fill after the colon
         
       columns_transcripts:
         ID:   # <- the gene symbols column name, fill after the colon
         values:   # <- the numeric column name,  fill after the colon
         
       compartment:
         en    
        ```
   (`**`) each metabolomics comparison must match with a DEG file.
   
   If you need to run the integration using the differential mean enrichment, 
   a separate analysis configuration file must be created.

   </details><!--closes  section analysis config for the integration-->  
   
  </details><!--closes  section analysis config-->


  <details><!--section general config-->
  <summary><b> The general configuration </b></summary>    
   
   The *general configuration* file is located directly inside the `config` folder. It must declare the name of its corrresponding *analysis configuration* file. Place yourself in `defaults:` then `- analysis:`, see below:
     
     ```
     hydra:
       job:
         chdir: true
       run:
         dir: outputs/${now:%Y-%m-%d}/${now:%H-%M-%S}/${analysis.dataset.label}-${analysis.method.label}

     defaults:
       - analysis:   # <--  the 'analysis configuration' name, fill after the colon

     figure_path: figures
     table_path: tables    
     ```

   
   This *general configuration* file also provides the parameters for generating the output folder names; keep these 
   output parameters unmodified, so DIMet will generate the full output folders and files names
   automatically.
     
   </details><!--closes section general config--> 
   
   Note that across all the types of configuration files, any referenced file name must be written without the extension.
     

</details><!--closes section  organize the configuration files--> 


 Examples of configuration files, with their respective datasets are provided 
in  [Zenodo (manuscript_data)](https://sandbox.zenodo.org/record/1234735). 
Also complementary minimal datasets with their configuration files are provided in
[Zenodo (minimal_examples)](https://sandbox.zenodo.org/record/1234736). 


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
