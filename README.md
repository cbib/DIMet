DIMet: Differential analysis of Isotope-labeled targeted Metabolomics data
===

# Introduction

DIMet is a bioinformatics pipeline for differential and time-course analysis of targeted isotope-labeled data.

DIMet supports the analysis of full metabolite abundances and isotopologue contributions, and allows to perform it either in the differential comparison mode or as a time-series analysis. As input, the DIMet accepts three types of measures: a) isotopologues’ contributions, b) fractional contributions (also known as mean enrichment), c) full metabolites’ abundances. Specific functions process each of the three types of measures separately.

_Note_: DIMet is intended for downstream analysis of tracer metabolomics data that has been corrected for the presence of natural isotopologues. 

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

### Data files

The user has to provide his data files within the data folder. Original data files themselves
have to be placed in the `raw` subfolder along with a metadata file that contains 
the experimental setup corresponding to the data. 

#### - Metadata file
The structure of the **metadata file** has to contain 6 columns named 
<code>name_to_plot</code>, <code>timepoint</code>, 
<code>timenum</code>, <code>condition</code>, 
<code>compartment</code>, <code>original_name</code>. 
Here is the semantics of the columns:

- <code>name_to_plot</code> is the string that will appear on the figures produced by DIMet
- <code>condition</code> is the experimental condition
- <code>timepoint</code> is the sampling time as it is defined in your experimental setup
(it is an arbitary string that can contain non numerical characters)
- <code>timenum</code> is the numerical encoding of the <code>timepoint</code>
- <code>compartment</code> is the name of the cellular compartment for which the measuring
has been done (e.g. "endo", "endocellular", "cyto", etc)
- <code>original_name</code> contains the column names that are provided in the quantification files

Example:

| name_to_plot | condition | timepoint | timenum | compartment | original_name |
|--------------|-----------|-----------|---------|-------------|---------------|
| Cond1 T0     | cond1     | T0        | 0       | comp_name   | T0_cond_1     |
| Cond1 T24    | cond1     | T24       | 24      | comp_name   | T24_cond_1    |
| Cond2 T0     | cond2     | T0        | 0       | comp_name   | T0_cond_2     |
| Cond3 T24    | cond2     | T24       | 24      | comp_name   | T24_cond_2    |



#### - Quantification files

The original data that the user provides is a set of the quantification files. 
Each quantification file is expected to correspond to one type of measure:

1. Isotopologue absolute values
2. total metabolites' Abundances
3. Mean enrichment (also called Fractional contribution) 
4. Isotopologue proportions

The user must organize his quantification files, in such a way that
the rows represent the molecules, whereas the columns (except the first one) represent the samples:
        
- The first column, named `ID`, contains the isotopologues or the metabolites identifiers.
- The second column and beyond, whose names must match with the column `original_name` in the **metadata**,
  contain the measures in numeric format (no letters or symbols allowed in the cells, only numbers).
- The isotopologues' `ID` must follow the convention: `MetaboliteID_m+X` (for example: `AMP_m+4`, `cit_m+0`, `cit_m+1`)


<details><!--format quantification-->

 <summary>Examples of the expected format of the quantification files
 <sub><sup>(click to show/hide)</sup></sub>
 </summary>
  
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


Both quantification and metadata files must be given as in .csv 'tab' delimited format.
Note that if you provide at least one type of measure, you can still 
run some of the analyses, by making sure that the data you provide is suitable for the analysis that you choose.


#### -  Data files for the omics integration (optional)

If the user has acquired both the metabolome and the transcriptome under the same
experimental conditions, she/he can perform the integration of both omics. DIMet performs
such integration at the metabolic pathway level, taking as input the differential metrics
and the pathways files, to generate the **Metabologram(s)**. 

<details>
  <summary>Details about the required data files for the omics integration<sup><sub> (click to show/hide)</sub></sup>
  </summary>
  The differential metrics for the omics integration are of two types:

   1. the results of the differential analysis of the metabolomics data, performed by DIMet automatically when the integration job is launched
   2. the results, provided by the user, of the differential analysis of the transcriptome data. External tools such as DESeq2, EdgeR, Limma, exist for this purpose. 
    
  The <code>raw</code> folder previously seen provides the necessary for (1). 
  For (2), another new folder  named <code>integration_files</code>
  must be created at the same level:
  
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

  This new folder will contain the
  differential transcriptomics or DEGs file(s), as well as the pathways files:
   
  * The *DEGs* (Differentially Expressed Genes) file must be given in tab delimited .csv format.
   The user can provide >1 DEGs files if available, with clear informative
   filenames to specify a coherent integration.

    
| ensembl          | name                           | FC        | log2FoldChange    | padj     | gene_symbol |
|------------------|--------------------------------|-----------|-------------------|----------|-------------|
| ENSG00000105220  | glucose-6-phosphate isomerase  | 0.0000136 | -16.1660338229612 | 1.00E-10 | GPI         |
| ENSG00000156515  | hexokinase 1                   | 10        | 3.32192809488736  | 1.00E-03 | HK1         |
| ENSG00000153574  | ribose 5-phosphate isomerase A | 5         | 2.32192809488736  | 0.001    | RPIA        |
| ENSG00000141959  | phosphofructokinase, liver     | 1.75      | 0.807354922057604 | 0.05     | PFKL        |

    
   
  * The *metabolites per pathway* file: the column names must be the names of the pathways. The values below each column name correspond to the
   metabolites present in each pathway. It is allowed that a same metabolite appears in several pathways.
   The metabolites' names or identifiers must match with those appearing in the metabolite total
   abundances -or the mean enrichment-  file.  Only one file of this type is accepted.

| GLYCOLYSIS | PENTOSE_PHOSPHATE | ... |
|------------|-------------------|-----|
| Glucose_6P | Ribose_5P         | ... |
| Pyruvate   | Xylulose_5P       | ... |
| PEP        | Glucose_6P        | ... |
| ...        | ...               | ... |


  * The *genes or transcripts per pathway* file: the column names must be the names of the pathways. The values below each column correspond to the gene
   symbols present in each pathway. It is allowed that a same gene symbol appears in several pathways. The
   gene symbols must match with those appearing in the DEGs file.    Only one file of this type is accepted

    
    
| GLYCOLYSIS | PENTOSE_PHOSPHATE | ... |
|------------|-------------------|-----|
| GPI        | RPIA              | ... |
| HK1        | PGD               | ... |
| PKFL       | RBKS              | ... |
| ...        | ...               | ... |

  
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
   <summary>The integration dataset configuration (for <b>omics integration</b>)</summary>
  
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
       - method: differential_analysis
     
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
         values:   # <- the numeric column name,  fillafter the colon
         
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
     
</details><!--closes  2. section -->



 Examples of configuration files, with their respective datasets are provided 
in  [Zenodo (manuscript_data)](https://sandbox.zenodo.org/record/todo:updaate). 
Also complementary minimal datasets with their configuration files are provided in
[Zenodo (minimal_examples)](https://sandbox.zenodo.org/record/todo:update). **todo:update**  

## Provided datasets

Datasets, configuration and bash scripts corresponding to the results presented in the manuscript "DIMet: An open-source tool for Differential analysis of Isotope-labeled targeted Metabolomics data" by J. Galvis *et al*. 
are available at [Zenodo (manuscript_data)](https://sandbox.zenodo.org/record/todo:update). 

Download and uncompress the file `datasets_manuscript_DIMet.zip`.

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

------------------------------------------------------------

<details><!--1.section-->

<summary><b>1. The <code>data</code> folder </b> <sub><sup>(click to show/hide)</sup></sub></summary>

  <details>
  <summary>The structure of the <code>data</code> folder is: <sub><sup>(click to show/hide)</sup></sub> </summary>

```
data
├── Cycloserine_data
│   └── raw
│       ├── CorrectedIsotopologues.csv
│       ├── FracContribution_C.csv
│       ├── metadata_cycloser.csv
│       └── rawAbundances.csv
└── LDHAB-Control_data
    ├── integration_files
    │   ├── DEG_Control_LDHAB.csv
    │   ├── pathways_kegg_metabolites.csv
    │   ├── pathways_kegg_transcripts.csv
    │   └── readme.txt
    └── raw
        ├── AbundanceCorrected.csv
        ├── IsotopologuesAbs.csv
        ├── IsotopologuesProp.csv
        ├── MeanEnrichment13C.csv
        └── metadata_endo_ldh.csv
```
   </details>

By zooming into the content of any of the **`raw`** subfolders we can understand the file formatting that is required for using DIMet, both for **quantification** file(s) and the **metadata** file:

   <details><!--1.1 section-->
   <summary><b>1.1. Quantification files</b></summary>
   
   Each quantification .csv file must have a first column named `ID` (containing the isotopologues or the metabolites identifiers), being the names of the rest of the columns the original names of the samples. 
   
   Importantly, the isotopologue identifier must follow the convention: `MetaboliteID_m+X` (for example: `AMP_m+4`, `cit_m+0`, `cit_m+1`)
   
   <details><!--the section showing the first lines of the quantification tables-->
    
   <summary>As an illustration, see here the first lines of the quantification files <sub><sup>(click to show/hide)</sup></sub></summary>

   The first lines and columns of the Isotopologue absolute values file (`IsotopologuesAbs.csv`) which is inside the `raw` subfolder of `LDHAB-Control_data`:
      
     ```
     TAB
     LE
     ```
     
 The first lines and columns of the total metabolite abundances (`AbundanceCorrected.csv`) in the same subfolder:
      
     ```
     TAB
     LE
     ```

  The first lines and columns of the MeanEnrichment13C (`MeanEnrichment13C.csv`) in the same subfolder:
      
     ```
     TAB
     LE
     ```

  The first lines and columns of the Isotopologue proportions file (`IsotopologuesAbs.csv`) in the same subfolder:
      
     ```
     TAB
     LE
     ```
     
  </details><!--closes the section showing the first lines of the quantification tables-->

  
  </details><!--closes the 1.1 section-->


  <details><!--1.2.section-->
   
  <summary><b>1.2. The metadata file</b></summary>

  For describing the samples that are present in the quantification file(s), you must provide a unique metadata file with the columns:
   - dd
   - h
   - i
   - c
   - e
   - f
   - h
   
   
   <details>
   <summary>
    As an illustration, these are the first lines of the <code>metadata_endo_ldh.csv</code>file:
   </summary> 
   
   ```
   TABLEEEE
   METADATAAAAA
   ```
   </details>
  
  
  </details><!--closes 1.2. section-->  

By zooming into the content of the `integration_files` subfolder in the `LDHAB-Control_data` folder, we see the requirements for performing the pathway-based **omics integration** of the labeled targeted metabolomics data and transcriptomics data (note that if you aquired the metabolomes but not the transcriptomes, you can omit this part, and just provide the `raw` subfolders):


   <details><!--1.3.section-->
   <summary>
   <b>1.3. The files for performing the omics integration </b><sup><sub>(click to show/hide)</sub></sup>
   </summary> 

   * The file(s) of *DEGs* (Differentially Expressed Genes):
     
     The LDHAB vs Control (but not Cycloserine) experiment had the two types of omics acquired
     under the same experimental conditions. The transcriptome differential analysis 
     (using DESeq2), produced the file `DEG_Control_LDHAB.csv`:
     ```
 
     ```
     For your data and analysis, you can add more DEGs files if available, with clear informative
     filenames to specify coherent integrations.
   
   * The file of *metabolites per pathway* (`pathways_kegg_metabolites.csv`):
     
     It has as column names the names of the pathways. The values below each column name correspond to the
     metabolites present in each pathway. It is allowed that a same metabolite appears in several pathways.
     The metabolites' names or identifiers must match with those appearing in the metabolite total
     abundances -or the mean enrichment- file.

     Only one file of this type is accepted.

   * The file of *genes or transcripts per pathway* (`pathways_kegg_transcripts.csv`):

     It has as column names the names of the pathways. The values below each column correspond to the gene
     symbols present in each pathway. It is allowed that a same gene symbol appears in several pathways. The
     gene symbols must match with those appearing in the DEGs file.
   
     Only one file of this type is accepted.

   </details><!--closes 1.3.section-->
  

    
</details><!--closes 1. section-->

-----------------------------------------------

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

## Available analyses

- _pca_analysis_ computes the PCA and outputs tables with principal components and explained variances
- _pca_plot_ generated classical PCA plots
- _abundance_plot_ plots with bars of total metabolite abundances
- _enrichment_plot_ generates lineplots of mean enrichment
- _isotopologues_plot_ generates stacked bars of isotopologue proportions
- _differential_analysis_ runs differential analysis and computes the corresponding statistics
- _differential_multigroup_analysis_ same as before, but for > 2 groups
- _timecourse_analysis_ runs differential analysis for time-course experiments in pairwise fashion for successif time points
- _metabologram_ network integration between SIRM and trascriptomic data, resulting in metabologram plots

  
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
