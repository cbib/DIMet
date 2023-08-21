DIMet: Differential analysis of Isotope-labeled targeted Metabolomics data
===

# Introduction

DIMet is a bioinformatics pipeline for differential and time-course analysis of targeted isotope-labelled data.

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

The user has to provide his data files within the data folder. Original data files temselves have to be placed in the <code>raw</code> subfolder along with a metadata file that contains the the experimental setup corresponding to the data. 

The structure of the <bf>metadata file</bf> has to contain 6 columns named <code>name_to_plot</code>, <code>timepoint</code>, <code>timenum</code>, <code>condition</code>, <code>compartment</code>, <code>original_name<code>. Here is the semantics of the columns:

- <code>name_to_plot</code> is the string that will appear on the figurs produced by DIMet
- <code>timepoint</code> is the sampling time as it is provided in the quantification files (it is an arbitary string that can contain non numerical characters)
- <code>timenum</code> is the numerical encoding of the <code>timepoint</code>
- <code>condition</code> is the experimental condition
- <code>compartment</code> is the name of the cellulaur compartment for which the measuring has been done (e.g. "endo", "endocellular", "cyto", etc)
- <code>original_name<code> contains the column names that is provided in the quantification files

<br>Example:</bf>

| name_to_plot | timepoint | timenum | condition | compartment | original_name |
|--------------|-----------|---------|-----------|-------------|---------------|
| Cond1 T0     | T0        | 0       | Cond1     | comp_name   | T0_Cond_1     |
| Cond1 T24    | T24       | 24      | Cond1     | comp_name   | T24_Cond_1    |
| Cond2 T0     | T0        | 0       | Cond2     | comp_name   | T0_Cond_2     |
| Cond1 T24    | T24       | 24      | Cond2     | comp_name   | T24_Cond_2    |


### Configuration files

As it is shown in the folder structure above, there are three types of configuration files, all of them in yaml format: 

1. `dataset configuration` files:  provide the quantifications and metadata file names that are themselves present in the `data` directory.

2. `analysis configuration` files: indicates which analysis (e.g. "differential analysis") is to be run for which data file and with which parameters.

3. `general configuration` files: configures the analysis to be run by both pointing to the `analysis configuration` and providing general information such as e.g. output folder names etc.

**Deeper explanations about the configuration files** will be given through the datasets used in our manuscript -datasets available at [Zenodo (manuscript_data)](https://sandbox.zenodo.org/record/1224020)-, in the section [Provided datasets](#provided-datasets), subsection 
<a href="#config_folder">The <code>config</code> folder</a>. 

Further examples of configuration files, with their respective minimal datasets are provided on [Zenodo (minimal_examples)](https://sandbox.zenodo.org/record/1224340).

## Provided datasets

Datasets, configuration and bash scripts corresponding to the results presented in the manuscript "DIMet: An open-source tool for Differential analysis of Isotope-labeled targeted Metabolomics data" by J. Galvis *et al*. are available at [Zenodo (manuscript_data)](https://sandbox.zenodo.org/record/1224020). 

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


<details><!--2. section-->
 
<summary id="config_folder"><b>2. The <code>config</code> folder</b> <sub><sup>(click to show/hide)</sup></sub></summary>

  <details>
  <summary>The structure of the  <code>config</code>  folder is:  
   <sub><sup>(click to show/hide)</sup></sub>
  </summary>

```
config
├── analysis
│   ├── abundance_plot_Cycloserine.yaml
│   ├── abundance_plot_LDHAB-Control.yaml
│   ├── dataset
│   │   ├── Cycloserine_data.yaml
│   │   ├── LDHAB-Control_data_integrate.yaml
│   │   └── LDHAB-Control_data.yaml
│   ├── differential_analysis_pairwise_LDHAB-Control.yaml
│   ├── enrichment_lineplot_Cycloserine.yaml
│   ├── isotopologues_plot_Cycloserine.yaml
│   ├── isotopologues_plot_LDHAB-Control.yaml
│   ├── metabologram_abundance_LDHAB-Control.yaml
│   ├── metabologram_enrichment_LDHAB-Control.yaml
│   ├── pca_plot_LDHAB-Control.yaml
│   ├── pca_tables_Cycloserine.yaml
│   └── timecourse_analysis_Cycloserine.yaml
├── general_config_abundance_plot_Cycloserine.yaml
├── general_config_abundance_plot_LDHAB-Control.yaml
├── general_config_differential_analysis_LDHAB-Control.yaml
├── general_config_enrichment_lineplot_Cycloserine.yaml
├── general_config_isotopologues_plot_Cycloserine.yaml
├── general_config_isotopologues_plot_LDHAB-Control.yaml
├── general_config_metabologram_abundance_LDHAB-Control.yaml
├── general_config_metabologram_enrichment_LDHAB-Control.yaml
├── general_config_pca_plot_LDHAB-Control.yaml
├── general_config_pca_tables_Cycloserine.yaml
└── general_config_timecourse_analysis_Cycloserine.yaml

```
   </details> 
   
We can appreciate that there are 3 different configuration files: 
* the *dataset configuration* files (such as `Cycloserine_data.yaml`), located inside the `config/analysis/dataset/` folder
* the *analysis configuration* files (such as `abundance_plot_Cycloserine.yaml`), located inside the `config/analysis/` folder
* the *general configuration* files (such as `general_config_pca_tables_Cycloserine.yaml`), located directly inside the `config` folder


The next sections explain each type of these files. Furthermore, when analyzing your own data, you can use these configuration files as templates and modify them coherently, for this purpose we will guide you across the next sections, to make it easier for you.

 <details><!--2.1. section -->
 <summary id="dataset-config"><b>2.1. The dataset configuration </b></summary>  
  
   For each given dataset, there must exist a corresponding *dataset configuration* file, which is located inside the `config/analysis/dataset/` folder.
   This file must specify:
  - the name of the metadata file
  - the names of the quantification files
  - the "order" of the conditions, where the *control* is the first to be listed (note that the name assigned to the control condition can be different than the one we used in our data).

    Below there is the `LDHAB-Control_data.yaml` file,  where the `# <-` comment indicates the
    parts that you can modify as needed when preparing your own analysis.
    
    Note that the parameter `_target_` must **not** be changed.

    If you do not have all the four types of quantification that are shown below, just delete
    the entire line making reference to the type of quantification that is absent.
  
  ```
  _target_: dimet.data.DatasetConfig

  label: LDHAB-Control  # <- change after the colon
  name: experiment in hypoxia, LDHAB vs Control  # <- change after the colon
  subfolder: LDHAB-Control_data  # <- subfolder name in the data folder, change after the colon
    
  # ALWAYS WITHIN THE RAW SUBFOLDER
  metadata: metadata_endo_ldh   # <- change after the colon
  abundances: AbundanceCorrected   # <- change after the colon
  mean_enrichment: MeanEnrichment13C   # <- change after the colon
  isotopologue_proportions: IsotopologuesProp   # <- change after the colon
  isotopologues: IsotopologuesAbs   # <- change after the colon
    
  conditions :
   - Control # <- first must be control
   - sgLDHAB # <- the rest of the conditions vertically listed
  ```
  This dataset configuration will be used by all the types of analysis we offer except the omics integration, see the section *The integration dataset configuration* for that purpose.

   <details><!--2.1.1. section-->
   <summary>2.1.1.The integration dataset configuration (for <b>omics integration</b>)</summary>
  
   When a data folder for performing the **omics integration** exists, there must also exist a *dataset configuration* file.
   
   Below there is the `LDHAB-Control_data_integrate.yaml` file,  where the `# <-` comment indicates the
   parts that you can modify as needed when preparing your own analysis. 
    
   Note that the parameter `_target_` must **not** be changed.
    
    
    ```
    _target_: dimet.data.DataIntegrationConfig  

    label: integrate_DAM_DEG_LDHAB-Control   # <- change after the colon
    name: Integration of two omics LDHAB vs Control   # <- change after the colon
    subfolder: LDHAB-Control_data  # <- subfolder name in the data folder, change after the colon
    
    conditions :
     - Control # <- first must be control
     - sgLDHAB # <- the rest of the conditions vertically listed
    
    # ALWAYS WITHIN THE RAW SUBFOLDER
    metadata: metadata_endo_ldh   # <- change after the colon
    abundances: AbundanceCorrected   # <- change after the colon
    mean_enrichment: MeanEnrichment13C   # <- change after the colon
    
    # WITHIN THE INTEGRATION_FILES SUBFOLDER
    transcripts:  
      - DEG_Control_LDHAB   # <- file name
    
    pathways: 
      metabolites: pathways_kegg_metabolites  # <- file name, change after the colon
      transcripts: pathways_kegg_transcripts   # <- file name, change after the colon

    ```

   </details><!--closes 2.1.1. section-->

  End of the dataset configuration section.


 
  </details><!--closes 2.1 section-->

  <details><!--2.2. section-->
  <summary><b>2.2. The analysis configuration </b></summary>   
   
   The *analysis configuration* is located inside the `config/analysis/` folder. For each analysis to be performed, one analysis configuration file must exist. It indicates which is the type of analysis we want to run, on which dataset this analysis will be applied, and the parameters that are specific to that analysis.

   Below the file `differential_analysis_pairwise_LDHAB-Control.yaml` illustrates a pairwise differential analysis configuration. The `# <-` comment indicates the parts that you can modify as needed when preparing your own analysis.
   
     ```
     label: differential-analysis-LDHAB-Control   # <- change after the colon

     defaults:
       - dataset: LDHAB-Control_data     # <- change after the colon
       - method: differential_analysis
     
     comparisons :
       - [[sgLDHAB, T48], [Control, T48]]   # <-      
     
     statistical_test:
       abundances: prm-scipy   # <- see statistic test options, change after the colon
       mean_enrichment: prm-scipy    # <- see options, change after the colon
       isotopologues: disfit    # <- see options, change after the colon
       isotopologue_proportions: ranksum    # <- see options, change after the colon

     ```
   
   <!-- COMPLETE: Currently the statistic test options are: etc etc -->
   
   
   Note that the parameters in the *analysis configuration* must be also coherent with our dataset.
   Feel free to reuse all our analysis configuration files, most of them are really short, simple and easy to understand and modify.
  
        
   <details><!--2.2.1. section-->
   <summary>2.1.1. The analysis configuration for the <b>omics integration</b> </summary> 

   The integration performed by DIMet displays the 'differential omics' in the form of Metabologram(s). The analysis configuration file `metabologram_abundance_LDHAB-Control.yaml` uses the differential metabolites' abundances (the `# <-` comment indicates the parts that you can modify as needed when preparing your own analysis):
    
       ```
       label: metabologram-using-abundance-LDHAB-Control  # <- change after the colon
   
       defaults:
         - dataset: LDHAB-Control_data_integrate  # <- integration dataset config
         - method: metabologram_integration
       
       comparisons :
         - [[sgLDHAB, T48], [Control, T48]]   # <-  see note (**) 
       
       # running for total abundances
       statistical_test:
         abundances: prm-scipy   # <- change after the colon
       
       columns_metabolites:
         ID : metabolite
         values : log2FC   # <- or FC, change after the colon
         
       columns_transcripts:
         ID: external_gene_name  # <- the gene symbols column, change after the colon
         values: log2FoldChange  # <- change after the colon
         
       compartment:
         en    
        ```
   (`**`) each metabolomics comparison must match with a DEG file. In this `metabologram_abundance_LDHAB-Control.yaml` there is only 1 comparison. When > 1 comparisons, they must be listed the same strict order as the DEG files listed in the *integration dataset configuration* (check inside section <a href="#dataset-config">dataset configuration</a>). 
   
   If you need to run the integration using the differential mean enrichment, a separate analysis configuration file must be created, we provide `metabologram_enrichment_LDHAB-Control.yaml` that you can use as template.    
        
   </details><!--closes 2.2.1. section-->  
   
  </details><!--closes 2.2. section-->


  <details><!--2.3. section-->
  <summary><b>2.3. The general configuration </b></summary>    
   
   The *general configuration* file is located directly inside the `config` folder. It must declare the name of its corrresponding *analysis configuration* file. Place yourself in `defaults:` then `- analysis:`, see below:
     
     ```
     hydra:
       job:
         chdir: true
       run:
         dir: outputs/${now:%Y-%m-%d}/${now:%H-%M-%S}/${analysis.dataset.label}-${analysis.method.label}

     defaults:
       - analysis: differential_analysis_example2  # <-- HERE the 'analysis configuration' name

     figure_path: figures
     table_path: tables    
     ```
   When preparing your own analysis, use the above as template and simply replace `differential_analysis_example2` by the name of your analysis.
   
   This *general configuration* file also provides the parameters for generating the output folder names; keep these 
   output parameters unmodified, so DIMet will generate the full output folders and files names
   automatically.
     
   </details><!--closes 2.3. section--> 
   
   Note that across all the types of configuration files, any referenced file name must be written without the extension.
     
</details><!--closes  2. section -->



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
