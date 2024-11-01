This repository includes the code used in the *Recursive clustering of cellular diversity in scRNA-seq* paper.

Benchmark dataset files are made available for download [here](https://drive.google.com/file/d/1Gbm7U6pvKWmEv3ZotuJRxA4oj4tIWSVS/view?usp=sharing).
Alternatively, benchmark dataset files can be reconstructed from scratch by following the Snakemake pipelines available [here](https://github.com/satijalab/azimuth-references) within the human_pbmc, human_adipose, human_tonsil, and human_fetus folders.  

Before running code for the benchmark analysis, place the folder azimuth-references (as obtained in one of the two linked locations above) which contains the folders human_pbmc, human_adipose, human_tonsil, and human_fetus, inside the directory for which the R code is run, keeping the default file structure.  If the azimuth-references folder was downloaded from the Azimuth references Github, the Snakemake pipelines must be run for each reference, according to the documentation provided in the Azimuth references Github.

The Crohn's disease data can be downloaded [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202052).  Download GSE202052_RAW.tar, and unzip it to to a folder named "GSE202052_RAW" in the directory of your R session.

Code is in a script format, written in R.

To reproduce the exact results in the paper, use R version 4.0.5, and run **Install Libraries.R** to install the exact library versions used.
For Fig. S13 through S15, you will also need to use Python version 3.10, and first run **Install Packages.py** to install the packages used.

A windows 10 environment was used to run **main.R** when producing the clustering results visualized in Fig. 1 through 6, Fig. S5 through S12, Fig. S16 and Fig. S17, and to run all other code besides **main.R**.
A Ubuntu V24 environment was used to run **main.R** when producing the clustering results visualized in Fig. S1 through S4, and Fig. S13 through S15.

To reproduce Fig. 1 and Fig. 2:
  - Run file **main.R** with reference = 'PBMC', c_resolution = 0.025, algorithm = 'Louvain', HVGs = 2000, downsample = 0
  - Run file **TreePlot.R**.  This will produce Fig. 1.
  - Run file **Generate PBMC Example.R**.  This will produce Fig. 2, and print out the relative mean cell type purity score between recursive and single-pass.

To reproduce Fig. 3, Fig. 4, Fig. S7 through S12:
  - Run file **main.R** with reference equal to each of the benchmark dataset values ('PBMC', 'Adipose', 'Tonsil', 'Fetus'), with c_resolution equal to each value of the default resolution range (0.015 to 0.15 by 0.005 increments), algorithm = 'Louvain', HVGs = 2000, downsample = 0
  - Run file **metrics.R** with algorithm = 'Louvain', HVGs = 2000, is_res_range_small = FALSE.
  - Run file **balanced_metrics.py**
  - Run file **Metric Plot.R** with algorithm = 'Louvain', HVGs = 2000, is_res_range_smalls = FALSE, evaluation_types = 'cell type'.  This will produce Fig. 3, Fig. S7 through Fig. S12
  - Run file **Relative Cell Type Performance.R**.  This will produce Fig. 4, and will also print out the frequently better captured analysis.

To reproduce Fig. 5 and Fig. 6:
  - Run file **main.R** with reference = 'Helmsley, with c_resolution = 0.1, algorithm = 'Louvain', HVGs = 2000, downsample = 0, entropy_cutoff = 0.5
  - Run file **TreePlot.R**
  - Run Crohn's R files 1 through 7 (**Crohn's 6 - Shared Variables Tree.R** produces Fig. 5, and **Crohn's 7 - Create combined tree plot.R** produces Fig. 6)

To reproduce Fig. S1 through S4:
  - With:
      - algorithm = 'Leiden', HVGs = 2000
      - algorithm = 'SLM', HVGs = 2000
      - algorithm = 'Louvain', HVGs = 1000
      - algorithm = 'Louvain', HVGs = 3000
  - Change the above algorithm and HVGs setting for each of the the following and:
      - Run file **main.R** with reference equal to each of the benchmark dataset values ('PBMC', 'Adipose', 'Tonsil', 'Fetus'), with c_resolution equal to each value of the default resolution range (0.015 to 0.15 by 0.005 increments), downsample = 0
      - Run file **metrics.R** with is_res_range_small = FALSE
      - Run file **Metric Plot.R** with is_res_range_smalls = FALSE, evaluation_types = 'cell type'.  This will produce Fig. S1 through S4.
   
  To reproduce Fig. S5 and S6:
  - Run file **main.R** with reference equal to each of the benchmark dataset values ('PBMC', 'Adipose', 'Tonsil', 'Fetus'), with c_resolution equal to 0.005 to 0.05 by increments of 0.002 for the PBMC, adipose, and tonsil references, and c_resolution equal to 0.0025 to 0.01 by increments of 0.0002 for the fetus reference, algorithm = 'Louvain', HVGs = 2000, downsample = 0 (assuming **main.R** has already been run for the default range at some point).
  - Run file **metrics.R** with algorithm = 'Louvain', HVGs = 2000, is_res_range_small = TRUE.
  - Run file **Metric Plot.R** with algorithm = 'Louvain', HVGs = 2000, is_res_range_smalls = c(TRUE, TRUE), evaluation_types = 'cell type'.  This will produce Fig. S5 and S6.

  To reproduce Fig. S13 through S15:
  - Run file **main.R** with reference equal to each of the benchmark dataset values ('PBMC', 'Adipose', 'Tonsil', 'Fetus'), with c_resolution equal to each value of the default resolution range (0.015 to 0.15 by 0.005 increments), algorithm = 'Louvain', HVGs = 2000, for each value of downsample equal to 1 through 5
  - Run file **consistency metrics.R**
  - Run file **shared cluster members.R**
  - Run file **Metric Plot.R** with algorithm = 'Louvain', HVGs = 2000, is_res_range_smalls = FALSE, evaluation_types = 'consistency'.  This will produce Fig. S13 through S15.

To reproduce Fig. S16 and S17:
  - Run file **main.R** with reference = 'Helmsley, with c_resolution = 0.1, algorithm = 'Louvain', HVGs = 2000, downsample = 0, entropy_cutoff = 0
  - Run file **TreePlot.R**
  - Run **Crohn's 6 - Shared Variables Tree.R** to produce Fig. S16
  - Run **UMAP.R** to produce Fig. S17


