This repository includes the code used in the *Recursive clustering of cellular diversity in scRNA-seq* paper.

Code is in a script format, written in R, where:
  - To peform rercusive and single-pass equivalent clustering analysis, file **main.R** is run for a variety of resolutions for each dataset, then:
      - **Benchmark Analysis.R**: Generates the benchmark mean cell type cluster purity metric plots for the four scRNA-seq benchmark datasets (Figure 3).
      - **Relative Cell Type Performance.R**: Generates scatter plots of cell type cluster purity for recursive verus single-pass method results (Figure 4).
      - **Generate PBMC Example.R**: Generates heatmap plots for recursive and seurat equivalent cluster purity results (Figure 2 when run after main.R has been run with c_resolution = 0.02 and reference_name = "human_PBMC/").
      - **TreePlot.R**:  Generates hierarchical clustering tree plots (Figure 1 when run after main.R has been run with c_resolution = 0.02 and reference_name = "human_PBMC/").
  - To analyze the Crohn's disease dataset, main.R was run with is_helmsley=TRUE and c_resolution = 0.1, then:
      - **TreePlot.R** was run.
      - Crohn's R files 1-7 were run.
  - The analysis was performed using Seurat V4 and R version 4.0.5.
      - **renv.lock** is a renv file containing the exact versions of all installed libraries at the time of analysis.
      - **Install Libraries.R** may be a useful script for getting matching R libraries installed.
   
Benchmark dataset files can be reconstructed from scratch by following the snakemake pipelines available [here](https://github.com/satijalab/azimuth-references) within the human_pbmc, human_adipose, human_tonsil, and human_fetus folders.
Benchmark dataset files are also made available for download [here](https://drive.google.com/file/d/1Gbm7U6pvKWmEv3ZotuJRxA4oj4tIWSVS/view?usp=sharing).

The Crohn's disease dataset is currently not publicly available due to instituntional policy.
