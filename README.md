# scAM

[![Build Status](https://github.com/murti-abhishek/scAM.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/murti-abhishek/scAM.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://travis-ci.com/murti-abhishek/scAM.jl.svg?branch=main)](https://travis-ci.com/murti-abhishek/scAM.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/murti-abhishek/scAM.jl?svg=true)](https://ci.appveyor.com/project/murti-abhishek/scAM-jl)
[![Coverage](https://codecov.io/gh/murti-abhishek/scAM.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/murti-abhishek/scAM.jl)

This is a Julia package for analysing single cell and single nuclei RNA sequencing data.

If you're familiar with the Seurat (R) and Scanpy (Python), this package has similar functionality.

| Function                 | Usage                                                                                                                      |
| ------------------------ | ----------------------------------------------------------------------------------------------------------------------------- |
| create_scAMobj           | Creates the a custom object that this package relies on; Just provide the path to the folder where your dge files (gene expression matrices) are stored |
| merge_scAMobjs           | Merge scAM objects                                                                                                                                      |
| percentage_set_feature   | Function for adding pct feature set [e.g. percent mito genes]                                                                                           |
| MetricsPlot              | Function to plot counts, features and prefix_mito                                                                                                       |
| modify_scAMobj           | Modify the object after QC                                                                                                                              |
| Normalize                | Normalize the raw counts                                                                                                                                |
| Scale                    | Scale the counts                                                                                                                                        |
| Calculate_UMAP_reduction | Calulate_UMAP_reduction                                                                                                                                 |
| Cluster                  | Run the clustering algorithm                                                                                                                            |
| plot_UMAP                | Plot the UMAP with the clusters                                                                                                                         |
| FeaturePlot              | Feature plot specific genes                                                                                                                             |
| FindAllDGE               | Find the differentially expressed genes for each cluster                                                                                                |
| FindDGE                  | Find the differentially expressed genes for a specific cluster                                                                                          |
| plot_proportions         | Plot proportions by any two metadata identity                                                                                                           |
| ViolinPlot               | Make violin plots of specific genes by any metadata identity                                                                                            |
