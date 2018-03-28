# Paper Repository

This repository contains all R code to replicate the results and the figures in our paper: [AdaPT: An interactive procedure for multiple testing with side information](https://arxiv.org/abs/1609.06035). 

## Introduction
All real experiments (Section 5.1 and Section 5.3), except the gene/drug example with random orderings, can be done in a laptop. The simulations, as well as the gene/drug example with random orderings, should be ran on a cluster. All R files are listed in the folder `R/`. The bash files to submit jobs to the cluster are listed in the folder `bash/` (note that this depends on your cluster and the bash file might need to be changed correspondingly. See the section "cluster jobs" below for details). The outputs and the plots have already been contained in the folder `data/` and the folder `figs/`, respectively. 

## Installing the package
We are developing an R package [`adaptMT`](https://github.com/lihualei71/adaptMT). The file `adaptMR_0.2.0.9000.tar.gz` is the version that generates the results and the plots in the published paper. To replicate the results, we recommend installing the source package directly because later we might update the package. 

```
git clone https://github.com/lihualei71/adaptPaper.git
cd adaptPaper
R CMD INSTALL adaptMT_0.2.0.9000.tar.gz
```

## R files
The folder `R/` contains all R files:

- `test.R` gives a toy test that checks whether the functions work;
- `exprs_laptop.R` implements all experiments (Section 5.1 and Section 5.3) that can be done on a laptop;
- `estrogen_random_cluster.R`, `simul1_cluster.R` and `simul2_cluster.R` implement the tasks (Section 5.1 and Section 5.2) that must be done on a cluster;
- `aggregate.R` aggregates the outputs from the tasks ran on clusters;
- `adapt_exprs_plot.R` reproduces all figures in Section 5;
- all other files are helpers.

### A toy test
Our experimental results take several hours. We put a test file `test.R` that contains a toy task to check whether the functions work.
We recommend running it before replicating the results. If no error is thrown to the console, the other files should be implemented without problems. The test takes around 2 minutes. If any error is encountered, please report the issue [here](https://github.com/lihualei71/adaptPaper/issues) or contact me by email. 

```
cd R
R CMD BATCH test.R
cd ..
```
### Experiments on the laptop
`exprs_laptop.R` generates the results for 

- gene/drug data, Section 5.1, with highly informative ordering (`data/estrogen_high_res.RData`), with moderately informative ordering (`data/estrogen_mod_res.RData`),
- Bottomly data, Section 5.3 (`data/Bottomly_res.RData`), 
- airway data, Section 5.3 (`data/airway_res.RData`), 
- pasilla data, Section 5.3 (`data/pasilla_res.RData`),
- SILAC data, Section 5.3 (`data/proteomics_res.RData`). 

The task takes around several hours. Please be patient.

```
cd R
Rscript exprs_laptop.R 
cd ..
```
### Experiments on the cluster
`estrogen_random_cluster.R`, `simul1_cluster.R` and `simul2_cluster.R` should be implemented on a cluster. Each file reads two environmental variables: `times` and `seed`. `times` specifies the number of repeated experiments and `seed` specifies the random seed. In the folder `bash/`, there are three `.sh` files corresponding to each task. Each `.sh` files generates 50 jobs with `times = 2` and `seed` ranging from 0 to 49. The outputs should contain 50 `.RData` files for each job, which are not listed in this repository. Finally, `aggregate.R` merges all outputs and turn them into a single `.RData` file for each task, namely `estrogen_random.RData`, `simul1.RData` and `simul2.RData`.

### Figures 
`adapt_exprs_plot.R` generates all figures in Section 5.

```
cd R
R CMD BATCH adapt_exprs_plot.R 
cd ..
```
