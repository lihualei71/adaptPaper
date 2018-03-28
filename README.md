# 'adaptPaper' Repository

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

## A toy test
Our experimental results take several hours. We put a test file `test.R` that contains a toy task to check whether the functions work.
We recommend running it before replicating the results. If any error is encountered, please report the issue [here](https://github.com/lihualei71/adaptPaper/issues) or contact me by email.

```
cd R
Rscript test.R
```
