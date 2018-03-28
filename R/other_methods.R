#---------------------------------------------------------------
# Helpers: Setup all other methods
#---------------------------------------------------------------

## Accumulation Tests (Li & Barber, 2016a)
## Code from 'https://www.stat.uchicago.edu/~rina/accumulationtests/accumulation_test_functions.R'
source("accumulation_test_functions.R")

## Auxiliary functions of SABHA (Li & Barber, 2016b)
## Code from 'http://www.stat.uchicago.edu/~rina/sabha/All_q_est_functions.R'
source("All_q_est_functions.R")

## Independent Hypothesis Weighting (Ignatiadis et al. 2016)
if (!requireNamespace("IHW", quietly = TRUE)){
    ## try http:// if https:// URLs are not supported
    source("https://bioconductor.org/biocLite.R")
    biocLite("IHW")
}
library("IHW")

## BH procedure (Benjamini & Hochberg, 1995)
BH_method <- function(pvals, alpha){
    n <- length(pvals)
    khat <- max(c(0, which(sort(pvals) <= alpha * (1:n) / n)))
    which(pvals <= alpha * khat / n)
}

## Storey's BH Procedure (Storey et al. 2004)
Storey_method <- function(pvals,alpha,thr){
    n <- length(pvals)
    pi0 <- min(1, mean(pvals > thr) / (1 - thr))
    pvals[pvals > thr] <-  Inf
    khat <- max(c(0, which(sort(pvals) <= alpha * (1:n) / n / pi0)))
    which(pvals <= alpha * khat / n / pi0)
}

## Barber-Candes procedure (Barber and Candes, 2015)
BC_method <- function(pvals,alpha){
    sorted_mask_pvals <- sort(pmin(pvals, 1 - pvals))
    fdphat <- sapply(sorted_mask_pvals, function(thresh){
        (1 + sum(pvals >= 1 - thresh))/max(1, sum(pvals <= thresh))
    })
    khat <- which(fdphat <= alpha)
    if (length(khat) == 0){
        return(khat)
    } else {
        khat <- max(khat)
        phat <- sorted_mask_pvals[khat]
        return(which(pvals <= phat))
    }
}

## SABHA (Li & Barber, 2016b)
Solve_q_ordered_simple <- function(pvals, tau, eps, ADMM_params){
    target_num <- 5000
    n <- length(pvals)
    if (n <= target_num){
        qhat <- Solve_q_ordered(pvals, tau, eps, ADMM_params) 
        return(qhat)
    }
    num_reps <- ceiling(n / target_num)
    new_pvals <- sapply(1:target_num, function(i){
        ind <- pmin((i - 1) * num_reps + (1:num_reps), n)
        pvals[ind[1]]
    })
    qhat <- Solve_q_ordered(new_pvals, tau, eps, ADMM_params)
    qhat <- rep(qhat, each = num_reps)[1:n]
    return(qhat)
}

SABHA_method <- function(pvals, qhat, alpha, tau){
  # Use the original, or estimated q as input
    n <- length(pvals)
    pvals[pvals > tau] <- Inf
    khat <- max(c(0, which(sort(qhat * pvals) <= alpha * (1:n) / n)))
    which(qhat * pvals <= alpha * khat / n)
}

## Adaptive Seqstep (Lei & Fithian, 2016)
Adaptive_SeqStep_method <- function(pvals, alpha, thr1, thr2){ # Lei & Fithian 2016's method
  # thr1 & thr2 correspond to s & lambda (with s<=lambda) in their paper
  fdphat <- thr1 / (1 - thr2) * (1 + cumsum(pvals > thr2)) / pmax(1, cumsum(pvals <= thr1))
  if(any(fdphat <= alpha)){
    khat <- max(which(fdphat <= alpha))
    return(which(pvals[1:khat] <= thr1))
}else{
    return(NULL)
  }
}

## Independent Filtering (with oracle threshold)
if (!requireNamespace("genefilter", quietly = TRUE)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("genefilter")
}
library("genefilter")

IF_oracle <- function(pvals, x, alpha){
    theta_list <- seq(0, 1, 0.05)
    R_list <- sapply(theta_list, function(theta){
        R1 <- filtered_R(filter = x, test = pvals,
                         theta = theta, method = "BH",
                         alpha = alpha)
        R2 <- filtered_R(filter = rev(x), test = rev(pvals),
                         theta = theta, method = "BH",
                         alpha = alpha)
    })
    return(max(R_list))
}

## IHW (with oracle nbins)
ihw_oracle <- function(pvals, x, alpha){
    ## nrejs <- sapply(1:15, function(nbins){
    ##     rejections(ihw(pvals, x, alpha, nbins = nbins))
    ## })
    ## max(nrejs)
    rejections(ihw(pvals, x, alpha, nbins = 15))
}

