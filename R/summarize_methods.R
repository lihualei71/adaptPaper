#---------------------------------------------------------------
# Helpers: summarize all other methods and compute FDR and power
#---------------------------------------------------------------

## BH procedure (Benjamini & Hochberg, 1995)
summary_BH <- function(pvals, H0,
                       alphas = seq(0.01, 0.3, 0.01)){
    n <- length(pvals)
    results <- sapply(alphas, function(alpha){
        khat <- max(c(0, which(sort(pvals) <= alpha * (1:n) / n)))
        alpha <- alpha * khat / n
        nfrej <- sum(pvals[!H0] < alpha, na.rm = TRUE)
        ntrej <- sum(pvals[H0] < alpha, na.rm = TRUE)
        return(c(nfrej, ntrej))
    })
    nfrej <- as.numeric(results[1, ])
    ntrej <- as.numeric(results[2, ])
    nrej <- nfrej + ntrej
    FDP <- nfrej / pmax(nrej, 1)
    power <- ntrej / max(sum(H0), 1)
    df <- data.frame(nrej = nrej, FDP = FDP, power = power)
    return(df)
}

## Storey's BH Procedure (Storey et al. 2004)
summary_storey <- function(pvals, H0, thr = 0.5,
                           alphas = seq(0.01, 0.3, 0.01)){
    pi0 <- min(1, mean(pvals > thr) / (1 - thr))
    pvals[pvals > thr] = Inf
    n <- length(pvals)
    results <- sapply(alphas, function(alpha){
        khat <- max(c(0, which(sort(pvals) <= alpha / pi0 * (1:n) / n)))
        alpha <- alpha * khat / n
        nfrej <- sum(pvals[!H0] < alpha, na.rm = TRUE)
        ntrej <- sum(pvals[H0] < alpha, na.rm = TRUE)
        return(c(nfrej, ntrej))        
    })
    nfrej <- as.numeric(results[1, ])
    ntrej <- as.numeric(results[2, ])
    nrej <- nfrej + ntrej
    FDP <- nfrej / pmax(nrej, 1)
    power <- ntrej / max(sum(H0), 1)
    df <- data.frame(nrej = nrej, FDP = FDP, power = power)
    return(df)
}

## Barber-Candes Procedure (Barber and Candes, 2016)
summary_BC <- function(pvals, H0,
                       alphas = seq(0.01, 0.3, 0.01)){
    sorted.mask.pvals <- sort(pmin(pvals, 1 - pvals))
    fdphat <- sapply(sorted.mask.pvals, function(thresh){
        (1 + sum(pvals >= 1 - thresh)) / max(1, sum(pvals <= thresh))
    })
    results <- sapply(alphas, function(alpha){
        khat <- which(fdphat <= alpha)
        if (length(khat) == 0){
            return(c(0, 0))
        } else {
            khat <- max(khat)
            phat <- sorted.mask.pvals[khat]
            nfrej <- sum(pvals[!H0] <= phat, na.rm = TRUE)
            ntrej <- sum(pvals[H0] <= phat, na.rm = TRUE)
            return(c(nfrej, ntrej))
        }
    })
    nfrej <- as.numeric(results[1, ])
    ntrej <- as.numeric(results[2, ])
    nrej <- nfrej + ntrej
    FDP <- nfrej / pmax(nrej, 1)
    power <- ntrej / max(sum(H0), 1)
    df <- data.frame(nrej = nrej, FDP = FDP, power = power)
    return(df)
}

## IHW (Ignatiadis and Huber, 2016)
library("IHW")
summary_IHW <- function(pvals, H0,
                        alphas = seq(0.01, 0.3, 0.01)){
    results <- sapply(alphas, function(alpha){
        rej <- rejected_hypotheses(ihw(pvals, 1:length(pvals), alpha))
        nfrej <- sum(rej[!H0], na.rm = TRUE)
        ntrej <- sum(rej[H0], na.rm = TRUE)
        return(c(nfrej, ntrej))
    })
    nfrej <- as.numeric(results[1, ])
    ntrej <- as.numeric(results[2, ])
    nrej <- nfrej + ntrej
    FDP <- nfrej / pmax(nrej, 1)
    power <- ntrej / max(sum(H0), 1)
    df <- data.frame(nrej = nrej, FDP = FDP, power = power)
    return(df)    
}

## AdaPT (Lei and Fithian, 2016)
summary_adapt <- function(adapt, H0, pvals){
    results <- apply(adapt$s, 2, function(s){
        tmp <- (pvals <= s)
        nfrej <- sum(tmp[!H0], na.rm = TRUE)
        ntrej <- sum(tmp[H0], na.rm = TRUE)
        return(c(nfrej, ntrej))
    })
    nfrej <- as.numeric(results[1, ])
    ntrej <- as.numeric(results[2, ])
    nrej <- nfrej + ntrej
    FDP <- nfrej / pmax(nrej, 1)
    power <- ntrej / max(sum(H0),1)
    df <- data.frame(nrej = nrej, FDP = FDP, power = power)
    return(df)
}

## SABHA (Li & Barber, 2016b)
summary_SABHA <- function(pvals, H0, qhat,
                          alphas = seq(0.01, 0.3, 0.01)){
    n <- length(pvals)
    results <- sapply(alphas, function(alpha){
        khat <- max(c(0, which(sort(qhat * pvals) <= alpha * (1:n) / n)))
        alpha <- alpha * khat / n
        nfrej <- sum(qhat[!H0] * pvals[!H0] < alpha, na.rm = TRUE)
        ntrej <- sum(qhat[H0] * pvals[H0] < alpha, na.rm = TRUE)
        return(c(nfrej, ntrej))
    })
    nfrej <- as.numeric(results[1, ])
    ntrej <- as.numeric(results[2, ])
    nrej <- nfrej + ntrej
    FDP <- nfrej / pmax(nrej, 1)
    power <- ntrej / max(sum(H0), 1)
    df <- data.frame(nrej = nrej, FDP = FDP, power = power)
    return(df)
}
