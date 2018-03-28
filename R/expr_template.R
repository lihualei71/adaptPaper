if (!requireNamespace("adaptMT", quietly = TRUE)){
    devtools::install_github("lihualei71/adaptMT")
}
library("adaptMT")
library("splines")

### set up all methods
## Accumulation Tests (Li & Barber, 2016a)
## Code from 'https://www.stat.uchicago.edu/~rina/accumulationtests/accumulation_test_functions.R'
source("accumulation_test_functions.R")

## Auxiliary functions of SABHA (Li & Barber, 2016b)
## Code from 'http://www.stat.uchicago.edu/~rina/sabha/All_q_est_functions.R'
source("All_q_est_functions.R")
source("other_methods.R")

adapt_expr <- function(x, pvals){
    oldw <- options(warn = -1)
    
    alphas <- seq(0.01, 0.3, by=0.01) # target FDR level
    tau <- 0.5; eps <- 0.1 # parameters for SABHA
    thr <- 0.5 # parameter for Storey-BH
    thr1 <- 0.1; thr2 <- 0.5 # parameters for adaptive SeqStep
    ADMM_params <- c(10^2, 10^3, 2, 5000, 1e-3) # alpha_ADMM,beta,eta,max_iters,converge_thr
    n <- length(pvals)

    num_alpha <- length(alphas)
    max_alpha <- max(alphas)
    ## gather results
    NumRej <- matrix(0, nrow = 13, ncol = num_alpha)

    ## methods: 1 SeqStep, 2 HingeExp, 3 ForwardStop, 4 Adaptive SeqStep, 5 BH, 6 Storey-BH, 7 Barber-Candes, 8 SABHA (step), 9 SABHA (ordered), 10 IHW, 11 IHW (oracle), 12 IF (oracle), 13 AdaPT
    qhat_step <- Solve_q_step(pvals, tau, eps)
    qhat_ordered <- Solve_q_ordered_simple(pvals, tau, eps, ADMM_params)
    for(i in 1:num_alpha){
        alpha <- alphas[i]
        NumRej[1, i] <- SeqStep(pvals, alpha = alpha, C = 2)
        NumRej[2, i] <- HingeExp(pvals, alpha = alpha)
        NumRej[3, i] <- ForwardStop(pvals, alpha = alpha)
        NumRej[4, i] <- length(Adaptive_SeqStep_method(pvals, alpha, alpha, tau))
        NumRej[5, i] <- length(BH_method(pvals, alpha))
        NumRej[6, i] <- length(Storey_method(pvals, alpha, thr))
        NumRej[7, i] <- length(BC_method(pvals, alpha))
        NumRej[8, i] <- length(SABHA_method(pvals, qhat_step, alpha, tau))
        NumRej[9, i] <- length(SABHA_method(pvals, qhat_ordered, alpha, tau))
        NumRej[10, i] <- rejections(ihw(pvals, x[, 1], alpha))
        NumRej[11, i] <- ihw_oracle(pvals, 1:n, alpha)
        NumRej[12, i] <- IF_oracle(pvals, x[, 1], alpha)
        print(paste0("Step ", i, " finished!"))
    }

    print("AdaPT starts!")

    pi.formulas <- paste0("ns(x, df = ", 6:10, ")")
    mu.formulas <- paste0("ns(x, df = ", 6:10, ")")
    formulas <- expand.grid(pi.formulas, mu.formulas)
    pi_formulas <- as.character(formulas[, 1])
    mu_formulas <- as.character(formulas[, 1])

    res_adapt <- adapt_glm(x, pvals,
                           dist = beta_family(),
                           pi_formulas = pi_formulas,
                           mu_formulas = mu_formulas,
                           nfits = 50)
    nrejs <- res_adapt$nrejs
    NumRej[13,] <- nrejs[1:num_alpha]

    result <- list(NumRej = NumRej, adapt = res_adapt)
    
    options(oldw)
    return(result)
}
