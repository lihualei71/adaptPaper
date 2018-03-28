#---------------------------------------------------------------
# Generate the results for Simulation 1.
# Need to be ran on clusters. See README.
#---------------------------------------------------------------

source("other_methods.R")
source("summarize_methods.R")
if (!requireNamespace("adaptMT", quietly = TRUE)){
    devtools::install_github("lihualei71/adaptMT")
}
library("adaptMT")
library("mgcv")

repeat_times <- as.numeric(Sys.getenv("times"))
seed <- as.numeric(Sys.getenv("seed"))
set.seed(seed)

simul1 <- function(x, mu, H0, 
                   pi_formula, mu_formula, 
                   alphas, repeat_times,
                   ...){
    n <- length(mu)
    n1 <- n2 <- floor(sqrt(n))
    m <- length(alphas)
    summary_FDP <- list()
    summary_power <- list()
    for (j in 1:7){
        summary_FDP[[j]] <- matrix(rep(0, repeat_times * m),
                                   ncol = repeat_times)
        summary_power[[j]] <- matrix(rep(0, repeat_times * m),
                                     ncol = repeat_times)
    }

    ## Settings for SABHA. See Li & Barber (2016)
    rhoG <- 0.9037629 # max_k ||(D_G+)_k||_2
    TV_bd <- 10
    tau <- 0.5
    eps <- 0.1
    ADMM_params <- c(10^2, 10^3, 2, 5000, 1e-3) 
    
    for (i in 1:repeat_times){
        z <- rnorm(n) + mu
        pvals <- 1 - pnorm(z)
        
        BH_result <- summary_BH(pvals, H0, alphas)
        summary_FDP[[1]][, i] <- BH_result[, 2]
        summary_power[[1]][, i] <- BH_result[, 3]

        storey_result <- summary_storey(pvals, H0,
                                        alphas = alphas)
        summary_FDP[[2]][, i] <- storey_result[, 2]
        summary_power[[2]][, i] <- storey_result[, 3]

        BC_result <- summary_BC(pvals, H0, alphas)
        summary_FDP[[3]][, i] <- BC_result[, 2]
        summary_power[[3]][, i] <- BC_result[, 3]

        IHW_result <- summary_IHW(pvals, H0, alphas)
        summary_FDP[[4]][, i] <- IHW_result[, 2]
        summary_power[[4]][, i] <- IHW_result[, 3]

        qhat <- Solve_q_TV_2dim(matrix(pvals, n1, n2), tau, eps, TV_bd, ADMM_params)
        SABHA_result <- summary_SABHA(pvals, H0, qhat, alphas)
        summary_FDP[[5]][, i] <- SABHA_result[, 2]
        summary_power[[5]][, i] <- SABHA_result[, 3]

        SABHA.factor <- 1 + 1 / eps / (1 - tau) / sqrt(n) +
            2 * rhoG * TV_bd * sqrt(log(n)) / eps^2 / n
        SABHA_result2 <- summary_SABHA(pvals, H0, qhat,
                                       alphas / SABHA.factor)
        summary_FDP[[6]][, i] <- SABHA_result2[, 2]
        summary_power[[6]][, i] <- SABHA_result2[, 3]
        
        res <- try(adapt_gam(x, pvals,
                             pi_formula = pi_formula,
                             mu_formula = mu_formula,
                             alphas = alphas))
        if (class(res) != "try-error"){
            adapt_result <- summary_adapt(res, H0, pvals)
            summary_FDP[[7]][, i] <- adapt_result[, 2]
            summary_power[[7]][, i] <- adapt_result[, 3]
        }

        print(paste0(i, "-th step finishes!"))
    }
    avg_FDP <- lapply(summary_FDP, function(FDP){
        apply(FDP, 1, function(x){mean(x, na.rm = TRUE)})
    })
    avg_power <- lapply(summary_power, function(power){
        apply(power, 1, function(x){mean(x, na.rm = TRUE)})
    })
    return(list(FDP = avg_FDP, power = avg_power))
}

####### Generate x
n <- 2500
x1 <- seq(-100, 100, length.out = 50)
x2 <- seq(-100, 100, length.out = 50)
x <- expand.grid(x1, x2)
colnames(x) <- c("x1", "x2")
pi_formula <- mu_formula <- "s(x1, x2)"
alphas <- seq(0.01, 0.3, 0.01)
output_filename <- paste0("../data/simul1_seed_", seed, ".RData")

## Case 1: a circle in the center
H0 <- apply(x, 1, function(coord){sum(coord^2) < 900})
mu <- ifelse(H0, 2, 0)

result1 <- simul1(x, mu, H0, 
                  pi_formula = pi_formula,
                  mu_formula = mu_formula,
                  alphas = alphas,
                  repeat_times = repeat_times)
result <- list(result1)
save(file = output_filename, result)

## Case 2: a circle in the corner
H0 <- apply(x, 1, function(coord){sum((coord - 65)^2) < 900})
mu <- ifelse(H0, 2, 0)

result2 <- simul1(x, mu, H0, 
                  pi_formula = pi_formula,
                  mu_formula = mu_formula,
                  alphas = alphas,
                  repeat_times = repeat_times)
result <- list(result1, result2)
save(file = output_filename, result)


## Case 3: a thin ellipsoid
shape_fun <- function(coord){
    transform_coord <- c(coord[1] + coord[2], coord[2] - coord[1])/sqrt(2)
    transform_coord[1]^2 / 100^2 + transform_coord[2]^2 / 15^2 < 1
}
H0 <- apply(x, 1, shape_fun)
mu <- ifelse(H0, 2, 0)

result3 <- simul1(x, mu, H0,
                  pi_formula = pi_formula,
                  mu_formula = mu_formula,
                  alphas = alphas,
                  repeat_times = repeat_times)

## Save data
result <- list(result1, result2, result3)
save(file = output_filename, result)
