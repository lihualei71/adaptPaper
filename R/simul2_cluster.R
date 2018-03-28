source("../helpers/other_methods.R")
source("../helpers/summarize_methods.R")
if (!requireNamespace("adaptMT", quietly = TRUE)){
    devtools::install_github("lihualei71/adaptMT")
}
library("adaptMT")
library("mgcv")

repeat_times <- as.numeric(Sys.getenv("times"))
seed <- as.numeric(Sys.getenv("seed"))
set.seed(seed)

simul2 <- function(x, mu, pi,
                   alphas,
                   repeat_times){
    n <- length(mu)
    m <- length(alphas)
    summary_FDP <- list()
    summary_power <- list()
    for (j in 1:5){
        summary_FDP[[j]] <- matrix(rep(0, repeat_times * m),
                                   ncol = repeat_times)
        summary_power[[j]] <- matrix(rep(0, repeat_times * m),
                                     ncol = repeat_times)
    }
    xx <- data.frame(x[, 1:2])
    names(xx) <- c("x1", "x2")
    
    for (i in 1:repeat_times){
        H0 <- as.logical(ifelse(runif(n) < pi, 1, 0))
        y <- ifelse(H0, rexp(n, 1/mu), rexp(n, 1))
        pvals <- exp(-y)
        
        BH_result <- summary_BH(pvals, H0, alphas)
        summary_FDP[[1]][, i] <- BH_result[, 2]
        summary_power[[1]][, i] <- BH_result[, 3]

        storey_result <- summary_storey(pvals, H0, thr = 0.5,
                                        alphas)
        summary_FDP[[2]][, i] <- storey_result[, 2]
        summary_power[[2]][, i] <- storey_result[, 3]

        BC_result <- summary_BC(pvals, H0, alphas)
        summary_FDP[[3]][, i] <- BC_result[, 2]
        summary_power[[3]][, i] <- BC_result[, 3]

        res <- try(adapt_glmnet(x, pvals, alphas = alphas))
        if (class(res) != "try-error"){
            adapt_result <- summary_adapt(res, H0, pvals)
            summary_FDP[[4]][, i] <- adapt_result[, 2]
            summary_power[[4]][, i] <- adapt_result[, 3]
        }

        res_oracle <- try(
            adapt_glm(xx, pvals, 
                      pi_formulas = "x1 + x2",
                      mu_formulas = "x1 + x2",
                      alphas = alphas))
        if (class(res_oracle) != "try-error"){
            adapt_oracle_result <- summary_adapt(res_oracle, H0, pvals)
            summary_FDP[[5]][, i] <- adapt_oracle_result[, 2]
            summary_power[[5]][, i] <- adapt_oracle_result[, 3]
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

#### Simulation 2
m <- 100
n <- 2000
alphas <- seq(0.01, 0.3, 0.01)
output_filename <- paste0("../../data/simul2_seed_", seed, ".RData")

x <- matrix(runif(n * m), n, m)
pi1 <- 0.3

inv_logit <- function(x) {exp(x) / (1 + exp(x))}
beta_pi <- c(3, 3, rep(0, m-2))
beta0_pi <- uniroot(function(b){
    mean(inv_logit(x %*% beta_pi + b)) - pi1
}, c(-100, 100))$root
pi <- inv_logit(x %*% beta_pi + beta0_pi)
beta_mu <- c(2, 2, rep(0, m-2))
beta0_mu <- 0
mu <- pmax(1, x %*% beta_mu + beta0_mu)

result <- simul2(x, mu, pi, alphas, repeat_times)

save(file = output_filename, result)
