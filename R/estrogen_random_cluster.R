#---------------------------------------------------------------
# Generate the results for the gene/drug example with random orderings) in Section 5.1.
# Need to be ran on clusters. See README.
#---------------------------------------------------------------

source("../helpers/other_methods.R")
source("expr_template.R")

repeat_times <- as.numeric(Sys.getenv("times"))
seed <- as.numeric(Sys.getenv("seed"))
output_filename <- paste0("../data/estrogen_random_", seed, ".RData")
set.seed(seed)

NumRej_list <- list()
data(estrogen)
pvals <- estrogen$pvals
n <- length(pvals)
x <- data.frame(x = 1:n)

for (i in 1:repeat_times){
    ind <- sample(n)
    pvals <- pvals[ind]

    res <- adapt_expr(x, pvals)
    NumRej_list[[i]] <- res$NumRej
    save(file = output_filename, NumRej_list)
}
