source("other_methods.R")
source("expr_template.R")

output_filename <- paste0("../data/estrogen_random_ihw.RData")

NumRej_list <- list()
data(estrogen)
pvals <- estrogen$pvals
n <- length(pvals)
x <- data.frame(x = 1:n)

## for (seed in 0:49){
##     set.seed(seed)
##     print("################################################")
##     print(paste0("seed: ", seed))
##     print("################################################")    
##     for (i in 1:2){
##         ind <- sample(n)
##         pvals <- pvals[ind]

##         res <- adapt_expr(x, pvals)
##         NumRej_list[[2 * seed + i]] <- res$NumRej
##         save(file = output_filename, NumRej_list)
##     }
## }

NumRej <- matrix(0, 13, 30)
load("../data/estrogen_random_ihw.RData")
for (i in 1:100){
    NumRej <- NumRej + NumRej_list[[i]]
}
NumRej <- NumRej / 100
