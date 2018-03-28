#---------------------------------------------------------------
# Test whether the functions work
#---------------------------------------------------------------

if (!requireNamespace("IHWpaper", quietly = TRUE)){
    ## try http:// if https:// URLs are not supported
    source("https://bioconductor.org/biocLite.R")
    biocLite("IHWpaper")
}
library("IHWpaper")    
source("other_methods.R")
source("expr_template.R")
set.seed(1)

## Microarray data
data(estrogen)
dataset <- "estrogen_high"
x <- estrogen$ord_high
pvals <- estrogen$pvals

reorder <- order(x)
x <- x[reorder]
pvals <- pvals[reorder]
x <- data.frame(x = x)

x <- x[1:1000,,drop=FALSE]
pvals <- pvals[1:1000]

result <- adapt_expr(x, pvals)

