#---------------------------------------------------------------
# Generate all results (which can be done using a laptop)
# This includes Section 5.1 (except the gene/drug example with random orderings) and Section 5.3
#---------------------------------------------------------------

if (!requireNamespace("adaptMT", quietly = TRUE)){
    ## try http:// if https:// URLs are not supported
    devtools::install_github("lihualei71/adaptMT")
}
library("adaptMT")

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
for (i in 1:2){
    if (i == 1){
        dataset <- "estrogen_high"
        x <- estrogen$ord_high
    } else {
        dataset <- "estrogen_mod"
        x <- estrogen$ord_mod
    }
    pvals <- estrogen$pvals

    reorder <- order(x)
    x <- x[reorder]
    pvals <- pvals[reorder]
    x <- data.frame(x = x)

    result <- adapt_expr(x, pvals)
    filename <- paste0("../data/", dataset, "_res.RData")
    save(file = filename, result)
}

## RNA-seq data
for (dataset in c("bottomly", "airway", "pasilla")){
    data <- analyze_dataset(dataset)

    pvals <- data$pvalue
    ind <- which(!is.na(pvals))
    pvals <- pvals[ind]
    x <- log(data$baseMean + 1)
    x <- x[ind]

    reorder <- rev(order(x))
    x <- x[reorder]
    pvals <- pvals[reorder]
    x <- data.frame(x = x)

    result <- adapt_expr(x, pvals)
    filename <- paste0("../data/", dataset, "_res.RData")
    save(file = filename, result)
}

## Proteomics data
dataset <- "proteomics"

proteomics_file <- system.file(
    "extdata/real_data",
    "science_signaling.csv",
    package = "IHWpaper"
    )

proteomics_df <- read.csv(proteomics_file, stringsAsFactors = F)

pvals <- rank(
    proteomics_df$p1,
    ties.method="first"
    ) * proteomics_df$p1 / nrow(proteomics_df) 
x <- log(proteomics_df$X..peptides)

reorder <- rev(order(x))
x <- x[reorder]
pvals <- pvals[reorder]
x <- data.frame(x = x)

result <- adapt_expr(x, pvals)
filename <- paste0("../data/", dataset, "_res.RData")
save(file = filename, result)
