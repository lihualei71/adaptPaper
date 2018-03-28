#---------------------------------------------------------------
# Generate all plots in Section 5.
#---------------------------------------------------------------

if (!requireNamespace("adaptMT", quietly = TRUE)){
    devtools::install_github("lihualei71/adaptMT")
}
library("adaptMT")
library("splines")
source("expr_plot_fun.R")
library("latex2exp")

#### Real data examples
methods <- c('SeqStep', 'HingeExp', 'ForwardStop', 'Ada. SeqStep', 'BH', 'Storey-BH', 'Barber-Candes', 'SABHA (step)', 'SABHA (ordered)', 'IHW', 'IHW (oracle)', 'IF (oracle)', 'AdaPT')
## Green for non-adaptive ordered testing procedures: SeqStep, Accumulation Test, Forward Stop and Adaptive SeqStep; Blue for non-adaptive BH-type methods: BH, Storey, BC; Orange for IHW and IF; Black for SABHA; Red for AdaPT
cols <- c('green', 'green', 'green', 'green', 'blue', 'blue','blue', 'black', 'black', 'orange', 'orange', 'orange', 'red')
ltys <- c(4,3,2,1,1,2,3,2,1,3,2,1,1)
pchs <- c(4,3,2,1,2,3,1,1,4,2,3,5,1)


#### gene/drug response example (Section 5.1)
load("../data/estrogen_high_res.RData")
NumRej2 <- result$NumRej
obj2 <- result$adapt

load("../data/estrogen_mod_res.RData")
NumRej1 <- result$NumRej
obj1 <- result$adapt

load("../data/estrogen_random.RData")
NumRej0 <- NumRej

## load("../data/estrogen_random_ihw.RData")
## NumRej <- matrix(0, 13, 30)
## load("../data/estrogen_random_ihw.RData")
## for (i in 1:100){
##     NumRej <- NumRej + NumRej_list[[i]]
## }
## NumRej <- NumRej / 100
## NumRej0[10:11, ] <- NumRej[10:11, ]
## NumRej <- NumRej0
## save(file = "../data/estrogen_random.RData", NumRej)

## Figure 2
pdf("../figs/estrogen_rejs.pdf", width = 10, height = 4)
par(mfrow = c(1, 3), oma = c(0, 0, 3, 0), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
plot_results(NumRej0, methods, "Random ordering", cols, ltys, pchs, c(0, 1500), "# of Rejections", TRUE, 1.2)
plot_results(NumRej1, methods, "Moderately informative ordering", cols, ltys, pchs, c(0, 4000), "# of Rejections", FALSE)
plot_results(NumRej2, methods, "Highly informative Ordering", cols, ltys, pchs, c(0, 4000), "", FALSE)
mtext("Number of Rejections (Gene/Drug Response)",outer=TRUE,cex=1.8,font=2)
dev.off()

## Figure 3 (left)
pdf("../figs/estrogen_moderate_info_05.pdf")
plot_thresh_lfdr(obj1, 0.05, "alpha = 0.05",
                 xlab = "x = rank", 
                 disp_ymax = 0.1,
                 legend_pos = "topright")
dev.off()

## Figure 3 (right)
pdf("../figs/estrogen_moderate_info_10.pdf")
plot_thresh_lfdr(obj1, 0.1, "alpha = 0.1",
                 xlab = "x = rank",
                 disp_ymax = 0.1,
                 legend_pos = "topright")
dev.off()

## Figure 4 (left)
pdf("../figs/estrogen_high_info_05.pdf")
plot_thresh_lfdr(obj2, 0.05, "alpha = 0.05",
                 xlab = "x = rank",
                 disp_ymax = 0.3, 
                 legend_pos = "topright")
dev.off()

## Figure 4 (right)
pdf("../figs/estrogen_high_info_10.pdf")
plot_thresh_lfdr(obj2, 0.1, "alpha = 0.1",
                 xlab = "x = rank",
                 disp_ymax = 0.3,
                 legend_pos = "topright")
dev.off()

corr_estrogen_mod <- corr_lfdr(obj1)
corr_estrogen_high <- corr_lfdr(obj2)

corr_estrogen <- list(corr_estrogen_mod,
                      corr_estrogen_high)

## Figure 5
pdf("../figs/estrogen_corr.pdf", width = 7, height = 4.5)
plot_corr(corr_estrogen, "Correlation of Estimated Local FDR", ylim = c(0.9, 1), cols = c("red", "blue"), ltys = 1:2, legend = c("moderately informative ordering", "highly informative ordering"), cex.lab = 1.3, cex.legend = 1.3, cex.main = 1.6)
dev.off()

#### Bottomly experiment
load("../data/Bottomly_res.RData")
NumRej <- result$NumRej
obj <- result$adapt

## Figure 11, left
pdf("../figs/Bottomly_rejs.pdf")
par(mfrow = c(1, 1), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
plot_results(NumRej, methods, "Number of Rejections (Bottomly)", cols, ltys, pchs, c(0, 5300), "# of Rejections")
dev.off()

## Figure 11, right
pdf("../figs/Bottomly_info_10.pdf")
plot_thresh_lfdr(obj, 0.1, "Bottomly, alpha = 0.1",
                 xlab = "x = log(baseMean)",
                 xlim = c(log(2), max(obj$data[["x"]]$x)),
                 disp_ymax = 0.1, legend_pos = "topleft")
dev.off()

## Figure 11, middle
corr_Bottomly <- corr_lfdr(obj)

pdf("../figs/Bottomly_corr.pdf")
plot_corr(corr_Bottomly, "Correlation of Estimated Local FDR (Bottomly)", ylim = c(0.98, 1), cex.main = 1.7, cex.lab = 1.7)
dev.off()

#### airway experiment
load("../data/airway_res.RData")
NumRej <- result$NumRej
obj <- result$adapt

## Figure 12, left
pdf("../figs/airway_rejs.pdf")
par(mfrow = c(1, 1), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
plot_results(NumRej, methods, "Number of Rejections (airway)", cols, ltys, pchs, c(0, 12000), "# of Rejections")
dev.off()

## Figure 12, right
pdf("../figs/airway_info_10.pdf")
plot_thresh_lfdr(obj, 0.1, "airway, alpha = 0.1",
                 xlab = "x = log(baseMean)",
                 disp_ymax = 0.1, legend_pos = "topleft")
dev.off()

## Figure 12, middle
corr_airway <- corr_lfdr(obj)

pdf("../figs/airway_corr.pdf")
plot_corr(corr_airway, "Correlation of Estimated Local FDR (airway)", ylim = c(0.98, 1), cex.main = 1.7, cex.lab = 1.7)
dev.off()

#### pasilla experiment
load("../data/pasilla_res.RData")
NumRej <- result$NumRej
obj <- result$adapt

## Figure 13, left
pdf("../figs/pasilla_rejs.pdf")
par(mfrow = c(1, 1), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
plot_results(NumRej, methods, "Number of Rejections (pasilla)", cols, ltys, pchs, c(0, 2000), "# of Rejections")
dev.off()

## Figure 13, right
pdf("../figs/pasilla_info_10.pdf")
plot_thresh_lfdr(obj, 0.1, "pasilla, alpha = 0.1",
                 xlab = "x = log(baseMean)",
                 xlim = c(log(2), max(obj$data[["x"]]$x)),
                 disp_ymax = 0.1, legend_pos = "topleft")
dev.off()

## Figure 13, middle
corr_pasilla <- corr_lfdr(obj)

pdf("../figs/pasilla_corr.pdf")
plot_corr(corr_pasilla, "Correlation of Estimated Local FDR (pasilla)", ylim = c(0.98, 1), cex.main = 1.7, cex.lab = 1.7)
dev.off()

#### proteomics experiment
load("../data/proteomics_res.RData")
NumRej <- result$NumRej
obj <- result$adapt

## Figure 14, left
pdf("../figs/proteomics_rejs.pdf")
par(mfrow = c(1, 1), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
plot_results(NumRej, methods, "Number of Rejections (SILAC)", cols, ltys, pchs, c(0, 1200), "# of Rejections")
dev.off()

## Figure 14, right
pdf("../figs/proteomics_info_10.pdf")
plot_thresh_lfdr(obj, 0.1, "SILAC, alpha = 0.1",
                 xlab = "x = log(number of peptides)",
                 disp_ymax = 0.1, legend_pos = "topleft")
dev.off()

## Figure 14, middle
corr_proteomics <- corr_lfdr(obj)

pdf("../figs/proteomics_corr.pdf")
plot_corr(corr_proteomics, "Correlation of Estimated Local FDR (SILAC)", ylim = c(0.98, 1), cex.main = 1.7, cex.lab = 1.7)
dev.off()

#### Simulation 1 meta info
set.seed(1)
n <- 2500
x1 <- seq(-100, 100, length.out = 50)
x2 <- seq(-100, 100, length.out = 50)
x <- expand.grid(x1, x2)
colnames(x) <- c("x1", "x2")
pi.formula <- mu.formula <- "s(x1, x2)"
alpha.list <- seq(0.01, 0.3, 0.01)

## Figure 6
pdf("../figs/simul1_truth.pdf", width = 5, height = 1.75)
par(mfrow = c(1, 3), mar = c(2, 2, 2, 2))
n <- 2500
x1 <- seq(-100, 100, length.out = 50)
x2 <- seq(-100, 100, length.out = 50)
x <- expand.grid(x1, x2)
colnames(x) <- c("x1", "x2")
mains <- c("Circle in the middle",
           "Circle in the corner",
           "Thin ellipse")
for (i in 1:3){
    if (i == 1){
        H0 <- apply(x, 1, function(coord){sum(coord^2) < 900})
    } else if (i == 2){
        H0 <- apply(x, 1, function(coord){sum((coord - 65)^2) < 900})
    } else {
        shape.fun <- function(coord){
            transform.coord <- c(coord[1] + coord[2], coord[2] - coord[1])/sqrt(2)
            transform.coord[1]^2 / 100^2 + transform.coord[2]^2 / 15^2 < 1
        }
        H0 <- apply(x, 1, shape.fun)
    }
    plot_mask_2d(x, H0, cex = 0.5,
                 col_bg = "#A9A9A9", col_fg = "#000000",
                 main = mains[i], xaxt = "n", yaxt = "n")
    axis(side = 1, at = c(-100, 0, 100))
    axis(side = 2, at = c(-100, 0, 100))
}
dev.off()

## Figure 8
H0 <- apply(x, 1, function(coord){sum(coord^2) < 900})
mu <- ifelse(H0, 2, 0)
zvals <- rnorm(n) + mu
pvals <- 1 - pnorm(zvals)
library("mgcv")
formula <- "s(x1, x2)"

res <- adapt_gam(x, pvals, formula, formula)
res_corr <- corr_lfdr(res, niter_oracle = 1)

pdf("../figs/simul1_fit.pdf", width = 10, height = 3.7)
par(mfrow = c(1, 3))
plot_lfdr_2d(x, 1-res_corr$lfdr[, 15],
             main = "Estimated local FDR (alpha = 0.5)",
             cex = 1.5, cex.main = 1.7,
             xaxt = "n", yaxt = "n", pch = 20)
plot_lfdr_2d(x, 1-res_corr$lfdr[, 16],
             main = "Estimated local FDR (alpha = 0.3)",
             cex = 1.5, cex.main = 1.7,
             xaxt = "n", yaxt = "n", pch = 20)
plot_lfdr_2d(x, 1-res_corr$lfdr[, 17],
             main = "Estimated local FDR (alpha = 0.1)",
             cex = 1.5, cex.main = 1.7,
             xaxt = "n", yaxt = "n", pch = 20)
dev.off() 


## Plots for simulation 1
methods <- c('BH', 'Storey-BH', 'Barber-Candes', 'SABHA', 'IHW', 'AdaPT')
## Blue for non-adaptive BH-type methods: BH, Storey, BC; Orange for IHW; Black for SABHA; Red for AdaPT
cols <- c('blue', 'blue','blue', 'black', 'orange', 'red')
ltys <- c(1,2,3,2,3,1)
pchs <- c(2,3,1,1,2,1)   
inds <- c(1:5,7)

load("../data/simul1.RData")
titles <- c("Circle in the middle", "Circle in the corner",
            "Thin ellipse")
FDR.filename <- "../figs/simul1_FDR.pdf"
ylim <- c(0, 0.35)
pdf(FDR.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
    FDP <- result[[k]]$FDP[inds,]
    legend <- (k == 1)
    plot_results(FDP, methods, titles[k], cols, ltys, pchs,
                 ylim = ylim, ylab = "FDR",
                 legend = legend, cex.legend = 1.1)
}
dev.off()

power.filename <- "../figs/simul1_power.pdf"
pdf(power.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
    power <- result[[k]]$power[inds,]
    legend <- (k == 1)
    plot_results(power, methods, titles[k], cols, ltys, pchs,
                 ylim = c(0, 1.05), ylab = "power",
                 legend = legend)
}
dev.off()

#### Simulation 2 meta info
m <- 100
n <- 2000

inv_logit <- function(x){exp(x) / (1 + exp(x))}
x <- matrix(runif(n * m), n, m)
pi1 <- 0.3
beta_pi <- c(3, 3, rep(0, m-2))
beta0_pi <- uniroot(function(b){
    mean(inv_logit(x %*% beta_pi + b)) - pi1
}, c(-100, 100))$root
pi <- inv_logit(x %*% beta_pi + beta0_pi)

beta_mu <- c(2, 2, rep(0, m-2))
beta0_mu <- 0
mu <- 1 / pmax(1, x %*% beta_mu + beta0_mu)

pdf("../figs/simul2_hist.pdf", width = 10, height = 3.7)
par(mfrow = c(1, 2), cex = 1.5, mar = c(4, 5, 1, 3))
hist(pi, 20, xlab = TeX("$\\pi_{1i}$"), main = "")
hist(1 / mu, 20, xlab = TeX("$\\mu_{i}$"), main = "")
dev.off()

## Plots for simulation 2
methods <- c('BH', 'Storey-BH', 'Barber-Candes', 'AdaPT', 'AdaPT (oracle)')
## Blue for non-adaptive BH-type methods: BH, Storey, BC; Orange for IHW; Black for SABHA; Red for AdaPT
cols <- c('blue', 'blue','blue', 'red', 'red')
ltys <- c(1,2,3,2,1)
pchs <- c(2,3,1,2,1)   

load("../data/simul2.RData")

filename <- "../figs/simul2_FDR_power.pdf"
pdf(filename, width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(5, 5, 1, 3), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
FDP <- result$FDP
plot_results(FDP, methods, '', cols, ltys, pchs,
             ylim = c(0, 0.35), ylab = "FDR",
             cex.legend = 1.2)
power <- result$power
plot_results(power, methods, '', cols, ltys, pchs,
             ylim = c(0, 0.45), ylab = "power",
             legend = FALSE)
dev.off()
