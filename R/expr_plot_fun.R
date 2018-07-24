#---------------------------------------------------------------
# Helpers: functions to generate plots in Section 5
#---------------------------------------------------------------

plot_results <- function(vals, methods, title,
                         cols, ltys, pchs,
                         ylim, ylab,
                         legend = TRUE,
                         cex.legend = 1.2){
### Plot results for each method
    alphalist <- seq(0.01, 0.3, 0.01)
    plot(0:1, 0:1, type = 'n',
         xlim = range(alphalist), ylim = ylim,
         xlab = expression(paste('Target FDR level ',alpha)),
         ylab = ylab,
         main = title, axes = FALSE)
    axis(side = 1, at = c(0, 0.1, 0.2, 0.3))
    axis(side = 2)
    alpha_pt = c(10, 20, 30)
    for (i in 1:length(methods)){
        points(alphalist, vals[i, ],
               type = 'l', col = cols[i], lty = ltys[i])
        points(alphalist[alpha_pt], vals[i, alpha_pt],
               col = cols[i], pch = pchs[i])
    }
    if (legend){
        legend("topleft", methods,
               col = cols, lty = ltys, pch = pchs,
               seg.len = 3, cex = cex.legend, bty = "n")
    }
}

plot_corr <- function(objs, title,
                      ylim = c(0.7, 1),
                      cols = 1,
                      ltys = 1,
                      legend = NULL,
                      legend.position = "bottomright",
                      cex.legend = 1.3,
                      ...){
    par(mfrow = c(1, 1), ...)
    plot(0:1, 0:1, xlim = c(0.5, 0.01), ylim = ylim, type = 'n',
         xlab = expression(paste("Target FDR ", alpha, " (large --> small)")),
         ylab = 'Correlation', main = title, axes=FALSE)
    axis(side = 1, at = c(0.5, 0.4, 0.3, 0.2, 0.1, 0))
    axis(side = 2)
    if (is.null(names(objs))){
        for (i in 1:length(objs)){
            obj <- objs[[i]]
            points(obj$alphas, obj$corr, type = 'l',
                   lwd = 2, col = cols[i], lty = ltys[i])
        }
        if (!is.null(legend)){
            legend(legend.position, legend,
                   col = cols, lty = ltys,
                   bty = "n", cex = cex.legend)
        }
    } else {
        points(objs$alphas, objs$corr, type = 'l',
               lwd = 2, col = cols, lty = ltys)
    }
}

plot_thresh_lfdr <- function(obj, x, pvals,
                             alpha, title, xlab, legend_pos,
                             disp_ymax = 0.2, xlim = NULL){
    if (!requireNamespace("adaptMT", quietly = TRUE)){
        stop("\'adaptMT\' package not found. Please install.")
    }

    par(mfrow = c(2, 1), mar=c(4.1, 4.1, 2, 0.15),
        cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
    title_thresh <- paste0("Rejection threshold, (", title, ")")
    adaptMT::plot_1d_thresh(obj, x, pvals,
                            alpha, title_thresh,
                            xlab = xlab, xlim = xlim,
                            disp_ymax = disp_ymax)
    title_lfdr <- paste0("Estimated local FDR, (", title, ")")
    adaptMT::plot_1d_lfdr(obj, x, pvals,
                          alpha, title_lfdr,
                          xlab = xlab, xlim = xlim,
                          disp_ymax = disp_ymax,
                          legend_pos = legend_pos)
}

plot_mask_2d <- function(x, mask, main, cex,
                         col_bg = "#FFB6C1",
                         col_fg = "#800000",
                         ...){
    par(...)
    color <- ifelse(mask, col_fg, col_bg)
    plot(x[, 1], x[, 2], type = "n", xlab = "", ylab = "",
         main = main)
    points(x[, 1], x[, 2], col = color, cex = cex)
}

plot_lfdr_2d <- function(x, lfdr, main, cex,
                         col_fg = "#000080",
                         col_bg = "#ADD8E6",
                         ...){
    par(...)
    rbPal <- colorRampPalette(c(col_bg, col_fg))
    color <- rbPal(5)[as.numeric(cut(lfdr, breaks=4))]
    plot(x[, 1], x[, 2], type = "n", xlab = "", ylab = "",
         main = main)
    points(x[, 1], x[, 2], col = color, cex = cex)
}
