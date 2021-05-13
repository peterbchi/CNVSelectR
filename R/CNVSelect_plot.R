#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Function to create summary plot of data and confidence intervals
#'
#' This function takes output from CNVSelect_test and creates a summary plot
#'
#'
#' @param test_out object that contains output from CNVSelect_test
#' @keywords phylogeny, CNV, neutral model
#' @export
#' @examples
#' CNVSelect_summary(test_out)

CNVSelect_plot <- function(test_out){
  y_max <- max(c(test_out$freqs, test_out$CIupper)) * 1.3
  x_max <- max(test_out$dS) * 1.3
  plot(test_out$freqs ~ as.vector(test_out$dS),
       xlim=c(0,x_max),
       ylim=c(0,y_max),
       xlab="dS",
       ylab="frequency",
       main="Critical Value cutoffs and data points",
       pch=19)

  # add confidence intervals
  n.seq <- length(test_out$freqs)
  CI_lw <- x_max/50

  for(i in 1:n.seq){
    lines(x=c(test_out$dS[i], test_out$dS[i]), y=c(test_out$crit_lower[i], test_out$crit_upper[i]))
    lines(x=c(test_out$dS[i]-CI_lw, test_out$dS[i]+CI_lw), y=c(test_out$crit_lower[i], test_out$crit_lower[i]))
    lines(x=c(test_out$dS[i]-CI_lw, test_out$dS[i]+CI_lw), y=c(test_out$crit_upper[i], test_out$crit_upper[i]))
  }
}
