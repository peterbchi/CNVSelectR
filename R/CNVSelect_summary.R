#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Function to write summary table of p-values and confidence intervals
#'
#' This function takes output from CNVSelect_test and writes a summary table
#'
#'
#' @param test_out object that contains output from CNVSelect_test
#' @param make_kable make kable object for Markdown file if TRUE; otherwise return matrix
#' @param sig_dig number of significant figures to include in summary table
#' @keywords phylogeny, CNV, neutral model
#' @export
#' @examples
#' CNVSelect_summary(test_out)

CNVSelect_summary <- function(test_out, make_kable=TRUE, sig_dig=3){
  # get number of pairs being tested
  n.pairs <- length(test_out$freqs)

  # start output matrix
  output <- matrix(NA, ncol=4, nrow=n.pairs)
  rownames(output) <- colnames(test_out$dS)

  # fill in easy numeric values
  output[,1] <- signif(test_out$dS, digits=sig_dig)
  output[,2] <- test_out$freqs
  output[,4] <- signif(test_out$p_val, digits=sig_dig)

  # format and fill in confidence intervals
  for(i in 1:n.pairs){
    output[i,3] <- paste("(", test_out$crit_lower[i], ", ", test_out$crit_upper[i], ")", sep="")
  }

  colnames(output) <- c("dS", "frequency", "0.05-level Crit. Vals", "p-value")

  # make a kable or return a matrix
  if(make_kable){
    kable(output)
  } else {
    return(output)
  }

}
