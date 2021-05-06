#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Read inputs from csv file
#'
#' This function reads inputs from a file corresponding to duplicated sequences
#'
#'
#' @param file file path
#' @keywords phylogeny, CNV, neutral model
#' @export
#' @examples
#' read_DAF(file)

read_DAF <- function(file){
  DAF <- read_csv(file, col_names = FALSE)

  first.col <- dim(DAF)[1]

  output <- list()
  output$N <- as.numeric(DAF[1,2])   # N
  output$ploidy <- as.numeric(DAF[2,2])   # ploidy
  output$freqs <- as.numeric(DAF[5:first.col,1]$X1)

  return(output)
}
