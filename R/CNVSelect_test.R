#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Read inputs from csv file
#'
#' This function reads inputs from a file corresponding to duplicated sequences,
#' to be used in Genedupdip_neutralgenerator and mexpv_dipneut functions
#'
#'
#' @param input_file csv file path with frequencies and stuff
#' @param fasta_file alignment file in fasta format
#' @keywords phylogeny, CNV, neutral model
#' @export
#' @examples
#' CNVSelect_test(input_file, fasta_file)

CNVSelect_test <- function(input_file, fasta_file){
  output <- list()
  inputs <- read_DAF(input_file)
  dS <- get_dS(fasta_file)

  output$freqs <- inputs$freqs
  output$dS <- dS

  for(i in 1:length(dS)){
    t <- (100/3e-2)*dS[i]
    up <- 35*dS[i]

    out1 <- Genedupdip_neutralgenerator(inputs$N, up)
    out2 <- mexpv_dipneut(t=t, A=t(out1[[1]]), v=out1[[2]], N=inputs$N, Pos=out1[[3]])

    output$CIlower[i] <- out2[[2]][length(out2[[2]])]
    output$CIupper[i] <- out2[[3]][length(out2[[3]])]

    lower.tail <- out2[[5]][inputs$freqs[i]*2*inputs$N]
    if(lower.tail < 0.5){
      output$p_val[i] <- lower.tail * 2
    } else{
      output$p_val[i] <- (1-lower.tail) * 2
    }
  }
  return(output)
}
