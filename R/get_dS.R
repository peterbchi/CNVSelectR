#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Read inputs from fasta file
#'
#' This function reads duplicated sequences from a fasta file and
#' calculates pairwise dS.
#'
#'
#' @param fasta_file file path
#' @keywords phylogeny, CNV, neutral model
#' @export
#' @examples
#' get_dS(fasta_file)

get_dS <- function(fasta_file){
  x <- read.alignment(fasta_file,format="fasta")
  n.seq <- length(x$seq)

  # Check if all sequences are of lengths that are multiples of 3
  # and also for internal stop codons
  for(i in 1:n.seq){
    tmp <- unlist(strsplit(x$seq[[i]], ""))
    if(length(tmp)/3 != round(length(tmp)/3, 0)){
      stop("All sequences must be of lengths that are a multiple of 3")
    }
    codon.starts <- seq(1, length(tmp), by=3)
    stop.codons <- c("taa", "tag", "tga", "TAA", "TAG", "TGA")
    for(j in codon.starts){
      codon <- paste(tmp[c(j, j+1, j+2)], collapse="")
      if(is.element(codon, stop.codons)){
        stop("Sequences must not have internal stop codons")
      }
    }
  }

  dS <- NULL
  col_names <- NULL
  if(n.seq > 2){
    seqs <- seq(1,n.seq, by=2)
    for(i in seqs){
      pair <- x
      pair$seq <- pair$seq[c(i, (i+1))]
      pair$nam <- pair$nam[c(i, (i+1))]
      pair$nb <- 2
      dS <- c(dS, kaks(pair)$ks)
      col_names <- c(col_names, paste(pair$nam[1], pair$nam[2], sep = "/"))
    }
  }
  dS <- matrix(dS, nrow=1)
  colnames(dS) <- col_names
  return(dS)
}
