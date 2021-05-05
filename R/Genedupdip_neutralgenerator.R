#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Neutral model generator
#'
#' This function creates the infinitesimal generator for the model
#'
#'
#' @param N population size
#' @param up Poisson rate of pseudogenization
#' @keywords phylogeny, CNV, neutral model
#' @export
#' @examples
#' Genedupdip_neutralgenerator(N, up)

Genedupdip_neutralgenerator <- function(N,up){
  #Generator for the model without neofunctionalization
  Pos <- matrix(0, ncol=N+1, nrow=N+1)
  count <- 0


  #Create indexing matrix Pos such that Q(Pos(i,j),:) is row of Q
  #corresponding to state (i,j), etc.
  for (i in 1:(N+1)){
    for (j in 1:(N+2-i)){
      #if ~(i==1 && j ==1) #&& ~(i== 1 && j == N+1)
      count <- count + 1
      Pos[i,j] <- count
    }
  }


  # Evaluate number of non-zero entries of Q
  nonzerolength <- 7*(N-2)*(N-1)/2 + (N-1)*15 + 5


  #Declare index vectors ii, jj and rate vector vv for construction of sparse
  #Q matrix (ii(n) = i, jj(n) = j, vv(n) = q_ij <-> Q(i,j) = q_ij)
  ii <- rep(0, nonzerolength)
  jj <- rep(0, nonzerolength)
  vv <- rep(0, nonzerolength)

  #i is number of AAAA
  #j is number of AAA-
  #k = N-i-j is number of AA--

  #First consider 'middle transitions' when nothing is 0
  count1 <- 1
  for (i in 2:N){#1:N-1, +1 for indexing
    for (j in 2:(N-i+1)){#similar
      if((N-i+1) > 1){
        k <- N-(i-1)-(j-1)#k not indexing so no need to modify
        # changed pbi, pbj, pbk to a single pb vector
        pb <- pbirth(i-1,j-1,k,N)

        pdi <- (i-1)/N
        pdj <- (j-1)/N
        pdk <- k/N

        ii[count1:(count1+6)] <- Pos[i,j]
        jj[count1:(count1+6)] <- c(Pos[i+1,j],Pos[i+1,j-1],Pos[i-1,j+1],Pos[i,j+1],Pos[i-1,j],Pos[i,j-1],Pos[i,j])
        vv[count1]   <- pdk*pb[1]
        vv[count1+1] <- pdj*pb[1]
        vv[count1+2] <- pdi*pb[2]+2*(i-1)*up
        vv[count1+3] <- pdk*pb[2]
        vv[count1+4] <- pdi*pb[3]
        vv[count1+5] <- pdj*pb[3]+(j-1)*up
        vv[count1+6] <- -sum(vv[count1:(count1+5)])
        count1 <- count1+7
      }
    }
    #Now transitions where k = 0
    j <- N-i+2
    k <- N-(i-1)-(j-1)#=0
    if (k != 0){
      stop("k should be 0 but isn't")
    }
    pb <- pbirth(i-1,j-1,k,N)
    pdi <- (i-1)/N
    pdj <- (j-1)/N
    pdk <- k/N

    ii[count1:(count1+4)] <- Pos[i,j]
    jj[count1:(count1+4)] <- c(Pos[i+1,j-1],Pos[i-1,j+1],Pos[i,j-1],Pos[i-1,j],Pos[i,j])
    vv[count1] <- pdj*pb[1]
    vv[count1+1] <- pdi*pb[2]+2*(i-1)*up
    vv[count1+2] <- pdj*pb[3]+(j-1)*up
    vv[count1+3] <- pdi*pb[3]
    vv[count1+4] <- -sum(vv[count1:(count1+3)])
    count1 <- count1+5
  }
  #Now transitions where i = 0
  i <- 1#indexing
  for (j in 2:N){
    k <- N-(i-1)-(j-1)
    pb <- pbirth(i-1,j-1,k,N)

    pdi <- (i-1)/N
    pdj <- (j-1)/N
    pdk <- k/N

    ii[count1:(count1+4)] <- Pos[i,j]
    jj[count1:(count1+4)] <- c(Pos[i,j+1],Pos[i,j-1],Pos[i+1,j-1],Pos[i+1,j],Pos[i,j])
    vv[count1] <- pdk*pb[2]
    vv[count1+1] <- pdj*pb[3]+(j-1)*up
    vv[count1+2] <- pdj*pb[1]
    vv[count1+3] <- pdk*pb[1]
    vv[count1+4] <- -sum(vv[count1:(count1+3)])
    count1 <- count1 + 5
  }
  #Now transitions where i = N;
  i <- N+1
  j <- 1
  k <- 0
  #only pseudogenization here
  ii[count1:(count1+1)] <- Pos[i,j]
  jj[count1:(count1+1)] <- c(Pos[i-1,j+1],Pos[i,j])
  vv[count1] <- 2*(i-1)*up
  vv[count1+1] <- -vv[count1]
  count1 <- count1+2

  #Now transitions where j = 0;
  j <- 1
  for (i in 2:N){
    k <- N-(i-1)-(j-1)
    pb <- pbirth(i-1,j-1,k,N)
    pdi <- (i-1)/N
    pdj <- (j-1)/N
    pdk <- k/N

    ii[count1:(count1+4)] <- Pos[i,j]
    jj[count1:(count1+4)] <- c(Pos[i+1,j],Pos[i-1,j],Pos[i,j+1],Pos[i-1,j+1],Pos[i,j])
    vv[count1] <- pdk*pb[1]
    vv[count1+1] <- pdi*pb[3]
    vv[count1+2] <- pdk*pb[2]
    vv[count1+3] <- pdi*pb[2]+2*(i-1)*up
    vv[count1+4] <- -sum(vv[count1:(count1+3)])
    count1 <- count1+5
  }

  #Now transitions where j = N;
  j <- N+1
  i <- 1
  k <- 0
  pb <- pbirth(i-1,j-1,k,N)#j's can give birth to both other types so we need these, obviously pdj = 1 so that's omitted.
  ii[count1:(count1+2)] <- Pos[i,j]
  jj[count1:(count1+2)] <- c(Pos[i+1,j-1],Pos[i,j-1],Pos[i,j])
  vv[count1] <- pb[1]
  vv[count1+1] <- pb[3] +(j-1)*up
  vv[count1+2] <- -vv[count1]-vv[count1+1]
  count1 <- count1+2

  #Finished getting rates

  # Q is the generator matrix
  Q <- sparseMatrix(i=ii,j=jj,x=vv,dims=c(count,count), symmetric=FALSE)
  #e1 is initial distribution
  e1 <- sparseMatrix(i=1,j=Pos[2,1],x=1,dims=c(1,count), symmetric=FALSE)
  return(list(Q, e1, Pos))
#  return(list(ii, jj, vv, count))

}

pbirth <- function(i,j,k,N){
    #This function gives the probability of birth of type i, j, k
    #individuals given state of population
    pbi <- (i/N)*(i-1)/(N-1) + (1/4)*(j/N)*(j-1)/(N-1) + (j/N)*i/(N-1)
    pbj <- (j/N)*(k/(N-1)) + (j/N)*(i/(N-1)) + (1/2)*(j/N)*((j-1)/(N-1)) +2*(i/N)*(k/(N-1))
    pbk <- (k/N)*((k-1)/(N-1)) + (j/N)*(k/(N-1)) + (1/4)*(j/N)*((j-1)/(N-1))
    # eps <- 1e-15 I don't know wtf this is for
    return(c(pbi, pbj, pbk))
}

