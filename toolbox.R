#R-code of "Structural changes within trophic levels are constrained by within-family assembly rules at lower trophic levels" by:
#Chuliang Song, Florian Altermatt, Ian Pearse, Serguei Saavedra
#published in: Ecology Letters


#function that computes the structural stability from the herbivore interaction matrix
#inputs: alpha = herbivore interaction matrix at given time
#output: out = structural stability of community persistence at given time
library(mvtnorm) 
Omega <- function(alpha){
  S <- nrow(alpha)
  omega <- function(S, Sigma){
    m <- matrix(0,S,1)
    a <- matrix(0,S,1)
    b <- matrix(Inf,S,1)  
    d <- pmvnorm(lower = rep(0,S), upper = rep(Inf,S), mean = rep(0,S), sigma = Sigma)
    out <- 2*(2*d[1])^(1/(S-1))
    return(out)
  }
  if(length(which(diag(alpha)==0))==0){
    Sigma <- chol2inv(alpha, size = NCOL(alpha), LINPACK = FALSE)
    return(omega(S,Sigma))
  }
  else{
    f <- function(m) class(try(solve(t(m)%*%m),silent=T))=="matrix"
    if(f(alpha)==FALSE){
      return(0)
    }
    else{
      Sigma <- solve(t(alpha) %*% alpha)
      return(omega(S,Sigma))
    }
  }
}

#function that infers the herbivore interaction matrixfrom the (binary) herbivore-plant interaction matrix
#inputs: Y = herbivore-plant interaction matrix at given time
#output: beta = herbivore interaction matrix at given time
get_interaction_matrix <- function(A){
  gram <- t(A) %*% A
  for(i in 1:ncol(gram)){
    gram[,i] <- gram[,i]/sum(gram[,i])
  }
  diag(gram) <- 1
  return(gram)
}

extract_binary_Inte <- function(present, matrix){
  ### extracting the corresponding network M_present
  M_present <- matrix[present$species,]
  w_interactions <- which(colSums(M_present)>0)
  M_present <- M_present[,w_interactions] 
  M_present
}

sample_all <- function(a){
  if(length(a)==1)
    a
  else
    sample(a)
}