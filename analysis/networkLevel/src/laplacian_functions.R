## keeps only spp with interactions
sppPresent<- function(M) M[rowSums(abs(M)) != 0, colSums(abs(M)) != 0]

## turns a adjacency matrix into a Laplacian matrix
toLap<-function(intMat){
  ## removing spp that do not interact
  intMat<-sppPresent(intMat)
  ## creating a square matrix
  lap.mat<- matrix(0,dim(intMat)[1]+dim(intMat)[2],dim(intMat)[1]+dim(intMat)[2]) 
  ## filling the top corner
  lap.mat[1:dim(intMat)[1], dim(intMat)[1]+1:dim(intMat)[2]]=-intMat
  ## filling the left low corner
  lap.mat[dim(intMat)[1]+1:dim(intMat)[2],1:dim(intMat)[1]]=-t(intMat)		
  ## fill the diagonals with the degrees
  diag(lap.mat)=-apply(lap.mat,1,sum)	
  return(lap.mat)
}

## calculates the number of compartments and the biggest eigenvalue of
## the Laplacian (which is the algebraic connectivity)
algCone<-function(int.mat){
  lapmat<-toLap(int.mat)
  ## getting the biggest non zero eigenvalue
  ## calculates all eigenvalues
  eigens=eigen(lapmat, only.values=TRUE)$values## gets rid of the
  ## number of zeroes and returns the number of components (= the
  ## number of zeroes)
  comps=round(sum(eigens<=0.00001),1) ## conta numero de autovalores
  ## iguais a zero
  re.eigen<-Re(round(eigens, 4)) ## geting the real parts and rounding
  alg.conn<-max(re.eigen)##  Algebraic connectivity
  return(list=c(comps,alg.conn))
}
