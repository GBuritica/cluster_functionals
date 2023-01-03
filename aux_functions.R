#######################################################################
#######################################################################
#######################################################################
#######################################################################
##### Gloria Buriticá
#### Auxiliar functions
####  
#### 
#######################################################################
#######################################################################
#######################################################################



clusterfunctional <- function(path,index,bn,p=alpha,knt,cfun){
  return( mean( sapply( index , function(k)  cfun(path[(((k-1)*bn) +1):(k*bn)],p,knt) ) ) )
}

cfun1 <- function(x0,p=alpha,knt){
  return( (max(x0)^p)/sum(x0^p))
}

cfun11 <- function(x0,p=alpha,knt){
  return(  sum( sapply( x0, function(k) (k > knt)*1 ) )  )
}

cfun11 <- function(x0,p=alpha,knt){
  return( sum(x0^p)/(max(x0)^p) )
}

cfun2 <- function(x0,p=1){
  return( sum(x0^alpha)/sum(x0)^alpha )
}

cfun22 <- function(x0,p=alpha){
  return( sum(x0)^alpha/sum(x0^alpha) )
}


