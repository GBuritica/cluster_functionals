#######################################################################
#######################################################################
#######################################################################
#######################################################################
##### Gloria Buriticá
####  SS - cluster functionals - e.g. extremal index - c(1)
####  k vs. block length choice
#### 
#######################################################################
##### Asymptotic normality cluster estimator
#######################################################################
hill <- function(ordered,klim0){
    hill0 <- sapply(klim0, function(k) 
    mean( sapply(1:max(1,(k-1)), function(l) log(ordered[l]/ordered[k]) ) ) )
  return(hill0)
}
#### PATHS
## AX + B     with log A ~ N-0.5, B ~ Unif(0,1)
#path <- ARCHm(n)
#plot.ts(path)
#theta <-  0.2792
#ordered <- sort(path,decreasing=T)
#xi <-  hill(ordered, 2:floor(n^0.7) )
#xi2 <- alphaestimator(path, k1 = 2:floor(n^0.7))
#plot(2:floor(n^0.7),xi)
#points(2:floor(n^0.7),xi2$xi, col = 2)
#abline(h=1)
#alpha <- 1
ARCHm <- function(n0){
  x0  <- runif(1,0,1)
  nor <- rnorm(2*n0)
  uni <- runif(2*n0,0,1)
  for(i in 1:(2*n0) ) x0  <- c(x0, ( (exp(nor[i] - 0.5)*(x0[i])) + uni[i] ) )
  return(x0[(n0+1):(2*n0)])
}

## ARCHm Robert
## X = (eta + lambda X) Z^2, Z iid Gaussian.
#path <- ARCH2m(n)
#plot.ts(path)
#theta <-  0.727
#alpha <- 1
ARCH2m <- function(n0){
  eta <- 2*10^{-5}
  lambda <- .5
  nor   <- rnorm(2*n0)
  x0    <- eta
  for(i in 1:(2*n0) ) x0  <- c(x0, ( (eta+ (lambda*x0[i]) )*nor[i]^2 ) )
  return(x0[(n0+1):(2*n0)])
}

## AR model
#path <- ARm(n,0.7)
#plot.ts(path)
#theta <- 1-0.5^1
ARm  <- function(n0,par0=0.7){
  x0 <- arima.sim(n = n0, list(ar=par0, ma=0), rand.gen=function(n) abs( rt(n,df=1) ) ) 
  return(x0)
}

#path <- sARCHm(n,0.5)
#plot.ts(path)
#theta <- 0.727
sARCHm <- function(n0,lambda0){
  x0           <- vector(mode = "double" , 3*n0)
  norm         <- rnorm(3*n0)
  x0[1]        <- (2*10^(-5))
  for( i in 2:(3*n0)) x0[i] <-  ( (2*10^(-5)) + lambda0*x0[(i-1)] )*(norm[i])^2
  return(x0[(2*n0+1):(3*n0)])
}

#path <- ARCHm2(n,0.5)
#plot.ts(path)
#theta <- 0.835
ARCHm2 <- function(n0,lambda0){
  x0           <- vector(mode = "double" , 3*n0)
  norm           <- rnorm(3*n0,0,1)
  x0[1]        <- 0#2*10^{-5}   
  for( i in 2:(3*n0)) x0[i] <-  (norm[i])*sqrt( (2*10^(-5)) + lambda0*(x0[(i-1)])^2 )
  return(x0[(2*n0+1):(3*n0)]) 
}
#######################################################################
#######################################################################
#######################################################################
#######################################################################
### Extremal Index Cluser process function 
eiCP   <- function(path0,alpha0,klim0,n0){
  b          <- floor(sqrt(n0/(klim0)))
  eireturn   <- vector(mode="double", length=length(klim0) )
  eivariance <- vector(mode="double", length=length(klim0) )
  for(k in 1:length(klim0)){
    b0 <- b[k]
    ## computes paths
    sumaalpha <- sapply( 1:floor(n0/b0),
                         function(l) sum(path0[((l-1)*b0+1):(l*b0)]^alpha0) )
    
    maxaalpha <- sapply( 1:floor(n0/b0),
                         function(l) max(path0[((l-1)*b0+1):(l*b0)]^alpha0) )
    
    ordered   <- order(sumaalpha, decreasing=T) ## returns de order (1), (2), (3)...
    ind       <- ordered[1:klim0[k]]                   ## I chose the ones less than k
    estimate  <- mean( maxaalpha[ind]/sumaalpha[ind] ) ## among these I compute the mean
    variance  <- mean( maxaalpha[ind]^2/sumaalpha[ind]^2)
    ### moves b forward
    eireturn[k] <- estimate
    eivariance[k] <- variance
  }
  return( data.frame('b'=b,'k'=klim0, 'extremal_index'=eireturn,'ei_variance'=eivariance))
}


piCP   <- function(path0,alpha0,klim0,n0){
  b   <- floor(sqrt(n0/(klim0)))
  estimate <- NULL
  variance <- NULL
  for(k in 1:length(klim0)){
    b0 <- b[k]
    ## computes paths
    sumaalpha <- sapply( 1:floor(n0/b0),
                         function(l) sum(path0[((l-1)*b0+1):(l*b0)]^alpha0) )
    
    maxaalpha <- sapply( 1:floor(n0/b0), function(l)
      sort(path0[((l-1)*b0+1):(l*b0)]^alpha0, decreasing = T )[1:5])
    
    ordered       <- order(sumaalpha, decreasing=T) ## returns de order (1), (2), (3)...
    ind           <- ordered[1:klim0[k]]           ## I chose the ones less than k
    estimate      <- rbind(estimate,  c( mean( maxaalpha[1,ind]/sumaalpha[ind]),
                                         sapply(1:4, function(k)
                                           mean((maxaalpha[k,ind]-maxaalpha[(k+1),ind])/sumaalpha[ind] ) ) ) )  ## among these I compute the mean
    colnames(estimate)  <- c('Theta', 'Pi1', 'Pi2', 'Pi3', 'Pi4')
    estimate        <- as.data.frame(estimate)
    variance      <- rbind(variance, c( mean( maxaalpha[1,ind]^2/sumaalpha[ind]^2),
                                        sapply(1:4, function(k)
                                          mean( (maxaalpha[k,ind]-maxaalpha[(k+1),ind])^2/sumaalpha[ind]^2 )  )))
    colnames(variance)  <- c('Theta', 'Pi1', 'Pi2', 'Pi3', 'Pi4')
    variance          <- as.data.frame(variance)
  }
  return( list('estimate' = estimate, 'variance' = variance, 'k' = klim0, 'b' = b) )
}
