#######################################################################
#######################################################################
#######################################################################
#######################################################################
##### Gloria Buriticá
####  SS - cluster functionals - e.g. extremal index - c(1)
####  k vs. block length choice
#### 
#######################################################################
#######################################################################
source("IndexofRV.R")
#######################################################################
#######################################################################
## Parameters
n         <- 5000
alpha0    <- 1.3
phi       <- 0.6 
es       <- NULL
esmax    <- NULL
kchoice  <- NULL

bn0    <-   sapply(2:6, function(k) 2^k)
knmax  <-   floor( n*0.04) 
thet   <-   1-phi^alpha#0.2792 # 
estimate    <- matrix(NA, nrow = length(bn0), ncol = knmax )
estimate2   <- matrix(NA, nrow=length(bn0), ncol = knmax )
estimate3   <- matrix(NA, nrow = length(bn0), ncol = knmax )


for(N in 1:2){
  print(N)
  path    <-   abs(arima.sim(n=n, list(ar=phi, ma=0), rand.gen=function(n) rt(n,df=alpha0) ))#abs(ARCHm(n))#
  rand_path <- sample(path) 
  
  alpha  <- (1/alphaestimator(path,k1=floor(n^0.8))$xi);print(alpha)#
  p      <- alpha 
  
  for(j in  1:length(bn0) ){
    bn  <- bn0[j]
    mn <- min(floor(n/bn),knmax)

    samplepsum  <- sapply( 1:mn, function(k)   sum( path[( ((k-1)*bn) +1):(k*bn)]^p ) )  ## |X_[1,bn]|_p^p sample.
    samplepmax  <- sapply( 1:mn, function(k)   max( path[( ((k-1)*bn) +1):(k*bn)] ) ) 
    samplepsum_rand  <- sapply( 1:mn, function(k) sum( rand_path[( ((k-1)*bn) +1):(k*bn)]^p ) )
    
    sortsum        <- sort(samplepsum,  partial = 1:mn )       ## places the kth largest at the n-kn+1 place 
    sortmax        <- sort(samplepmax,  partial = 1:mn )
    sortsum_rand   <- sort(samplepsum_rand,  partial = 1:mn )
    
    
    for(kn  in 1:(mn-1) ){
      
      kth    <- sortsum[(mn-kn)]       ## places the kth largest at the n-kn+1 place 
      kth2   <- sortsum_rand[(mn-kn)]
      kth3   <- sortmax[(mn-kn)]
        
      index <-  which( samplepsum > kth )                               ## gives  index. of the k largest. 
      index2<-  which( samplepsum_rand > kth2 ) 
      index3<-  which( samplepmax > kth3 ) 
      
      ## Computes estimates
      estimate[j,kn] <-  clusterfunctional(path,index,bn,p,kth,cfun1)              ## computes kn estimates
      estimate2[j,kn] <- clusterfunctional(rand_path,index2,bn,p,kth2,cfun1)              ## computes kn estimates
      estimate3[j,kn] <- 1/clusterfunctional(path,index3,bn,p,kth3,cfun11)
      
      
      
    }
    
  }
  
  
  est <- estimate2
  for(k in 1:(mn-1)) est[,k] <- abs(estimate2[ ,k] -1)
  
  #plot.ts(est[1,], ylim = c(0,5))
  #for(j in 1:length(bn0) ) lines(est[j,], col = j)
  
  
  k0 <- sapply(1:length(bn0), function(j) which.min(est[j, ]  ) ) ;# k0
  kchoice <- rbind(kchoice, k0)
  es      <-   rbind( es, sapply(1:length(bn0), function(j) estimate[j,k0[j]] ) )   
  esmax   <- rbind( esmax, sapply(1:length(bn0), function(j) estimate3[j,k0[j]] ) )   
  
}



par(mfrow=c(1,1))
estim <- NULL
for(j in 1:length(bn0)) estim <- cbind(estim, es[,j],esmax[ ,j] )
boxplot(estim)
abline(h=1-phi^alpha0)

par(mfrow = c(1,2))
boxplot(es, ylim = c(0,1))
abline(h=1-phi^alpha0)
boxplot(esmax, ylim = c(0,1))
abline(h=1-phi^alpha0)

summary(kchoice)
boxplot(kchoice)

