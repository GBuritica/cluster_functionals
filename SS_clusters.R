#######################################################################
#######################################################################
#######################################################################
#######################################################################
##### Gloria Buriticá
#### SS - Extremal index - c(1)
####  
#### 
#### Boxplot of estimates of the cluster index of sums and the extremal index
#### Comparison: alpha-cluster-based vs. maxima-cluster-based estimator
#######################################################################
source("IndexofRV.R")
source("Plots_code.R")
library(ggplot2)
library(latex2exp)
library(gridExtra)


### Parameters
n      <- 20000
thr2   <- thr <- c(2)
bu     <- sapply(1:6, function(k) 2^k)
phi0    <- c(0.8,0.5); 
infor1 <- data.frame( "DATA" = NULL, "TYPE"=NULL, "zn" =NULL, "BL" = NULL, "TH" = NULL, "MODE" = NULL )
infor2 <- data.frame( "DATA" = NULL, "TYPE"=NULL, "zn" =NULL, "BL" = NULL, "TH" = NULL, "MODE" = NULL )
alpha0 <- 0.7#1.3#
thet <-   sapply(1:length(phi0), function(k) 1-phi0[k]^alpha0)
c1   <-   1/sapply(1:length(phi0), function(k) (1-phi0[k])^alpha0/thet[k])
title<-   TeX('$\\widehat{c}(1)$')# TeX('$\\widehat{\\theta}_{| \\mathbf{X}|}$') #
## functions theta 

for(N in 1:10){
  print(N)
  for(j in 1:length(phi0)){
    phi    <- phi0[j]
    path   <- abs(arima.sim(n=n, list(ar=phi, ma=0), rand.gen=function(n) rt(n,df=alpha0) ))
    n      <- length(path)
    alpha  <- (1/alphaestimator(path,k1=floor(n^0.8))$xi)#;print(alpha)
    #alpha <- alpha0
    sorted <- sort(path,decreasing = T)
    p      <- 1#alpha#
    for(i in 1:length(bu)){
      
      mn  <- floor(n/bu[i])
      
      maxi  <- sapply( 1:mn, function(k) max(   path[((k-1)*bu[i] + 1):(k*bu[i])] )  )
      sumaalpha <- sapply( 1:mn, function(k) sum(   path[((k-1)*bu[i] + 1):(k*bu[i])]^alpha )  )
      
      #Sorted Paths
      n2 <- mn
      #sortedmaxi <- sort(maxi, decreasing = TRUE); #print(n2)
      #sortedsumaalpha <- sort(sumaalpha, decreasing=T)
      
      #Thresholds
      kn      <-  max(1,floor( (n)/bu[i]^thr2[1] ) )+1
      th3     <- sort(maxi, partial = mn - kn +1 )[(mn-kn+1)]   #sortedmaxi[ max(1,floor( (n)/bu[i]^thr2[1])) +1 ]                ## supremom norm th.
      th4     <- sort(sumaalpha, partial = mn - kn +1 )[(mn-kn+1)] #sortedsumaalpha[ max(1,floor( (n)/bu[i]^thr2[1]))+1 ] 
      
      #Index
      ind2     <- which(maxi > th3 )                                         ## Careful not to use >=. 
      ind4     <-  which(sumaalpha > th4)
      
      if(p != alpha){
        suma  <- sapply( 1:floor(n/bu[i]), function(k) sum(   path[((k-1)*bu[i] + 1):(k*bu[i])]^p )  )
        
        #sortedsuma <- sort(suma, decreasing = TRUE)
        ## thresholds
        #th1     <- th11 <- sorted[ max(2,floor(n^thr[1])) ]                        ## empirical th.
        th2     <-  sort(suma, partial = mn-kn+1)[(mn-kn+1)]    #sortedsuma[ max(1,floor( (n)/bu[i]^thr2[1]))+1 ]                ## p-norm th.          ## a-norm th 
        ## index
        #ind1     <- which(suma > th1^(alpha) )
        ind1     <- which(suma > th2 )                                            ## choosing clusters
        # functions c(1)
        
        g1 <- function(vector, p)      (sum(vector^alpha)/(sum(vector^p)^alpha))          ##1/c(p)
        g2 <- function(vector,th3)     sum( sum(vector) > th3 )        
        g3 <- function(vector,th3)     sum( vector > th3 )
        g4 <- function(vector, p)      (sum(vector^p)^alpha/(sum(vector^alpha)))            ## c(p)
        
        #estim1       <-  mean( sapply(ind1,  function(k)    g1(path[((k-1)*bu[i] +1):(k*bu[i])] , p)  )) 
        estim2       <-  sum(  sapply(1:n2,  function(k)    g2(path[((k-1)*bu[i] +1):(k*bu[i])] , th3 )  )) 
        estim3       <-  sum(  sapply(1:n2,  function(k)    g3(path[((k-1)*bu[i] +1):(k*bu[i])] , th3 )  )) 
        estim4       <-  mean( sapply(ind4, function(k)     g4(path[((k-1)*bu[i] +1):(k*bu[i])] , p)  )) 
      
        estim2       <-   estim2/estim3 
        
      }
      else{
        # functions theta
        g1 <- function(vector, p)      (max(vector^p)/sum(vector^p))
        g2 <- function(vector,th3)     sum( vector > th3)#sum(vector^p)/max(vector^p)  #
        estim4       <-  mean( sapply(ind4, function(k)     g1(path[((k-1)*bu[i] +1):(k*bu[i])] , p)  )) 
        estim2      <-   mean(sapply(ind2, function(k)      g2(path[((k-1)*bu[i] +1):(k*bu[i])] , th3 )  )) 
        
        estim2 <- 1/estim2 
      }

        
        if(j==1){
          #if(p!=alpha)  infor1       <- rbind( infor1, c(estim1, 3, th2, bu[i], thr2[1], mode=j))  ##lp inference
          infor1       <- rbind( infor1, c(estim2, 1, th3, bu[i], thr2[1], mode=j))
          infor1       <- rbind( infor1, c(estim4, 2, th4, bu[i], thr2[1], mode=j))  ##la-inference
        }
        else{
          #if(p!=alpha)  infor2       <- rbind( infor2, c(estim1, 3, th2, bu[i], thr2[1], mode=j))
          infor2       <- rbind( infor2, c(estim2, 1, th3, bu[i], thr2[1], mode=j))
          infor2       <- rbind( infor2, c(estim4, 2, th4, bu[i], thr2[1], mode=j))
        } 
      }
  }
}

names(infor2) <- names(infor1) <-  c("DATA", "TYPE", "zn", "BL", "TH", "MODE")
g1 <- niceplot(infor1,estim=c1[1],title)
g2 <- niceplot(infor2,estim=c1[2],title)
grid.arrange(g1,g2,nrow=1)


