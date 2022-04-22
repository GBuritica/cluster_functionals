##### Simulation extremal index
##### Asymptotic normality cluster estimator
library(ggplot2)
library(plotly)
library(plot3D)
library(latex2exp)

#path    <- ARm(n)
#ordered <- sort(path, decreasing=T)
#alpha   <- rep(1,length(klim1))#hill(ordered,klim1)

#for(k1 in 1:length(klim1)){  ## for different alpha values
#  estimateei   <- eiCP(path,alpha[k1],klim1,n)
#  for( l in 1:length(klim1) ) results <- rbind( results, c(estimateei[l,2],klim1[k1],klim1[l]) ) 
  #print(k1)
#}

#results <- as.data.frame(results)
#colnames(results) <- c("estimate","k1","k2")

#ggplot(data = results, aes(y=as.factor(k1), x=as.factor(k2), fill=abs(estimate-theta) )) + 
#  geom_tile()


#### produces results
#### SS
n       <- 8000
klim    <- floor(n^0.7)
b       <- unique(floor(sqrt(n/1:klim)))
klim1   <- floor(n/b^2)                    ## getting index
klim1   <- klim1[-1]

results  <- matrix(data=0, nrow = length(klim1), ncol = length(klim1) )
results2 <- matrix(data=0, nrow = length(klim1), ncol = length(klim1) )
N0      <- 500
for(N in 1:N0){
  print(N)
  path    <- ARm(n,0.7)
  ordered <- sort(path, decreasing=T)
  alpha   <- hill(ordered,klim1)
  
  for(k1 in 1:length(klim1)){  ## for different alpha values
    estimateei   <- eiCP(path,alpha[k1],klim1,n)
    for( l in 1:length(klim1) ){
      mean          <- results[k1,l]
      var           <- results2[k1,l]
      mean2         <- mean +  ((estimateei[l,2] - mean)/N)

      results2[k1,l] <- var  + mean^2 - mean2^2 + ((estimateei[l,2]^2 - var - mean^2)/N)
      results[k1,l]  <- mean2
      }
    }
}


#load(results,N0,theta,file = "2610ARCHm.Rdata")
#load(results,N0,theta,file = "2610ARm.7.Rdata")

#load(results,results2,N0,theta,klim1,b,n,file = "2311ARCHm.Rdata")
#save(results,results2,N0,theta,klim1,b,n,file = "2311AR.7.Rdata")
#load("2311AR.7.Rdata")


### setting up data frame 
load("2311AR.7.Rdata")
results3 <- NULL
for( k1 in 1:length(klim1)){
  for(k2 in 1:length(klim1)){
    results3<- rbind(results3, c(results[k1,k2],klim1[k1],klim1[k2]))
  }
}

plot_ly(z = ~sqrt(results/N))
N <- N0
### plot in b
fig <- plot_ly(x = sqrt(n/klim1), y = sqrt(n/klim1), z =sqrt(results/N) )
fig <- fig %>% add_surface()
fig

par(mfrow=c(2,3))
## variance
image2D( x = klim1, y = sqrt(n/klim1), z =(sqrt(results2)),
         ylab = TeX('b2 : $\\widehat{\\theta}^{scp}(b2,\\widehat{\\alpha})$'), 
         xlab = TeX('k1 : $\\widehat{\\alpha}(k1)$'),cex.axis=0.8,
         shade = 0.2, resfac = 1,#breaks = c(seq(0.2,1.2,.2),2,4,8),
         #contour = list(col = "white", labcex = 0.5, lwd = 2, alpha = 0.3),
         main = "squared variance")
# bias
image2D( x = klim1, y = sqrt(n/klim1), z =abs(results-theta),
         ylab = TeX('b2 : $\\widehat{\\theta}^{scp}(b2,\\widehat{\\alpha})$'), 
         xlab = TeX('k1 : $\\widehat{\\alpha}(k1)$'),cex.axis=0.8,
         shade = 0.2, resfac = 1,#breaks = c(seq(0.2,1.2,.2),2,4,8),
         #contour = list(col = "white", labcex = 0.5, lwd = 2, alpha = 0.3),
         main = "absolute bias")
## mean squared errors
image2D(x = klim1, y = sqrt(n/klim1), z =sqrt(results2 +abs(results-theta)^2),
        shade = 0.2, cex.axis=0.8,
        ylab = TeX('b2 : $\\widehat{\\theta}^{scp}(b2,\\widehat{\\alpha})$'), 
        xlab = TeX('k1 : $\\widehat{\\alpha}(k1)$'),cex.axis=0.8,
        #contour = list(col = "white", labcex = 0.5, lwd = 3, alpha = 0.3)
        main = "MSE")
contour2D(x = klim1, y = sqrt(n/klim1), z =(sqrt(results2)),
          ylab = TeX('b2 : $\\widehat{\\theta}^{scp}(b2,\\widehat{\\alpha})$'), 
          xlab = TeX('k1 : $\\widehat{\\alpha}(k1)$'),cex.axis=0.8,
          #contour = TRUE, shade = 0.5, lphi = 0, 
          main = "squared variance")
## bias
contour2D(x = klim1, y = sqrt(n/klim1), z =abs(results-theta), 
          ylab = TeX('b2 : $\\widehat{\\theta}^{scp}(b2,\\widehat{\\alpha})$'), 
          xlab = TeX('k1 : $\\widehat{\\alpha}(k1)$'),cex.axis=0.8,
          #contour = TRUE, shade = 0.5, lphi = 0, 
          main = "absolute bias")
contour2D(x = klim1, y = sqrt(n/klim1), z =sqrt(results2 +abs(results-theta)^2), 
          ylab = TeX('b2 : $\\widehat{\\theta}^{scp}(b2,\\widehat{\\alpha})$'), 
          xlab = TeX('k1 : $\\widehat{\\alpha}(k1)$')
          ,contour = TRUE, shade = 0.5, lphi = 0, main = "MSE",cex.axis=0.8)

### plot in k = o(n/b) = n/b^2
fig <- plot_ly(x = (klim1), y = (klim1), z =sqrt(results/N) )
fig <- fig %>% add_surface()
fig

par(mfrow= c(2,3))
## variance
image2D( x = klim1, y = klim1, z =(sqrt(results2)),
         ylab = TeX('k2 : $\\widehat{\\theta}^{scp}(k2,\\widehat{\\alpha})$'), 
         xlab = TeX('k1 : $\\widehat{\\alpha}(k1)$'),cex.axis=0.8,
         shade = 0.2, resfac = 1,#breaks = c(seq(0.2,1.2,.2),2,4,8),
         #contour = list(col = "white", labcex = 0.5, lwd = 2, alpha = 0.3),
         main = "squared variance")
## bias
image2D( x = klim1, y = klim1, z =abs(results-theta),
         ylab = TeX('k2 : $\\widehat{\\theta}^{scp}(k2,\\widehat{\\alpha})$'), 
         xlab = TeX('k1 : $\\widehat{\\alpha}(k1)$'),cex.axis=0.8,
         shade = 0.2, resfac = 1,#breaks = c(seq(0.2,1.2,.2),2,4,8),
         #contour = list(col = "white", labcex = 0.5, lwd = 2, alpha = 0.3),
         main = "absolute bias")
## mean squared errors
image2D( x = klim1, y = klim1, z = sqrt(results2 +abs(results-theta)^2)  ,
         ylab = TeX('k2 : $\\widehat{\\theta}^{scp}(k2,\\widehat{\\alpha})$'), 
         xlab = TeX('k1 : $\\widehat{\\alpha}(k1)$'),cex.axis=0.8,
         shade = 0.2, resfac = 1,#breaks = c(seq(0.2,1.2,.2),2,4,8),
         #contour = list(col = "white", labcex = 0.5, lwd = 2, alpha = 0.3),
         main = "mean squared errors")
contour2D(x = klim1, y = klim1, z =(sqrt(results2)),ylim = c(0,200),
          ylab = TeX('k2 : $\\widehat{\\theta}^{scp}(k2,\\widehat{\\alpha})$'), 
          xlab = TeX('k1 : $\\widehat{\\alpha}(k1)$'),cex.axis=0.8,
          #contour = TRUE, shade = 0.5, lphi = 0, 
          main = "squared variance")
contour2D(x = klim1, y = klim1, z =abs(results-theta), ylim = c(0,200),
          ylab = TeX('k2 : $\\widehat{\\theta}^{scp}(k2,\\widehat{\\alpha})$'), 
          xlab = TeX('k1 : $\\widehat{\\alpha}(k1)$'),cex.axis=0.8,
          #contour = TRUE, shade = 0.5, lphi = 0, 
          main = "absolute bias")
contour2D(x = klim1, y = klim1, z =sqrt(results2 +abs(results-theta)^2), ylim = c(0,200),
          ylab = TeX('k2 : $\\widehat{\\theta}^{scp}(k2,\\widehat{\\alpha})$'), 
          xlab = TeX('k1 : $\\widehat{\\alpha}(k1)$'),cex.axis=0.8,
          #contour = TRUE, shade = 0.5, lphi = 0, 
          main = "mean squared errors")


par(mfrow= c(2,3))
## variance
image2D( x = sqrt(n/klim1), y = sqrt(n/klim1), z =(sqrt(results2)),
         ylab = TeX('b2 : $\\widehat{\\theta}^{scp}(b2,\\widehat{\\alpha})$'), 
         xlab = TeX('b1 : $\\widehat{\\alpha}(n/b1^2)$'),cex.axis=0.8,
         shade = 0.2, resfac = 1,#breaks = c(seq(0.2,1.2,.2),2,4,8),
         #contour = list(col = "white", labcex = 0.5, lwd = 2, alpha = 0.3),
         main = "squared variance")
## bias
image2D( x = sqrt(n/klim1), y = sqrt(n/klim1), z =abs(results-theta),
         ylab = TeX('b2 : $\\widehat{\\theta}^{scp}(b2,\\widehat{\\alpha})$'), 
         xlab = TeX('b1 : $\\widehat{\\alpha}(n/b1^2)$'),cex.axis=0.8,
         shade = 0.2, resfac = 1,#breaks = c(seq(0.2,1.2,.2),2,4,8),
         #contour = list(col = "white", labcex = 0.5, lwd = 2, alpha = 0.3),
         main = "absolute bias")
## mean squared errors
image2D( x = sqrt(n/klim1), y = sqrt(n/klim1), z = sqrt(results2 +abs(results-theta)^2)  ,
         ylab = TeX('b2 : $\\widehat{\\theta}^{scp}(b2,\\widehat{\\alpha})$'), 
         xlab = TeX('b1 : $\\widehat{\\alpha}(n/b1^2)$'),cex.axis=0.8,
         shade = 0.2, resfac = 1,#breaks = c(seq(0.2,1.2,.2),2,4,8),
         #contour = list(col = "white", labcex = 0.5, lwd = 2, alpha = 0.3),
         main = "mean squared errors")
contour2D(x = sqrt(n/klim1), y = sqrt(n/klim1), z =(sqrt(results2)),
          ylab = TeX('b2 : $\\widehat{\\theta}^{scp}(b2,\\widehat{\\alpha})$'), 
          xlab = TeX('b1 : $\\widehat{\\alpha}(n/b1^2)$'),cex.axis=0.8,
          #contour = TRUE, shade = 0.5, lphi = 0, 
          main = "squared variance")
contour2D(x = sqrt(n/klim1), y = sqrt(n/klim1), z =abs(results-theta),
          ylab = TeX('b2 : $\\widehat{\\theta}^{scp}(b2,\\widehat{\\alpha})$'), 
          xlab = TeX('b1 : $\\widehat{\\alpha}(n/b1^2)$'),cex.axis=0.8,
          #contour = TRUE, shade = 0.5, lphi = 0, 
          main = "absolute bias")
contour2D(x = sqrt(n/klim1), y = sqrt(n/klim1), z =sqrt(results2 +abs(results-theta)^2), 
          ylab = TeX('b2 : $\\widehat{\\theta}^{scp}(b2,\\widehat{\\alpha})$'), 
          xlab = TeX('b1 : $\\widehat{\\alpha}(n/b1^2)$'),cex.axis=0.8,
          #contour = TRUE, shade = 0.5, lphi = 0, 
          main = "mean squared errors")

#### CLUSTER
n       <- 8000
klim    <- floor(n^0.7)

b       <- unique(floor(sqrt(n/1:klim)))
klim1   <- floor(n/b^2)                    ## getting index

res    <-  eiCP(path,alpha,klim1,n)
res    <-  eiCP2(path,alpha,klim,n)
par(mfrow=c(1,2))
plot(res, ylim = c(0,1), type="l");points(res, ylim = c(0,1), type="p");abline(h=theta)
plot.ts(res[,2]);abline(h=theta)

eiCP  <- function(path0,alpha0,klim0,n0){
  b <- floor(sqrt(n0/(klim0)))
  eireturn <- vector(mode="double", length=length(klim0) )
  for(k in 1:length(klim0)){
    b0 <- b[k]
    #print(b0)
    ## fixes b
    
    ## computes paths
    sumaalpha <- sapply( 1:floor(n0/b0), 
                    function(l) sum(path0[((l-1)*b0+1):(l*b0)]^alpha0) )
    maxaalpha <- sapply( 1:floor(n0/b0), 
                    function(l) max(path0[((l-1)*b0+1):(l*b0)]^alpha0) )
    ordered   <- order(sumaalpha, decreasing=T) ## returns de order (1), (2), (3)...
    ind       <- ordered[1:klim0[k]]                   ## I chose the ones less than k
    estimate  <- mean( maxaalpha[ind]/sumaalpha[ind] ) ## among these I compute the mean
    ### moves b forward
    eireturn[k] <- estimate

    
  }
  return( cbind(b,eireturn))
} 

eiCP2 <- function(path0,alpha0,klim0,n0){
  b <- floor(sqrt(n0/(1:klim0)))
  k <- 1
  eireturn <- vector(mode="double", length=klim0)
  while(k <= klim0){
    b0 <- b[k]
    #print(b0)
    ## fixes b
    
    ## computes paths
    sumaalpha <- sapply( 1:floor(n0/b0), 
                         function(l) sum(path0[((l-1)*b0+1):(l*b0)]^alpha0) )
    maxaalpha <- sapply( 1:floor(n0/b0), 
                         function(l) max(path0[((l-1)*b0+1):(l*b0)]^alpha0) )
    ordered   <- order(sumaalpha, decreasing=T) ## returns de order (1), (2), (3)...
    ind       <- ordered[1:k]                   ## I chose the ones less than k
    estimate  <- mean( maxaalpha[ind]/sumaalpha[ind] ) ## among these I compute the mean
    ### moves b forward
    while(b[k]==b0 && k <= klim0){
      eireturn[k] <- estimate
      k <- k+1
    }
    
  }
  return( cbind(b,eireturn))
} 

#### HILL
path    <- ARCHm(n)
ordered <- sort(path, decreasing=T)
plot(hill(ordered,2:n^0.8))
abline(h=1)

hill <- function(ordered,klim0){
  hill0 <- sapply(klim0, function(k) 
              mean( sapply(1:max(1,(k-1)), function(l) log(ordered[l]/ordered[k]) ) ) )
  return(hill0)
}

#### PATHS

## AX + B     with log A ~ N-0.5, B ~ Unif(0,1)
path <- ARCHm(n)
plot.ts(path)
theta <-  0.2792
alpha <- 1
ARCHm <- function(n0){
  x0  <- runif(1,0,1)
  nor <- rnorm(2*n0)
  uni <- runif(2*n0,0,1)
  for(i in 1:(2*n0) ) x0  <- c(x0, ( (exp(nor[i] - 0.5)*(x0[i])) + uni[i] ) )
  return(x0[(n0+1):(2*n0)])
}

## AR model
path <- ARm(n,0.7)
plot.ts(path)
theta <- 1-0.7^1
ARm  <- function(n0,par0=0.7){
  x0 <- arima.sim(n = n0, list(ar=par0, ma=0), rand.gen=function(n) abs( rt(n,df=1) ) ) 
  return(x0)
}

path <- sARCHm(n,0.5)
plot.ts(path)
theta <- 0.727
sARCHm <- function(n0,lambda0){
  x0           <- vector(mode = "double" , 3*n0)
  norm         <- rnorm(3*n0)
  x0[1]        <- (2*10^(-5))
  for( i in 2:(3*n0)) x0[i] <-  ( (2*10^(-5)) + lambda0*x0[(i-1)] )*(norm[i])^2
  return(x0[(2*n0+1):(3*n0)])
}

path <- ARCHm2(n,0.5)
plot.ts(path)
theta <- 0.835
ARCHm2 <- function(n0,lambda0){
  x0           <- vector(mode = "double" , 3*n0)
  norm           <- rnorm(3*n0,0,1)
  x0[1]        <- 0#2*10^{-5}   
  for( i in 2:(3*n0)) x0[i] <-  (norm[i])*sqrt( (2*10^(-5)) + lambda0*(x0[(i-1)])^2 )
  return(x0[(2*n0+1):(3*n0)]) 
}
#######################################################################



