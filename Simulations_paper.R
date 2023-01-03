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
source("IndexofRV.R")
source("aux2_functions.R")
library(ggplot2)
library(plotly)
library(plot3D)
library(latex2exp)



########################################################################
########################################################################
########################################################################
#### produces SS
#### Comparison k1 vs k2 for the extremal index
n       <- 12000
klim    <- floor(n^0.75)
b       <- unique(floor(sqrt(n/1:klim)))
klim1   <- floor(n/b^2)                    ## getting index
klim1   <- klim1[-1]
phi     <- 0.5
theta  <-  1-phi
results  <- matrix(data=0, nrow = length(klim1), ncol = length(klim1) ) ## saves variance
results2 <- matrix(data=0, nrow = length(klim1), ncol = length(klim1) ) ## saves mean
N0       <- 500
for(N in 1:2){
  print(N)
  path    <- ARCHm(n)#ARm(n,phi)#
  ordered <- sort(path, decreasing=T)
  alpha   <- 1/alphaestimator(ordered,k1=klim1)$xi 
  
  
  for(k1 in 1:length(klim1)){  ## for different alpha values
    al <- max(0.01,alpha[k1])
    al <- min(al,10)
    estimateei      <- eiCP(path,al,klim1,n)
    for( l in 1:length(klim1) ){
      mean          <- results[k1,l]
      var           <- results2[k1,l]
      mean2         <- mean +  ((estimateei[l,2] - mean)/N)

      results2[k1,l] <- var  + mean^2 - mean2^2 + ((estimateei[l,2]^2 - var - mean^2)/N)
      results[k1,l]  <- mean2
      }
  }
}


#######################################################################
#######################################################################
#######################################################################
##### Plot Heat maps
lim <- 40
MSE <- sqrt((results2 +abs(results-theta)^2))
fig1 <-plot_ly(x = klim1[1:lim], y = klim1[1:lim],
        z = ~MSE, type = "contour",
        autocontour = F,
        contours = list(
          start = 0,  end = 0.4, size = 0.02, showlabels = TRUE
        ))%>% hide_colorbar() %>% 
        layout(title  = 'Mean squared errors', xaxis = list(title = 'k1'),
       yaxis = list(title = 'k2'))

fig1
Var <- sqrt(results2)
fig2 <- plot_ly(x = klim1[1:lim], y = klim1[1:lim],
        z = ~Var, type = "contour",
        autocontour = F,
        contours = list(
          start = 0,  end = 0.4, size = 0.01, showlabels = TRUE
        ))%>% hide_colorbar()%>% 
      layout(title  = 'Standard deviations', xaxis = list(title = 'k1'),
             yaxis = list(title = 'k2'))
fig2
bias <- sqrt((abs(results-theta)^2)) 
fig1 <-plot_ly(x = klim1[1:lim], y = klim1[1:lim],
               z = ~bias, type = "contour",
               autocontour = F,
               contours = list(
                 start = 0,  end = 0.4, size = 0.02, showlabels = TRUE
               ))%>% hide_colorbar() %>% 
  layout(title  = 'Bias', xaxis = list(title = 'k1'),
         yaxis = list(title = 'k2'))

fig1


#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#### Cluster sizes
#### SS
#### Checking the Gaussian profile of estimates
n            <- 12000
klim         <- floor(n^0.75)
b            <- unique(floor(sqrt(n/1:klim)))
klim1        <- floor(n/b^2) 
phi          <- 0.5
path         <- ARCHm(n)
alphae       <- 1
new          <- eiCP2(path,alphae,klim1,n) 
pi0          <- new$estimate[,1] 
pi0var       <- new$variance[,1] 
pi1          <- new$estimate[,2] 
pi1var       <- new$variance[,2] 
pi2          <- new$estimate[,3] 
pi2var       <- new$variance[,3] 
pi3          <- new$estimate[,4] 
pi3var       <- new$variance[,4] 
pi5          <- new$estimate[,5] 
pi5var       <- new$variance[,5] 
for(i in 1:10){
  print(i)
  path         <- ARCHm(n)
  alpha        <- 1
  new          <- eiCP2(path,alpha,klim1,n) 
  
  ## Save data
  pi0          <- rbind(pi0, new$estimate[,1] )
  pi0var       <- rbind(pi0var, new$variance[,1] )
  pi1          <- rbind(pi1, new$estimate[,2] )
  pi1var       <- rbind(pi1var, new$variance[,2] )
  pi2          <- rbind(pi2, new$estimate[,3] )
  pi2var       <- rbind(pi2var, new$variance[,3] )
  pi3          <- rbind(pi3, new$estimate[,4] )
  pi3var       <- rbind(pi3var, new$variance[,4] )
  pi5          <- rbind(pi5, new$estimate[,5] )
  pi5var       <- rbind(pi5var, new$variance[,5] )

  alphae       <- c(alphae,alpha)
}

#######################################################################
#######################################################################
#######################################################################
#### Plot histogram
k     <- 2
ki    <- klim1[k]
pi    <- pi0[,k]
pivar <- pi0var[,k]
theta <- 0.2792# 
median(pi)

hist(pi,prob=T,  breaks = 80, ylim = c(0,200),
     density=20, xlab = " ", main = "Cluster Size 2")
abline( v = theta , lty = 1 , col = "red")

sd <- sqrt( mean(pivar - pi^2) )
sd2 <- sqrt(mean( (pi-mean(pi))^2 ) )


curve(dnorm(x, mean=median(pi), sd=sd/sqrt(ki) ), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n" )

curve(dnorm(x, mean=median(pi), sd=sd(pi) ), 
      col="darkblue", lwd=1, add=TRUE, yaxt="n" ,lty=2 )

#######################################################################
#### Plot boxplot
par(mfrow = c(1,1))
boxplot(pi,freq=F, ylim = c(0,1), ylab = "Extremal index"); abline( h = theta , lty = 2 , col = "red")
boxplot(pivar,freq=F, ylim = c(0,1))


#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
## Estimate as a function of k with confidence intervals 
## for one single trajectory

n <- 12000
klim         <- floor(n^0.7)
b            <- unique(floor(sqrt(n/1:klim)))
klim1        <- floor(n/b^2) 
path         <- ARCHm(n)
alphae       <- 1/alphaestimator(path,k1=1200)$xi 
estimate     <- eiCP2(path,alphae,klim1,n) 
ks           <- estimate$k
plot(ks, estimate$estimate[,1],type = "l", xlim = c(0,200) , ylim = c(0,1), 
     xlab = "k", ylab = " ")
lines(ks, estimate$estimate[,1] + qnorm(0.975)*sqrt(abs(estimate$variance[,1] - estimate$estimate[,1]^2)/ks) , lty = 3, col = 'black' )
lines(ks, estimate$estimate[,1] - qnorm(0.975)*sqrt(abs(estimate$variance[,1] - estimate$estimate[,1]^2)/ks) , lty = 3, col = 'black' )
points(ks,estimate$estimate[,1], pch = 16, cex = 0.5)
abline(h=theta, col = 'red', lty = 2)




