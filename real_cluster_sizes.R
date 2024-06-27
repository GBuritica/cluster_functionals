#######################################################################
#######################################################################
#######################################################################
#######################################################################
##### Gloria Buriticá
####  Computation of cluster lengths for SRE
####  
#### 
#### Boxplot of estimates of the cluster index of sums and the extremal index
#### Comparison: alpha-cluster-based vs. maxima-cluster-based estimator
#######################################################################

real_cluster_sizes <- function(N=10000){
  res <- NULL
  for(i in 1:N){
    n <- 50000
    E <- rexp(1)
    St <- cumsum(c(rnorm(n) - 0.5))
    s  <- sum(St > -E)
    res <- c(res, s)
    #print(i)
  }
  ### pi_j
  theta <- mean(res==0)
  ci    <- qnorm(0.975)*sd((res==0))/sqrt(length(res))
  theta <- c(theta-ci,theta,theta+ci)
  #print(theta)
  for(j in 0:4){
    ci  <- qnorm(0.975)*sd( (res==j) - (res ==(j+1) ) )/sqrt(length(res))
    pij <- (mean(res==j) - mean(res==(j+1) )) 
    pij <- c(pij- ci, pij,pij+ci) 
    theta <- cbind(theta,pij)
  }
  colnames(theta) <- c('theta', 'pi0', 'pi1', 'pi2', 'pi3', 'pi4' )
  rownames(theta) <- c( 'lbound', 'estimate' , 'ubound' )
  as.data.frame(theta)  
}
#real_cluster_sizes(10)

