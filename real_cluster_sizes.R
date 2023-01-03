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
res <- NULL
for(i in 1:10000){
  n <- 50000
  E <- rexp(1)
  St <- cumsum(c(rnorm(n) - 0.5))
  s  <- sum(St > -E)
  res <- c(res, s)
  print(i)
}

### pi_j
j <- 0
ci  <- qnorm(0.975)*sd( (res==j) - (res ==(j+1) ) )/sqrt(length(res))
pij <- (mean(res==j) - mean(res==(j+1) )) 
pij <- c(pij- ci, pij,pij+ci) 
print(pij)

### Extremal index
theta <- mean(res==0)
ci    <- qnorm(0.975)*sd((res==0))/sqrt(length(res))
theta <- c(theta-ci,theta,theta+ci)
print(theta)

