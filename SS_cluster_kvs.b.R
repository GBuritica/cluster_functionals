#######################################################################
#######################################################################
#######################################################################
#######################################################################
##### Gloria Buriticá
#### SS - cluster functionals - e.g. extremal index - c(1)
#### k vs. block length
#### 
#######################################################################
source("/Users/Buritica/Dropbox/Thèse/git/index_regular_variation/IndexofRV.R")
library(ggplot2)
library(plotly)
#library(latex2exp)
library(gridExtra)



N0 <- 100
n  <- 8000


phi0    <- c(0.8,0.6,0.1)
alpha00  <-  c(0.8,1.3)
thet    <-  sapply(1:length(phi0), function(k) 1-phi0[k]^alpha00)
c1      <-  sapply(1:length(phi0), function(k) (1-phi0[k])^alpha00/thet[k])


a <- 2
b <- 2
estim <- thet[a,b]#0.2792#
phi    <- phi0[b]
alpha0 <- alpha00[a]#11#
name0  <- "theta"
kmax  <-  min(floor(n/4),floor(n*0.1))              ## values in k
bn0   <-  seq(2,sqrt(n),1) #sapply(2:6, function(k) 2^k)                ## values in bn

##MSE
info  <-  matrix(data=0, nrow= length(bn0), ncol = kmax, dimnames=list( bn0 , 1:kmax ) )  #data.frame( matrix(ncol=3,nrow=0, dimnames=list(NULL, c("bn", "kn",name0))))
info2 <-  matrix(data=0, nrow= length(bn0), ncol = kmax, dimnames=list( bn0 , 1:kmax ) )

## VARIANCE
infov  <-  matrix(data=0, nrow= length(bn0), ncol = kmax, dimnames=list( bn0 , 1:kmax ) )  #data.frame( matrix(ncol=3,nrow=0, dimnames=list(NULL, c("bn", "kn",name0))))
infov2 <-  matrix(data=0, nrow= length(bn0), ncol = kmax, dimnames=list( bn0 , 1:kmax ) )

## mean
infom  <-  matrix(data=0, nrow= length(bn0), ncol = kmax, dimnames=list( bn0 , 1:kmax ) )  #data.frame( matrix(ncol=3,nrow=0, dimnames=list(NULL, c("bn", "kn",name0))))
infom2 <-  matrix(data=0, nrow= length(bn0), ncol = kmax, dimnames=list( bn0 , 1:kmax ) )

for(N in 1:1000){
  print(N)
  #path   <-   abs(ARCHm(n))#
  path   <-   abs(arima.sim(n=n, list(ar=phi, ma=0), rand.gen=function(n) rt(n,df=alpha0) ))
  alpha  <-   alpha0#(1/alphaestimator(path,k1=floor(n^0.8))$xi);print(alpha)
  p      <- alpha 
  for(j in  1:length(bn0) ){
    bn  <- bn0[j]
    mn <- min(floor(n/bn),floor( n*0.1) )
    samplepsum  <- sapply( 1:mn, function(k)   sum( path[( ((k-1)*bn) +1):(k*bn)]^p ) )  ## |X_[1,bn]|_p^p sample.
    samplemaxsum  <- sapply( 1:mn, function(k) max( path[( ((k-1)*bn) +1):(k*bn)] ) )
    
    for(kn  in 1:(mn-1) ){
      kth   <- sort(samplepsum,  partial = (mn-kn) )[(mn-kn)]       ## places the kth largest at the n-kn+1 place 
      kth2   <- sort(samplemaxsum,  partial = (mn-kn) )[(mn-kn)]
      index <-  which( samplepsum > kth )                               ## gives  index. of the k largest. 
      index2<-  which( samplemaxsum > kth2 ) 
      
      ## Computes estimates
      estimate <-  clusterfunctional(path,index,bn,p,kth,cfun1)              ## computes kn estimates
      estimate2 <- 1/clusterfunctional(path,index2,bn,p,kth2,cfun11)              ## computes kn estimates
        
      mean <-  ( (N-1)*infom[j,kn] + estimate)/N 
      mean2 <- ( (N-1)*infom2[j,kn] + estimate2)/N 
      ##MSE
      info[j,kn] <- ( (N-1)*info[j,kn] + (estimate-estim)^2)/N                                ## info into matrix
      info2[j,kn]<- ( (N-1)*info2[j,kn]+ (estimate2-estim)^2)/N     
      
      ##absolute squared error
      infov[j,kn] <- infov[j,kn] + (infom[j,kn])^2 - mean^2 + 
                            ( (estimate^2-infov[j,kn]-infom[j,kn]^2)/N)                              ## info into matrix
      infov2[j,kn] <-infov2[j,kn] + (infom2[j,kn])^2 - mean2^2 + 
                            ( (estimate2^2-infov2[j,kn]-infom2[j,kn]^2)/N) 
    
      ## Updates info
      ##Mean
      infom[j,kn]  <- mean                       ## info into matrix
      infom2[j,kn] <- mean2
      
    }
  }
}

kmatrix <- info
for(i in 1:kmax) kmatrix[ ,i] <- rep(i,length(bn0))


#save(info,info2,infov,infov2,infom,infom2,alpha,phi,n,estim,kmax,bn0, file='2704AR.61000.RData')
#save(info,info2,infov,infov2,infom,infom2,alpha,phi,n,estim,kmax,bn0, file='2704ARCH3000.RData')
#load( '2704AR.81000RData')

bias1 <- (info*(infom-estim)^2/info)
bias2 <- (info*(infom2-estim)^2/info)
mse1 <-  sqrt(info*info/info)
mse2 <-  sqrt(info2*info2/info2)
var1 <-  (infov*infov/infov)
var2 <-  (infov2*infov2/infov2)

min(bias1, na.rm=T) < min(bias2, na.rm=T)
(min( sqrt(bias1), na.rm=T)-min( sqrt(bias2), na.rm=T))*100/min( sqrt(bias2), na.rm=T)

which(bias1 == min(bias1, na.rm=T), arr.ind = TRUE)
which(bias2 == min(bias2, na.rm=T), arr.ind = TRUE)

min(mse1,na.rm = T) < min(mse2,na.rm = T)
(min(mse1,na.rm = T)- min(mse2,na.rm = T))*100/min(mse2,na.rm = T)


which(mse1 == min(mse1, na.rm=T), arr.ind = TRUE)
which(mse2 == min(mse2, na.rm=T), arr.ind = TRUE)



subtit <- TeX('\\alpha = 1.3, \\varphi = .6, \\theta = .4852$')


#fig11 <-plot_ly(x=1:kmax, y= bn0, z = ~bias1, type = "heatmap", 
#                coloraxis = 'coloraxis') 


#fig22 <- plot_ly(x=1:kmax, y= bn0, z = ~bias2, type = "heatmap",
#                 coloraxis = 'coloraxis' ) 


#fig1 <- subplot(fig11, fig22)
#fig1 <- fig1 %>% layout(coloraxis=list(colorscale='Viridis',cmin=0 ))
#fig1 <- fig1 %>% layout(title="Bias" ); fig1



fig33 <- plot_ly(x=1:kmax,y= bn0, z = ~(mse1), type = "heatmap", coloraxis = 'coloraxis'   ) 
fig33 <- fig33 %>% layout(xaxis=list(), yaxis=list(title="b - block lengths"))

fig44 <- plot_ly(x=1:kmax,y= bn0, z = ~(mse2), type = "heatmap" ,coloraxis = 'coloraxis')
fig44 <- fig44 %>% layout(xaxis= list(title="k - order statistics",showgrid = F), yaxis=list(title="b - block lengths",showgrid = F))

fig2 <- subplot(fig33, fig44, titleY = TRUE, titleX = TRUE, margin = 0.01 , nrows=2)

fig2 <- fig2 %>% layout(coloraxis=list(cmin=0, cmax=max(mse1,na.rm=T)), colorscale='Reds',
                        title = 'Root mean squared errors',
                        xaxis = list( 
                          zerolinecolor = '#ffff', 
                          zerolinewidth = 2, 
                          gridcolor = 'ffff'), 
                        yaxis = list( 
                          zerolinecolor = '#ffff', 
                          zerolinewidth = 2, 
                          gridcolor = 'ffff'),
                        annotations = list( 
                          list( x = 0.5,  
                                y = .995,  
                                text = subtit, 
                                xref = "paper",  
                                yref = "paper",  
                                xanchor = "center",  
                                yanchor = "bottom",  
                                showarrow = FALSE ), 
                          list( x = 0.4,  
                                y = 0.9,  
                                text = TeX('$\\alpha-\\text{inference}$') , 
                                xref = "paper",  
                                yref = "paper",  
                                xanchor = "center",  
                                yanchor = "bottom",  
                                showarrow = FALSE ), 
                          list( x = 0.4,  
                                y = 0.4,  
                                text = TeX('$\\infty-\\text{inference}$'),  
                                xref = "paper",  
                                yref = "paper",  
                                xanchor = "center",  
                                yanchor = "bottom",  
                                showarrow = FALSE )) )
fig2 <- fig2 %>% config(mathjax = 'cdn') ## For latex to Appear
fig2




fig55 <- plot_ly(x=1:kmax,y= bn0, z = ~sqrt(var1), type = "heatmap",coloraxis = 'coloraxis') 
fig55 <- fig55 %>% layout(yaxis=list(title="b - block lengths"))

fig66 <- plot_ly(x=1:kmax,y= bn0, z = ~sqrt(var2), type = "heatmap",coloraxis = 'coloraxis' ) 
fig66 <- fig66 %>% layout( xaxis= list(title="k - order statistics", showgrid=F), yaxis=list(title="b - block lengths", showgrid=F))


fig3 <- subplot(fig55, fig66,  nrows=2, titleY = TRUE, titleX = TRUE, margin = 0.02 )
fig3 <- fig3 %>% layout( coloraxis=list(cmin=0, cmax=max(sqrt(var1),na.rm = T) ), 
                      title = TeX('$\\text{Root variance}'),
                      xaxis = list( 
                        zerolinecolor = '#ffff', 
                        zerolinewidth = 2, 
                        gridcolor = 'ffff'), 
                      yaxis = list( 
                        zerolinecolor = '#ffff', 
                        zerolinewidth = 2, 
                        gridcolor = 'ffff'),
                      annotations = list( 
                        list( x = 0.54,  
                              y = .99,  
                              text = subtit, 
                              xref = "paper",  
                              yref = "paper",  
                              xanchor = "center",  
                              yanchor = "bottom",  
                              showarrow = FALSE ), 
                        list( x = 0.4,  
                              y = 0.9,  
                              text = TeX('$\\alpha-\\text{inference}$') , 
                              xref = "paper",  
                              yref = "paper",  
                              xanchor = "center",  
                              yanchor = "bottom",  
                              showarrow = FALSE ), 
                        list( x = 0.4,  
                              y = 0.4,  
                              text = TeX('$\\infty-\\text{inference}$'),  
                              xref = "paper",  
                              yref = "paper",  
                              xanchor = "center",  
                              yanchor = "bottom",  
                              showarrow = FALSE )) )
fig3 <- fig3 %>% config(mathjax = 'cdn'); fig3


fig11 <- plot_ly(x=1:kmax, y= bn0, z = ~sqrt((bias1)), type = "heatmap",  coloraxis = 'coloraxis' ) 
fig11 <- fig11 %>% layout( yaxis=list(title="b - block lengths"))

fig22 <- plot_ly(x=1:kmax, y= bn0, z = ~sqrt((bias2)), type = "heatmap", coloraxis = 'coloraxis' ) 
fig22 <- fig22 %>% layout(xaxis= list(title="k - order statistics",showgrid=F), yaxis=list(title="b - block lengths",showgrid=F))


fig3 <- subplot(fig11, fig22, nrows=2, titleY = TRUE, titleX = TRUE, margin = 0.02 )
fig3 <- fig3 %>% layout( coloraxis=list(cmin=0,cmax=max(sqrt(bias1),na.rm=T)) , 
                         title = TeX('$\\text{Root bias}'),
                         xaxis = list( 
                           zerolinecolor = '#ffff', 
                           zerolinewidth = 2, 
                           gridcolor = 'ffff'), 
                         yaxis = list( 
                           zerolinecolor = '#ffff', 
                           zerolinewidth = 2, 
                           gridcolor = 'ffff'),
                         annotations = list( 
                           list( x = 0.53,  
                                 y = .99,  
                                 text = subtit, 
                                 xref = "paper",  
                                 yref = "paper",  
                                 xanchor = "center",  
                                 yanchor = "bottom",  
                                 showarrow = FALSE ), 
                           list( x = 0.4,  
                                 y = 0.9,  
                                 text = TeX('$\\alpha-\\text{inference}$') , 
                                 xref = "paper",  
                                 yref = "paper",  
                                 xanchor = "center",  
                                 yanchor = "bottom",  
                                 showarrow = FALSE ), 
                           list( x = 0.4,  
                                 y = 0.4,  
                                 text = TeX('$\\infty-\\text{inference}$'),  
                                 xref = "paper",  
                                 yref = "paper",  
                                 xanchor = "center",  
                                 yanchor = "bottom",  
                                 showarrow = FALSE )) )
fig3 <- fig3 %>% config(mathjax = 'cdn'); fig3


############## Relative gain

rmse <- ( (sqrt(info)-sqrt(info2))*100/sqrt(info2))
rbias <-  info*( sqrt(abs(infom-estim))-sqrt(abs(infom2-estim)) )*100/(info*sqrt(abs(infom2-estim)) )


fig1 <- plot_ly(y=bn0,z = ~rmse ,
                type = "contour",  
                colorscale='Reds',
                showscale=F, contours = 
                  list(showlabels = TRUE,
                       start=-100,
                       end=100,
                       size=20)) 

fig1 <- fig1 %>% layout(  
                          title = ('Relative percentage change RMSE'),
                          xaxis= list(title="k - order statistics",showgrid=F), 
                          yaxis=list(title="b - block lengths",showgrid=F),
                          annotations = list( list( x = 0.5,  
                          y = .995,  
                         text = subtit , 
                          xref = "paper",  
                          yref = "paper",  
                           xanchor = "center",  
                           yanchor = "bottom",  
                           showarrow = FALSE ))) 
fig1<- fig1 %>% config(mathjax = 'cdn'); fig1



#fig2 <- plot_ly(y=bn0,z = ~rbias,
#                type = "contour",contours = list(showlabels = TRUE,start=-1000,
#                                                 end=1000,
#                                                 size=500) ) 


#subplot(fig1, fig2)




par(mfrow=c(1,3))
plot.ts((infom[1,]-estim)^2, ylim=c(0,0.02), xlim = c(0,250))
sapply(2:4, function(k) lines((infom[k,]-estim)^2 , col = k))
sapply(1:4, function(k) lines((infom2[k,]-estim)^2, col = k, lty=2) )


plot.ts(info[1,], ylim=c(0,0.02), xlim = c(0,250))
sapply(2:4, function(k) lines(info[k,] , col = k))
sapply(1:4, function(k) lines(info2[k,], col = k, lty=2) )

plot.ts(infov[1,], ylim=c(0,0.02), xlim = c(0,250))
sapply(2:4, function(k) lines(infov[k,] , col = k))
sapply(1:4, function(k) lines(infov2[k,], col = k, lty=2) )

plot.ts(infom[1,], ylim=c(0,1))
sapply(2:4, function(k) lines(infom[k,] , col = k))
sapply(1:4, function(k) lines(infom2[k,], col = k, lty=2) )
abline(h=estim)


#######################################################################
#######################################################################
########### Auxiliar functions


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





