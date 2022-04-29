#######################################################################
#######################################################################
#######################################################################
#######################################################################
##### Gloria Buriticá
#### SS - Extremal index - c(1)
####  
#### 
#######################################################################
source("/Users/Buritica/Dropbox/Thèse/git/index_regular_variation/IndexofRV.R")
library(ggplot2)
library(latex2exp)
library(gridExtra)

n <-8000
thr2   <- thr <- c(2)
bu     <- sapply(1:6, function(k) 2^k)
phi0    <- c(0.8,0.6); 
infor1 <- data.frame( "DATA" = NULL, "TYPE"=NULL, "zn" =NULL, "BL" = NULL, "TH" = NULL, "MODE" = NULL )
infor2 <- data.frame( "DATA" = NULL, "TYPE"=NULL, "zn" =NULL, "BL" = NULL, "TH" = NULL, "MODE" = NULL )
alpha0 <- 1
thet <-   sapply(1:length(phi0), function(k) 1-phi0[k]^alpha0)
c1   <-   1/sapply(1:length(phi0), function(k) (1-phi0[k])^alpha0/thet[k])
title<-   TeX('$\\widehat{\\theta}_{| \\mathbf{X}|}$') #TeX('$\\frac{1}{\\widehat{c}(1)$')# 
## functions theta 
for(N in 11:500){
  for(j in 1:length(phi0)){
    phi    <- phi0[j]
    path   <- abs(arima.sim(n=n, list(ar=phi, ma=0), rand.gen=function(n) rt(n,df=alpha0) ))
    n      <- length(path)
    alpha  <- (1/alphaestimator(path,k1=floor(n^0.8))$xi);print(alpha)
    #alpha <- alpha0
    sorted <- sort(path,decreasing = T)
    p      <- 1#alpha#
    for(i in 1:length(bu)){
      suma       <- sapply( 1:floor(n/bu[i]), function(k) sum(   path[((k-1)*bu[i] + 1):(k*bu[i])]^p )  )
      if(p != alpha) sumaalpha  <- sapply( 1:floor(n/bu[i]), function(k) sum(   path[((k-1)*bu[i] + 1):(k*bu[i])]^alpha )  )
      maxi       <- sapply( 1:floor(n/bu[i]), function(k) max(   path[((k-1)*bu[i] + 1):(k*bu[i])] )  )
    
      sortedsuma      <- sort(suma, decreasing = TRUE); n2 <- length(suma)
      if(p != alpha) sortedsumaalpha <- sort(sumaalpha, decreasing = TRUE);
      sortedmaxi      <- sort(maxi, decreasing = TRUE); #print(n2)
    
        ## thresholds
        #th1     <- th11 <- sorted[ max(2,floor(n^thr[1])) ]                            ## empirical th.
        if(p != alpha) th1     <- sortedsumaalpha[ max(1,floor( (n)/bu[i]^thr2[1]))+1 ]                ## p-norm th.
        th2     <- sortedsuma[ max(1,floor( (n)/bu[i]^thr2[1]))+1 ]                     ## p-norm th.
        th3     <- sortedmaxi[ max(1,floor( (n)/bu[i]^thr2[1]))+1 ]                     ## supremom norm th.
        ## index
        #ind1     <- which(suma > th1^(alpha) )
        if(p != alpha) ind11    <- which(sumaalpha > th1)
        ind1     <- which(suma > th2 )                                         ## choosing clusters
        ind2     <- which(maxi > th3 )                                         ## Careful not to use >=. 
        # functions c(1)
        
        g1  <- function(vector)      (sum(vector^alpha)/(sum(vector^p)^alpha))
        g3  <- function(vector, p)      (sum(vector^p)^alpha/(sum(vector^alpha)))
        g2  <- function(vector,th3)     ( sum(vector) > th3 )*1

        estim1       <-  mean( sapply(ind1, function(k)    g1(path[((k-1)*bu[i] +1):(k*bu[i])])  )) 
        estim2       <-  sum(  sapply(1:n2, function(k)    g2(path[((k-1)*bu[i] +1):(k*bu[i])] , th3 )  )) 
        estim3       <-  sum(  sapply(1:n2, function(k)    g3(path[((k-1)*bu[i] +1):(k*bu[i])] , p )  )) 
        
        # functions theta
        #g1 <- function(vector, p)      (max(vector^p)/sum(vector^p))
        #g2 <- function(vector,th3)     sum( vector > th3)
        #estim1       <-  mean( sapply(ind1, function(k)    g1(path[((k-1)*bu[i] +1):(k*bu[i])] , p)  )) 
        #estim2      <-   mean(  sapply(ind2, function(k)    g2(path[((k-1)*bu[i] +1):(k*bu[i])] , th3 )  )) 
        
      
        estim2       <-   estim3/estim2 # 1/estim2 #
        if(j==1){
          infor1       <- rbind( infor1, c(estim1, 1, th2, bu[i], thr2[1], mode=j))
          infor1       <- rbind( infor1, c(estim2, 2, th3, bu[i], thr2[1], mode=j))
        }
        else{
          infor2       <- rbind( infor2, c(estim1, 1, th2, bu[i], thr2[1], mode=j))
          infor2       <- rbind( infor2, c(estim2, 2, th3, bu[i], thr2[1], mode=j))
        } 
      }
  }
}

names(infor2) <- names(infor1) <-  c("DATA", "TYPE", "zn", "BL", "TH", "MODE")


g1 <- niceplot(infor1,estim=c1[1],title)
g2 <- niceplot(infor2,estim=c1[2],title)
grid.arrange(g1,g2,nrow=1)
#ggplot(infor1, aes( y = DATA,  type=as.factor(BL), fill = as.factor(TYPE) ))+ geom_boxplot() + geom_hline(yintercept =c1) +ylim(0,5)

#load(infor1,infor2,title,phi0,alpha0,bu,n, file="0112C18000.Rdata")
#load(infor1,infor2,title,phi0,alpha0,bu,n, file="0112C13000.Rdata")
#load(infor1,infor2,title,phi0,alpha0,bu,n, file="0112C11000.Rdata")

#load(infor1,infor2,title,phi0,alpha0,bu,n, file="0112thet8000.Rdata")
#load(infor1,infor2,title,phi0,alpha0,bu,n, file="0112thet3000.Rdata")
#save(infor1,infor2,title,phi0,alpha0,bu,n, file="0112thet1000.Rdata")

#thet <- c1
niceplot <- function(infor,estim, title){
  gg<- ggplot(infor, aes( y = 1/DATA,  type=as.factor(BL), fill = as.factor(TYPE) ))+ 
    #fill=as.factor(MODE), colour = as.factor(MODE)) +
    geom_boxplot(#"#4271AE",
      colour = "#4271AE", alpha = 1,outlier.colour = "azure2", outlier.shape = 1, outlier.size = 0.1, outlier.alpha = 1)+ #+geom_point(alpha=0.1,pch=16,size=0.5)+
      scale_fill_manual(values = c("lightsteelblue3",'white'))
    #geom_hline(yintercept=thet, color = "darkgrey", alpha = 0.5,lty=2) 
  gg <-  addingtitles(gg,"block length",title,estim)
  return(gg)
}
addingtitles <- function(gg,title,title2=" ",estim){
  gg <- gg +
    geom_hline(yintercept=estim,lty=2, col="blue", alpha=1) +
    scale_y_continuous(name =" " , labels = c(0,round(estim,3),1), breaks = c(0,round(estim,3),1), limits=c(0,5)
                       )+ 
    scale_x_continuous(name = title, labels=sapply(1:6, function(k) 2^k), breaks=seq(-0.32,0.32,0.125) )+
    theme_bw()+              
    ggtitle(title2)+
    theme(#panel.grid.major = element_line(colour = "#d3d3d3"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size = 10 ,face = "bold", hjust=0),
      text=element_text(),
      axis.title = element_text(size=50),
      axis.title.y= element_text(size=24,angle=0),
      axis.title.x= element_text(size=10,angle=0),
      axis.text.x = element_text(colour="black", size = 10, angle = 15 ),
      axis.text.y = element_text(colour="black", size = 10),
      axis.line =   element_line(size=0.1, colour = "black"),
      legend.position="none")
  
  return(gg)
}  

