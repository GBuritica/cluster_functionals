#######################################################################
#######################################################################
#######################################################################
#######################################################################
##### Gloria Buriticá
####  Code for producting plots
####  
#### 
#######################################################################
niceplot <- function(infor,estim, title){
  gg<- ggplot(infor, aes( y = DATA,  type=as.factor(BL), fill = as.factor(TYPE) ))+ 
    #fill=as.factor(MODE), colour = as.factor(MODE)) +
    geom_boxplot(#"#4271AE",
      colour = "#4271AE", alpha = 1,outlier.colour = "azure2", outlier.shape = 1, outlier.size = 0.1, outlier.alpha = 1)+ #+geom_point(alpha=0.1,pch=16,size=0.5)+
    scale_fill_manual(values = c('white',"lightsteelblue3","gray"))
  #geom_hline(yintercept=thet, color = "darkgrey", alpha = 0.5,lty=2) 
  gg <-  addingtitles(gg,"block length",title,estim)
  return(gg)
}
addingtitles <- function(gg,title,title2=" ",estim){
  gg <- gg +
    geom_hline(yintercept=estim,lty=2, col="blue", alpha=1) +
    scale_y_continuous(name =" " , labels = c(0,round(estim,3),1), breaks = c(0,round(estim,3),1), limits=c(0,3.5)
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
