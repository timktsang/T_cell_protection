
setwd("/Users/timtsang/SPH Dropbox/Tim Tsang/")

data1 <- read.csv("Shared kiddivax CMI/protection/data/2020_08_30_data.csv")
data1 <- data1[data1$CMI==1&data1$vac==0,]

## season1 to 6, the difference of pbs, sh1, ph1 sh3


## compute the bootstrap mean
bootCI <- function(input){
  input <- input[!is.na(input)]
  rec <- rep(NA,1000)
  for (i in 1:1000){
    rec[i] <- mean(input[sample(1:length(input),replace = T)])
  }  
  return(c(mean(input),quantile(rec,c(0.025,0.975))))
}


plotdata <- matrix(NA,12,24)
pvalue <- matrix(NA,12,4)
pvalue2 <- matrix(NA,12,4)
for (i in 1:6){
  for (j in 0:3){
    # CD4 infected
    if (sum(!is.na(data1[data1$season==i-1&data1$inf==1,7+2*j]))>0){
      a1 <- t.test(data1[data1$season==i-1&data1$inf==1,7+2*j])
      a1 <- bootCI(data1[data1$season==i-1&data1$inf==1,7+2*j])
      plotdata[i,1:3+6*j] <- a1
    }
    # CD4 not infected
    if (sum(!is.na(data1[data1$season==i-1&data1$inf==0,7+2*j]))>0){
      a1 <- t.test(data1[data1$season==i-1&data1$inf==0,7+2*j])
      a1 <- bootCI(data1[data1$season==i-1&data1$inf==0,7+2*j])
      plotdata[i,1:3+3+6*j] <-  a1
    }
    # pvalue
    if (sum(!is.na(data1[data1$season==i-1&data1$inf==1,7+2*j]))>0&sum(!is.na(data1[data1$season==i-1&data1$inf==0,7+2*j]))>0){
      a1 <- t.test(data1[data1$season==i-1&data1$inf==0,7+2*j],data1[data1$season==i-1&data1$inf==1,7+2*j]) 
      pvalue[i,j+1] <- a1$p.value
      a1 <- wilcox.test(data1[data1$season==i-1&data1$inf==0,7+2*j],data1[data1$season==i-1&data1$inf==1,7+2*j]) 
      pvalue2[i,j+1] <- a1$p.value
    }
    
    # CD8 infected
    if (sum(!is.na(data1[data1$season==i-1&data1$inf==1,8+2*j]))>0){
      a1 <- t.test(data1[data1$season==i-1&data1$inf==1,8+2*j])
      a1 <- bootCI(data1[data1$season==i-1&data1$inf==1,8+2*j])
      plotdata[i+6,1:3+6*j] <- a1
    }
    # CD8 not infected
    if (sum(!is.na(data1[data1$season==i-1&data1$inf==0,8+2*j]))>0){
      a1 <- t.test(data1[data1$season==i-1&data1$inf==0,8+2*j])
      a1 <- bootCI(data1[data1$season==i-1&data1$inf==0,8+2*j])
      plotdata[i+6,1:3+3+6*j] <- a1
    }
    # pvalue
    if (sum(!is.na(data1[data1$season==i-1&data1$inf==1,8+2*j]))>0&sum(!is.na(data1[data1$season==i-1&data1$inf==0,8+2*j]))>0){
      a1 <- t.test(data1[data1$season==i-1&data1$inf==0,8+2*j],data1[data1$season==i-1&data1$inf==1,8+2*j]) 
      pvalue[i+6,j+1] <- a1$p.value
      a1 <- wilcox.test(data1[data1$season==i-1&data1$inf==0,8+2*j],data1[data1$season==i-1&data1$inf==1,8+2*j]) 
      pvalue2[i+6,j+1] <- a1$p.value
    }
    
  }  
}


plotfunction <- function(data1,titlename,pa,whichrow,trun){
  
  dif1 <- -0.2
  dif2 <- 0.2
  
  alldata <- data1[data1$season==whichrow-7,8+2*(0:3)]  
  
  get <- max(alldata,na.rm=T)
  
  ############################################################################################################################################
  ## panel A
  par(mar=c(3,4,2,1))
  
  plot(0,0,xlab="",ylab="", main="", axes=F, xlim=c(0,4), ylim=c(0, trun*1.05),type="n")
  
  axis(1,at=c(-1,0.5,1.5,2.5,3.5,5),labels=c(NA,"PBS","sH1","pH1","H3",NA),cex.axis=1)
  axis(2,at=0:5*trun/5, las=1, pos=-0.05)
  
  
  points(1:4-0.5+dif1,plotdata[whichrow,0:3*6+1],pch=16,col="black")
  points(1:4-0.5+dif2,plotdata[whichrow,0:3*6+4],pch=17,col="black")
  for (i in 1:4){
    lines(rep(i-0.5+dif1,2),plotdata[whichrow,(i-1)*6+2:3]) 
    lines(rep(i-0.5+dif2,2),plotdata[whichrow,(i-1)*6+3+2:3])  
  }
  
  jpara1 <- 0.15
  jpara <- 0.25
  sd <- 0.1
  for (j in 1:4){
    sd2 <- sd
    if (j==1){
      sd2 <- 0.03
    }  
    plotvec <- pmin(trun*1.00,data1[data1$season==whichrow-7&data1$inf==1,8+2*(j-1)])
    points( jitter(rep(j-0.5+dif1,length(plotvec)),amount=jpara1), plotvec,col=rgb(1,0,0,sd2),pch=16)
    plotvec <- pmin(trun*1.00,data1[data1$season==whichrow-7&data1$inf==0,8+2*(j-1)])
    points( jitter(rep(j-0.5+dif2,length(plotvec)),amount=jpara1), plotvec,col=rgb(0,0,1,sd2),pch=17)
  
    if (pvalue2[whichrow,j]<0.05&!is.na(pvalue2[whichrow,j])){
     
      text(j-0.5,trun*1.08,"*",cex=2.3)
      lines(c(j-0.5+dif1,j-0.5+dif2),rep(trun*1.05,2),lty=1)
      lines(rep(j-0.5+dif1,2),trun*c(1.05,1.02),lty=1)
      lines(rep(j-0.5+dif2,2),trun*c(1.05,1.02),lty=1)
      
    }
    }
  
  lines(c(-1,5),rep(trun,2),lty=2)
  
  mtext(expression(paste("CD69"^"+","IFN-",gamma^"+"," in CD8 (%)",sep="")),side=2,line=2.25)
  
  
  mtext(titlename,side=3,line=0.5)
  title(main=pa, adj=0)

  if (whichrow==6+6){
    legend(2,1.9,c("Infection","Non-infection"),pch=16:17,lty=1,cex=1,bty="n",col=c("red","blue"))
  }  
}

pdf("Shared kiddivax CMI/protection/figureS3.pdf",width=9, height=6)
layout(matrix( 1:6, nrow=2,byrow=T))
plotfunction(data1,"Epidemic 1: pH1N1","A",7,2)
plotfunction(data1,"Epidemic 2: H3N2","B",8,2)
plotfunction(data1,"Epidemic 3: pH1N1","C",9,2)
plotfunction(data1,"Epidemic 4: H3N2","D",10,5)
plotfunction(data1,"Epidemic 5: pH1N1","E",11,2)
plotfunction(data1,"Epidemic 6: H3N2","F",12,2)

dev.off()

## here also the correlation between CD8 and strain
library(ppcor)

a1 <- lm(data1$CD8.sH1~as.factor(data1$season))
tt1 <- a1$residuals
a2 <- lm(data1$CD8.pH1~as.factor(data1$season))

cor(data1$CD8.sH1,data1$CD8.pH1,use="pairwise.complete.obs")
cor(data1$CD8.sH3,data1$CD8.pH1,use="pairwise.complete.obs")

a1 <- lm(data1$CD8.sH1~data1$CD8.sH3 + data1$CD8.pH1 + as.factor(data1$season))
summary(a1)

