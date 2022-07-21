
setwd("/Users/timtsang/SPH Dropbox/Tim Tsang/")

data1 <- read.csv("Shared kiddivax CMI/protection/data/2020_08_30_data.csv")
data1 <- data1[data1$CMI==1&data1$vac==0,]



## compute the bootstrap mean
bootCI <- function(temp2,input1,input2){

  rec <- rep(NA,1000)
  for (i in 1:1000){
    temp <- temp2[sample(1:nrow(temp2),replace=T),]
    rec[i] <-  cor(temp[,input1],temp[,input2],method="spearman",use="pairwise")
  }  
  return(c(mean(rec),quantile(rec,c(0.025,0.975))))
}


# pH1 and h3
# spearman correlation

output <- matrix(NA,3,12)

for (i in 1:2){
temp <- data1[data1$type==i-1,]

output[1,1:3+3*(i-1)] <- bootCI(temp,"AT1","CD4.sH1")
output[2,1:3+3*(i-1)] <- bootCI(temp,"AT1","CD4.pH1")
output[3,1:3+3*(i-1)] <- bootCI(temp,"AT1","CD4.sH3")

output[1,6+1:3+3*(i-1)] <- bootCI(temp,"AT1","CD8.sH1")
output[2,6+1:3+3*(i-1)] <- bootCI(temp,"AT1","CD8.pH1")
output[3,6+1:3+3*(i-1)] <- bootCI(temp,"AT1","CD8.sH3")

}







rowname <- c("HAI titer against pH1N1"," CD4/CD8 against sH1N1"," CD4/CD8 against pH1N1"," CD4/CD8 against H3N2",
             "HAI titer against H3N2"," CD4/CD8 against sH1N1"," CD4/CD8 against pH1N1"," CD4/CD8 against H3N2")

NAvec <- rep(NA,6)
a2 <- rbind(NAvec,cbind(output[,1:3],output[,4:6]),NAvec,cbind(output[,1:3+6],output[,4:6+6]))

a2round <- format(round(a2,digits=2),nsmall=2)
a1 <- cbind(paste(a2round[,1]," (",a2round[,2],", ",a2round[,3],")",sep=""),paste(a2round[,4]," (",a2round[,5],", ",a2round[,6],")",sep=""))
a1[grepl("NA",a1)] <- NA

pdf("Shared kiddivax CMI/protection/figureS4.pdf",width=10, height=5)


par(mar=c(0,0,0,0))

plot(0,0,xlab="",ylab="", main="", axes=F, xlim=c(-5.5,10.5), ylim=c(3.5,15.5),type="n")


cc <- 12 - (1:8-0.5)
left <- 15
panelb <- 6.5
textadj <- 1  

for (i in 5:45){
  if (i%%2==0){
    adjfac <- 0.08
    #   polygon(c(-18,-18,6.5,6.5),c(i+1+ adjfac,i+ adjfac,i+ adjfac,i+1+ adjfac),col=rgb(0,0,0,0.05),border=F)  
  }
}

for (i in 1:8){
  boldind <- c(1,5)
  if (1*(i%in%boldind)){
    text(-21+left,cc[i]+0.08,rowname[i],las=1,pos=4,font=2) 
  }
  else{
    text(-21+left,cc[i]+0.08,rowname[i],las=1,pos=4,font=1) 
  }
}


for (sym in 0:1){
  
#  for (u in 1:10){
#    arrows(-0.3+sym*panelb,2.5,-3+sym*panelb,2.5,col="black",length=0.05)
#    arrows(0.3+sym*panelb,2.5,3+sym*panelb,2.5,col="black",length=0.05)
#  }
  
  #text(-1.8+sym*panelb,1.5,"Higher protection") 
  #text(1.7+sym*panelb,1.5,"Lower protection")  
  
  
  axis(1,at=-1:1+sym*panelb,labels=c(-0.25,0,0.25),cex.axis=1,pos=4.3)
  #axis(2,at=c(-1.5,cc,17),labels=NA,las=1, pos=-2.2)
  
  plotvector <- a2[,1+3*sym]*4
  #plotvector[plotvector< -2] <- NA
  points(plotvector+sym*panelb,cc+0.08,pch=16,col="black",cex=0.8)
  
  
  for (i in 1:8){
    #+sym*panelb
    lines(c(max(-2,a2[i,2+3*sym]*4+sym*panelb),min(1,a2[i,3+3*sym]*4)+sym*panelb),rep(cc[i]+0.08,2),col="black")  
  
    text(4+sym*panelb-textadj,cc[i]+0.08,a1[i,1+sym],las=1)
    
  }
  
  lines(c(0,0)+sym*panelb,c(4.3,11.5),lty=2)
  
  
  text(0+sym*panelb,14.5,"Spearman correlation")
  text(0+sym*panelb,13.5,"(95% CI)")
  text(0,15.5,"CD4",cex=1.2)
  text(0+panelb,15.5,"CD8",cex=1.2)
  text(-21+left,12.5,"Correlation between pre-epidemic HAI titers and CD4/CD8 T cell response",las=1,pos=4,font=1) 
  
  text(-2.5,15.5,"A",cex=1.5)
  text(5,15.5,"B",cex=1.5)
  
}


dev.off()

