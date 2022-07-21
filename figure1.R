
## read the data here to get the date range for each range  
library(network)
library(chron)

setwd("/Users/timtsang/SPH Dropbox/Tim Tsang/")
## use Jul 1, 2008 as baseline

###########################################
## read the raw data

sero_y0 <- read.csv("Shared kiddivax analysis/data/cumulative incidence data/2009 Kiddivax Pilot/2012_01_09_sero.csv",as.is = T)

sero_y1 <- read.csv("Shared kiddivax analysis/data/cumulative incidence data/2010 Kiddivax Main/2012_12_04_serology.csv",as.is=T)

sero_y2 <- read.csv("Shared kiddivax analysis/data/cumulative incidence data/2011 Kiddivax Follow Up/2012_09_20_serology.csv",as.is=T)

sero_y3 <- read.csv("Shared kiddivax analysis/data/cumulative incidence data/2012 Kiddivax Follow Up/2013_08_23_serology.csv",as.is=T)

sero_y4 <- read.csv("Shared kiddivax analysis/data/cumulative incidence data/2014 Kiddivax Follow Up/2016_07_12_kiddivax_yr4to5_serology.csv",as.is=T)

#R0
(min(dates(sero_y0$date.mids,format="d/m/y"),na.rm=T) - dates("12/31/2008"))/7
(max(dates(sero_y0$date.mids,format="d/m/y"),na.rm=T) - dates("12/31/2008"))/7

#R1
(min(dates(sero_y1$prevax,format="d/m/y"),na.rm=T) - dates("12/31/2008"))/7
(max(dates(sero_y1$prevax,format="d/m/y"),na.rm=T) - dates("12/31/2008"))/7

#R2
(min(dates(sero_y1$mid.season,format="d/m/y"),na.rm=T) - dates("12/31/2008"))/7
(max(dates(sero_y1$mid.season,format="d/m/y"),na.rm=T) - dates("12/31/2008"))/7

#R3
(min(dates(sero_y1$post.season,format="d/m/y"),na.rm=T) - dates("12/31/2008"))/7
(max(dates(sero_y1$post.season,format="d/m/y"),na.rm=T) - dates("12/31/2008"))/7

#R4
(min(dates(sero_y2$mid.season,format="d/m/y"),na.rm=T) - dates("12/31/2008"))/7
(max(dates(sero_y2$mid.season,format="d/m/y"),na.rm=T) - dates("12/31/2008"))/7 

#R5
(min(dates(sero_y2$post.season,format="d/m/y"),na.rm=T) - dates("12/31/2008"))/7
(max(dates(sero_y2$post.season,format="d/m/y"),na.rm=T) - dates("12/31/2008"))/7

#R6
(min(dates(sero_y3$mid.season,format="d/m/y"),na.rm=T) - dates("12/31/2008"))/7
(max(dates(sero_y3$mid.season,format="d/m/y"),na.rm=T) - dates("12/31/2008"))/7

#R7
(min(dates(sero_y3$post.season,format="d/m/y"),na.rm=T) - dates("12/31/2008"))/7
(max(dates(sero_y3$post.season,format="d/m/y"),na.rm=T) - dates("12/31/2008"))/7

#R8
(min(dates(sero_y4$date.y4.mid,format="d/m/y"),na.rm=T) - dates("12/31/2008"))/7
(max(dates(sero_y4$date.y4.mid,format="d/m/y"),na.rm=T) - dates("12/31/2008"))/7


#R9
(min(dates(sero_y4$date.y5.pre,format="d/m/y"),na.rm=T) - dates("12/31/2008"))/7
(max(dates(sero_y4$date.y5.pre,format="d/m/y"),na.rm=T) - dates("12/31/2008"))/7


proxy <- read.csv("Shared kiddivax CMI/protection/data/ILILAB_wei_2019.csv")

proxy$GPrate <- proxy$GPratePer1000/1000
proxy$H1N1.prox <- proxy$A.H1N1/proxy$Specimens.tested*proxy$GPrate
proxy$H3N2.prox <- proxy$A.H3N2/proxy$Specimens.tested*proxy$GPrate
proxy$H1N1pdm.prox <- proxy$A.H1N1pdm/proxy$Specimens.tested*proxy$GPrate
proxy$FluB.prox <- proxy$B/proxy$Specimens.tested*proxy$GPrate
proxy$index <- 1:nrow(proxy)

proxy <- proxy[1:313,]  # from May 2010
#proxy$H3N2.prox[proxy$index>=261] <- 0  # set to zero to avoid confusion
#proxy$H1N1pdm.prox[proxy$index>=261] <- 0
#proxy$FluB.prox[proxy$index>=261] <- 0

#
# plot
#


pdf("Shared kiddivax CMI/protection/figure1_new.pdf",width=10, height=4)
layout(matrix(1))
par(mar=c(2,4,0,2))
plot(NA,xlim=c(1,313),ylim=c(0,0.096),ylab="",xlab="",axes=F)

mtext("Influenza activity",side=2,at=0.035,line = 3)
polygon(proxy$index[c(1,1:nrow(proxy),nrow(proxy),1)],c(0,proxy$H3N2.prox,0,0),
        col=rgb(0,0,0.9,0.7),border=F)  
polygon(proxy$index[c(1,1:nrow(proxy),nrow(proxy),1)],
        c(0,proxy$H1N1pdm.prox,0,0),col=rgb(1,0,0,0.7),border=F)
#polygon(proxy$index[c(1,1:nrow(proxy),nrow(proxy),1)],
#        c(0,proxy$FluB.prox,0,0),col=rgb(0,0.9,0,0.5),border=F)

for (i in 1:10){
lines(proxy$index[1:nrow(proxy)],proxy$FluB.prox,col="green")
}
axis(1,at=c(1,cumsum(table(proxy$Year))),labels=NA,tck=0)
axis(1,at=c(1,cumsum(table(proxy$Year)[1:6])),labels=NA)
#mtext(2009,side=1,at=34,line=0.5)
for(i in 1:6){mtext(2008+i,side=1,at=cumsum(table(proxy$Year))[i]-26,line=0.5)}

axis(2,at=0:7*0.01,las=1)

legend(265,0.07,legend=c("H3N2","pH1N1"),fill=c(rgb(0,0,0.9,0.5),rgb(1,0,0,0.5)),border=NA,bty="n")
legend(265,0.055,legend=c("HAI only","HAI+CD4/CD8"),fill=c(rgb(0,0,1,0.1),rgb(1,0,1,0.1)),border=NA,bty="n")



epidlist <- list(c(13,17),c(34,59),c(67,71),c(81,101),c(119,124),c(143,153),c(173,176),c(196,204),c(225,228),c(250,257))
roundpos2 <- 0.07
for (i in 1:10){
        paravalue <- epidlist[[i]]
        colset <- 0
        if (i %in% c(1,3,5,7)){
        colset <- 1        
        }
        polygon(c(paravalue[1],paravalue[1],paravalue[2],paravalue[2]),
                c(roundpos2,0.0,0.0,0.07),col=rgb(colset,0,1,0.1),border=F)      
        text(sum(paravalue)/2, 0.065,paste("R",i,sep=""),cex=1)
     #  text(sum(paravalue)/2, 0.06,"HAI",cex=0.8)
        if (colset==1){
    #    text(sum(paravalue)/2, 0.056,"+",cex=0.8)
    #    text(sum(paravalue)/2, 0.052,"CD4/CD8",cex=0.8)
        }
}


adj <- 0.005
rep1 <- 0.078
rep3 <- 0.086
rep2 <- 0.095
adj3 <- 0.003
adj2 <- 0.012
lines(c(13,13,17,17),c(0.0785,0.08,0.08,0.0785)-adj)
text(sum(13+17)/2, rep1,"Baseline CD4/CD8",cex=0.8)
lines(c(15,15,54,54),c(0.0785,0.08,0.08,0.0785)+adj2)
text(sum(15+54)/2, rep2,"Protection in epidemic 1",cex=0.8)

lines(c(67,67,71,71),c(0.0785,0.08,0.08,0.0785)-adj)
text(sum(67+71)/2, rep1,"Baseline CD4/CD8",cex=0.8)
lines(c(69,69,115,115),c(0.0785,0.08,0.08,0.0785)+adj3)
text(sum(69+115)/2, rep3,"Protection in epidemic 2/3",cex=0.8)

lines(c(119,119,124,124),c(0.0785,0.08,0.08,0.0785)-adj)
text(sum(119+124)/2, rep1,"Baseline CD4/CD8",cex=0.8)
lines(c(121.5,121.5,185,185),c(0.0785,0.08,0.08,0.0785)+adj2)
text(sum(121.5+185)/2, rep2,"Protection in epidemic 4",cex=0.8)

lines(c(173,173,176,176),c(0.0785,0.08,0.08,0.0785)-adj)
text(sum(173+176)/2, rep1,"Baseline CD4/CD8",cex=0.8)
lines(c(174.5,174.5,242,242),c(0.0785,0.08,0.08,0.0785)+adj3)
text(sum(174.5+242)/2, rep3,"Protection in epidemic 5/6",cex=0.8)


points(c(39,90,109,178,219,244),c(0.071,0.023,0.023,0.028,0.007,0.0075),pch=1,cex=2.5)  # circle
points(c(39,90,109,178,219,244),c(0.071,0.023,0.023,0.028,0.007,0.0075),pch=c("1","2","3","4","5","6"),
       cex=1)  # number

dev.off()
#
# End of script
#