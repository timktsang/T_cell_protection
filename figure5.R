library(lme4)
library(geepack)

setwd("/Users/timtsang/SPH Dropbox/Tim Tsang/")

data1 <- read.csv("Shared kiddivax CMI/protection/data/2020_08_30_data.csv")
data1$uid <- as.factor(paste(data1$hhID,data1$member,sep="-"))
data1 <- data1[order(data1$uid),]

data1 <- data1[data1$CMI==1&data1$vac==0,]
data1$agegp <- 1*(data1$age>=18)

data1$season <- as.factor(data1$season)
data1$HAIpos <- 1*(data1$AT1>=3)

## change to log scale
for (i in 9:14){
  mimv <- min(data1[data1[,i]>0,i],na.rm = T)
  print(  mimv )
  data1[,i] <- log2(pmax(data1[,i],mimv/2))
}

for (i in 0:5){
  for (j in 9:14){
#  temp1 <- data1[data1$season==i,j] 
#  temp1sd <- (temp1-mean(temp1,na.rm = T))/sd(temp1,na.rm = T)  
#  data1[data1$season==i,j] <- temp1sd 
}
}

# by epidemic
plotdata <- matrix(NA,12,9)
multipler <- 1


for (i in 1:2){
tempdata <- data1[data1$type==i-1,]

# multi

if (sum(tempdata$AT1>=3)>0){
tempdata1 <- tempdata[!is.na(tempdata$CD4.sH1)&!is.na(tempdata$CD8.sH1),]
tempdata1$season2 <- as.factor(as.character(tempdata1$season))
  
#  a2 <- glm(inf ~ I(CD4.sH1*multipler) + I(CD8.sH1*multipler) + agegp + season + HAIpos  ,data=tempdata,family="binomial")
  
#  a1 <- glmer(as.factor(inf) ~ (1|uid) + I(CD4.sH1*multipler) + I(CD8.sH1*multipler) + as.factor(agegp) + as.factor(season) + as.factor(AT1>=3)  ,data=tempdata,family="binomial")
  a1 <- geeglm(inf ~ I(CD4.sH1*multipler) + I(CD8.sH1*multipler) + agegp + season2 + HAIpos  ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
  
  a1 <- geeglm(inf ~ I(CD4.sH1*multipler) + I(CD8.sH1*multipler) + agegp + season2 + HAIpos  ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
  
  }
else{
  tempdata1 <- tempdata[!is.na(tempdata$CD4.sH1)&!is.na(tempdata$CD8.sH1),]
  tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
#  a1 <- glm(as.factor(inf) ~ I(CD4.sH1*multipler) + I(CD8.sH1*multipler) + as.factor(agegp) + as.factor(season)  ,data=tempdata,family="binomial")  
  a1 <- geeglm(inf ~ I(CD4.sH1*multipler) + I(CD8.sH1*multipler) + agegp + season2  ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
  
  }
a2 <- coef(summary(a1))
plotdata[i+4,1:3] <- exp(a2[2,1] + c(0,-1.96,1.96)*a2[2,2])
plotdata[i+6,1:3] <- exp(a2[3,1] + c(0,-1.96,1.96)*a2[3,2])

if (sum(!is.na(tempdata$CD4.pH1))>0){
  if (sum(tempdata$AT1>=3)>0&sum(!is.na(tempdata$CD4.pH1))>0){
    
    tempdata1 <- tempdata[!is.na(tempdata$CD4.pH1)&!is.na(tempdata$CD8.pH1),]
    tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
    
   # a1 <- glm(as.factor(inf) ~ I(CD4.pH1*multipler) + I(CD8.pH1*multipler) + as.factor(agegp) + as.factor(season) + as.factor(AT1>=3)  ,data=tempdata,family="binomial")
  
    a1 <- geeglm(inf ~ I(CD4.pH1*multipler) + I(CD8.pH1*multipler) + agegp + season2 + HAIpos ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
    
    }
  if (sum(tempdata$AT1>=3)==0&sum(!is.na(tempdata$CD4.pH1))>0){
    
    tempdata1 <- tempdata[!is.na(tempdata$CD4.pH1)&!is.na(tempdata$CD8.pH1),]
    tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
    
    
  #  a1 <- glm(as.factor(inf) ~ I(CD4.pH1*multipler) + I(CD8.pH1*multipler) + as.factor(agegp) + as.factor(season)  ,data=tempdata,family="binomial")  
  
    a1 <- geeglm(inf ~ I(CD4.pH1*multipler) + I(CD8.pH1*multipler) + agegp + season2  ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
    
    }
  a2 <- coef(summary(a1))
  plotdata[i+4,1:3+3] <- exp(a2[2,1] + c(0,-1.96,1.96)*a2[2,2])
  plotdata[i+6,1:3+3] <- exp(a2[3,1] + c(0,-1.96,1.96)*a2[3,2])
}


if (sum(tempdata$AT1>=3)>0){
  
  tempdata1 <- tempdata[!is.na(tempdata$CD4.sH3)&!is.na(tempdata$CD8.sH3),]
  tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
  
  
  
  #a1 <- glm(as.factor(inf) ~ I(CD4.sH3*multipler) + I(CD8.sH3*multipler) + as.factor(agegp) + as.factor(season) + as.factor(AT1>=3)  ,data=tempdata,family="binomial")

  a1 <- geeglm(inf ~ I(CD4.sH3*multipler) + I(CD8.sH3*multipler) + agegp + season2 + HAIpos ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
  
  }
else{
  
  
  tempdata1 <- tempdata[!is.na(tempdata$CD4.sH3)&!is.na(tempdata$CD8.sH3),]
  tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
  
  
  #a1 <- glm(as.factor(inf) ~ I(CD4.sH3*multipler) + I(CD8.sH3*multipler) + as.factor(agegp) + as.factor(season)  ,data=tempdata,family="binomial")  

  a1 <- geeglm(inf ~ I(CD4.sH3*multipler) + I(CD8.sH3*multipler) + agegp + season2  ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
  
  
  }
a2 <- coef(summary(a1))
plotdata[i+4,1:3+6] <- exp(a2[2,1] + c(0,-1.96,1.96)*a2[2,2])
plotdata[i+6,1:3+6] <- exp(a2[3,1] + c(0,-1.96,1.96)*a2[3,2])

## uni
if (sum(tempdata$AT1>=3)>0){
  
  tempdata1 <- tempdata[!is.na(tempdata$CD4.sH1),]
  tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
  
  
#a1 <- glm(as.factor(inf) ~ I(CD4.sH1*multipler) + as.factor(agegp) + as.factor(season) + as.factor(AT1>=3)  ,data=tempdata,family="binomial")

#a2 <- glmer(as.factor(inf) ~ (1|uid) + I(CD4.sH1*multipler)  + as.factor(agegp) + as.factor(season) + as.factor(AT1>=3)  ,data=tempdata,family="binomial")

a1 <- geeglm(inf ~ I(CD4.sH1*multipler) + agegp + season2 + HAIpos ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")



}
else{
#a1 <- glm(as.factor(inf) ~ I(CD4.sH1*multipler) + as.factor(agegp) + as.factor(season)  ,data=tempdata,family="binomial")  

tempdata1 <- tempdata[!is.na(tempdata$CD4.sH1),]
tempdata1$season2 <- as.factor(as.character(tempdata1$season))  

a1 <- geeglm(inf ~ I(CD4.sH1*multipler) + agegp + season2  ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")


}
a2 <- coef(summary(a1))
plotdata[i,1:3] <- exp(a2[2,1] + c(0,-1.96,1.96)*a2[2,2])

if (sum(!is.na(tempdata$CD4.pH1))>0){
if (sum(tempdata$AT1>=3)>0&sum(!is.na(tempdata$CD4.pH1))>0){
  
  #a1 <- glm(as.factor(inf) ~ I(CD4.pH1*multipler) + as.factor(agegp) + as.factor(season) + as.factor(AT1>=3)  ,data=tempdata,family="binomial")

  
  tempdata1 <- tempdata[!is.na(tempdata$CD4.pH1),]
  tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
  
  a1 <- geeglm(inf ~ I(CD4.pH1*multipler) + agegp + season2 + HAIpos ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
  
  
  }
if (sum(tempdata$AT1>=3)==0&sum(!is.na(tempdata$CD4.pH1))>0){
  
  
  tempdata1 <- tempdata[!is.na(tempdata$CD4.pH1),]
  tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
  
  a1 <- geeglm(inf ~ I(CD4.pH1*multipler) + agegp + season2  ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
  
  
  #a1 <- glm(as.factor(inf) ~ I(CD4.pH1*multipler) + as.factor(agegp) + as.factor(season)  ,data=tempdata,family="binomial")  

  }
a2 <- coef(summary(a1))
plotdata[i,1:3+3] <- exp(a2[2,1] + c(0,-1.96,1.96)*a2[2,2])
}


if (sum(tempdata$AT1>=3)>0){
  
  tempdata1 <- tempdata[!is.na(tempdata$CD4.sH3),]
  tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
  
  a1 <- geeglm(inf ~ I(CD4.sH3*multipler) + agegp + season2 + HAIpos ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
  
  
  #a1 <- glm(as.factor(inf) ~ I(CD4.sH3*multipler) + as.factor(agegp) + as.factor(season) + as.factor(AT1>=3)  ,data=tempdata,family="binomial")

  }
else{
  
  tempdata1 <- tempdata[!is.na(tempdata$CD4.sH3),]
  tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
  
  a1 <- geeglm(inf ~ I(CD4.sH3*multipler) + agegp + season2  ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
  
  
 # a1 <- glm(as.factor(inf) ~ I(CD4.sH3*multipler) + as.factor(agegp) + as.factor(season)  ,data=tempdata,family="binomial")  

  }
a2 <- coef(summary(a1))
plotdata[i,1:3+6] <- exp(a2[2,1] + c(0,-1.96,1.96)*a2[2,2])

## CD8
if (sum(tempdata$AT1>=3)>0){
  
 # a1 <- glm(as.factor(inf) ~ I(CD8.sH1*multipler) + as.factor(agegp) + as.factor(season) + as.factor(AT1>=3)  ,data=tempdata,family="binomial")

  tempdata1 <- tempdata[!is.na(tempdata$CD8.sH1),]
  tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
  
  a1 <- geeglm(inf ~ I(CD8.sH1*multipler) + agegp + season2 + HAIpos  ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
  
  
  }
else{
  
  tempdata1 <- tempdata[!is.na(tempdata$CD8.sH1),]
  tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
  
  a1 <- geeglm(inf ~ I(CD8.sH1*multipler) + agegp + season2   ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
  
  
  
 # a1 <- glm(as.factor(inf) ~ I(CD8.sH1*multipler) + as.factor(agegp) + as.factor(season)  ,data=tempdata,family="binomial")  

  }
a2 <- coef(summary(a1))
plotdata[i+2,1:3] <- exp(a2[2,1] + c(0,-1.96,1.96)*a2[2,2])

if (sum(!is.na(tempdata$CD4.pH1))>0){
  if (sum(tempdata$AT1>=3)>0&sum(!is.na(tempdata$CD4.pH1))>0){
    
    #a1 <- glm(as.factor(inf) ~ I(CD8.pH1*multipler) + as.factor(agegp) + as.factor(season) + as.factor(AT1>=3)  ,data=tempdata,family="binomial")
    
    tempdata1 <- tempdata[!is.na(tempdata$CD8.pH1),]
    tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
    
    a1 <- geeglm(inf ~ I(CD8.pH1*multipler) + agegp + season2 + HAIpos  ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
    
    
    }
  if (sum(tempdata$AT1>=3)==0&sum(!is.na(tempdata$CD4.pH1))>0){
    
    
    tempdata1 <- tempdata[!is.na(tempdata$CD8.pH1),]
    tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
    
    a1 <- geeglm(inf ~ I(CD8.pH1*multipler) + agegp + season2   ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
    
    
    
    #a1 <- glm(as.factor(inf) ~ I(CD8.pH1*multipler) + as.factor(agegp) + as.factor(season)  ,data=tempdata,family="binomial")  
  
    }
  a2 <- coef(summary(a1))
  plotdata[i+2,1:3+3] <- exp(a2[2,1] + c(0,-1.96,1.96)*a2[2,2])
}


if (sum(tempdata$AT1>=3)>0){
  
 # a1 <- glm(as.factor(inf) ~ I(CD8.sH3*multipler) + as.factor(agegp) + as.factor(season) + as.factor(AT1>=3)  ,data=tempdata,family="binomial")

  
  tempdata1 <- tempdata[!is.na(tempdata$CD8.sH3),]
  tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
  
  a1 <- geeglm(inf ~ I(CD8.sH3*multipler) + agegp + season2 + HAIpos   ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
  
  
  
  }
else{
  
  #a1 <- glm(as.factor(inf) ~ I(CD8.sH3*multipler) + as.factor(agegp) + as.factor(season)  ,data=tempdata,family="binomial")  

  
  tempdata1 <- tempdata[!is.na(tempdata$CD8.sH3),]
  tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
  
  a1 <- geeglm(inf ~ I(CD8.sH3*multipler) + agegp + season2   ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
  
  
  }
a2 <- coef(summary(a1))
plotdata[i+2,1:3+6] <- exp(a2[2,1] + c(0,-1.96,1.96)*a2[2,2])
}

plotdataA <- plotdata


library(lme4)
library(geepack)

data1 <- read.csv("Shared kiddivax CMI/protection/data/2020_08_30_data.csv")
data1$uid <- as.factor(paste(data1$hhID,data1$member,sep="-"))
data1 <- data1[order(data1$uid),]

data1 <- data1[data1$CMI==1&data1$vac==0&data1$AT1<3,]
data1$agegp <- 1*(data1$age>=18)

data1$season <- as.factor(data1$season)
data1$HAIpos <- 1*(data1$AT1>=3)

## change to log scale
for (i in 9:14){
  mimv <- min(data1[data1[,i]>0,i],na.rm = T)
  print(  mimv )
  data1[,i] <- log2(pmax(data1[,i],mimv/2))
}

for (i in 0:5){
  for (j in 9:14){
    #  temp1 <- data1[data1$season==i,j] 
    #  temp1sd <- (temp1-mean(temp1,na.rm = T))/sd(temp1,na.rm = T)  
    #  data1[data1$season==i,j] <- temp1sd 
  }
}

# by epidemic
plotdata <- matrix(NA,12,9)
multipler <- 1


for (i in 1:2){
  tempdata <- data1[data1$type==i-1,]
  
  # multi
  
  if (sum(tempdata$AT1>=3)>0){
    tempdata1 <- tempdata[!is.na(tempdata$CD4.sH1)&!is.na(tempdata$CD8.sH1),]
    tempdata1$season2 <- as.factor(as.character(tempdata1$season))
    
    #  a2 <- glm(inf ~ I(CD4.sH1*multipler) + I(CD8.sH1*multipler) + agegp + season + HAIpos  ,data=tempdata,family="binomial")
    
    #  a1 <- glmer(as.factor(inf) ~ (1|uid) + I(CD4.sH1*multipler) + I(CD8.sH1*multipler) + as.factor(agegp) + as.factor(season) + as.factor(AT1>=3)  ,data=tempdata,family="binomial")
    a1 <- geeglm(inf ~ I(CD4.sH1*multipler) + I(CD8.sH1*multipler) + agegp + season2 + HAIpos  ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
    
  }
  else{
    tempdata1 <- tempdata[!is.na(tempdata$CD4.sH1)&!is.na(tempdata$CD8.sH1),]
    tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
    #  a1 <- glm(as.factor(inf) ~ I(CD4.sH1*multipler) + I(CD8.sH1*multipler) + as.factor(agegp) + as.factor(season)  ,data=tempdata,family="binomial")  
    a1 <- geeglm(inf ~ I(CD4.sH1*multipler) + I(CD8.sH1*multipler) + agegp + season2  ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
    
  }
  a2 <- coef(summary(a1))
  plotdata[i+4,1:3] <- exp(a2[2,1] + c(0,-1.96,1.96)*a2[2,2])
  plotdata[i+6,1:3] <- exp(a2[3,1] + c(0,-1.96,1.96)*a2[3,2])
  
  if (sum(!is.na(tempdata$CD4.pH1))>0){
    if (sum(tempdata$AT1>=3)>0&sum(!is.na(tempdata$CD4.pH1))>0){
      
      tempdata1 <- tempdata[!is.na(tempdata$CD4.pH1)&!is.na(tempdata$CD8.pH1),]
      tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
      
      # a1 <- glm(as.factor(inf) ~ I(CD4.pH1*multipler) + I(CD8.pH1*multipler) + as.factor(agegp) + as.factor(season) + as.factor(AT1>=3)  ,data=tempdata,family="binomial")
      
      a1 <- geeglm(inf ~ I(CD4.pH1*multipler) + I(CD8.pH1*multipler) + agegp + season2 + HAIpos ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
      
    }
    if (sum(tempdata$AT1>=3)==0&sum(!is.na(tempdata$CD4.pH1))>0){
      
      tempdata1 <- tempdata[!is.na(tempdata$CD4.pH1)&!is.na(tempdata$CD8.pH1),]
      tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
      
      
      #  a1 <- glm(as.factor(inf) ~ I(CD4.pH1*multipler) + I(CD8.pH1*multipler) + as.factor(agegp) + as.factor(season)  ,data=tempdata,family="binomial")  
      
      a1 <- geeglm(inf ~ I(CD4.pH1*multipler) + I(CD8.pH1*multipler) + agegp + season2  ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
      
    }
    a2 <- coef(summary(a1))
    plotdata[i+4,1:3+3] <- exp(a2[2,1] + c(0,-1.96,1.96)*a2[2,2])
    plotdata[i+6,1:3+3] <- exp(a2[3,1] + c(0,-1.96,1.96)*a2[3,2])
  }
  
  
  if (sum(tempdata$AT1>=3)>0){
    
    tempdata1 <- tempdata[!is.na(tempdata$CD4.sH3)&!is.na(tempdata$CD8.sH3),]
    tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
    
    
    
    #a1 <- glm(as.factor(inf) ~ I(CD4.sH3*multipler) + I(CD8.sH3*multipler) + as.factor(agegp) + as.factor(season) + as.factor(AT1>=3)  ,data=tempdata,family="binomial")
    
    a1 <- geeglm(inf ~ I(CD4.sH3*multipler) + I(CD8.sH3*multipler) + agegp + season2 + HAIpos ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
    
  }
  else{
    
    
    tempdata1 <- tempdata[!is.na(tempdata$CD4.sH3)&!is.na(tempdata$CD8.sH3),]
    tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
    
    
    #a1 <- glm(as.factor(inf) ~ I(CD4.sH3*multipler) + I(CD8.sH3*multipler) + as.factor(agegp) + as.factor(season)  ,data=tempdata,family="binomial")  
    
    a1 <- geeglm(inf ~ I(CD4.sH3*multipler) + I(CD8.sH3*multipler) + agegp + season2  ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
    
    
  }
  a2 <- coef(summary(a1))
  plotdata[i+4,1:3+6] <- exp(a2[2,1] + c(0,-1.96,1.96)*a2[2,2])
  plotdata[i+6,1:3+6] <- exp(a2[3,1] + c(0,-1.96,1.96)*a2[3,2])
  
  ## uni
  if (sum(tempdata$AT1>=3)>0){
    
    tempdata1 <- tempdata[!is.na(tempdata$CD4.sH1),]
    tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
    
    
    #a1 <- glm(as.factor(inf) ~ I(CD4.sH1*multipler) + as.factor(agegp) + as.factor(season) + as.factor(AT1>=3)  ,data=tempdata,family="binomial")
    
    #a2 <- glmer(as.factor(inf) ~ (1|uid) + I(CD4.sH1*multipler)  + as.factor(agegp) + as.factor(season) + as.factor(AT1>=3)  ,data=tempdata,family="binomial")
    
    a1 <- geeglm(inf ~ I(CD4.sH1*multipler) + agegp + season2 + HAIpos ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
    
    
    
  }
  else{
    #a1 <- glm(as.factor(inf) ~ I(CD4.sH1*multipler) + as.factor(agegp) + as.factor(season)  ,data=tempdata,family="binomial")  
    
    tempdata1 <- tempdata[!is.na(tempdata$CD4.sH1),]
    tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
    
    a1 <- geeglm(inf ~ I(CD4.sH1*multipler) + agegp + season2  ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
    
    
  }
  a2 <- coef(summary(a1))
  plotdata[i,1:3] <- exp(a2[2,1] + c(0,-1.96,1.96)*a2[2,2])
  
  if (sum(!is.na(tempdata$CD4.pH1))>0){
    if (sum(tempdata$AT1>=3)>0&sum(!is.na(tempdata$CD4.pH1))>0){
      
      #a1 <- glm(as.factor(inf) ~ I(CD4.pH1*multipler) + as.factor(agegp) + as.factor(season) + as.factor(AT1>=3)  ,data=tempdata,family="binomial")
      
      
      tempdata1 <- tempdata[!is.na(tempdata$CD4.pH1),]
      tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
      
      a1 <- geeglm(inf ~ I(CD4.pH1*multipler) + agegp + season2 + HAIpos ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
      
      
    }
    if (sum(tempdata$AT1>=3)==0&sum(!is.na(tempdata$CD4.pH1))>0){
      
      
      tempdata1 <- tempdata[!is.na(tempdata$CD4.pH1),]
      tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
      
      a1 <- geeglm(inf ~ I(CD4.pH1*multipler) + agegp + season2  ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
      
      
      #a1 <- glm(as.factor(inf) ~ I(CD4.pH1*multipler) + as.factor(agegp) + as.factor(season)  ,data=tempdata,family="binomial")  
      
    }
    a2 <- coef(summary(a1))
    plotdata[i,1:3+3] <- exp(a2[2,1] + c(0,-1.96,1.96)*a2[2,2])
  }
  
  
  if (sum(tempdata$AT1>=3)>0){
    
    tempdata1 <- tempdata[!is.na(tempdata$CD4.sH3),]
    tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
    
    a1 <- geeglm(inf ~ I(CD4.sH3*multipler) + agegp + season2 + HAIpos ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
    
    
    #a1 <- glm(as.factor(inf) ~ I(CD4.sH3*multipler) + as.factor(agegp) + as.factor(season) + as.factor(AT1>=3)  ,data=tempdata,family="binomial")
    
  }
  else{
    
    tempdata1 <- tempdata[!is.na(tempdata$CD4.sH3),]
    tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
    
    a1 <- geeglm(inf ~ I(CD4.sH3*multipler) + agegp + season2  ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
    
    
    # a1 <- glm(as.factor(inf) ~ I(CD4.sH3*multipler) + as.factor(agegp) + as.factor(season)  ,data=tempdata,family="binomial")  
    
  }
  a2 <- coef(summary(a1))
  plotdata[i,1:3+6] <- exp(a2[2,1] + c(0,-1.96,1.96)*a2[2,2])
  
  ## CD8
  if (sum(tempdata$AT1>=3)>0){
    
    # a1 <- glm(as.factor(inf) ~ I(CD8.sH1*multipler) + as.factor(agegp) + as.factor(season) + as.factor(AT1>=3)  ,data=tempdata,family="binomial")
    
    tempdata1 <- tempdata[!is.na(tempdata$CD8.sH1),]
    tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
    
    a1 <- geeglm(inf ~ I(CD8.sH1*multipler) + agegp + season2 + HAIpos  ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
    
    
  }
  else{
    
    tempdata1 <- tempdata[!is.na(tempdata$CD8.sH1),]
    tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
    
    a1 <- geeglm(inf ~ I(CD8.sH1*multipler) + agegp + season2   ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
    
    
    
    # a1 <- glm(as.factor(inf) ~ I(CD8.sH1*multipler) + as.factor(agegp) + as.factor(season)  ,data=tempdata,family="binomial")  
    
  }
  a2 <- coef(summary(a1))
  plotdata[i+2,1:3] <- exp(a2[2,1] + c(0,-1.96,1.96)*a2[2,2])
  
  if (sum(!is.na(tempdata$CD4.pH1))>0){
    if (sum(tempdata$AT1>=3)>0&sum(!is.na(tempdata$CD4.pH1))>0){
      
      #a1 <- glm(as.factor(inf) ~ I(CD8.pH1*multipler) + as.factor(agegp) + as.factor(season) + as.factor(AT1>=3)  ,data=tempdata,family="binomial")
      
      tempdata1 <- tempdata[!is.na(tempdata$CD8.pH1),]
      tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
      
      a1 <- geeglm(inf ~ I(CD8.pH1*multipler) + agegp + season2 + HAIpos  ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
      
      
    }
    if (sum(tempdata$AT1>=3)==0&sum(!is.na(tempdata$CD4.pH1))>0){
      
      
      tempdata1 <- tempdata[!is.na(tempdata$CD8.pH1),]
      tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
      
      a1 <- geeglm(inf ~ I(CD8.pH1*multipler) + agegp + season2   ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
      
      
      
      #a1 <- glm(as.factor(inf) ~ I(CD8.pH1*multipler) + as.factor(agegp) + as.factor(season)  ,data=tempdata,family="binomial")  
      
    }
    a2 <- coef(summary(a1))
    plotdata[i+2,1:3+3] <- exp(a2[2,1] + c(0,-1.96,1.96)*a2[2,2])
  }
  
  
  if (sum(tempdata$AT1>=3)>0){
    
    # a1 <- glm(as.factor(inf) ~ I(CD8.sH3*multipler) + as.factor(agegp) + as.factor(season) + as.factor(AT1>=3)  ,data=tempdata,family="binomial")
    
    
    tempdata1 <- tempdata[!is.na(tempdata$CD8.sH3),]
    tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
    
    a1 <- geeglm(inf ~ I(CD8.sH3*multipler) + agegp + season2 + HAIpos   ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
    
    
    
  }
  else{
    
    #a1 <- glm(as.factor(inf) ~ I(CD8.sH3*multipler) + as.factor(agegp) + as.factor(season)  ,data=tempdata,family="binomial")  
    
    
    tempdata1 <- tempdata[!is.na(tempdata$CD8.sH3),]
    tempdata1$season2 <- as.factor(as.character(tempdata1$season))  
    
    a1 <- geeglm(inf ~ I(CD8.sH3*multipler) + agegp + season2   ,data=tempdata1,id = hhID, waves = uid,family="binomial", corstr="exchangeable")
    
    
  }
  a2 <- coef(summary(a1))
  plotdata[i+2,1:3+6] <- exp(a2[2,1] + c(0,-1.96,1.96)*a2[2,2])
}


# ph1n1 , then H3N2
# cd4 cd8
# univariate
# sh1n1
# ph1n1
# sh3n2

plotdataB <- plotdata

# ph1n1 , then H3N2
# cd4 cd8
# univariate
# sh1n1
# ph1n1
# sh3n2

rowname <- c("All participants",
             " pH1N1 infections",
             "  sH1N1-specific CD4/CD8 (cross-strain)",
             "  pH1N1-specific CD4/CD8 (homosubtypic)",
             "  H3N2-specific CD4/CD8 (heterosubtypic)",
             " H3N2 infections",
             "  sH1N1-specific CD4/CD8 (heterosubtypic)",
             "  pH1N1-specific CD4/CD8 (heterosubtypic)",
             "  H3N2-specific CD4/CD8 (homosubtypic)",
             "Participants with HAI titer < 40",
             " pH1N1 infections",
             "  sH1N1-specific CD4/CD8 (cross-strain)",
             "  pH1N1-specific CD4/CD8 (homosubtypic)",
             "  H3N2-specific CD4/CD8 (heterosubtypic)",
             " H3N2 infections",
             "  sH1N1-specific CD4/CD8 (heterosubtypic)",
             "  pH1N1-specific CD4/CD8 (heterosubtypic)",
             "  H3N2-specific CD4/CD8 (homosubtypic)")

plotdata2A <- rbind(plotdataA[,1:3],plotdataA[,4:6],plotdataA[,7:9])
plotdata2B <- rbind(plotdataB[,1:3],plotdataB[,4:6],plotdataB[,7:9])

a2 <- rbind(cbind(plotdata2A[c(9,9,c(1,13,25),9,c(1,13,25)+1),],plotdata2A[c(9,9,c(1,13,25),9,c(1,13,25)+1)+2,]),
            cbind(plotdata2B[c(9,9,c(1,13,25),9,c(1,13,25)+1),],plotdata2B[c(9,9,c(1,13,25),9,c(1,13,25)+1)+2,]))
            
a2round <- format(round(a2,digits=2),nsmall=2)
a1 <- cbind(paste(a2round[,1]," (",a2round[,2],", ",a2round[,3],")",sep=""),paste(a2round[,4]," (",a2round[,5],", ",a2round[,6],")",sep=""))
a1[grepl("NA",a1)] <- NA

pdf("Shared kiddivax CMI/protection/figure2.pdf",width=12, height=10)

par(mar=c(0,0,0,0))


plot(0,0,xlab="",ylab="", main="", axes=F, xlim=c(-7.5,11.5), ylim=c(1.5,25.5),type="n")


cc <- 22 - (1:18-0.5)
left <- 13
panelb <- 7.5
textadj <- 1  

#for (i in 5:45){
#  if (i%%2==0){
#    adjfac <- 0.08
 #   polygon(c(-18,-18,6.5,6.5),c(i+1+ adjfac,i+ adjfac,i+ adjfac,i+1+ adjfac),col=rgb(0,0,0,0.05),border=F)  
#  }
#}

for (i in 1:18){
boldind <- c(1,2,6,10,11,15)
if (1*(i%in%boldind)){
  text(-21+left,cc[i]+0.08,rowname[i],las=1,pos=4,font=2) 
}
else{
  text(-21+left,cc[i]+0.08,rowname[i],las=1,pos=4,font=1) 
}
}


for (sym in 0:1){
  
for (u in 1:10){
  arrows(-0.3+sym*panelb,2.5,-3+sym*panelb,2.5,col="black",length=0.05)
  arrows(0.3+sym*panelb,2.5,3+sym*panelb,2.5,col="black",length=0.05)
}
  
  text(-1.8+sym*panelb,1.5,"Higher protection") 
  text(1.7+sym*panelb,1.5,"Lower protection")  


axis(1,at=-2:1+sym*panelb,labels=2^(-2:1),cex.axis=1,pos=4.3)
#axis(2,at=c(-1.5,cc,17),labels=NA,las=1, pos=-2.2)

plotvector <- log2(a2[,1+3*sym])
#plotvector[plotvector< -2] <- NA
points(plotvector+sym*panelb,cc+0.08,pch=16,col="black",cex=0.8)


for (i in 1:18){
  #+sym*panelb
  lines(c(max(-2,log2(a2[i,2+3*sym])+sym*panelb),min(1,log2(a2[i,3+3*sym]))+sym*panelb),rep(cc[i]+0.08,2),col="black")  
  testaddarrow <- max(-2+sym*panelb,log2(a2[i,2+3*sym]))==-2
  if (  testaddarrow & !is.na(  testaddarrow)){
    arrows(-1+sym*panelb,cc[i]+0.08,-2+sym*panelb,cc[i]+0.08,col="black",length=0.05)
  }
  
  testaddarrow <- min(1,log2(a2[i,3+3*sym])) == 1
  if (  testaddarrow & !is.na(  testaddarrow)){
    arrows(0+sym*panelb,cc[i]+0.08,1+sym*panelb,cc[i]+0.08,col="black",length=0.05)
  }
  text(4+sym*panelb-textadj,cc[i]+0.08,a1[i,1+sym],las=1)

}

lines(c(0,0)+sym*panelb,c(4.3,20.5),lty=2)


text(0+sym*panelb,24.5,"Adjusted Odds Ratio")
text(0+sym*panelb,23.5,"(95% CI)")
text(0,25.5,"CD4",cex=1.2)
text(0+panelb,25.5,"CD8",cex=1.2)
text(-21+left,22.5,"Every fold increase in CD4 or CD8 T cell response",las=1,pos=4,font=1) 

text(-2.5,25.5,"A",cex=1.5)
text(5,25.5,"B",cex=1.5)

}


dev.off()


#legend(13,2, cex=0.7, legend=c("A(H3N2)","B"),lty=1,pch=c(16,17),col=c("red","blue"),bty="n")

#title(main="A", adj=0)
#mtext("Step 1",side=3,cex=1,line=0, at=0)
#mtext("Step 2",side=3,cex=1,line=0, at=5)

#legend(5.5,14.5, cex=0.7, legend=c("A(H3N2)","B"),lty=1,pch=c(16,17),col=c("red","blue"),bty="n")

#mtext("A",side=3,cex=1,line=0,at=-2.5)
#mtext("B",side=3,cex=1,line=0,at=2.5)
