##########################################################################
##        Dicerandra christmanii - vital rate modeling                  ##
##        Date: 05/20/2020    - Version: 10 with plots                  ##
##        Pedro F. Quintana Ascencio & Federico Lopez Borghesi          ##
##########################################################################

rm(list=ls())
pkgs<-c("nlme","bbmle","lme4","lattice","MuMIn","MASS","dplyr",
        "optimx","rethinking","rstan")
### Check if package is installed, and if not, install them
for (i in 1: length(pkgs)){
  if(pkgs[i] %in% rownames(installed.packages()) == FALSE) {install.packages(pkgs[i])}
}

### load all packages
lapply(pkgs, library, character.only = TRUE)


##set working directory and data call
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

dat2 <- read.csv("Dicerandra working file.csv", header=T)
dat2 <- subset(dat2, coded_stage!=0 & coded_stage!=9 & coded_stage !=2)
dim(dat2)

##############################################
# Change in stems as natural logarithms
###############################################

dat2$log_ini_branches <- log(dat2$ini_branches)
dat2$log_fin_branches <- log(dat2$fin_branches)

##################################################################################
########  assigning mortality for plants sampled multiple times within a year

#dat2$ann_sur[ which(is.na(dat2$ann_sur) & dat2$ini_branches>0 & dat2$year != 2018)  ] <- 0

##################################################################################

#### MANAGEMENT TREATMENT DEFINITION  ############################################

dat2$treatment <- 1
dat2$treatment[dat2$sites >1 & dat2$sites <4] <- 2
dat2$treatment[dat2$sites == 6] <- 3
dat2$treatment[dat2$sites == 7] <- 4
dat2$ini_branches <- as.numeric(dat2$ini_branches)
dat2$fin_branches <- as.numeric(dat2$fin_branches)
dat2$treatment <- factor(dat2$treatment)

#### IDENTIFICTION OF STAGES  ######################################

dat2$stg <- dat2$coded_stage
dat2$stg[dat2$stg==1] <- 3 
dat2$stg <- factor(dat2$stg)

# 3 adults
# 5 yearlings

###### Re-classifying 1994 misclassified seedlings
dat2$stg[dat2$year==1994 & dat2$ini_branches >=4]  <-3

###### Correcting a plant with less than a branch
dat2$ini_branches[ which(dat2$ini_branches < 1) ] <- 1


#### STATUS FOR AUGMENTATION ##############################################
### For comparisson with preburn (below)
dat2$aug <- 0
dat2$aug[dat2$year>2009 & dat2$year <2012 & dat2$treatment==3] <- 1
dat2$aug1 <- 0
dat2$aug1[dat2$year==2010 & dat2$treatment==3] <- 1
dat2$aug2 <- 0
dat2$aug2[dat2$year==2011 & dat2$treatment==3] <- 1

## 0 no pulse augmentation
## 1 pulse augmentation

#### BURN STATUS FOR INTRODUCTION & AUGMENTATION ######################
#### THIS IS DIFFERENT FROM PREBURN AND BURN ##########################

dat2$burn_status[dat2$burn_status==2] <- 1

## 0 no burn
## 1 burn

####### PRETREATMENT DEFINITION FOR AUGMENTATION #####################

dat2$preburn <-0
dat2$preburn[dat2$trt_type==4] <- 1

## 0 no preburn
## 1 preburn


#### PULSE EFFECT FOR INTRODUCTION ####################################

dat2$intr <- 0
dat2$intr[dat2$year>2011 & dat2$year <2014 & dat2$treatment==4] <- 1
dat2$intr1 <- 0
dat2$intr1[dat2$year==2012 & dat2$treatment==4] <- 1
dat2$intr2 <- 0
dat2$intr2[dat2$year==2013 & dat2$treatment==4] <- 1

dat2$intryear <- 0
dat2$intryear [dat2$year==2013 & dat2$treatment==4 & dat2$stg==5] <- 1

## 0 no pulse of introduction
## 1 pulse of introduction 

### HABITAT CONDITION DEFINITION FOR INTRODUCTION ######################


dat2$habitat <-0
dat2$habitat[dat2$trt_type==6] <-1

#table(datoriginal$trt_type[ datoriginal$survival !=0 & datoriginal$survival !=9 & datoriginal$survival !=2])
## interior (0)
## edge (1)


### BURN CONDITION DEFINITION FOR INTRODUCTION ##########################

dat2$burn <- 0
dat2$burn[dat2$treatment==4 & dat2$burn_status==1] <- 1

###################################################################################
####   Dummy variables to test hypothesis
###################################################################################

yr<-c(1994:2018)
yrlab<-as.character(yr)

dat2$T2 <- 0
dat2$T2[dat2$treatment==2] <- 1
dat2$T3 <- 0
dat2$T3[dat2$treatment==3] <- 1
dat2$T4 <- 0
dat2$T4[dat2$treatment==4] <- 1
dat2$yn <- dat2$year-1994

dat2$py <- 0
k <-0
for (i in 1:length(yr)){
  for(j in 1:7){
    k <- k+1
    #print(c(i,j))
    dat2$py[dat2$sites==j & dat2$yn==i]  <- k
  }
}

# write.csv(dat2,"Dicerandra.work.file.csv") #### IF NECESSARY

###########################################################################
###########################################################################

yr<-c(1994:2018)
yrlab<-as.character(yr)

############################################################################

load("growth_adultsi.RData", .GlobalEnv)
load("growth_seedlingsni.RData", .GlobalEnv)
load("Preprodi.RData", .GlobalEnv)
load("Fecundity_prop.RData", .GlobalEnv)
load("Survival.RData", .GlobalEnv)
load("fruits.RData", .GlobalEnv)
load("Recruitment_sbk.RData", .GlobalEnv)
load("Seeds.RData", .GlobalEnv)


######################### END Fecundity ################################### 

###########################################################################
###########################################################################

##### SURVIVAL (bothe coded_stages in the same model)

###########################################################################
###########################################################################
#tiff(filename = "Vital rates all together.tiff", width = 1000, height = 1000,
#     pointsize =1, res=300) 
par(mfcol=c(2,3), mar=c(4,4,2,2))

x <- seq(log(1),log(800),length.out = 100)



trt <- array(0,c(5,4))
# plot conditions
trt[1,] <- c(0,0,0,0)
trt[2,] <- c(0,1,0,0)
trt[3,] <- c(0,0,1,0)
trt[4,] <- c(0,0,0,1)
trt[5,] <- c(0,0,0,1)
trtlab<-c("Wild","Road","Augmented","Introduced","Introduced")

clrs<-list(rgb(0,104,139,255,maxColorValue = 255),
           rgb(252,111,70,255,maxColorValue = 255),
           rgb(113,90,160,255,maxColorValue = 255),
           rgb(63,144,94,255,maxColorValue = 255),
           rgb(63,144,94,255,maxColorValue = 255))

clrs2<-list(rgb(0,104,139,50,maxColorValue = 255),
            rgb(252,111,70,50,maxColorValue = 255),
            rgb(113,90,160,50,maxColorValue = 255),
            rgb(63,144,94,50,maxColorValue = 255),
            rgb(63,144,94,50,maxColorValue = 255))

dat5 <-dat2 
dat5$survivor<-1
dat5$survivor[is.na(dat5$fin_branches)]<-0


dat5$T2 <- 0
dat5$T2[dat5$treatment==2] <- 1
dat5$T3 <- 0
dat5$T3[dat5$treatment==3] <- 1
dat5$T4 <- 0
dat5$T4[dat5$treatment==4] <- 1
dat5$yn <- dat5$year-1994
dat5 <- dat5[,c(-3,-7,-16,-18)]
dat5 <- dat5[-which(is.na(dat5$tag)),]
row.has.na <- apply(dat5, 1, function(x){any(is.na(x))})
sum(row.has.na)
dat5$py <- 0
k <-0
for (i in 1:109){
  for(j in 1:7){
    k <- k+1
    #print(c(i,j))
    dat5$py[dat5$sites==j & dat5$yn==i]  <- k
  }
}

dat5$stg_C[dat5$stg==3] <- 0
dat5$stg_C[dat5$stg==5] <- 1
dat5 <- subset(dat5,!is.na(yn) & !is.na(py))
row.has.na <- apply(dat5, 1, function(x){any(is.na(x))})
#sum(row.has.na)

dat5$treatment <- as.factor(dat5$treatment)
dat5$stg <- as.factor(dat5$stg)
dat5$log_ini_branches <- log(dat5$ini_branches)
dat5$logsquare <- dat5$log_ini_branches^2


dat5$burn[dat5$year!=2014 & dat5$treatment==4] <- 0
dat5$postburn <- 0
dat5$postburn[dat5$year> 2014 & dat5$treatment==4] <- 1
table(dat5$burn,dat5$year,dat5$treatment)
table(dat5$postburn,dat5$year,dat5$treatment)

######################   plot survival GLM


x <- seq(log(1),log(1000),0.1)



# plot conditions
trt2 <- c(0,1,0,0,0) 
trt3 <- c(0,0,1,0,0)
trt4 <- c(0,0,0,1,1)
wtrtlab<-c("Wild","Road","Augmented","Introduced","Introduced")

stg <- c(0,1) # 0-adult; 1:seedling
stgname<-c("3","5")
stsur<-c(3,5)
labstg <- c("Adults","Yearlings")
hab <- c(0,0,0,0,1)
lintyp<-c(1,1,1,2,1)


sum_model5 <- precis(model5,digits=3,depth=2)
sum_model5 <- precis(model5,digits=3,depth=1)


##################################### Plot Bayesian P(Survival) together ###
seq <- c(2,1)

pj.zeros <- matrix(0,1000,109)
yj.zeros <- matrix(0,1000,25)

for (j in seq){
  labelmain <- paste("Survival - ",labstg[j])
  if (j==1){x  <-seq(log(1),log(800),0.1)}
  if (j==2){x  <-seq(log(1),log(40),0.1)}
  plot(x,seq(0,1,length.out=length(x)), xlab="ln(Number of Branches)",
       ylim=c(0,1),main=labelmain,xlim=c(0, 8),type="n",ylab="Probability",cex.lab=1.2)
  #if (j==1){text(0.12,0.97,"B",cex=2)}
  #if (j==2){text(0.12,0.97,"A",cex=2)}
  #points(dat5$log_ini_branches[dat5$stg == stsur[j]],dat5$survivor[dat5$stg == stsur[j]],col=COL,pch=16,cex=0.6)
  #abline(v=mean(dat5$log_ini_branches[dat5$stg==stgname[j]]),lty=3)
  #abline(v=max(dat5$log_ini_branches[dat5$stg==stgname[j]]),lty=3,col="red") 
  
  for (i in 1:length(trt2)){
    v.trt2 <- trt2[i]
    v.trt3 <- trt3[i]
    v.trt4 <- trt4[i]
    new.data <-  data.frame(cbind(rep(1,length(x)),
                                  log_ini_branches = x,
                                  logsquare = x^2,
                                  T2 = rep(v.trt2,length(x)),
                                  T3 =rep(v.trt3,length(x)),
                                  T4 = rep(v.trt4,length(x)),
                                  stg_C= stg[j],
                                  habitat = rep(hab[i],length(x)),
                                  aug = rep(0,length(x)), 
                                  burn = rep(0,length(x)),
                                  preburn = rep(0,length(x)),
                                  postburn = rep(0,length(x)),
                                  py=1,
                                  yn=1))
    #change aug
    mu <- NULL
    mu <- link(model5,n=1000,data=new.data,
               replace=list(pp=pj.zeros,y=yj.zeros))
    mu.PI <- apply(mu,2,PI,prob=0.75)
    shade(mu.PI,x,col=clrs2[[i]])
    #lines(x,mu.PI[1,],col=clrs[[i]],lwd=2,lty=2)
    #lines(x,mu.PI[2,],col=clrs[[i]],lwd=2,lty=2)
  }
  for (i in 1:length(trt2)){
    v.trt2 <- trt2[i]
    v.trt3 <- trt3[i]
    v.trt4 <- trt4[i]
    new.data <-  data.frame(cbind(rep(1,length(x)),
                                  log_ini_branches = x,
                                  logsquare = x^2,
                                  T2 = rep(v.trt2,length(x)),
                                  T3 =rep(v.trt3,length(x)),
                                  T4 = rep(v.trt4,length(x)),
                                  stg_C= stg[j],
                                  habitat = rep(hab[i],length(x)),
                                  aug = rep(0,length(x)), 
                                  burn = rep(0,length(x)),
                                  preburn = rep(0,length(x)),
                                  postburn = rep(0,length(x)),
                                  py=1,
                                  yn=1))
    #change aug
    mu <- NULL
    mu <- link(model5,n=1000,data=new.data,
               replace=list(pp=pj.zeros,y=yj.zeros))
    mu.mean <- apply(mu,2,mean)
    lines(x,mu.mean,col=clrs[[i]],lwd=3,lty=lintyp[i])
  }
  
    }
   



###########################################################################
###########################################################################

# GROWTH logarithmic change in number of branches for YEARLINGS

###########################################################################
###########################################################################


dat3y <- subset(dat2, stg ==5 & !is.na(ini_branches)& !is.na(fin_branches))
dat3y$treatment <- as.factor(dat3y$treatment)
dat3y$log_ini_branches <- log(dat3y$ini_branches)
dat3y$log_fin_branches <- log(dat3y$fin_branches)


dat3y$burn_status[is.na(dat3y$burn_status)] <- 0
dat3y$gap[is.na(dat3y$gap)] <- 0
dat3y <- dat3y[,c(-4)]
dat3y <- na.omit(dat3y)
row.has.na <- apply(dat3y, 1, function(x){any(is.na(x))})
sum(row.has.na)


####################################################################### 

x <- seq(log(1),log(40),0.1)

post <- extract.samples(model2)
total_a <- sapply(1:40,function(beta1)post$a+post$p[,beta1])

#################################### Plot Bayesian

trt <- array(0,c(5,4))
# plot conditions
trt[1,] <- c(0,0,0,0)
trt[2,] <- c(0,1,0,0)
trt[3,] <- c(0,0,1,0)
trt[4,] <- c(0,0,0,1)
trt[5,] <- c(0,0,0,1)
trtlab<-c("Wild","Road","Augmented","Introduced")

COL <- adjustcolor(c("black", "red", "green", "blue")[dat4$treatment], alpha.f = 0.5)
colt<- c("black", "red", "green", "blue")

pj.zeros <- matrix(0,1000,97)
yj.zeros <- matrix(0,1000,24)


  plot(log(x),log(seq(0,1,length.out=length(x))), xlab="ln ( Number of Branches (t))",
       ylim=c(0,8),main="Growth Yearlings",xlim=c(0,8),type="n",ylab="ln( Number of Branches (t+1)",cex.lab=1.2)
  #text(0.12,7.9,"C",cex=2)
  #points(log(dat3y$ini_branches[dat3y$treatment==i]),log(dat3y$fin_branches[dat3y$treatment==i]),
  #       col=colt[i],pch=16,cex=0.4) 
  for (i in 1:dim(trt)[1]){
  v.trt2 <- trt[i,2]
  v.trt3 <- trt[i,3]
  v.trt4 <- trt[i,4]
  new.data <-  data.frame(cbind(rep(1,length(x)),
                                log_ini_branches = x,
                                T2 = rep(v.trt2,length(x)),
                                T3 =rep(v.trt3,length(x)),
                                T4 = rep(v.trt4,length(x)),
                                habitat = rep(hab[i],length(x)),
                                aug1 = rep(0,length(x)),
                                aug2 = rep(0,length(x)),
                                preburn = rep(0,length(x)),
                                py=1,
                                yn=1))
  
  mu <- NULL
  mu <- link(model2,n=1000,data=new.data,
             replace=list(pp=pj.zeros,y=yj.zeros))
  mu.PI <- apply(mu,2,PI,prob=0.75)
  shade(mu.PI,x,col=clrs2[[i]])
  
  }
  
  for (i in 1:dim(trt)[1]){
    v.trt2 <- trt[i,2]
    v.trt3 <- trt[i,3]
    v.trt4 <- trt[i,4]
    new.data <-  data.frame(cbind(rep(1,length(x)),
                                  log_ini_branches = x,
                                  T2 = rep(v.trt2,length(x)),
                                  T3 =rep(v.trt3,length(x)),
                                  T4 = rep(v.trt4,length(x)),
                                  habitat = rep(hab[i],length(x)),
                                  aug1 = rep(0,length(x)),
                                  aug2 = rep(0,length(x)),
                                  preburn = rep(0,length(x)),
                                  py=1,
                                  yn=1))
    
    mu <- NULL
    mu <- link(model2,n=1000,data=new.data,
               replace=list(pp=pj.zeros,y=yj.zeros))
    mu.mean <- apply(mu,2,mean)
    lines(x,mu.mean,col=clrs[[i]],lwd=3,lty=lintyp[i])
  }    
    lines(0:8,0:8, lty=2, col="darkgrey")
    
#legend(0.5,8,trtlab,col=c(1:i),lty = c(1, 1, 1, 1),cex=0.7,ncol=2,bty = "n")

########################## END YEARLING GROWTH ###########################
###########################################################################
###########################################################################
  
  # GROWTH logarithmic change in number of branches for ADULT PLANTS
  
###########################################################################
###########################################################################
  
  dat4 <- subset(dat2, stg ==3 & !is.na(ini_branches) & !is.na(fin_branches))
  table(dat4$treatment,dat4$year)
  table(dat4$treatment,dat4$sites)
  length( which(is.na(dat4$ann_sur))  )
  length( which(is.na(dat4$ann_sur) & dat4$ini_branches>0 & dat4$year != 2018)  )
  length( which(is.na(dat4$ann_sur) & dat4$ini_branches>0 & dat4$year == 2018)  )
  
  ###########################################################################
  #plot Growth
  x <- seq(log(1),log(800),length.out = 100)
  
  labelmain <- paste("Growth Adults")
  
  
##################################### Plot Bayesian
  
  #  load("growth_adultsi.RData", .GlobalEnv)
  
  
  pj.zeros <- matrix(0,1000,102)
  yj.zeros <- matrix(0,1000,24)
  
  
  plot(log(x),log(seq(0,1,length.out=length(x))), xlab="ln (Number of Branches (t))",
       ylim=c(0,8),main="Growth Adults",xlim=c(0,8),type="n",ylab="ln (Number of Branches (t+1)",
       cex.lab=1.2)
  #text(0.12,7.9,"D",cex=2)
  for (i in 1:dim(trt)[1]){
    #points(log(dat4$ini_branches[dat4$treatment==i]),log(dat4$fin_branches[dat4$treatment==i]),
    #       col=COL[i],pch=16,cex=0.3)   
    v.trt2 <- trt[i,2]
    v.trt3 <- trt[i,3]
    v.trt4 <- trt[i,4]
    new.data <-  data.frame(cbind(rep(1,length(x)),
                                  log_ini_branch = x,
                                  T2 = rep(v.trt2,length(x)),
                                  T3 =rep(v.trt3,length(x)),
                                  T4 = rep(v.trt4,length(x)),
                                  intr = rep(0,length(x)),
                                  habitat = rep(hab[i],length(x)),
                                  py=rep(1,length(x)),
                                  yn=rep(1,length(x))))
    mu <- NULL
    mu <- link(model1,n=1000,data=new.data,
               replace=list(p=pj.zeros,y=yj.zeros))
    mu.PI <- apply(mu,2,PI,prob=0.95)
    shade(mu.PI,x,col=clrs2[[i]])

  }
  
  for (i in 1:dim(trt)[1]){
    #points(log(dat4$ini_branches[dat4$treatment==i]),log(dat4$fin_branches[dat4$treatment==i]),
    #       col=COL[i],pch=16,cex=0.3)   
    v.trt2 <- trt[i,2]
    v.trt3 <- trt[i,3]
    v.trt4 <- trt[i,4]
    new.data <-  data.frame(cbind(rep(1,length(x)),
                                  log_ini_branch = x,
                                  T2 = rep(v.trt2,length(x)),
                                  T3 =rep(v.trt3,length(x)),
                                  T4 = rep(v.trt4,length(x)),
                                  intr = rep(0,length(x)),
                                  habitat = rep(hab[i],length(x)),
                                  py=rep(1,length(x)),
                                  yn=rep(1,length(x))))
    mu <- NULL
    mu <- link(model1,n=1000,data=new.data,
               replace=list(p=pj.zeros,y=yj.zeros))
    mu.mean <- apply(mu,2,mean)
    lines(x,mu.mean,col=clrs[[i]],lwd=3,lty=lintyp[i])
    
  }
  
  lines(0:8,0:8,lty=2,col="darkgrey")
  
########################### END ADULT GROWTH ##############################
  
###########################################################################
###########################################################################

##### PROBABILITY OF REPRODUCTION only adults

###########################################################################
###########################################################################

dat2p <- subset(dat2,stg==3)
dat2p$prep <- 0
dat2p$prep[dat2p$flw_branches>0] <- 1
dat2p$log_ini_branches <- log(dat2p$ini_branches)


######################################################### Bayesian P(rep)
dat2$stg_C[dat2$stg==3] <- 0
dat2$stg_C[dat2$stg==5] <- 1

dat2p <- subset(dat2p, !is.na(log_ini_branches))

stg <- c(0,1) # 0-adult; 1:yearling
stgname<-c("3","5")
stsur<-c(3,5)
labstg <- c("Adult","Yearling")
o <- order(x)

##################################### Bayesian P(Rep), average ###
# load("Preprodi.RData", .GlobalEnv)

pj.zeros <- matrix(0,1000,109)
yj.zeros <- matrix(0,1000,25)
x1 <- seq(log(1),log(800),length.out = 100)

plot(x[o],x[o], xlab="ln (Number of Branches)",
     ylim=c(0,1),main="Adult P( Reproduction )",xlim=c(0, 8),type="n",ylab="Probability",cex.lab=1.2)
#text(0.12,0.97,"E",cex=2)
for (j in 1:1){ 
  if (j==1){x1 <- seq(log(1),log(800),0.1)}
  if (j==2){x1 <- seq(log(1),log(40),0.1)}
  #par(mfrow=c(2,2))
  for (i in 1:dim(trt)[1]){
  #c.prep <- seq(min(dat2p$log_ini_branches[dat2p$stg == stsur[j] & dat2p$treatment==i]-0.001), 
  #              max(dat2p$log_ini_branches[dat2p$stg == stsur[j] & dat2p$treatment==i]),0.3)
  #v.prep <- cut(dat2p$log_ini_branches[dat2p$stg == stsur[j] & dat2p$treatment==i],c.prep)
  #t.prep <- table( dat2p$prep[dat2p$stg == stsur[j] & dat2p$treatment==i],v.prep)
  #if(sum(t.prep[1,]) != sum(colSums(t.prep))){
  #p.prep <- t.prep[2,]/colSums(t.prep)
  #points(dat2p$log_ini_branches[dat2p$stg == stsur[j] & dat2p$treatment==i],
  #       dat2p$prep[dat2p$stg == stsur[j] & dat2p$treatment==i],col=colt[i],pch=16,cex=0.6)
  #if (j==1){
  #  points(c.prep[-1],p.prep,col=colt[i],pch=1,cex=log(colSums(t.prep)))}
  #if (j==2){
  #  points(c.prep[-1],p.prep,col=colt[i],pch=1,cex=log(colSums(t.prep)))}
  #abline(v=mean(dat2p$log_ini_branches[dat2p$stg==stgname[j]]),lty=3)
  #abline(v=max(dat2p$log_ini_branches[dat2p$stg==stgname[j]]),lty=3,col="red") }
  #if(sum(t.prep[1,]) == sum(colSums(t.prep))){
  #  plot(x[o],x[o], xlab="ln (Number of Branches)",
  #       ylim=c(0,1),main=labelmain,xlim=c(0, 8),type="n",ylab="Probability")
  #  points(dat2$log_ini_branches[dat2p$stg == stsur[j] & dat2p$treatment==i],
  #         dat2$prep[dat2p$stg == stsur[j] & dat2p$treatment==i],col=colt[i],pch=16,cex=log(colSums(t.prep)))
  #  abline(v=mean(dat2p$log_ini_branches[dat2p$stg==stgname[j]]),lty=3)
  #  abline(v=max(dat2p$log_ini_branches[dat2p$stg==stgname[j]]),lty=3,col="red") }
  
  v.trt2 <- trt[i,2]
  v.trt3 <- trt[i,3]
  v.trt4 <- trt[i,4]
  new.data <-  data.frame(cbind(rep(1,length(x1)),
                                log_ini_branch = x1,
                                T2 = rep(v.trt2,length(x1)),
                                T3 =rep(v.trt3,length(x1)),
                                T4 = rep(v.trt4,length(x1)),
                                stg_C = rep(stg[j],length(x1)),
                                aug1 = rep(0,length(x1)),
                                aug2 = rep(0,length(x1)),
                                habitat = rep(hab[i],length(x1)),
                                preburn = rep(0,length(x1)),
                                burn = rep(0,length(x1)),
                                py=1,
                                yn=1))
  mu <- NULL
  mu <- link(model3,n=1000,data=new.data,
             replace=list(pp=pj.zeros,y=yj.zeros))
  mu.PI <- apply(mu,2,PI,prob=0.95)
  shade(mu.PI,x1,col=clrs2[[i]])
  }

  for (i in 1:dim(trt)[1]){
#    c.prep <- seq(min(dat2p$log_ini_branches[dat2p$stg == stsur[j] & dat2p$treatment==i]-0.001), 
#                  max(dat2p$log_ini_branches[dat2p$stg == stsur[j] & dat2p$treatment==i]),0.3)
#    v.prep <- cut(dat2p$log_ini_branches[dat2p$stg == stsur[j] & dat2p$treatment==i],c.prep)
#    t.prep <- table( dat2p$prep[dat2p$stg == stsur[j] & dat2p$treatment==i],v.prep)

    v.trt2 <- trt[i,2]
    v.trt3 <- trt[i,3]
    v.trt4 <- trt[i,4]
    new.data <-  data.frame(cbind(rep(1,length(x1)),
                                  log_ini_branch = x1,
                                  T2 = rep(v.trt2,length(x1)),
                                  T3 =rep(v.trt3,length(x1)),
                                  T4 = rep(v.trt4,length(x1)),
                                  stg_C = rep(stg[j],length(x1)),
                                  aug1 = rep(0,length(x1)),
                                  aug2 = rep(0,length(x1)),
                                  habitat = rep(hab[i],length(x1)),
                                  preburn = rep(0,length(x1)),
                                  burn = rep(0,length(x1)),
                                  py=1,
                                  yn=1))
    mu <- NULL
    mu <- link(model3,n=1000,data=new.data,
               replace=list(pp=pj.zeros,y=yj.zeros))
    mu.mean <- apply(mu,2,mean)
    lines(x1,mu.mean,col=clrs[[i]],lwd=3,lty=lintyp[i])
  }
}

#legend(2.2,0.2,trtlab,col=c(1:i),lty = c(1, 1, 1, 1),cex=0.7,ncol=2,bty = "n")


###########################################################################
###########################################################################

##### FECUNDITY ONLY FOR ADULTS (because small sample size)

###########################################################################
###########################################################################

dat3a <- subset(dat2,flw_branches>0 & stg == "3" & !is.na(ini_branches) & !is.na(flw_branches ))
dat3a$log_flw_branches <- log(dat3a$flw_branches)
dat3a$log_ini_branches <- log(dat3a$ini_branches)


#################### Bayesian Fecundity ########################## 
dat3a $T2 <- 0
dat3a $T2[dat3a $treatment==2] <- 1
dat3a $T3 <- 0
dat3a $T3[dat3a $treatment==3] <- 1
dat3a $T4 <- 0
dat3a $T4[dat3a $treatment==4] <- 1
dat3a $yn <- dat3a $year-2011
dat3a  <- dat3a [,c(-4,-5)]
#dat3a  <- na.omit(dat3a )
row.has.na <- apply(dat3a , 1, function(x){any(is.na(x))})
sum(row.has.na)
#dim(dat3a)

dat3a $py <- 0
dat3a$yn <- dat3a$year-1994

sites <- table(dat3a$sites,dat3a$year)
dsites <- dim(sites)  
k <-0
for (i in 1:dsites[2]){
  for(j in 1:dsites[1]){
    if (sites[j,i]>0){
    k <- k+1
    #print(c(i,j))
    dat3a $py[dat3a $sites==j & dat3a $yn==i]  <- k}
  }
}
#length(unique(dat3a$py))
#unique(dat3a$py[order(dat3a$py)])
#dim(dat3a)

dat3a  <- subset(dat3a ,!is.na(yn) & !is.na(py) & py>0)
row.has.na <- apply(dat3a , 1, function(x){any(is.na(x))})
#sum(row.has.na)
#dim(dat3a)

#unique(dat3a$py)[order(unique(dat3a$py))]

######################################################### Bayesian P(fec proportions)
##### USES STAN CODE IN FecundityR.R

#Outname <- paste("Fecundity_prop.RData",sep="")
#save(model4, file = Outname , ascii=T)


#  load("Fecundity_prop.RData", .GlobalEnv)
sum_model4 <- precis(model4,digits=3,depth=1)
#pairs(model4,pars=c("a","b1","b22","b23","b24","b32","b33",
#                    "b34","b4","sigmay","sigmap"))

postfec<- extract(model4)
#print(names(postfec))
#head(postfec$a)
#head(postfec$s)
#dim(postfec$y)
matrix_of_post <- as.matrix(model4)
#print(colnames(matrix_of_post))
df_of_draws <- as.data.frame(model4)
#print(colnames(df_of_draws))

# sum_model4 <- precis(model4,digits=3,depth=2)
sum_model4 <- precis(model4,digits=3,depth=1)

##################################### Bayesian P(fecundity, proportions) ###

calcpred <-function(matrixs,n,data) {
  preds<- array(0,c(n,nrow(data)))
  d <- dim(matrixs)
  chos <- sample(d[1],n)
  for(i in 1:1000){
    mui <- as.matrix(data)%*%matrixs[chos[i],1:13]
    preds[i,] <- mui
  }
  return(preds)
} 

##################################### Average both parts together####################
ag <- c(0,1)
pb <- c(0,1)
hab <- c(0,0,0,0,1)
br <- c(0,1)

x <- seq(log(1),log(800),length.out = 100)

  labelmain <- c("Adult Fecundity") 
plot(x,x, xlab="ln (Number of Branches)",
     ylim=c(0,8),main=labelmain,xlim=c(0, 8),type="n",ylab="ln (Flowering Branches)",cex.lab=1.2)
#points(dat3a$log_ini_branches[dat3a$treatment==i] ,dat3a$log_flw_branch[dat3a$treatment==i],
#       col=colt[i],pch=16,cex=0.3) 
#text(0.12,7.9,"F",cex=2)
for (i in 1:dim(trt)[1]){
  
  v.trt2 <- trt[i,2]
  v.trt3 <- trt[i,3]
  v.trt4 <- trt[i,4]
  new.data <-  data.frame(cbind(rep(1,length(x)),
                                x,
                                rep(v.trt2,length(x)),
                                rep(v.trt3,length(x)),
                                rep(v.trt4,length(x)),
                                x*v.trt2,
                                x*v.trt3,
                                x*v.trt4,
                                rep(0,length(x)),
                                rep(0,length(x)),
                                rep(hab[i],length(x)),
                                rep(0,length(x)),
                                rep(0,length(x)))
  )
  
  mu <- 1/(1 + (1/exp(calcpred(matrix_of_post,n=1000,new.data))))
  mu.PI <- apply(mu,2,PI,prob=0.95)
  shade(log(sweep(mu.PI, 2, exp(x), FUN = "*")),x,col=clrs2[[i]])

}

for (i in 1:dim(trt)[1]){

  v.trt2 <- trt[i,2]
  v.trt3 <- trt[i,3]
  v.trt4 <- trt[i,4]
  new.data <-  data.frame(cbind(rep(1,length(x)),
                                x,
                                rep(v.trt2,length(x)),
                                rep(v.trt3,length(x)),
                                rep(v.trt4,length(x)),
                                x*v.trt2,
                                x*v.trt3,
                                x*v.trt4,
                                rep(0,length(x)),
                                rep(0,length(x)),
                                rep(hab[i],length(x)),
                                rep(0,length(x)),
                                rep(0,length(x)))
                                )
  
  mu <- 1/(1 + (1/exp(calcpred(matrix_of_post,n=1000,new.data))))
  mu.mean <- apply(mu,2,mean)
  #mu.PI <- apply(mu,2,PI,prob=0.95)
  lines(x,log(mu.mean*exp(x)),col=clrs[[i]],lwd=2.5,lty=lintyp[i])
 
  }
  
lines(0:9,0:9,lty=2,col="darkgrey")


######################### END Fecundity ################################### 


# dev.off()
#legend(2,0.4,trtlab,col=c(1:i),lty = c(1, 1, 1, 1),cex=0.7,ncol=2,bty = "n")




###################


### APPENDIX GRAPHS

## SURVIVAL YEARLINGS

par(mfrow=c(2,3))

trt <- array(0,c(6,4))
# plot conditions
trt[1,] <- c(0,0,0,0)
trt[2,] <- c(0,1,0,0)
trt[3,] <- c(0,0,1,0)
trt[4,] <- c(0,0,1,0)
trt[5,] <- c(0,0,0,1)
trt[6,] <- c(0,0,0,1)

clrs<-list(rgb(0,104,139,255,maxColorValue = 255),
           rgb(252,111,70,255,maxColorValue = 255),
           rgb(113,90,160,255,maxColorValue = 255),
           rgb(113,90,160,255,maxColorValue = 255),
           rgb(63,144,94,255,maxColorValue = 255),
           rgb(63,144,94,255,maxColorValue = 255))

clrs2<-list(rgb(0,104,139,50,maxColorValue = 255),
            rgb(252,111,70,50,maxColorValue = 255),
            rgb(113,90,160,50,maxColorValue = 255),
            rgb(113,90,160,50,maxColorValue = 255),
            rgb(63,144,94,50,maxColorValue = 255),
            rgb(63,144,94,50,maxColorValue = 255))


pj.zeros <- matrix(0,1000,109)
yj.zeros <- matrix(0,1000,25)

labPopTyp<-c("Natural gap pop.","Natural roadside pop.","Augmentation unburned",
             "Augmentation pre-burned","Introduction in gap","Introduction in roadside")

trt2 <- c(0,1,0,0,0,0) 
trt3 <- c(0,0,1,1,0,0)
trt4 <- c(0,0,0,0,1,1)
brn<- c(0,0,0,1,0,0)
hab<-c(0,0,0,0,0,1)

x<-seq(log(1),log(40),0.1)
for (i in 1:6){
  plot(x,seq(0,1,length.out=length(x)), xlab="ln(Number of Branches)",
       ylim=c(0,1),main=labPopTyp[i],xlim=c(0, 4),type="n",ylab="Probability",cex.lab=1.2)

    siz<-dat5$log_ini_branches[dat5$stg==5 & dat5$trt_type==i & dat5$burn==0]
    sur<-dat5$survivor[dat5$stg==5 & dat5$trt_type==i & dat5$burn==0]
    
    b <- c(-0.1,0.2,0.8,1.6,4)
    pos<-c(0,0.5,1.2,3.1)
    z <- cut(siz, b)
    prebyden <- tapply(sur,z,sum)
    tab <- table(z)
    probs <-prebyden/tab
    probs <- as.vector(probs)
    points(pos,probs,pch=16,cex=(tab)^(2/7),col=rgb(0.3,0.3,0.3,0.5))

    
    
    v.trt2 <- trt2[i]
    v.trt3 <- trt3[i]
    v.trt4 <- trt4[i]
    new.data <-  data.frame(cbind(rep(1,length(x)),
                                  log_ini_branches = x,
                                  logsquare = x^2,
                                  T2 = rep(v.trt2,length(x)),
                                  T3 =rep(v.trt3,length(x)),
                                  T4 = rep(v.trt4,length(x)),
                                  stg_C= rep(1,length(x)),
                                  habitat = rep(hab[i],length(x)),
                                  aug = rep(0,length(x)), 
                                  burn = rep(0,length(x)),
                                  preburn = rep(brn[i],length(x)),
                                  postburn = rep(0,length(x)),
                                  py=1,
                                  yn=1))
    #change aug
    mu <- NULL
    mu <- link(model5,n=1000,data=new.data,
               replace=list(pp=pj.zeros,y=yj.zeros))
    mu.PI <- apply(mu,2,PI,prob=0.75)
    shade(mu.PI,x,col=clrs2[[i]])
    mu.mean <- apply(mu,2,mean)
    lines(x,mu.mean,col=clrs[[i]],lwd=3)
  }
  

## SURVIVAL ADULTS
x<-seq(log(1),log(800),0.1)
for (i in 1:6){
  plot(x,seq(0,1,length.out=length(x)), xlab="ln(Number of Branches)",
       ylim=c(0,1),main=labPopTyp[i],xlim=c(0,8),type="n",ylab="Probability",cex.lab=1.2)
  
  siz<-dat5$log_ini_branches[dat5$stg==3 & dat5$trt_type==i & dat5$burn==0]
  sur<-dat5$survivor[dat5$stg==3 & dat5$trt_type==i & dat5$burn==0]
  
  b <- c(-0.1,0.2,0.8,1.2,1.6,1.8,2,2.2,2.7,3,3.5,4,8)
  pos<-c(0,0.5,1,1.4,1.7,1.9,2.1,2.5,2.9,3.3,3.8,6)
  z <- cut(siz, b)
  prebyden <- tapply(sur,z,sum)
  tab <- table(z)
  probs <-prebyden/tab
  probs <- as.vector(probs)
  points(pos,probs,pch=16,cex=(tab)^(2/7),col=rgb(0.3,0.3,0.3,0.5))
  
  v.trt2 <- trt2[i]
  v.trt3 <- trt3[i]
  v.trt4 <- trt4[i]
  new.data <-  data.frame(cbind(rep(1,length(x)),
                                log_ini_branches = x,
                                logsquare = x^2,
                                T2 = rep(v.trt2,length(x)),
                                T3 =rep(v.trt3,length(x)),
                                T4 = rep(v.trt4,length(x)),
                                stg_C= rep(0,length(x)),
                                habitat = rep(hab[i],length(x)),
                                aug = rep(0,length(x)), 
                                burn = rep(0,length(x)),
                                preburn = rep(brn[i],length(x)),
                                postburn = rep(0,length(x)),
                                py=1,
                                yn=1))
  #change aug
  mu <- NULL
  mu <- link(model5,n=1000,data=new.data,
             replace=list(pp=pj.zeros,y=yj.zeros))
  mu.PI <- apply(mu,2,PI,prob=0.75)
  shade(mu.PI,x,col=clrs2[[i]])
  mu.mean <- apply(mu,2,mean)
  lines(x,mu.mean,col=clrs[[i]],lwd=3)
}


### GROWTH YEARLINGS
par(mfrow=c(2,3))
x <- seq(log(1),log(40),0.1)

post <- extract.samples(model2)
total_a <- sapply(1:40,function(beta1)post$a+post$p[,beta1])

for (i in 1:dim(trt)[1]){
  
  plot(log(x),log(seq(0,1,length.out=length(x))), xlab="ln ( Number of Branches (t))",
  ylim=c(0,8),main=labPopTyp[i],xlim=c(0,4),type="n",ylab="ln( Number of Branches (t+1)",cex.lab=1.2)
  
  iniS<-dat2$log_ini_branches[dat2$stg==5 & dat2$trt_type==i]
  finS<-dat2$log_fin_branches[dat2$stg==5 & dat2$trt_type==i]
  
  points(iniS,finS,col=rgb(0.3,0.3,0.3,0.2),pch=16)
  
  v.trt2 <- trt[i,2]
  v.trt3 <- trt[i,3]
  v.trt4 <- trt[i,4]
  new.data <-  data.frame(cbind(rep(1,length(x)),
                                log_ini_branches = x,
                                T2 = rep(v.trt2,length(x)),
                                T3 =rep(v.trt3,length(x)),
                                T4 = rep(v.trt4,length(x)),
                                habitat = rep(hab[i],length(x)),
                                aug1 = rep(0,length(x)),
                                aug2 = rep(0,length(x)),
                                preburn = rep(brn[i],length(x)),
                                py=1,
                                yn=1))
  
  mu <- NULL
  mu <- link(model2,n=1000,data=new.data,
             replace=list(pp=pj.zeros,y=yj.zeros))
  mu.PI <- apply(mu,2,PI,prob=0.75)
  shade(mu.PI,x,col=clrs2[[i]])
  mu.mean <- apply(mu,2,mean)
  lines(x,mu.mean,col=clrs[[i]],lwd=3)
  lines(0:8,0:8, lty=2, col="darkgrey")
}    

##### ADULT GROWTH
par(mfrow=c(2,3))

x <- seq(log(1),log(800),length.out = 100)

labelmain <- paste("Growth Adults")


##################################### Plot Bayesian

#  load("growth_adultsi.RData", .GlobalEnv)


pj.zeros <- matrix(0,1000,102)
yj.zeros <- matrix(0,1000,24)


#text(0.12,7.9,"D",cex=2)
for (i in 1:dim(trt)[1]){
  
  plot(log(x),log(seq(0,1,length.out=length(x))), xlab="ln (Number of Branches (t))",
       ylim=c(0,8),main=labPopTyp[i],xlim=c(0,8),type="n",ylab="ln (Number of Branches (t+1)",
       cex.lab=1.2)
  
  iniS<-dat4$log_ini_branches[dat4$stg==3 & dat4$trt_type==i]
  finS<-dat4$log_fin_branches[dat4$stg==3 & dat4$trt_type==i]
  
  points(iniS,finS,col=rgb(0.3,0.3,0.3,0.2),pch=16)
   
  v.trt2 <- trt[i,2]
  v.trt3 <- trt[i,3]
  v.trt4 <- trt[i,4]
  new.data <-  data.frame(cbind(rep(1,length(x)),
                                log_ini_branch = x,
                                T2 = rep(v.trt2,length(x)),
                                T3 =rep(v.trt3,length(x)),
                                T4 = rep(v.trt4,length(x)),
                                intr = rep(0,length(x)),
                                habitat = rep(hab[i],length(x)),
                                py=rep(1,length(x)),
                                yn=rep(1,length(x))))
  mu <- NULL
  mu <- link(model1,n=1000,data=new.data,
             replace=list(p=pj.zeros,y=yj.zeros))
  mu.PI <- apply(mu,2,PI,prob=0.95)
  shade(mu.PI,x,col=clrs2[[i]])
  mu.mean <- apply(mu,2,mean)
  lines(x,mu.mean,col=clrs[[i]],lwd=3)
  lines(0:8,0:8,lty=2,col="darkgrey")
}


########## P(reproductive) ############

par(mfrow=c(2,3))

pj.zeros <- matrix(0,1000,109)
yj.zeros <- matrix(0,1000,25)
x1 <- seq(log(1),log(800),length.out = 100)


  for (i in 1:dim(trt)[1]){

    plot(x[o],x[o], xlab="ln (Number of Branches)",
      ylim=c(0,1),main=labPopTyp[i],xlim=c(0, 8),type="n",ylab="Probability",cex.lab=1.2)
    
    iniS<-dat2p$log_ini_branches[dat2p$stg==3 & dat2p$trt_type==i & dat2p$burn==0]
    prepa<-dat2p$prep[dat2p$stg==3 & dat2p$trt_type==i & dat2p$burn==0]
    
    b <- c(-0.1,0.2,0.8,1.2,1.6,1.8,2,2.2,2.7,3,3.5,4,8)
    pos<-c(0,0.5,1,1.4,1.7,1.9,2.1,2.5,2.9,3.3,3.8,6)
    z <- cut(iniS, b)
    prebyden <- tapply(prepa,z,sum)
    tab <- table(z)
    probs <-prebyden/tab
    probs <- as.vector(probs)
    points(pos,probs,pch=16,cex=tab^(2/7),col=rgb(0.3,0.3,0.3,0.5))

    
    v.trt2 <- trt[i,2]
    v.trt3 <- trt[i,3]
    v.trt4 <- trt[i,4]
    new.data <-  data.frame(cbind(rep(1,length(x1)),
                                  log_ini_branch = x1,
                                  T2 = rep(v.trt2,length(x1)),
                                  T3 =rep(v.trt3,length(x1)),
                                  T4 = rep(v.trt4,length(x1)),
                                  stg_C = rep(0,length(x1)),
                                  aug1 = rep(0,length(x1)),
                                  aug2 = rep(0,length(x1)),
                                  habitat = rep(hab[i],length(x1)),
                                  preburn = rep(brn[i],length(x1)),
                                  burn = rep(0,length(x1)),
                                  py=1,
                                  yn=1))
    mu <- NULL
    mu <- link(model3,n=1000,data=new.data,
               replace=list(pp=pj.zeros,y=yj.zeros))
    mu.PI <- apply(mu,2,PI,prob=0.95)
    shade(mu.PI,x1,col=clrs2[[i]])
    mu.mean <- apply(mu,2,mean)
    lines(x1,mu.mean,col=clrs[[i]],lwd=3)
  }



#### FECUNDITY
par(mfrow=c(2,3))

ag <- c(0,1)
pb <- c(0,1)
br <- c(0,1)

x <- seq(log(1),log(800),length.out = 100)

for (i in 1:dim(trt)[1]){
  
  plot(x,x, xlab="ln (Number of Branches)",
       ylim=c(0,8),main=labPopTyp[i],xlim=c(0, 8),type="n",ylab="ln (Flowering Branches)",cex.lab=1.2)
  
  iniS<-dat3a$log_ini_branches[dat3a$stg==3 & dat3a$trt_type==i]
  flw<-dat3a$log_flw_branches[dat3a$stg==3 & dat3a$trt_type==i]
  
  
  points(iniS,flw,col=rgb(0.3,0.3,0.3,0.2),pch=16)

  
  v.trt2 <- trt[i,2]
  v.trt3 <- trt[i,3]
  v.trt4 <- trt[i,4]
  new.data <-  data.frame(cbind(rep(1,length(x)),
                                x,
                                rep(v.trt2,length(x)),
                                rep(v.trt3,length(x)),
                                rep(v.trt4,length(x)),
                                x*v.trt2,
                                x*v.trt3,
                                x*v.trt4,
                                rep(0,length(x)),
                                rep(0,length(x)),
                                rep(hab[i],length(x)),
                                rep(brn[i],length(x)),
                                rep(0,length(x)))
  )
  
  mu <- 1/(1 + (1/exp(calcpred(matrix_of_post,n=1000,new.data))))
  mu.PI <- apply(mu,2,PI,prob=0.95)
  shade(log(sweep(mu.PI, 2, exp(x), FUN = "*")),x,col=clrs2[[i]])
  mu.mean <- apply(mu,2,mean)
  #mu.PI <- apply(mu,2,PI,prob=0.95)
  lines(x,log(mu.mean*exp(x)),col=clrs[[i]],lwd=2.5)
  lines(0:9,0:9,lty=2,col="darkgrey") 
}




