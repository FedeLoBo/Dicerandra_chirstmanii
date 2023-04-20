##########################################################################
##        Dicerandra christmanii - vital rate modeling                  ##
##        Date: 06/11/2020    - Version: 19 dynamics for LTRE           ##
##        Pedro F. Quintana Ascencio & Federico Lopez Borghesi          ##
##########################################################################

rm(list=ls())
pkgs<-c("nlme","bbmle","lme4","lattice","MuMIn","MASS","dplyr",
        "optimx","rethinking","rstan","scales")
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


############################################################################
############################################################################
############################################### Main loop   ################

ttt <- c(1,2,2,1,1,3,4)
cells <- c(50,60,80,90,100,200,300)

hbtat<-c(1,1,1,1,1,1,1,1)
intrd<-c(1,1,0,0,0,0,0,0)
brn<-c(0,0,1,0,0,0,0,0)
pstbrn<-c(0,0,0,1,1,1,1,1)

iterations <- 200
years <- 100
yrint<-8

for (ki in 3){
  ns <- cells[ki]
  x <- seq(log(1),log(900),length.out = ns)  
  
  ### BURNED INTRODUCTION              
  mean.adlts.intb <- array(0,c(iterations ,length(x),yrint))
  mean.yrl.intb <- array(0,c(iterations ,length(x),yrint))
  mean.rep.intb <- array(0,c(iterations ,length(x),2,yrint))
  mean.fec.intb <- array(0,c(iterations ,length(x),yrint))
  mean.sur.intb <- array(0,c(iterations ,length(x),2,yrint))
  # NOT BURNED INTRODUCTION
  mean.adlts.intnb <- array(0,c(iterations ,length(x),yrint))
  mean.yrl.intnb <- array(0,c(iterations ,length(x),yrint))
  mean.rep.intnb <- array(0,c(iterations ,length(x),2,yrint))
  mean.fec.intnb <- array(0,c(iterations ,length(x),yrint))
  mean.sur.intnb <- array(0,c(iterations ,length(x),2,yrint))
  # GAP/ROAD
  mean.adults <- array(0,c(iterations ,length(x),years))
  mean.yearlings <- array(0,c(iterations ,length(x),years))
  mean.rep <- array(0,c(iterations ,length(x),2,years))
  mean.fec <- array(0,c(iterations ,length(x),years))
  mean.sur <- array(0,c(iterations ,length(x),2,years))
  
  trt <- array(0,c(4,4))
  # plot conditions
  trt[1,] <- c(0,0,0,0) #wild gaps  
  trt[2,] <- c(0,1,0,0)  # wild roads
  trt[3,] <- c(0,0,1,0) # augmentation
  trt[4,] <- c(0,0,0,1)  # introduction
  
  for(sqi in 1:yrint){
    
    print(sqi)
    
    ##########################  BURNED INTRODUCTION ##############################
    
    ########################### ADULT GROWTH ####################################
    
    # only for this subfile
    dat4 <- subset(dat2, stg ==3 & !is.na(ini_branches) & !is.na(fin_branches))
    dat4$burn_status[is.na(dat4$burn_status)] <- 0 
    dat4$gap[is.na(dat4$gap)] <- 0
    dat4$flw_branches[is.na(dat4$flw_branches)] <- 0 
    dat4$ann_sur[ which(is.na(dat4$ann_sur)) ] <- 0
    dat4 <- dat4[,c(-4)]
    row.has.na <- apply(dat4, 1, function(x1){ any(is.na(x1)) })
    sum(row.has.na)
    dat4$ID_num[ which(is.na(dat4$ann_sur)) ][order(dat4$ID_num[ which(is.na(dat4$ann_sur))])]
    dim(dat4)
    dat4$treatment <- as.factor(dat4$treatment)
    
    #############################################################################
    sum_model1 <- precis(model1,digits=3,depth=1)
    #############################################################################
    
    post1 <- extract.samples(model1)
    stoc <- rnorm(1,0,post1$sigmay[sqi])  # ask Pedro
    
    
    mat.stoch1 <- cbind(post1$a,post1$b1,post1$b22,post1$b23,post1$b24,post1$b32,post1$b33,
                        post1$b34,post1$b4,post1$b5)
    mat.stoch1 <- cbind(mat.stoch1[1:iterations,], rep(0,iterations),rep(stoc,iterations) )
    
    labelmain <- paste("Growth Adults")
    
    v.trt2 <- trt[ttt[7],2]
    v.trt3 <- trt[ttt[7],3]
    v.trt4 <- trt[ttt[7],4]
    
    new.data <-  data.frame(cbind(rep(1,length(x)),
                                  log_ini_branch = x,
                                  T2 = rep(v.trt2,length(x)),
                                  T3 =rep(v.trt3,length(x)),
                                  T4 = rep(v.trt4,length(x)),
                                  intr = rep(intrd[sqi],length(x)),
                                  habitat = rep(hbtat[sqi],length(x)),
                                  py=1,
                                  yn=1))
    
    mu <-array(0,c(iterations,cells[ki])) 
    new.data.t <- as.matrix(cbind(new.data[,1:5],new.data[,3:5]*new.data[,2],new.data[,6:7],array(1,c(cells[ki],2))))  
    for(wi in 1:iterations) {mu[wi,]<-as.numeric(new.data.t %*% mat.stoch1[wi,])}
    mean.adlts.intb[,,sqi] <-mu
    
    
    lammat<-link(model1,n=iterations) 
    pred.mean <- apply(lammat,2,mean)
    res.mean <- dat4$log_fin_branches-pred.mean
    m3i <- lme(log_fin_branches ~ log_ini_branches * treatment + intr, random= ~1|yn/sites, data=dat4)
    dat4$res2 <- resid(m3i)
    res.var <- summary(m3i)$sigma
    dat4$pz <-log(dat4$res2^2)
    dat4$pzB <- log(res.mean^2)
    grr <-lm(pzB~ log_ini_branches+treatment,data= dat4)
    g1 <- lm(pzB~1,data=dat4)
    
    cgra <-coef(grr)
    
    ########################### END ADULT GROWTH ####################################
    
    ########################### SEEDLING GROWTH #####################################
    
    
    dat3y <- subset(dat2, stg ==5 & !is.na(ini_branches)& !is.na(fin_branches))
    dat3y$treatment <- as.factor(dat3y$treatment)
    dat3y$log_ini_branches <- log(dat3y$ini_branches)
    dat3y$log_fin_branches <- log(dat3y$fin_branches)
    
    ## no seedlings in 7/2012 nor 6/2015
    
    dat3y$burn_status[is.na(dat3y$burn_status)] <- 0
    dat3y$gap[is.na(dat3y$gap)] <- 0
    dat3y <- dat3y[,c(-4)]
    dat3y <- na.omit(dat3y)
    row.has.na <- apply(dat3y, 1, function(x1){any(is.na(x1))})
    sum(row.has.na)
    
    #############################################################################
    sum_model2 <- precis(model2,digits=3,depth=1)
    #############################################################################
    
    
    
    #par(mfrow=c(1,1))
    labelmain <- paste("Growth Yearlings")
    
    post2 <- extract.samples(model2)
    stoc <- rnorm(1,0,post2$sigmay[sqi])
    
    
    mat.stoch2 <- cbind(post2$a,post2$b1,post2$b22,post2$b23,post2$b24,post2$b32,post2$b33,
                        post2$b34,post2$b4,post2$b5,post2$b6,post2$b7)
    mat.stoch2 <- cbind(mat.stoch2[1:iterations,], rep(0,iterations),rep(stoc,iterations) )
    
    
    v.trt2 <- trt[ttt[7],2]
    v.trt3 <- trt[ttt[7],3]
    v.trt4 <- trt[ttt[7],4]
    
    new.data <-  data.frame(cbind(rep(1,length(x)),
                                  log_ini_branches = x,
                                  T2 = rep(v.trt2,length(x)),
                                  T3 =rep(v.trt3,length(x)),
                                  T4 = rep(v.trt4,length(x)),
                                  habitat = rep(hbtat[sqi],length(x)),
                                  aug1 = rep(0,length(x)),
                                  aug2 = rep(0,length(x)),
                                  preburn = rep(0,length(x)),
                                  py=1,
                                  yn=1))
    
    mu <-array(0,c(iterations,cells[ki])) 
    new.data.t <- as.matrix(cbind(new.data[,1:5],new.data[,3:5]*new.data[,2],new.data[,6:9],array(1,c(cells[ki],2))))  
    for(wi in 1:iterations) {mu[wi,]<-as.numeric(new.data.t %*% mat.stoch2[wi,])}
    mean.yrl.intb[,,sqi] <-mu
    
    
    
    #############################################################################
    lammaty<-link(model2,n=iterations) 
    pred.meany <- apply(lammaty,2,mean)
    res.meany <- pred.meany-dat3y$log_fin_branches
    
    
    m3iy <- lme(log_fin_branches ~ log_ini_branches + treatment  , random= ~1|yn/sites, data=dat3y)
    dat3y$res2 <- resid(m3iy)
    
    
    res.vary <- summary(m3iy)$sigma
    dat3y$pz <-log(dat3y$res2^2)
    dat3y$pzBy <- log(res.meany^2)
    grry <-lm(pzBy~ log_ini_branches+treatment,data= dat3y)
    g1 <- lm(pzBy~1,data=dat3y)
    #AICctab(grry,g1)
    #summary(grry)
    cgry <-coef(grry)
    
    #par(mfrow=c(1,2))
    #plot(dat3sy$log_ini_branch,dat3sy$pz,pch=16,cex=0.5,ylim=c(-20,5))
    #plot(dat3sy$log_ini_branch,pzBy,pch=16,cex=0.5,ylim=c(-20,5))
    
    ########################### END YEARLING GROWTH ####################################
    
    ########################### PROBABILITY OF REPRODUCTION ############################
    
    
    dat2p <- subset(dat2,stg==3)
    dat2p$prep <- 0
    dat2p$prep[dat2p$flw_branches>0] <- 1
    dat2p$log_ini_branches <- log(dat2p$ini_branches)
    
    
    stg <- c(0,1) # 0-adult; 1:yearling
    stgname<-c("3","5")
    stsur<-c(3,5)
    labstg <- c("Adult","Yearling")
    o <- order(x)
    
    COL <- adjustcolor(c("black", "red", "green", "blue")[dat2p$treatment], alpha.f = 0.5)
    
    dat2$stg_C[dat2$stg==3] <- 0
    dat2$stg_C[dat2$stg==5] <- 1
    
    dat2p <- subset(dat2p, !is.na(log_ini_branches))
    
    #par(mfrow=c(2,2))
    
    #############################################################################
    sum_model3 <- precis(model3,digits=3,depth=1)
    #############################################################################
    
    post3 <- extract.samples(model3)
    stoc <- rnorm(1,0,post3$sigmay[sqi])
    
    
    mat.stoch3 <- cbind(post3$a,post3$b1,post3$b22,post3$b23,post3$b24,post3$b32,post3$b33,
                        post3$b34,post3$b4,post3$b5,post3$b6,post3$b7)
    mat.stoch3 <- cbind(mat.stoch3[1:iterations,], rep(0,iterations),rep(stoc,iterations) )
    
    
    for (j in 1:1){ 
      
      labelmain <- c("Probability of Reproduction")
      c.prep <- seq(min(dat2p$log_ini_branches[dat2p$stg == stsur[j] & dat2p$treatment==ttt[6]])-0.001, 
                    max(dat2p$log_ini_branches[dat2p$stg == stsur[j] & dat2p$treatment==ttt[6]]),0.3)
      v.prep <- cut(dat2p$log_ini_branches[dat2p$stg == stsur[j] & dat2p$treatment==ttt[6]],c.prep)
      t.prep <- table( dat2p$prep[dat2p$stg == stsur[j] & dat2p$treatment==ttt[6]],v.prep)
      if(sum(t.prep[1,]) != sum(colSums(t.prep))){
        p.prep <- t.prep[2,]/colSums(t.prep)}
      
      v.trt2 <- trt[ttt[7],2]
      v.trt3 <- trt[ttt[7],3]
      v.trt4 <- trt[ttt[7],4]
      new.data <-  data.frame(cbind(rep(1,length(x)),
                                    log_ini_branch = x,
                                    T2 = rep(v.trt2,length(x)),
                                    T3 =rep(v.trt3,length(x)),
                                    T4 = rep(v.trt4,length(x)),
                                    aug1 = rep(0,length(x)),
                                    aug2 = rep(0,length(x)),
                                    habitat = rep(hbtat[sqi],length(x)),
                                    preburn = rep(0,length(x)),
                                    py=1,
                                    yn=1))
      
      mu <-array(0,c(iterations,cells[ki])) 
      new.data.t <- as.matrix(cbind(new.data[,1:5],new.data[,3:5]*new.data[,2],new.data[,6:9],array(1,c(cells[ki],2))))  
      for(wi in 1:iterations) {mu[wi,]<-1/(1+1/exp(as.numeric(new.data.t %*% mat.stoch3[wi,])))}
      mean.rep.intb[,,j,sqi] <- mu
      
    }
    
    
    ########################### END PROBABILITY OF REPRODUCTION ########################
    
    ########################### FECUNDITY  #############################################
    # Number of flowering branches (Fecundity), Estimated originally as proportions
    ####################################################################################
    
    dat3a <- subset(dat2,flw_branches>0 & stg == "3" & !is.na(ini_branches) & !is.na(flw_branches ))
    dat3a$log_flw_branches <- log(dat3a$flw_branches)
    dat3a$log_ini_branches <- log(dat3a$ini_branches)
    
    dat3a $T2 <- 0
    dat3a $T2[dat3a $treatment==2] <- 1
    dat3a $T3 <- 0
    dat3a $T3[dat3a $treatment==3] <- 1
    dat3a $T4 <- 0
    dat3a $T4[dat3a $treatment==4] <- 1
    dat3a $yn <- dat3a $year-2011
    dat3a  <- dat3a [,c(-4,-5)]
    #dat3a  <- na.omit(dat3a )
    row.has.na <- apply(dat3a , 1, function(x1){any(is.na(x1))})
    sum(row.has.na)
    dim(dat3a)
    
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
    length(unique(dat3a$py))
    #unique(dat3a$py[order(dat3a$py)])
    dim(dat3a)
    
    dat3a  <- subset(dat3a ,!is.na(yn) & !is.na(py) & py>0)
    row.has.na <- apply(dat3a , 1, function(x1){any(is.na(x1))})
    sum(row.has.na)
    dim(dat3a)
    
    post4 <- extract.samples(model4)
    matrix_of_post <- as.matrix(model4)
    df_of_draws <- as.data.frame(model4)
    
    #############################################################################
    sum_model4 <- precis(model4,digits=3,depth=1)
    #############################################################################
    #####################################   Plot Model  ####################
    
    
    stoc <- rnorm(1,0,post4$sigmay[sqi])
    #COL <- adjustcolor(c("black", "red", "green", "blue")[dat4$treatment], alpha.f = 0.5)
    colt<- c("black", "red", "green", "blue")
    
    
    
    calcpred <-function(matrixs,iters,new.data,stoc) {
      preds<- array(0,c(iters,nrow(new.data)))
      for(i in 1:iters){
        mui <- as.matrix(new.data)%*%c(matrixs[i,c(1:13)],0, stoc)
        preds[i,] <- mui
      }
      return(preds)
    } 
    
    
    labelmain <- c("Fecundity") 
    #plot(x,x, xlab="ln (Number of Branches)",
    #     ylim=c(0,8),main=labelmain,xlim=c(0, 8),type="n",ylab="ln (Flowering Branches)",cex.lab=1.2)
    
    v.trt2 <- trt[ttt[7],2]
    v.trt3 <- trt[ttt[7],3]
    v.trt4 <- trt[ttt[7],4]
    new.data <-  data.frame(cbind(rep(1,length(x)),
                                  x,
                                  rep(v.trt2,length(x)), 
                                  rep(v.trt3,length(x)),
                                  rep(v.trt4,length(x)),
                                  x*v.trt2,
                                  x*v.trt3,
                                  x*v.trt4,
                                  rep(0,length(x)), # aug 1
                                  rep(0,length(x)), # aug 2
                                  rep(hbtat[sqi],length(x)), # habitat
                                  rep(0,length(x)), # preburn (aug)
                                  rep(brn[sqi],length(x)), # burn
                                  rep(1,length(x)),
                                  rep(1,length(x)))
    )
    
    mu <- 1/(1 + (1/exp( calcpred(matrix_of_post,iterations,new.data, stoc) ) ))
    mean.fec.intb[,,sqi] <-mu
    
    
    
    
    ########################### END FECUNDITY ##########################################
    
    ########################### SURVIVAL ###############################################
    
    
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
    row.has.na <- apply(dat5, 1, function(x1){any(is.na(x1))})
    sum(row.has.na)
    dat5$py <- 0
    k <-0
    
    sites <- table(dat5$sites,dat5$year)
    dsites <- dim(sites) 
    
    for (i in 1:25){
      for(j in 1:7){
        k <- k+1
        #print(c(i,j))
        dat5$py[dat5$sites==j & dat5$yn==i]  <- k
      }
    }
    
    dat5$stg_C[dat5$stg==3] <- 0
    dat5$stg_C[dat5$stg==5] <- 1
    dat5 <- subset(dat5,!is.na(yn) & !is.na(py))
    row.has.na <- apply(dat5, 1, function(x1){any(is.na(x1))})
    sum(row.has.na)
    
    dat5$treatment <- as.factor(dat5$treatment)
    dat5$stg <- as.factor(dat5$stg)
    dat5$log_ini_branches <- log(dat5$ini_branches)
    dat5$logsquare <- dat5$log_ini_branches^2
    dat5 <- subset(dat5, !is.na(log_ini_branches))
    
    dat5$burn[dat5$year!=2014 & dat5$treatment==4] <- 0
    dat5$postburn <- 0
    dat5$postburn[dat5$year> 2014 & dat5$treatment==4] <- 1
    
    #############################################################################
    sum_model5 <- precis(model5,digits=3,depth=1)
    #############################################################################
    
    #pj.zeros <- matrix(0,1000,12)
    #yj.zeros <- matrix(0,1000,25)
    
    
    post5 <- extract.samples(model5)
    stoc <- rnorm(1,0,post5$sigmay[sqi])
    
    
    mat.stoch5 <- cbind(post5$a,post5$b11,post5$b12,post5$b22,post5$b23,post5$b24,post5$b32,post5$b33,
                        post5$b34,post5$b41,post5$b42,post5$b43, post5$b52,post5$b53,post5$b54,
                        post5$b6,post5$b7,post5$b8,post5$b9,post5$b10,post5$b111  )
    mat.stoch5 <- cbind(mat.stoch5[1:iterations,], rep(0,iterations),rep(stoc,iterations) )
    
    
    for (j in 1:length(stg)){
      #par(mfrow=c(2,2))
      labelmain <- c("Survival",labstg[j])
      v.trt2 <- trt[ttt[7],2]
      v.trt3 <- trt[ttt[7],3]
      v.trt4 <- trt[ttt[7],4]
      new.data <-  data.frame(cbind(rep(1,length(x)),
                                    log_ini_branches = x,
                                    logsquare = x^2,
                                    T2 = rep(v.trt2,length(x)),
                                    T3 = rep(v.trt3,length(x)),
                                    T4 = rep(v.trt4,length(x)),
                                    stg_C= stg[j],
                                    habitat = rep(hbtat[sqi],length(x)),
                                    aug = rep(0,length(x)), 
                                    burn = rep(brn[sqi],length(x)),
                                    preburn = rep(0,length(x)),
                                    postburn = rep(pstbrn[sqi],length(x)),
                                    py=1,
                                    yn=1))
      #change aug 
      
      mu <-array(0,c(iterations,cells[ki])) 
      new.data.t <- as.matrix(cbind(new.data[,1:6],new.data[,4:6]*new.data[,2],new.data[,7],new.data[,7]*new.data[,2:3], 
                                    new.data[,4:6]*new.data[,3],new.data[,8:12],new.data[,8]*new.data[,10] ,array(1,c(cells[ki],2))))  
      for(wi in 1:iterations) {mu[wi,]<-1/(1+1/exp(as.numeric(new.data.t %*% mat.stoch5[wi,])))}
      mean.sur.intb[,,j,sqi] <- mu
      
      
    }
    
    
  } ##### sqi
  
  for(sqi in 1:yrint){
    
    print(sqi)
    
    
    
    ########################## NOT BURNED INTRODUCTION ##############################
    
    ########################### ADULT GROWTH ####################################
    
    # only for this subfile
    dat4 <- subset(dat2, stg ==3 & !is.na(ini_branches) & !is.na(fin_branches))
    dat4$burn_status[is.na(dat4$burn_status)] <- 0 
    dat4$gap[is.na(dat4$gap)] <- 0
    dat4$flw_branches[is.na(dat4$flw_branches)] <- 0 
    dat4$ann_sur[ which(is.na(dat4$ann_sur)) ] <- 0
    dat4 <- dat4[,c(-4)]
    row.has.na <- apply(dat4, 1, function(x1){ any(is.na(x1)) })
    sum(row.has.na)
    dat4$ID_num[ which(is.na(dat4$ann_sur)) ][order(dat4$ID_num[ which(is.na(dat4$ann_sur))])]
    dim(dat4)
    dat4$treatment <- as.factor(dat4$treatment)
    
    #############################################################################
    sum_model1 <- precis(model1,digits=3,depth=1)
    #############################################################################
    
    post1 <- extract.samples(model1)
    stoc <- rnorm(1,0,post1$sigmay[sqi])  # ask Pedro
    
    
    mat.stoch1 <- cbind(post1$a,post1$b1,post1$b22,post1$b23,post1$b24,post1$b32,post1$b33,
                        post1$b34,post1$b4,post1$b5)
    mat.stoch1 <- cbind(mat.stoch1[1:iterations,], rep(0,iterations),rep(stoc,iterations) )
    
    labelmain <- paste("Growth Adults")
    
    v.trt2 <- trt[ttt[7],2]
    v.trt3 <- trt[ttt[7],3]
    v.trt4 <- trt[ttt[7],4]
    
    new.data <-  data.frame(cbind(rep(1,length(x)),
                                  log_ini_branch = x,
                                  T2 = rep(v.trt2,length(x)),
                                  T3 =rep(v.trt3,length(x)),
                                  T4 = rep(v.trt4,length(x)),
                                  intr = rep(intrd[sqi],length(x)),
                                  habitat = rep(hbtat[sqi],length(x)),
                                  py=1,
                                  yn=1))
    
    mu <-array(0,c(iterations,cells[ki])) 
    new.data.t <- as.matrix(cbind(new.data[,1:5],new.data[,3:5]*new.data[,2],new.data[,6:7],array(1,c(cells[ki],2))))  
    for(wi in 1:iterations) {mu[wi,]<-as.numeric(new.data.t %*% mat.stoch1[wi,])}
    mean.adlts.intnb[,,sqi] <-mu
    
    
    lammat<-link(model1,n=iterations) 
    pred.mean <- apply(lammat,2,mean)
    res.mean <- dat4$log_fin_branches-pred.mean
    m3i <- lme(log_fin_branches ~ log_ini_branches * treatment + intr, random= ~1|yn/sites, data=dat4)
    dat4$res2 <- resid(m3i)
    res.var <- summary(m3i)$sigma
    dat4$pz <-log(dat4$res2^2)
    dat4$pzB <- log(res.mean^2)
    grr <-lm(pzB~ log_ini_branches+treatment,data= dat4)
    g1 <- lm(pzB~1,data=dat4)
    
    cgra <-coef(grr)
    
    ########################### END ADULT GROWTH ####################################
    
    ########################### SEEDLING GROWTH #####################################
    
    
    dat3y <- subset(dat2, stg ==5 & !is.na(ini_branches)& !is.na(fin_branches))
    dat3y$treatment <- as.factor(dat3y$treatment)
    dat3y$log_ini_branches <- log(dat3y$ini_branches)
    dat3y$log_fin_branches <- log(dat3y$fin_branches)
    
    ## no seedlings in 7/2012 nor 6/2015
    
    dat3y$burn_status[is.na(dat3y$burn_status)] <- 0
    dat3y$gap[is.na(dat3y$gap)] <- 0
    dat3y <- dat3y[,c(-4)]
    dat3y <- na.omit(dat3y)
    row.has.na <- apply(dat3y, 1, function(x1){any(is.na(x1))})
    sum(row.has.na)
    
    #############################################################################
    sum_model2 <- precis(model2,digits=3,depth=1)
    #############################################################################
    
    
    
    #par(mfrow=c(1,1))
    labelmain <- paste("Growth Yearlings")
    
    post2 <- extract.samples(model2)
    stoc <- rnorm(1,0,post2$sigmay[sqi])
    
    
    mat.stoch2 <- cbind(post2$a,post2$b1,post2$b22,post2$b23,post2$b24,post2$b32,post2$b33,
                        post2$b34,post2$b4,post2$b5,post2$b6,post2$b7)
    mat.stoch2 <- cbind(mat.stoch2[1:iterations,], rep(0,iterations),rep(stoc,iterations) )
    
    
    v.trt2 <- trt[ttt[7],2]
    v.trt3 <- trt[ttt[7],3]
    v.trt4 <- trt[ttt[7],4]
    
    new.data <-  data.frame(cbind(rep(1,length(x)),
                                  log_ini_branches = x,
                                  T2 = rep(v.trt2,length(x)),
                                  T3 =rep(v.trt3,length(x)),
                                  T4 = rep(v.trt4,length(x)),
                                  habitat = rep(hbtat[sqi],length(x)),
                                  aug1 = rep(0,length(x)),
                                  aug2 = rep(0,length(x)),
                                  preburn = rep(0,length(x)),
                                  py=1,
                                  yn=1))
    
    mu <-array(0,c(iterations,cells[ki])) 
    new.data.t <- as.matrix(cbind(new.data[,1:5],new.data[,3:5]*new.data[,2],new.data[,6:9],array(1,c(cells[ki],2))))  
    for(wi in 1:iterations) {mu[wi,]<-as.numeric(new.data.t %*% mat.stoch2[wi,])}
    mean.yrl.intnb[,,sqi] <-mu
    
    
    
    #############################################################################
    lammaty<-link(model2,n=iterations) 
    pred.meany <- apply(lammaty,2,mean)
    res.meany <- pred.meany-dat3y$log_fin_branches
    
    
    m3iy <- lme(log_fin_branches ~ log_ini_branches + treatment  , random= ~1|yn/sites, data=dat3y)
    dat3y$res2 <- resid(m3iy)
    
    
    res.vary <- summary(m3iy)$sigma
    dat3y$pz <-log(dat3y$res2^2)
    dat3y$pzBy <- log(res.meany^2)
    grry <-lm(pzBy~ log_ini_branches+treatment,data= dat3y)
    g1 <- lm(pzBy~1,data=dat3y)
    #AICctab(grry,g1)
    #summary(grry)
    cgry <-coef(grry)
    
    #par(mfrow=c(1,2))
    #plot(dat3sy$log_ini_branch,dat3sy$pz,pch=16,cex=0.5,ylim=c(-20,5))
    #plot(dat3sy$log_ini_branch,pzBy,pch=16,cex=0.5,ylim=c(-20,5))
    
    ########################### END YEARLING GROWTH ####################################
    
    ########################### PROBABILITY OF REPRODUCTION ############################
    
    
    dat2p <- subset(dat2,stg==3)
    dat2p$prep <- 0
    dat2p$prep[dat2p$flw_branches>0] <- 1
    dat2p$log_ini_branches <- log(dat2p$ini_branches)
    
    
    stg <- c(0,1) # 0-adult; 1:yearling
    stgname<-c("3","5")
    stsur<-c(3,5)
    labstg <- c("Adult","Yearling")
    o <- order(x)
    
    COL <- adjustcolor(c("black", "red", "green", "blue")[dat2p$treatment], alpha.f = 0.5)
    
    dat2$stg_C[dat2$stg==3] <- 0
    dat2$stg_C[dat2$stg==5] <- 1
    
    dat2p <- subset(dat2p, !is.na(log_ini_branches))
    
    #par(mfrow=c(2,2))
    
    #############################################################################
    sum_model3 <- precis(model3,digits=3,depth=1)
    #############################################################################
    
    post3 <- extract.samples(model3)
    stoc <- rnorm(1,0,post3$sigmay[sqi])
    
    
    mat.stoch3 <- cbind(post3$a,post3$b1,post3$b22,post3$b23,post3$b24,post3$b32,post3$b33,
                        post3$b34,post3$b4,post3$b5,post3$b6,post3$b7)
    mat.stoch3 <- cbind(mat.stoch3[1:iterations,], rep(0,iterations),rep(stoc,iterations) )
    
    
    for (j in 1:1){ 
      
      labelmain <- c("Probability of Reproduction")
      c.prep <- seq(min(dat2p$log_ini_branches[dat2p$stg == stsur[j] & dat2p$treatment==ttt[6]])-0.001, 
                    max(dat2p$log_ini_branches[dat2p$stg == stsur[j] & dat2p$treatment==ttt[6]]),0.3)
      v.prep <- cut(dat2p$log_ini_branches[dat2p$stg == stsur[j] & dat2p$treatment==ttt[6]],c.prep)
      t.prep <- table( dat2p$prep[dat2p$stg == stsur[j] & dat2p$treatment==ttt[6]],v.prep)
      if(sum(t.prep[1,]) != sum(colSums(t.prep))){
        p.prep <- t.prep[2,]/colSums(t.prep)}
      
      v.trt2 <- trt[ttt[7],2]
      v.trt3 <- trt[ttt[7],3]
      v.trt4 <- trt[ttt[7],4]
      new.data <-  data.frame(cbind(rep(1,length(x)),
                                    log_ini_branch = x,
                                    T2 = rep(v.trt2,length(x)),
                                    T3 =rep(v.trt3,length(x)),
                                    T4 = rep(v.trt4,length(x)),
                                    aug1 = rep(0,length(x)),
                                    aug2 = rep(0,length(x)),
                                    habitat = rep(hbtat[sqi],length(x)),
                                    preburn = rep(0,length(x)),
                                    py=1,
                                    yn=1))
      
      mu <-array(0,c(iterations,cells[ki])) 
      new.data.t <- as.matrix(cbind(new.data[,1:5],new.data[,3:5]*new.data[,2],new.data[,6:9],array(1,c(cells[ki],2))))  
      for(wi in 1:iterations) {mu[wi,]<-1/(1+1/exp(as.numeric(new.data.t %*% mat.stoch3[wi,])))}
      mean.rep.intnb[,,j,sqi] <- mu
      
    }
    
    
    ########################### END PROBABILITY OF REPRODUCTION ########################
    
    ########################### FECUNDITY  #############################################
    # Number of flowering branches (Fecundity), Estimated originally as proportions
    ####################################################################################
    
    dat3a <- subset(dat2,flw_branches>0 & stg == "3" & !is.na(ini_branches) & !is.na(flw_branches ))
    dat3a$log_flw_branches <- log(dat3a$flw_branches)
    dat3a$log_ini_branches <- log(dat3a$ini_branches)
    
    dat3a $T2 <- 0
    dat3a $T2[dat3a $treatment==2] <- 1
    dat3a $T3 <- 0
    dat3a $T3[dat3a $treatment==3] <- 1
    dat3a $T4 <- 0
    dat3a $T4[dat3a $treatment==4] <- 1
    dat3a $yn <- dat3a $year-2011
    dat3a  <- dat3a [,c(-4,-5)]
    #dat3a  <- na.omit(dat3a )
    row.has.na <- apply(dat3a , 1, function(x1){any(is.na(x1))})
    sum(row.has.na)
    dim(dat3a)
    
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
    length(unique(dat3a$py))
    #unique(dat3a$py[order(dat3a$py)])
    dim(dat3a)
    
    dat3a  <- subset(dat3a ,!is.na(yn) & !is.na(py) & py>0)
    row.has.na <- apply(dat3a , 1, function(x1){any(is.na(x1))})
    sum(row.has.na)
    dim(dat3a)
    
    post4 <- extract.samples(model4)
    matrix_of_post <- as.matrix(model4)
    df_of_draws <- as.data.frame(model4)
    
    #############################################################################
    sum_model4 <- precis(model4,digits=3,depth=1)
    #############################################################################
    #####################################   Plot Model  ####################
    
    
    stoc <- rnorm(1,0,post4$sigmay[sqi])
    #COL <- adjustcolor(c("black", "red", "green", "blue")[dat4$treatment], alpha.f = 0.5)
    colt<- c("black", "red", "green", "blue")
    
    
    
    calcpred <-function(matrixs,iters,new.data,stoc) {
      preds<- array(0,c(iters,nrow(new.data)))
      for(i in 1:iters){
        mui <- as.matrix(new.data)%*%c(matrixs[i,c(1:13)],0, stoc)
        preds[i,] <- mui
      }
      return(preds)
    } 
    
    
    labelmain <- c("Fecundity") 
    #plot(x,x, xlab="ln (Number of Branches)",
    #     ylim=c(0,8),main=labelmain,xlim=c(0, 8),type="n",ylab="ln (Flowering Branches)",cex.lab=1.2)
    
    v.trt2 <- trt[ttt[7],2]
    v.trt3 <- trt[ttt[7],3]
    v.trt4 <- trt[ttt[7],4]
    new.data <-  data.frame(cbind(rep(1,length(x)),
                                  x,
                                  rep(v.trt2,length(x)), 
                                  rep(v.trt3,length(x)),
                                  rep(v.trt4,length(x)),
                                  x*v.trt2,
                                  x*v.trt3,
                                  x*v.trt4,
                                  rep(0,length(x)), # aug 1
                                  rep(0,length(x)), # aug 2
                                  rep(hbtat[sqi],length(x)), # habitat
                                  rep(0,length(x)), # preburn (aug)
                                  rep(0,length(x)), # burn
                                  rep(1,length(x)),
                                  rep(1,length(x)))
    )
    
    mu <- 1/(1 + (1/exp( calcpred(matrix_of_post,iterations,new.data, stoc) ) ))
    mean.fec.intnb[,,sqi] <-mu
    
    
    
    
    ########################### END FECUNDITY ##########################################
    
    ########################### SURVIVAL ###############################################
    
    
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
    row.has.na <- apply(dat5, 1, function(x1){any(is.na(x1))})
    sum(row.has.na)
    dat5$py <- 0
    k <-0
    
    sites <- table(dat5$sites,dat5$year)
    dsites <- dim(sites) 
    
    for (i in 1:25){
      for(j in 1:7){
        k <- k+1
        #print(c(i,j))
        dat5$py[dat5$sites==j & dat5$yn==i]  <- k
      }
    }
    
    dat5$stg_C[dat5$stg==3] <- 0
    dat5$stg_C[dat5$stg==5] <- 1
    dat5 <- subset(dat5,!is.na(yn) & !is.na(py))
    row.has.na <- apply(dat5, 1, function(x1){any(is.na(x1))})
    sum(row.has.na)
    
    dat5$treatment <- as.factor(dat5$treatment)
    dat5$stg <- as.factor(dat5$stg)
    dat5$log_ini_branches <- log(dat5$ini_branches)
    dat5$logsquare <- dat5$log_ini_branches^2
    dat5 <- subset(dat5, !is.na(log_ini_branches))
    
    dat5$burn[dat5$year!=2014 & dat5$treatment==4] <- 0
    dat5$postburn <- 0
    dat5$postburn[dat5$year> 2014 & dat5$treatment==4] <- 1
    
    #############################################################################
    sum_model5 <- precis(model5,digits=3,depth=1)
    #############################################################################
    
    #pj.zeros <- matrix(0,1000,12)
    #yj.zeros <- matrix(0,1000,25)
    
    
    post5 <- extract.samples(model5)
    stoc <- rnorm(1,0,post5$sigmay[sqi])
    
    
    mat.stoch5 <- cbind(post5$a,post5$b11,post5$b12,post5$b22,post5$b23,post5$b24,post5$b32,post5$b33,
                        post5$b34,post5$b41,post5$b42,post5$b43, post5$b52,post5$b53,post5$b54,
                        post5$b6,post5$b7,post5$b8,post5$b9,post5$b10,post5$b111  )
    mat.stoch5 <- cbind(mat.stoch5[1:iterations,], rep(0,iterations),rep(stoc,iterations) )
    
    
    for (j in 1:length(stg)){
      #par(mfrow=c(2,2))
      labelmain <- c("Survival",labstg[j])
      v.trt2 <- trt[ttt[7],2]
      v.trt3 <- trt[ttt[7],3]
      v.trt4 <- trt[ttt[7],4]
      new.data <-  data.frame(cbind(rep(1,length(x)),
                                    log_ini_branches = x,
                                    logsquare = x^2,
                                    T2 = rep(v.trt2,length(x)),
                                    T3 = rep(v.trt3,length(x)),
                                    T4 = rep(v.trt4,length(x)),
                                    stg_C= stg[j],
                                    habitat = rep(hbtat[sqi],length(x)),
                                    aug = rep(0,length(x)), 
                                    burn = rep(0,length(x)),
                                    preburn = rep(0,length(x)),
                                    postburn = rep(0,length(x)),
                                    py=1,
                                    yn=1))
      #change aug 
      
      mu <-array(0,c(iterations,cells[ki])) 
      new.data.t <- as.matrix(cbind(new.data[,1:6],new.data[,4:6]*new.data[,2],new.data[,7],new.data[,7]*new.data[,2:3], 
                                    new.data[,4:6]*new.data[,3],new.data[,8:12],new.data[,8]*new.data[,10] ,array(1,c(cells[ki],2))))  
      for(wi in 1:iterations) {mu[wi,]<-1/(1+1/exp(as.numeric(new.data.t %*% mat.stoch5[wi,])))}
      mean.sur.intnb[,,j,sqi] <- mu
      
      
    }
    
    
  } ##### sqi
  
  for(sqi in 1:years){
    
    print(sqi)
    
    
    ##############################      ROAD      ###############################
    
    ########################### ADULT GROWTH ####################################
    
    # only for this subfile
    dat4 <- subset(dat2, stg ==3 & !is.na(ini_branches) & !is.na(fin_branches))
    dat4$burn_status[is.na(dat4$burn_status)] <- 0 
    dat4$gap[is.na(dat4$gap)] <- 0
    dat4$flw_branches[is.na(dat4$flw_branches)] <- 0 
    dat4$ann_sur[ which(is.na(dat4$ann_sur)) ] <- 0
    dat4 <- dat4[,c(-4)]
    row.has.na <- apply(dat4, 1, function(x1){ any(is.na(x1)) })
    sum(row.has.na)
    dat4$ID_num[ which(is.na(dat4$ann_sur)) ][order(dat4$ID_num[ which(is.na(dat4$ann_sur))])]
    dim(dat4)
    dat4$treatment <- as.factor(dat4$treatment)
    
    #############################################################################
    sum_model1 <- precis(model1,digits=3,depth=1)
    #############################################################################
    
    post1 <- extract.samples(model1)
    stoc <- rnorm(1,0,post1$sigmay[sqi])
    
    
    mat.stoch1 <- cbind(post1$a,post1$b1,post1$b22,post1$b23,post1$b24,post1$b32,post1$b33,
                        post1$b34,post1$b4,post1$b5)
    mat.stoch1 <- cbind(mat.stoch1[1:iterations,], rep(0,iterations),rep(stoc,iterations) )
    
    labelmain <- paste("Growth Adults")
    
    v.trt2 <- trt[ttt[2],2]
    v.trt3 <- trt[ttt[2],3]
    v.trt4 <- trt[ttt[2],4]
    new.data <-  data.frame(cbind(rep(1,length(x)),
                                  log_ini_branch = x,
                                  T2 = rep(v.trt2,length(x)),
                                  T3 =rep(v.trt3,length(x)),
                                  T4 = rep(v.trt4,length(x)),
                                  intr = rep(0,length(x)),
                                  habitat = rep(0,length(x)),
                                  py=1,
                                  yn=1))
    
    mu <-array(0,c(iterations,cells[ki])) 
    new.data.t <- as.matrix(cbind(new.data[,1:5],new.data[,3:5]*new.data[,2],new.data[,6:7],array(1,c(cells[ki],2))))  
    for(wi in 1:iterations) {mu[wi,]<-as.numeric(new.data.t %*% mat.stoch1[wi,])}
    mean.adults[,,sqi] <-mu
    
    
    lammat<-link(model1,n=iterations) 
    pred.mean <- apply(lammat,2,mean)
    res.mean <- dat4$log_fin_branches-pred.mean
    m3i <- lme(log_fin_branches ~ log_ini_branches * treatment + intr, random= ~1|yn/sites, data=dat4)
    dat4$res2 <- resid(m3i)
    res.var <- summary(m3i)$sigma
    dat4$pz <-log(dat4$res2^2)
    dat4$pzB <- log(res.mean^2)
    grr <-lm(pzB~ log_ini_branches+treatment,data= dat4)
    g1 <- lm(pzB~1,data=dat4)
    
    cgra <-coef(grr)
    
    ########################### END ADULT GROWTH ####################################
    
    ########################### SEEDLING GROWTH #####################################
    
    
    dat3y <- subset(dat2, stg ==5 & !is.na(ini_branches)& !is.na(fin_branches))
    dat3y$treatment <- as.factor(dat3y$treatment)
    dat3y$log_ini_branches <- log(dat3y$ini_branches)
    dat3y$log_fin_branches <- log(dat3y$fin_branches)
    
    ## no seedlings in 7/2012 nor 6/2015
    
    dat3y$burn_status[is.na(dat3y$burn_status)] <- 0
    dat3y$gap[is.na(dat3y$gap)] <- 0
    dat3y <- dat3y[,c(-4)]
    dat3y <- na.omit(dat3y)
    row.has.na <- apply(dat3y, 1, function(x1){any(is.na(x1))})
    sum(row.has.na)
    
    #############################################################################
    sum_model2 <- precis(model2,digits=3,depth=1)
    #############################################################################
    
    
    
    #par(mfrow=c(1,1))
    labelmain <- paste("Growth Yearlings")
    
    post2 <- extract.samples(model2)
    stoc <- rnorm(1,0,post2$sigmay[sqi])
    
    
    mat.stoch2 <- cbind(post2$a,post2$b1,post2$b22,post2$b23,post2$b24,post2$b32,post2$b33,
                        post2$b34,post2$b4,post2$b5,post2$b6,post2$b7)
    mat.stoch2 <- cbind(mat.stoch2[1:iterations,], rep(0,iterations),rep(stoc,iterations) )
    
    
    v.trt2 <- trt[ttt[2],2]
    v.trt3 <- trt[ttt[2],3]
    v.trt4 <- trt[ttt[2],4]
    new.data <-  data.frame(cbind(rep(1,length(x)),
                                  log_ini_branches = x,
                                  T2 = rep(v.trt2,length(x)),
                                  T3 =rep(v.trt3,length(x)),
                                  T4 = rep(v.trt4,length(x)),
                                  habitat = rep(0,length(x)),
                                  aug1 = rep(0,length(x)),
                                  aug2 = rep(0,length(x)),
                                  preburn = rep(0,length(x)),
                                  py=1,
                                  yn=1))
    
    mu <-array(0,c(iterations,cells[ki])) 
    new.data.t <- as.matrix(cbind(new.data[,1:5],new.data[,3:5]*new.data[,2],new.data[,6:9],array(1,c(cells[ki],2))))  
    for(wi in 1:iterations) {mu[wi,]<-as.numeric(new.data.t %*% mat.stoch2[wi,])}
    mean.yearlings[,,sqi] <-mu
    
    
    
    #############################################################################
    lammaty<-link(model2,n=iterations) 
    pred.meany <- apply(lammaty,2,mean)
    res.meany <- pred.meany-dat3y$log_fin_branches
    
    
    m3iy <- lme(log_fin_branches ~ log_ini_branches + treatment  , random= ~1|yn/sites, data=dat3y)
    dat3y$res2 <- resid(m3iy)
    
    
    res.vary <- summary(m3iy)$sigma
    dat3y$pz <-log(dat3y$res2^2)
    dat3y$pzBy <- log(res.meany^2)
    grry <-lm(pzBy~ log_ini_branches+treatment,data= dat3y)
    g1 <- lm(pzBy~1,data=dat3y)
    #AICctab(grry,g1)
    #summary(grry)
    cgry <-coef(grry)
    
    #par(mfrow=c(1,2))
    #plot(dat3sy$log_ini_branch,dat3sy$pz,pch=16,cex=0.5,ylim=c(-20,5))
    #plot(dat3sy$log_ini_branch,pzBy,pch=16,cex=0.5,ylim=c(-20,5))
    
    ########################### END YEARLING GROWTH ####################################
    
    ########################### PROBABILITY OF REPRODUCTION ############################
    
    
    dat2p <- subset(dat2,stg==3)
    dat2p$prep <- 0
    dat2p$prep[dat2p$flw_branches>0] <- 1
    dat2p$log_ini_branches <- log(dat2p$ini_branches)
    
    
    stg <- c(0,1) # 0-adult; 1:yearling
    stgname<-c("3","5")
    stsur<-c(3,5)
    labstg <- c("Adult","Yearling")
    o <- order(x)
    
    COL <- adjustcolor(c("black", "red", "green", "blue")[dat2p$treatment], alpha.f = 0.5)
    
    dat2$stg_C[dat2$stg==3] <- 0
    dat2$stg_C[dat2$stg==5] <- 1
    
    dat2p <- subset(dat2p, !is.na(log_ini_branches))
    
    #par(mfrow=c(2,2))
    
    #############################################################################
    sum_model3 <- precis(model3,digits=3,depth=1)
    #############################################################################
    
    #pj.zeros <- matrix(0,1000,109)
    #yj.zeros <- matrix(0,1000,25)
    
    
    post3 <- extract.samples(model3)
    stoc <- rnorm(1,0,post3$sigmay[sqi])
    
    
    mat.stoch3 <- cbind(post3$a,post3$b1,post3$b22,post3$b23,post3$b24,post3$b32,post3$b33,
                        post3$b34,post3$b4,post3$b5,post3$b6,post3$b7)
    mat.stoch3 <- cbind(mat.stoch3[1:iterations,], rep(0,iterations),rep(stoc,iterations) )
    
    
    for (j in 1:1){ 
      
      labelmain <- c("Probability of Reproduction")
      c.prep <- seq(min(dat2p$log_ini_branches[dat2p$stg == stsur[j] & dat2p$treatment==ttt[1]])-0.001, 
                    max(dat2p$log_ini_branches[dat2p$stg == stsur[j] & dat2p$treatment==ttt[1]]),0.3)
      v.prep <- cut(dat2p$log_ini_branches[dat2p$stg == stsur[j] & dat2p$treatment==ttt[1]],c.prep)
      t.prep <- table( dat2p$prep[dat2p$stg == stsur[j] & dat2p$treatment==ttt[1]],v.prep)
      if(sum(t.prep[1,]) != sum(colSums(t.prep))){
        p.prep <- t.prep[2,]/colSums(t.prep)}
      #plot(x[o],x[o], xlab="ln (Number of Branches)",
      #     ylim=c(0,1),main=labelmain,xlim=c(0, 8),type="n",ylab="Probability of reproduction",cex.lab=1.2)
      #abline(v=mean(dat2p$log_ini_branches[dat2p$stg==stgname[j]]),lty=3)
      #abline(v=max(dat2p$log_ini_branches[dat2p$stg==stgname[j]]),lty=3,col="red") }
      
      
      v.trt2 <- trt[ttt[2],2]
      v.trt3 <- trt[ttt[2],3]
      v.trt4 <- trt[ttt[2],4]
      new.data <-  data.frame(cbind(rep(1,length(x)),
                                    log_ini_branch = x,
                                    T2 = rep(v.trt2,length(x)),
                                    T3 =rep(v.trt3,length(x)),
                                    T4 = rep(v.trt4,length(x)),
                                    aug1 = rep(0,length(x)),
                                    aug2 = rep(0,length(x)),
                                    habitat = rep(0,length(x)),
                                    preburn = rep(0,length(x)),
                                    py=1,
                                    yn=1))
      
      mu <-array(0,c(iterations,cells[ki])) 
      new.data.t <- as.matrix(cbind(new.data[,1:5],new.data[,3:5]*new.data[,2],new.data[,6:9],array(1,c(cells[ki],2))))  
      for(wi in 1:iterations) {mu[wi,]<-1/(1+1/exp(as.numeric(new.data.t %*% mat.stoch3[wi,])))}
      mean.rep[,,j,sqi] <- mu
      
    }
    
    
    ########################### END PROBABILITY OF REPRODUCTION ########################
    
    ########################### FECUNDITY  #############################################
    # Number of flowering branches (Fecundity), Estimated originally as proportions
    ####################################################################################
    
    dat3a <- subset(dat2,flw_branches>0 & stg == "3" & !is.na(ini_branches) & !is.na(flw_branches ))
    dat3a$log_flw_branches <- log(dat3a$flw_branches)
    dat3a$log_ini_branches <- log(dat3a$ini_branches)
    
    dat3a $T2 <- 0
    dat3a $T2[dat3a $treatment==2] <- 1
    dat3a $T3 <- 0
    dat3a $T3[dat3a $treatment==3] <- 1
    dat3a $T4 <- 0
    dat3a $T4[dat3a $treatment==4] <- 1
    dat3a $yn <- dat3a $year-2011
    dat3a  <- dat3a [,c(-4,-5)]
    #dat3a  <- na.omit(dat3a )
    row.has.na <- apply(dat3a , 1, function(x1){any(is.na(x1))})
    sum(row.has.na)
    dim(dat3a)
    
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
    length(unique(dat3a$py))
    #unique(dat3a$py[order(dat3a$py)])
    dim(dat3a)
    
    dat3a  <- subset(dat3a ,!is.na(yn) & !is.na(py) & py>0)
    row.has.na <- apply(dat3a , 1, function(x1){any(is.na(x1))})
    sum(row.has.na)
    dim(dat3a)
    
    post4 <- extract(model4)
    matrix_of_post <- as.matrix(model4)
    df_of_draws <- as.data.frame(model4)
    
    #############################################################################
    sum_model4 <- precis(model4,digits=3,depth=1)
    #############################################################################
    #####################################   Plot Model  ####################
    
    
    stoc <- rnorm(1,0,post4$sigmay[sqi])
    #COL <- adjustcolor(c("black", "red", "green", "blue")[dat4$treatment], alpha.f = 0.5)
    colt<- c("black", "red", "green", "blue")
    
    
    
    calcpred <-function(matrixs,iters,new.data,stoc) {
      preds<- array(0,c(iters,nrow(new.data)))
      for(i in 1:iters){
        mui <- as.matrix(new.data)%*%c(matrixs[i,c(1:13)],0, stoc)
        preds[i,] <- mui
      }
      return(preds)
    } 
    
    
    labelmain <- c("Fecundity") 
    #plot(x,x, xlab="ln (Number of Branches)",
    #     ylim=c(0,8),main=labelmain,xlim=c(0, 8),type="n",ylab="ln (Flowering Branches)",cex.lab=1.2)
    
    v.trt2 <- trt[ttt[2],2]
    v.trt3 <- trt[ttt[2],3]
    v.trt4 <- trt[ttt[2],4]
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
                                  rep(0,length(x)),
                                  rep(0,length(x)),
                                  rep(0,length(x)),
                                  rep(1,length(x)),
                                  rep(1,length(x)))
    )
    
    mu <- 1/(1 + (1/exp( calcpred(matrix_of_post,iterations,new.data, stoc) ) ))
    mean.fec[,,sqi] <-mu
    
    
    
    
    ########################### END FECUNDITY ##########################################
    
    ########################### SURVIVAL ###############################################
    
    
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
    row.has.na <- apply(dat5, 1, function(x1){any(is.na(x1))})
    sum(row.has.na)
    dat5$py <- 0
    k <-0
    
    sites <- table(dat5$sites,dat5$year)
    dsites <- dim(sites) 
    
    for (i in 1:25){
      for(j in 1:7){
        k <- k+1
        #print(c(i,j))
        dat5$py[dat5$sites==j & dat5$yn==i]  <- k
      }
    }
    
    dat5$stg_C[dat5$stg==3] <- 0
    dat5$stg_C[dat5$stg==5] <- 1
    dat5 <- subset(dat5,!is.na(yn) & !is.na(py))
    row.has.na <- apply(dat5, 1, function(x1){any(is.na(x1))})
    sum(row.has.na)
    
    dat5$treatment <- as.factor(dat5$treatment)
    dat5$stg <- as.factor(dat5$stg)
    dat5$log_ini_branches <- log(dat5$ini_branches)
    dat5$logsquare <- dat5$log_ini_branches^2
    dat5 <- subset(dat5, !is.na(log_ini_branches))
    
    #############################################################################
    sum_model5 <- precis(model5,digits=3,depth=1)
    #############################################################################
    
    #pj.zeros <- matrix(0,1000,12)
    #yj.zeros <- matrix(0,1000,25)
    
    
    post5 <- extract.samples(model5)
    stoc <- rnorm(1,0,post5$sigmay[sqi])
    
    
    mat.stoch5 <- cbind(post5$a,post5$b11,post5$b12,post5$b22,post5$b23,post5$b24,post5$b32,post5$b33,
                        post5$b34,post5$b41,post5$b42,post5$b43, post5$b52,post5$b53,post5$b54,
                        post5$b6,post5$b7,post5$b8,post5$b9,post5$b10,post5$b111  )
    mat.stoch5 <- cbind(mat.stoch5[1:iterations,], rep(0,iterations),rep(stoc,iterations) )
    
    
    for (j in 1:length(stg)){
      #par(mfrow=c(2,2))
      labelmain <- c("Survival",labstg[j])
      v.trt2 <- trt[ttt[2],2]
      v.trt3 <- trt[ttt[2],3]
      v.trt4 <- trt[ttt[2],4]
      new.data <-  data.frame(cbind(rep(1,length(x)),
                                    log_ini_branches = x,
                                    logsquare = x^2,
                                    T2 = rep(v.trt2,length(x)),
                                    T3 = rep(v.trt3,length(x)),
                                    T4 = rep(v.trt4,length(x)),
                                    stg_C= stg[j],
                                    habitat = rep(0,length(x)),
                                    aug = rep(0,length(x)), 
                                    burn = rep(0,length(x)),
                                    preburn = rep(0,length(x)),
                                    postburn = rep(0,length(x)),
                                    py=1,
                                    yn=1))
      #change aug 
      
      mu <-array(0,c(iterations,cells[ki])) 
      new.data.t <- as.matrix(cbind(new.data[,1:6],new.data[,4:6]*new.data[,2],new.data[,7],new.data[,7]*new.data[,2:3], 
                                    new.data[,4:6]*new.data[,3],new.data[,8:12],new.data[,8]*new.data[,10] ,array(1,c(cells[ki],2))))  
      for(wi in 1:iterations) {mu[wi,]<-1/(1+1/exp(as.numeric(new.data.t %*% mat.stoch5[wi,])))}
      mean.sur[,,j,sqi] <- mu
      
      
    }
    
    
  } ##### sqi
} ###### end ki



#############################################################################
##########################   FUNCTIONS    ###################################
#############################################################################

###############################################  Function adult growth ######

gxya <-function(x1,y1,trtmnt,cofa,mean.g,indxi) {
  
  mu.mean <- mean.g[indxi]
  
  sigmax2 <- exp(cofa[1] + cofa[2]*x1 + cofa[3:5]%*%trtmnt[2:4])
  
  fac1<-sqrt(2*pi*sigmax2)
  fac2<-((y1-mu.mean)^2)/(2*sigmax2)
  mux <- exp(-(fac2))/fac1
  return(mux)
  
} 

###############################################  Function yearling growth ######

gxyy <-function(x1,y1,trtmnt,cofy,mean.g,indxi) {
  
  mu.mean <- mean.g[indxi]
  
  sigmax2 <- exp(cofy[1] + cofy[2]*x1 + cofy[3:5]%*%trtmnt[2:4])
  fac1<-sqrt(2*pi*sigmax2)
  fac2<-((y1-mu.mean)^2)/(2*sigmax2)
  muxy <- exp(-(fac2))/fac1
  return(muxy)
} 


#########################################  Function P(reproduction)  ###############


px <- function(x1,stage, mean.repr,indxi) {
  
  mu.mean <- mean.repr[indxi,stage+1]
  return(mu.mean)
  
} 

#px(x[50],stg[1],mean.rep[1,,,1],50)

################################################  Function Fecundity ###############

fxp<-function(x1,meanfecu,indxi) {
  
  mu.mean <- meanfecu[indxi]
  return(mu.mean)
  
} 

# fxp(x[20],mean.fec[1,,1],20)
###############################################   Function P(Survival) #############

sx<-function(x1,stage,meansurv,indxi) {
  
  mu.mean <- meansurv[indxi,stage+1]
  return(0.95*mu.mean)
} 

# sx(x[ni],stg[1],mean.sur[1,,,1],50) 

########################    Function that combines survival and growth #############

pxy <- function(x1,y1,stage,trtmnt,cofa,mean.g,meansurv,indxi)
{return(sx(x1,stage,meansurv,indxi)*gxya(x1,y1,trtmnt,cofa,mean.g,indxi)) }


###########################################  Function Fruits per branch #############

fts<-function(tt,ni) {
  post <- extract.samples(modelf)
  f_aug <- exp(post$aug)
  f_intro <- exp(post$aug+post$B1)
  f_wild <- exp(post$aug+post$B2)
  f_road <- exp(post$aug+post$B3)
  if (tt==1){mean.fruits <- f_wild[ni] }
  if (tt==2){mean.fruits <- f_road[ni] }
  if (tt==3){mean.fruits <- f_aug[ni] }
  if (tt==4){mean.fruits <- f_intro[ni] }
  return(mean.fruits)
}
#fts(matsite1[yri,6],200) 
######################################## Function Seeds  per fruit #########

post <- extract.samples(model7)
semillas <- exp(post$seeds)


##############################   END vital ratesFunctions  ###################

##############################################################################  

####################  Defining matrices and Kernel functions   ###############

##############################################################################

######################################## Plant size range  ###################

minsize <- 0
maxsize <- 6.8   

#################### THE KERNEL K(y,x) for survival and growth  ##############

Kyx <-function(x1,y1,stage,trtmnt,cofa,mean.g,meansurv,indxi,indx) {
  xeval<-max(x1,minsize)
  xeval<-min(xeval,maxsize)
  yeval<-max(y1,minsize)
  yeval<-min(yeval,maxsize)
  return(pxy(xeval,yeval,stage,trtmnt,cofa,mean.g,meansurv,indxi))
}

bigmatrix<-function(ns,trtmnt,stage,cofa,mean.g,meansurv) {
  # upper and lower integration limits
  L<-minsize; U<-1*maxsize
  
  # boundary points b and mesh points y
  b<-L+c(0:ns)*(U-L)/ns
  y<-0.5*(b[1:ns]+b[2:(ns+1)])
  h <- y[2]-y[1]
  
  # loop to construct the matrix
  M<-matrix(0,ns,ns)
  for (i in 1:ns){
    for(j in 1:ns){
      M[i,j]<-Kyx(y[i],y[j],stage,trtmnt,cofa,mean.g,meansurv,i)
    }
  }
  M<-M*h
  return(list(matrix=M,meshpts=y))
}

          
          ############################################################################
          #### Introduction Burned ###################################################
          ############################################################################
          
          lambda.intb <- array(0,c(iterations,yrint))
          mat.year.intb <- array(0,c(96,96,iterations,yrint))
          
          for (ni in 1:iterations){
            for (yri in 1:yrint) {
              
              ############################################# Plant contribution ############
              
              datplants <- read.csv("yearling_recruitment.csv", header=T)
              ysec  <- 1994:2018
              y0 <- order(ysec)
              
              plant_contribution <- trunc(rnorm(1,mean(datplants$adults[datplants$treatment==4]),sqrt(var(datplants$adults[datplants$treatment==4]))))
              if(plant_contribution < 1) {plant_contribution  <- mean(datplants$adults[datplants$treatment==4])}
              plant_contribution_f <- 1/plant_contribution
              plant_contribution
              
              
              ################################################  Seeds per plant ############
              
              seeds_by_plant <-function(ns,stage,meanrep,meanfec,ni) {
                
                # upper and lower integration limits
                L<-minsize; U<-1*maxsize
                
                # boundary points b and mesh points y
                b <-L+c(0:ns)*(U-L)/ns
                y <-0.5*(b[1:ns]+b[2:(ns+1)])
                h <- y[2]-y[1]
                Spp <- array(0,c(ns))  
                for (z in 1:ns){
                  Spp[z] <- px(y[z],stage,meanrep,z) *
                    fxp(y[z],meanfec,z)*exp(y[z])*plant_contribution_f*
                    fts(ttt[7],ni)*semillas[ni]}
                return(Spp)
              }
              
              
              Sppi <- seeds_by_plant(ns,stg[1],mean.rep.intb[ni,,,yri],mean.fec.intb[ni,,yri],ni)
              
              #################################   Function Recruitment and Dormancy ##########
              posts <- extract.samples(model9)
              
              Ger <- posts$G[ni]
              Dor <- posts$D[ni]
              bc <- rep(0,4)
              a <- posts$a[ni]
              bc[1] <- 0
              bc[2] <- posts$b1[ni]
              bc[3] <- posts$b2[ni]
              bc[4] <- posts$b3[ni]
              Rec_seeds <- array(0,dim(Sppi))
              
              Recruits <- a+bc[ttt[7]]+log(Sppi*Ger) 
              
              
              Rec_seeds <- Dor*Sppi-exp(Recruits)
              
              ##################################  Yearling size distribution (proportion) 
              #The mean and variance of the gamma are E(X) = a*s and Var(X) = a*s^2.
              cfa <- rep(0.0001)
              Soffs  <- OffspringDistrib <- array(0,c(ns))
              
              par(mfrow=c(3,2))
              
              bks <-c(x,7)
              
              offmean <- mean(dat3y$log_ini_branch[dat3y$treatment==ttt[7]])
              offvar <- var(dat3y$log_ini_branch[dat3y$treatment==ttt[7]])
              hoff <- hist(dat3y$log_ini_branch[dat3y$treatment==ttt[7]],breaks=bks,plot=FALSE)#}
              
              OffspringDistrib <- dgamma(x+0.00001,(offmean^2)/offvar, offvar/offmean)
              Soffs <- OffspringDistrib/sum(OffspringDistrib)
              
              ############################################################################## 
              
              #################   Goldman matrix components ##############################
              
              # sum(Soffs[,4])
              
              MATYEARL <- array(0,c(ns,ns))
              
              FROMSBNK <- c(Soffs[1:14],sum(Soffs[15:80]),rep(0,ns))
              SBNK<-c(Dor,rep(0,15),Rec_seeds)
              length(SBNK)
              for (i in 1:ns){
                MATYEARL[i,] <- exp(Recruits[i]) * Soffs
              }
              KADULTS <- bigmatrix(ns,trt[ttt[7],],stg[1],cgra,mean.adlts.intb[ni,,yri],mean.sur.intb[ni,,,yri])
              
              KYEARLINGS <- bigmatrix(ns,trt[4,],stg[2],cgry,mean.yrl.intb[ni,,yri],mean.sur.intb[ni,,,yri])
              dim(KADULTS$matrix)
              YEARNULL <- matrix(rep(min(KADULTS$matrix),ns*ns),ns,ns)
              dim( YEARNULL )
              
              ###################################
              ####   OVERALL MATRIX   ###########
              ###################################
              
              ##############################################################################
              overalladults_s <-rbind(t(MATYEARL[1:80,1:15]),t(KADULTS$matrix ))
              # dim(overalladults)
              overallyearlings_s <-rbind(YEARNULL[1:15,1:15],t(KYEARLINGS$matrix[1:15,1:80]))
              # dim(overallyearlings)
              overall_no_sbnk_s <- cbind(FROMSBNK[1:95],overallyearlings_s,overalladults_s)
              # dim(overall_no_sbnk)
              overall <- rbind(SBNK[1:96], overall_no_sbnk_s )
              # dim(overall)
              ####################################################################################
              
              ########### Hc_kyx ###########################
              ##############################################
              mat.year.intb[,,ni,yri] <-overall
              lambda.intb[ni,yri] <- eigen(overall)$values[1]
              print(as.numeric(eigen(overall)$values[1]) )
              print(c(ni,yri))
            } ###### end yri
          } ####### end ni
          
          
          #####################################################
          
          ############################################################################
          #### Introduction Not Burned ###################################################
          ############################################################################
          
          lambda.intnb <- array(0,c(iterations,yrint))
          mat.year.intnb <- array(0,c(96,96,iterations,yrint))
          
          for (ni in 1:iterations){
            for (yri in 1:yrint) {
              
              ############################################# Plant contribution ############
              
              datplants <- read.csv("yearling_recruitment.csv", header=T)
              ysec  <- 1994:2018
              y0 <- order(ysec)
              
              plant_contribution <- trunc(rnorm(1,mean(datplants$adults[datplants$treatment==4]),sqrt(var(datplants$adults[datplants$treatment==4]))))
              if(plant_contribution < 1) {plant_contribution  <- mean(datplants$adults[datplants$treatment==4])}
              plant_contribution_f <- 1/plant_contribution
              plant_contribution
              
              
              ################################################  Seeds per plant ############
              
              seeds_by_plant <-function(ns,stage,meanrep,meanfec,ni) {
                
                # upper and lower integration limits
                L<-minsize; U<-1*maxsize
                
                # boundary points b and mesh points y
                b <-L+c(0:ns)*(U-L)/ns
                y <-0.5*(b[1:ns]+b[2:(ns+1)])
                h <- y[2]-y[1]
                Spp <- array(0,c(ns))  
                for (z in 1:ns){
                  Spp[z] <- px(y[z],stage,meanrep,z) *
                    fxp(y[z],meanfec,z)*exp(y[z])*plant_contribution_f*
                    fts(ttt[7],ni)*semillas[ni]}
                return(Spp)
              }
              
              
              Sppi <- seeds_by_plant(ns,stg[1],mean.rep.intnb[ni,,,yri],mean.fec.intnb[ni,,yri],ni)
              
              #################################   Function Recruitment and Dormancy ##########
              posts <- extract.samples(model9)
              
              Ger <- posts$G[ni]
              Dor <- posts$D[ni]
              bc <- rep(0,4)
              a <- posts$a[ni]
              bc[1] <- 0
              bc[2] <- posts$b1[ni]
              bc[3] <- posts$b2[ni]
              bc[4] <- posts$b3[ni]
              Rec_seeds <- array(0,dim(Sppi))
              
              Recruits <- a+bc[ttt[7]]+log(Sppi*Ger) 
              
              
              Rec_seeds <- Dor*Sppi-exp(Recruits)
              
              ##################################  Yearling size distribution (proportion) 
              #The mean and variance of the gamma are E(X) = a*s and Var(X) = a*s^2.
              cfa <- rep(0.0001)
              Soffs  <- OffspringDistrib <- array(0,c(ns))
              
              par(mfrow=c(3,2))
              
              bks <-c(x,7)
              
              offmean <- mean(dat3y$log_ini_branch[dat3y$treatment==ttt[7]])
              offvar <- var(dat3y$log_ini_branch[dat3y$treatment==ttt[7]])
              hoff <- hist(dat3y$log_ini_branch[dat3y$treatment==ttt[7]],breaks=bks,plot=FALSE)#}
              
              OffspringDistrib <- dgamma(x+0.00001,(offmean^2)/offvar, offvar/offmean)
              Soffs <- OffspringDistrib/sum(OffspringDistrib)
              
              ############################################################################## 
              
              #################   Goldman matrix components ##############################
              
              # sum(Soffs[,4])
              
              MATYEARL <- array(0,c(ns,ns))
              
              FROMSBNK <- c(Soffs[1:14],sum(Soffs[15:80]),rep(0,ns))
              SBNK<-c(Dor,rep(0,15),Rec_seeds)
              length(SBNK)
              for (i in 1:ns){
                MATYEARL[i,] <- exp(Recruits[i]) * Soffs
              }
              KADULTS <- bigmatrix(ns,trt[ttt[7],],stg[1],cgra,mean.adlts.intnb[ni,,yri],mean.sur.intnb[ni,,,yri])
              
              KYEARLINGS <- bigmatrix(ns,trt[4,],stg[2],cgry,mean.yrl.intnb[ni,,yri],mean.sur.intnb[ni,,,yri])
              dim(KADULTS$matrix)
              YEARNULL <- matrix(rep(min(KADULTS$matrix),ns*ns),ns,ns)
              dim( YEARNULL )
              
              ###################################
              ####   OVERALL MATRIX   ###########
              ###################################
              
              ##############################################################################
              overalladults_s <-rbind(t(MATYEARL[1:80,1:15]),t(KADULTS$matrix ))
              # dim(overalladults)
              overallyearlings_s <-rbind(YEARNULL[1:15,1:15],t(KYEARLINGS$matrix[1:15,1:80]))
              # dim(overallyearlings)
              overall_no_sbnk_s <- cbind(FROMSBNK[1:95],overallyearlings_s,overalladults_s)
              # dim(overall_no_sbnk)
              overall <- rbind(SBNK[1:96], overall_no_sbnk_s )
              # dim(overall)
              ####################################################################################
              
              ########### Hc_kyx ###########################
              ##############################################
              mat.year.intnb[,,ni,yri] <-overall
              lambda.intnb[ni,yri] <- eigen(overall)$values[1]
              print(as.numeric(eigen(overall)$values[1]) )
              print(c(ni,yri))
            } ###### end yri
          } ####### end ni
          
          
          #####################################################
          
          
          ############################################################################
          #### PULSES INTO WILD ######################################################
          ############################################################################
          
          lambda <- array(0,c(iterations,years))
          mat.year <- array(0,c(96,96,iterations,years))
          
          for (ni in 1:iterations){
            for (yri in 1:years) {
              
              ############################################# Plant contribution ############
              
              #plot(dat3a$log_ini_branches,dat3a$log_flw_branches,pch=16,cex=0.2)
              datplants <- read.csv("yearling_recruitment.csv", header=T)
              ysec  <- 1994:2018
              y0 <- order(ysec)
              
              plant_contribution <- trunc(rnorm(1,mean(datplants$adults[datplants$treatment==2]),sqrt(var(datplants$adults[datplants$treatment==2]))))
              if(plant_contribution < 1) {plant_contribution  <- mean(datplants$adults[datplants$treatment==2])}
              plant_contribution_f <- 1/plant_contribution
              plant_contribution
              
              
              ################################################  Seeds per plant ############
              
              seeds_by_plant <-function(ns,stage,meanrep,meanfec,ni) {
                
                # upper and lower integration limits
                L<-minsize; U<-1*maxsize
                
                # boundary points b and mesh points y
                b <-L+c(0:ns)*(U-L)/ns
                y <-0.5*(b[1:ns]+b[2:(ns+1)])
                h <- y[2]-y[1]
                Spp <- array(0,c(ns))  
                for (z in 1:ns){
                  Spp[z] <- px(y[z],stage,meanrep,z) *
                    fxp(y[z],meanfec,z)*exp(y[z])*plant_contribution_f*
                    fts(ttt[1],ni)*semillas[ni]}
                return(Spp)
              }
              
              Sppi <- seeds_by_plant(ns,stg[1],mean.rep[ni,,,yri],mean.fec[ni,,yri],ni)
              
              #################################   Function Recruitment and Dormancy ##########
              posts <- extract.samples(model9)
              
              Ger <- posts$G[ni]
              Dor <- posts$D[ni]
              bc <- rep(0,4)
              a <- posts$a[ni]
              bc[1] <- 0
              bc[2] <- posts$b1[ni]
              bc[3] <- posts$b2[ni]
              bc[4] <- posts$b3[ni]
              Rec_seeds <- array(0,dim(Sppi))
              
              Recruits <- a+bc[ttt[1]]+log(Sppi*Ger) 
              
              
              Rec_seeds <- Dor*Sppi-exp(Recruits)
              
              ##################################  Yearling size distribution (proportion) 
              #The mean and variance of the gamma are E(X) = a*s and Var(X) = a*s^2.
              cfa <- rep(0.0001)
              Soffs  <- OffspringDistrib <- array(0,c(ns))
              
              par(mfrow=c(3,2))
              
              bks <-c(x,7)
              
              if(ttt[2] < 3){
                offmean <- mean(dat3y$log_ini_branch[dat3y$treatment==ttt[2]])
                offvar <- var(dat3y$log_ini_branch[dat3y$treatment==ttt[2]])
                hoff <- hist(dat3y$log_ini_branch[dat3y$treatment==ttt[2]],breaks=bks,plot=FALSE)}
              
              OffspringDistrib <- dgamma(x+0.00001,(offmean^2)/offvar, offvar/offmean)
              Soffs <- OffspringDistrib/sum(OffspringDistrib)
              
              ############################################################################## 
              
              #################   Goldman matrix components ##############################
              
              # sum(Soffs[,4])
              
              MATYEARL <- array(0,c(ns,ns))
              
              FROMSBNK <- c(Soffs[1:15],rep(0,ns))
              SBNK<-c(Dor,rep(0,15),Rec_seeds)
              length(SBNK)
              for (i in 1:ns){
                MATYEARL[i,] <- exp(Recruits[i]) * Soffs
              }
              KADULTS <- bigmatrix(ns,trt[ttt[2],],stg[1],cgra,mean.adults[ni,,yri],mean.sur[ni,,,yri])
              
              KYEARLINGS <- bigmatrix(ns,trt[2,],stg[2],cgry,mean.yearlings[ni,,yri],mean.sur[ni,,,yri])
              dim(KADULTS$matrix)
              YEARNULL <- matrix(rep(min(KADULTS$matrix),ns*ns),ns,ns)
              dim( YEARNULL )
              if(length(x)==301){
                tiff(filename = "Matrices.tiff", width = 700, height = 700,
                     pointsize =5, res=300)
                par(mar=c(3,3,1,1), mgp=c(1.5,0.5,0), tck=-0.02)
                par(mfrow=c(2,2))
                image(KYEARLINGS$matrix^(1/8))
                image(KADULTS$matrix^(1/8))
                image(YEARNULL)
                image(MATYEARL^(1/8))
                dev.off()
              }
              ###################################
              ####   OVERALL MATRIX   ###########
              ###################################
              
              ##############################################################################
              overalladults_s <-rbind(t(MATYEARL[1:80,1:15]),t(KADULTS$matrix ))
              # dim(overalladults)
              overallyearlings_s <-rbind(YEARNULL[1:15,1:15],t(KYEARLINGS$matrix[1:15,1:80]))
              # dim(overallyearlings)
              overall_no_sbnk_s <- cbind(FROMSBNK[1:95],overallyearlings_s,overalladults_s)
              # dim(overall_no_sbnk)
              overall <- rbind(SBNK[1:96], overall_no_sbnk_s )
              # dim(overall)
              
              ########### Hc_kyx ###########################
              ##############################################
              mat.year[,,ni,yri] <-overall
              lambda[ni,yri] <- eigen(overall)$values[1]
              print(as.numeric(eigen(overall)$values[1]) )
              print(c(ni,yri))
            } ###### end yri
          } ####### end ni

########################
          
### ROAD DYNAMICS:

  rtracker <- array(0,c(iterations,years))
  ## initial population vector
  ssd <- read.csv("SSD.csv", header=T)
  n0 <- round(ssd[,2],03)  * 300
  sum(n0[16:96])
        
  for (k in 1: iterations){
      #n0 <- rep(1/(96),length=96)
      n0 <- round(ssd[,2],03)  * 300
      for(g in 1:years){
              
        K.ipm <- mat.year[,,k,g]
        n0<- K.ipm %*% n0
        N <- sum(n0)
        rtracker[k,g]<-log(N)
        #n0<-n0/N
      }}
          
    burnin <- round(years*0.1) ## drop the first 10% of the time series
    #rtracker <- rtracker[,-c(1:burnin)]
    rtracker <- rtracker[,1:years]

    write.csv(rtracker,file="rtracker_RoadSurv95.csv")

#### ROAD INTRODUCTION DYNAMICS
    
## Introduction as performed
    ## Introduction as performed
    
    intrplantsa<-log(dat2$ini_branches[dat2$year==2012 & dat2$trt==4 & dat2$stg ==3])
    
    brnSampl<-sample(1:208,round(0.61*208))
    brninit<-intrplantsa[brnSampl]
    hist(brninit)
    
    nbrnSampl<-sample(1:208,round(0.39*208))
    nbrninit<-intrplantsa[nbrnSampl]
    hist(nbrninit)
    
    
    brks <-0+c(0:ns)*(6.8)/80
    
    #mtp.fctr<-5 ### how much bigger is the pulse compared to what was actually done
    
    yaib<-hist(brninit,breaks=brks,plot=F)
    brninitAdts<-yaib$counts
    
    yainb<-hist(nbrninit,breaks=brks,plot=F)
    nbrninitAdts<-yainb$counts
    
    #yy<-hist(augmplantsy,breaks=brks,plot=F)
    #augmyearlings<-yy$counts*mtp.fctr
    inisds<-0
    
    ### FIRST EIGTH YEARS AS BURNED INTRODUCTION DYNAMICS
    
    introBdyn<-array(0,c(yrint,96,iterations))
    
    rtrackerIB <- array(0,c(iterations,years))
    
    for (k in 1: iterations){
      n0 <- c(inisds,rep(0,15),brninitAdts)
      
      for(g in 1:yrint){
        
        K.ipm <- mat.year.intb[,,k,g]
        n0<- K.ipm %*% n0
        introBdyn[g,,k]<- n0
        N <- sum(n0)
        rtrackerIB[k,g]<-N
        
        
      }
    }
    
    ### FIRST EIGTH YEARS AS NOT BURNED INTRODUCTION DYNAMICS
    
    introNBdyn<-array(0,c(yrint,96,iterations))
    
    rtrackerINB <- array(0,c(iterations,years))
    
    for (k in 1: iterations){
      n0 <- c(inisds,rep(0,15),nbrninitAdts)
      
      for(g in 1:yrint){
        
        K.ipm <- mat.year.intnb[,,k,g]
        n0<- K.ipm %*% n0
        introNBdyn[g,,k]<- n0
        N <- sum(n0)
        rtrackerINB[k,g]<-N
        
        
      }
    }
    
    
    ##### COMBINING BURNED AND NOT BURNED INTRODUCTIONS
    
    rtrackerRI<-rtrackerIB+rtrackerINB
    rtrackerRI
    
    introdyn<-introBdyn+introNBdyn
    
    #### SWITCH TO GAP DYNAMICS
    
    
    for (k in 1: iterations){
      n0 <- introdyn[8,,k]
      
      for(g in 9:years){
        
        K.ipm <- mat.year[,,k,g-8]
        n0<- K.ipm %*% n0
        N <- sum(n0)
        rtrackerRI[k,g]<-N
        #n0<-n0/N
      }}
          
    rtrackerRI<-log(rtrackerRI)
    
    burnin <- round(years*0.1) ## drop the first 10% of the time series
    #rtracker <- rtracker[,-c(1:burnin)]
    rtrackerRI <- rtrackerRI[,1:years]


    write.csv(rtrackerRI,file="rtracker_IntroRoadSurvr95.csv")

         