### ---------------
###
### Create: Dongjie Chen
### Date: 2023-03-02 21:59:01
### Email: chen_dj@sjtu.edu.cn
### Pancreatic Disease Center, Ruijin Hospital, SHSMU, Shanghai, China. 
###
### ---------------

rm (list=ls()); gc()

###########################################################################################################
####=======================================1. loading data=============================================####
###########################################################################################################

library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(future)
availableCores()  
suppressWarnings(suppressMessages(future::plan("multiprocess",workers = 40)))
nbrOfWorkers() 

load('cohorts_dat.rds')  # load 11 cohorts (1 training/10 testing set) data
load('SE_prog_gene.rds') # load 38 prognostic HM-SE genes 


###########################################################################################################
####=====================================2. data processing============================================####
###########################################################################################################

# z-score normalization of the expression matrix was applied across all dataset z-score normalization
mm <- cohorts_dat

mm <- lapply(mm,function(x){
  x <- cbind(x[,1:3],x[,SE_prog_gene])
})

mm <- lapply(mm,function(x){
  x[,-c(1:3)] <- scale(x[,-c(1:3)])
  return(x)})

est_data <- mm$AU_Array  # ICGC-AU-Array dataset (training set)
val_data_list <- mm
pre_var <- colnames(est_data)[-c(1:3)]
est_dd <- est_data[,c('OS.time','OS',pre_var)]
val_dd_list <- lapply(val_data_list,function(x){x[,c('OS.time','OS',pre_var)]})

rf_nodesize <- 5
seed <- 1

funlists <- c(StepCox,Enet,CoxBoost,RSF,plsRcox,Survivalsvm,GBM,superpc)
names <- c("StepCox","Enet","CoxBoost","RSF","plsRcox","Survivalsvm","GBM","superpc")
dirnames <- c("both","backward","forward")
alphalist <- c(seq(0,1,0.1))

###########################################################################################################
####===============================3. ML combinations==================================================####
###########################################################################################################

# RSF
RSF <- function(est_dd,val_dd_list){
  fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd,
               ntree = 1000,nodesize = rf_nodesize,##该值建议多调整  
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T,
               seed = seed)
  RS=predict(fit,newdata = val_dd_list)$predicted
  rid <- names(fit$importance[fit$importance>0])
  est_dd2 <- est_dd[,c('OS.time','OS',rid)]
  val_dd_list2 <- val_dd_list[,c('OS.time','OS',rid)]
  return(list(fit,RS,est_dd2,val_dd_list2))
}

# Enet,Lasso (Enet alpha=1), Ridge (Enet alpha=0)
Enet <- function(est_dd,val_dd_list,alpha){
  fit = cv.glmnet(as.matrix(est_dd[,3:ncol(est_dd)]),
                  as.matrix(Surv(est_dd$OS.time,est_dd$OS)),
                  family = "cox",
                  alpha = alpha,
                  nfolds = 10)
  RS = as.numeric(predict(fit,type = 'link', newx = as.matrix(val_dd_list[,-c(1, 2)]),s = fit$lambda.min))
  rid <- rownames(coef(fit,s = "lambda.min"))[which(coef(fit,s = "lambda.min")!=0)]
  est_dd2 <- est_dd[,c('OS.time','OS',rid)]
  val_dd_list2 <- val_dd_list[,c('OS.time','OS',rid)]
  return(list(fit,RS,est_dd2,val_dd_list2))
}

# StepCox
StepCox <- function(est_dd,val_dd_list,direction){
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = direction)
  RS = predict(fit, type = 'risk', newdata = val_dd_list)
  rid <- names(coef(fit))
  est_dd2 <- est_dd[,c('OS.time','OS',rid)]
  val_dd_list2 <- val_dd_list[,c('OS.time','OS',rid)]
  return(list(fit,RS,est_dd2,val_dd_list2))
}

# CoxBoost
CoxBoost <- function(est_dd,val_dd_list){
  pen <- optimCoxBoostPenalty(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                              trace=TRUE,start.penalty=500,parallel = T)
  
  cv.res <- cv.CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                        maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
  
  fit <- CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                  stepno=cv.res$optimal.step,penalty=pen$penalty)
  RS=as.numeric(predict(fit,newdata=val_dd_list[,-c(1,2)], newtime=val_dd_list[,1], newstatus=val_dd_list[,2], type="lp"))
  rid<-names(coef(fit)[coef(fit)!=0])
  est_dd2 <- est_dd[,c('OS.time','OS',rid)]
  val_dd_list2 <- val_dd_list[,c('OS.time','OS',rid)]
  return(list(fit,RS,est_dd2,val_dd_list2))
}

# plsRcox
plsRcox <- function(est_dd,val_dd_list){
  cv.plsRcox.res=cv.plsRcox(list(x=est_dd[,3:ncol(est_dd)],time=est_dd$OS.time,status=est_dd$OS),nt=10,verbose = FALSE)
  fit <- plsRcox(est_dd[,3:ncol(est_dd)],time=est_dd$OS.time,event=est_dd$OS,nt=as.numeric(cv.plsRcox.res[5]))
  RS=as.numeric(predict(fit,type="lp",newdata=val_dd_list[,-c(1,2)]))
  return(list(fit,RS))
}

# superpc
superpc <- function(est_dd,val_dd_list){
  data <- list(x=t(est_dd[,-c(1,2)]),y=est_dd$OS.time,censoring.status=est_dd$OS,
               featurenames=colnames(est_dd)[-c(1,2)])
  fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
  cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                       n.fold = 10,
                       n.components=3,
                       min.features=5,
                       max.features=nrow(data$x),
                       compute.fullcv= TRUE,
                       compute.preval=TRUE)
  test <-
    list(
      x = t(val_dd_list[, -c(1, 2)]),
      y = val_dd_list$OS.time,
      censoring.status = val_dd_list$OS,
      featurenames = colnames(val_dd_list)[-c(1, 2)]
    )
  ff <-
    superpc.predict(fit,
                    data,
                    test,
                    threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, ])],
                    n.components = 1)
  
  RS <- as.numeric(ff$v.pred)
  return(list(fit,RS))
}

# Survivalsvm
Survivalsvm <- function(est_dd,val_dd_list){
  fit = survivalsvm(Surv(OS.time,OS)~., data= est_dd, gamma.mu = 1)
  RS=as.numeric(predict(fit, val_dd_list)$predicted)
  return(list(fit,RS))
}

###########################################################################################################
####===============================4. Model Training with SEMet classifier=============================####
###########################################################################################################

# Enet(alpha=0.1)+Survivalsvm (reached max mean C-index)

#Enet(alpha=0.1)
seed <- 39
set.seed(seed)
cv.fit <- cv.glmnet(as.matrix(est_dd[,3:ncol(est_dd)]),
                    as.matrix(Surv(est_dd$OS.time,est_dd$OS)),
                    family = "cox",
                    alpha = 0.1,
                    nfolds = 10)
rid <- rownames(coef(cv.fit,s = "lambda.min"))[which(coef(cv.fit,s = "lambda.min")!=0)]

#Survivalsvm
est_dd2 <- est_dd[,c('OS.time','OS',rid)]
fit <- survivalsvm(Surv(OS.time,OS)~., data= est_dd2, gamma.mu = 1)

###########################################################################################################
####===============================5. Calculation of SEMet score=======================================####
###########################################################################################################

SEMet_score <- list()
for (i in names(cohorts_dat)) {
  val_data <- cohorts_dat[[i]] %>% column_to_rownames('ID')
  val_dd <- val_data[,c('OS.time','OS',rid)]
  SEMet_score[[i]] <- cbind(val_dd[,1:2],SEMet_score=as.numeric(predict(fit, val_dd)$predicted))
}
