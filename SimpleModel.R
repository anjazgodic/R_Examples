#clean the workspace
rm(list=ls())

libraries(“plyr”, “lattice”, “latticeExtra”, “grid”, “zoo”, “lubridate”, “mondate”, “reshape2”, “splines”, “rstan”, “survival”, “glmnet”, “gbm”, “gamlr”, “fitdistrplus”, “snowfall”, “snow”, “data.table”, “pROC”)

################################################################################################################
################################################## LOAD DATA ###################################################
################################################################################################################

load("/Volumes/Analytics/Projects/Data/lb2.RData")
data=lb2
data <- data[order(data$PId, data$CID),]

for(i in c(search.field("diag.", data), “x1”, “x2”)){
  data[,i][is.na(data[,i])] <- 0
}

colnames(quants) <- c("Bottom", "Top")

#################################################
#################### MODELING ###################
#################################################

ff = status ~ x3+x4+x5+x6+x7+x8+x9+
  x10+x11+x12+x13+x14+x15+x16+
  x17+x18+x19+x20+x21+x22+   
  x23+x24+x25+x26+x27+
  x28+x29+x30+x31+
  x32+x33+x34+year+x35+x36 + x37 

wd = model.matrix(ff,data)
y = (data$status) 

set.seed(22)

lasso1 = cv.glmnet(wd[,-1],y,nfolds=30,family='binomial')
plot(lasso1)
print(which(lasso1$lambda == lasso1$lambda.min))
print(lasso1$glmnet.fit$beta[,which(lasso1$lambda == lasso1$lambda.min)])
cc=(lasso1$glmnet.fit$beta[,which(lasso1$lambda == lasso1$lambda.min)])
prednames = names(cc)[as.numeric(abs(cc))>0]

