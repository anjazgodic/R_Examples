#READ DATA IN
s <- read.csv("schiz.csv")

#REMOVE UNNECESSARY VARIABLES
s <- s[, -c(1, 3:5)]
s$Schizo <- as.factor(s$Schizo+1)

#SEPARATE GROUPS
s0 <- s[s$Schizo == 1, ]
s1 <- s[s$Schizo == 2, ]
s2 <- s[s$Schizo == 3, ]

#DO Hotelling's T2 TO ENSURE THAT GROUPS MEANS ARE DIFFERENT AND IT IS WORTH DOING LDA
library(Hotelling)
hotelling.test(s0[,-1], s1[,-1])$pval ##the means are different
hotelling.test(s0[,-1], s2[,-1])$pval ##the means are different
hotelling.test(s1[,-1], s2[,-1])$pval ##the means are not different


#CHECK EQUAL COVARIANCE ASSUMPTION USING GRAPHICAL DISPLAY
pdf("Graphical_Display_EqualCov.pdf")
#pairwise scatterplots showing groups
library(psych)
pairs.panels(s[,-1], smooth = TRUE, scale = TRUE, density=TRUE,ellipses=F,digits = 2,method="pearson", pch = 21,lm=FALSE, cor=TRUE,jiggle=FALSE,factor=2,hist.col="cyan",show.points=TRUE,rug=TRUE, bg=c("green","red","blue")[s$Schizo], main="Pairwise Correlation, Histogram and Scatterplot Matrix")
dev.off()

pdf("Scatterplots_Colored_Groups.pdf")
#scatterplot of each variable, showing the three groups
par(mfrow = c(2, 5), mar=c(5,4,1,2)+0.1, oma = c(0, 0, 5, 0))
for (i in 2:length(s)) {
  plot(s0[, i], col = "Blue", pch = 20, xlab = "Patients", ylab = colnames(s0)[i])
  points(s1[, i], col = "Red", pch = 20)
  points(s2[, i], col = "Green", pch = 20)
}
title("Scatterplots of Clinical Scores Displaying Groups: \nDiagnostic 1 in Blue \nDiagnostic 2 in Red and \nDiagnostic 3 in Green", outer = TRUE)
dev.off()

#OTHER USEFUL PLOTS
library(ade4)
par(mar = c(0, 0, 0, 0))
pan1 <- function(x, y, ...) {
  xy <- cbind.data.frame(x, y)
  s.class(xy, as.factor(s$Schizo), include.ori = F, add.p = T, clab = 1.5,
          col = c("blue", "black", "red"), cpoi = 2, csta = 0.5)
}
pairs(s[, 2:3], panel = pan1) #Hs, D
pairs(s[, 5:6], panel = pan1) #Pd, Mf
pairs(s[, 8:9], panel = pan1) #Pt, Sc
pairs(s[, 10:11], panel = pan1) #Ma, Si
pairs(s[, 2:11], panel = pan1) #All,


#CHECK EQUAL COVARIANCE ASSUMPTION USING STATISTICAL TEST, USING SAS
#see SAS output with Box-M Test for equality of covariance matrices
libname aa "C:\Users\azgodic\Desktop";
data schiz;
set aa.schiz;
run;
title 'Test of Homogeneity of covariance matrices';
proc discrim data = schiz pool=test;
class Schizo;
var Hs D Hy Pd Mf Pa Pt Sc Ma Si;
run;


#CONDUCT LDA AND GET INFO FROM SUMMARY
library(MASS)
library(DiscriMiner)

fit <- lda(s[, -1], grouping = (s$Schizo))
pred <- predict(fit, s[, -1])
mylda = linDA(s[,-1], group=(s$Schizo))

#PLOT LDA
plot(mylda$functions[,1], mylda$functions[,2], xlim=c(-5,5), ylim=c(-5,5))

x <- -4:5
mycol <- as.numeric(s$Schizo) + 7
plot(fit, dimen=3, col = mycol, xlim = c(-4, 5), main = "Class Discrimination") #plot each observation in the space of the first 2 linear discriminant functions

plot(pred$x, ##LD1 e LD2 locations
     type="n", xlab="LD1", ylab="LD2", 
     main="Class Discrimination", xlim=c(-4, 5), ylim=c(-4,4)) # does not show points because type = "n"
text(pred$x,as.character(pred$class), ##use text to plot the LDA classification codes
     col=mycol) # color of true diagnosis, symbol of classification
abline(h=0)
abline(v=0) 

plot(fit, dimen = 1, type = "both")

log(pi1/pi2)-(mu1+mu2)%*%solve(S)%*%t(mu1-mu2)+(as.matrix(x)%*%solve(S)%*%t(mu1-mu2))[,1]
as.numeric(log(pi1/pi2)-(mu1+mu2)%*%solve(S)%*%t(mu1-mu2)+(as.matrix(x)%*%solve(S)%*%t(mu1-mu2))[,1]<0)+1

log(pi1/pi3)-(mu1+mu3)%*%solve(S)%*%t(mu1-mu3)+(as.matrix(x)%*%solve(S)%*%t(mu1-mu3))[,1]
as.numeric(log(pi1/pi3)-(mu1+mu3)%*%solve(S)%*%t(mu1-mu3)+(as.matrix(x)%*%solve(S)%*%t(mu1-mu3))[,1]<0)+1

log(pi2/pi3)-(mu2+mu3)%*%solve(S)%*%t(mu2-mu3)+(as.matrix(x)%*%solve(S)%*%t(mu2-mu3))[,1]
as.numeric(log(pi2/pi3)-(mu2+mu3)%*%solve(S)%*%t(mu2-mu3)+(as.matrix(x)%*%solve(S)%*%t(mu2-mu3))[,1]<0)+1

#CHECK ERROR RATE AND CLASSIFICATION COUNTS
accuracy <- table(s$Schizo, pred$class)
colnames(accuracy) <- c("Allocated to Group 1", "Allocated to Group 2", "Allocated to Group 3")
rownames(accuracy) <- c("Is Group 1", "Is Group 2", "Is Group 3")
write.csv(accuracy, "LDAcounts.csv")

#CONDUCT QDA AND GET INFO FROM SUMMARY
fit2 <- qda(s[, -1], grouping = s$Schizo)
myqda = quaDA(s[,-1], group=(s$Schizo))
summary(myqda)

#CHECK ERROR RATE AND CLASSIFICATION COUNTS
pred2 <- predict(fit2, s[, -1])
accuracy2 <- table(s$Schizo, pred2$class)
colnames(accuracy2) <- c("Allocated to Schizo=0", "Allocated to Schizo=1", "Allocated to Schizo=2")
rownames(accuracy2) <- c("Is Schizo=0", "Is Schizo=1", "Is Schizo=2")

#CHECKING WITH SAS
libname aa "C:\Users\azgodic\Desktop";
data schiz;
set aa.schiz;
run;
/*QDA FOR FUN*/
  title 'QDA JUST TO CHECK';
proc discrim data = schiz pool=no;
class Schizo;
var Hs D Hy Pd Mf Pa Pt Sc Ma Si;
run;

#EXPLORE DECISION TREES
library(rpart)

#CHECK GINI INDEX VS ENTROPY
fitgini<-rpart(Schizo~Hs+D+Hy+Pd+Mf+Pa+Pt+Sc+Ma+Si, data=schiz, method="class", parms=list(split='gini'), control=rpart.control(xval=0))
fitinfo<-rpart(Schizo~Hs+D+Hy+Pd+Mf+Pa+Pt+Sc+Ma+Si, data=schiz, method="class", parms=list(split='information'), control=rpart.control(xval=0))

#PLOT THE TREE
library(rpart.plot)
par(mfrow=c(1,2))
prp(fitgini); 
prp(fitinfo); 

par(mfrow=c(1,1))
printcp(fitgini)
printcp(fitinfo)

table(schiz$Schizo, predict(fitgini, type = "class"))
table(schiz$Schizo, predict(fitinfo, type = "class"))

#DO CROSS-VALIDATION
cvgini<-rpart(Schizo~Hs+D+Hy+Pd+Mf+Pa+Pt+Sc+Ma+Si, data=schiz, method="class", parms=list(split='gini'), control=rpart.control(xval=402))
cvinfo<-rpart(Schizo~Hs+D+Hy+Pd+Mf+Pa+Pt+Sc+Ma+Si, data=schiz, method="class", parms=list(split='information'), control=rpart.control(xval=402))

#CHECK PLOTS
par(mfrow=c(1,2))
prp(fitgini, split.col="blue", Margin=0.05, main="GINI", extra=104);
prp(fitinfo, split.col="darkgreen", main="Entropy", extra=104); 

#CHECK CPTABLE
cvgini$cptable
cvinfo$cptable

#PLOT CROSS-VALIDATION ERROR (XERROR) VS CP
par(mfrow=c(1,1), oma=rep(2,4))
plotcp(cvgini, cex.lab=1, upper="splits", col=2) 
title("Complexity Penalty and CV Relative Error: GINI", outer=T)
cp3 = which(cvgini$cptable[, 2] == 3)

#PRUNING TREES
prgini<-prune(cvgini, cvgini$cptable[cp3, 1])


#PLOT PRUNED TREES
table(schiz$Schizo, predict(prgini, type = "class"))




