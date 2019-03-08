#Used for manuscript of greenhouse experiment#

###############################
##########Data#################
###############################
library(tidyverse)
both<-cbind(loc[,-c(6:7)],tri[,-1],env[,-1])
str(both)
both2<-sapply(both[,-c(1:2)],function(x){
  ntest <- shapiro.test(x)
  if (ntest$p.value<0.05){
    x <- log(x)
  } else{
    x <- x
  }
})%>%scale
both2<-as.data.frame(both2)
str(both2)

a<-as.matrix(both2[,c("Shoot.number","Shoot.diameter","Greatest.height","Average.height","SER","SLA","Leaf.DM","Stem.DM","Aboveground.RGR",
                      "Rhizome.length","Rhizome.diameter","SRL","Belowground.DM","Belowground.RGR")])
b<-cor(a)
abs(b)>0.6

a2<-as.matrix(both2[,c("Lat","Long","bio36","bio37","bio38","bio39","bio40")])
b2<-cor(a2)
b2


##############################
#########corrplot#############
##############################
#####select pearson or spearman
corPlot <- function(data, cor_method){
  # data is a matrix or dataframe
  # cor_method indicates which kind of correlation coefficient is to be used
  # if cor_method is missing, then "pearson" coefficient will be used
  require(corrplot)
  require(grid)
  if (missing(cor_method)) {
    cor_method <- "spearman"
    exact <- NULL
  } else {
    cor_method <- cor_method
    exact <- FALSE
  }
  # mat : is a matrix of data
  # ... : further arguments to pass to the native R cor.test function
  cor.mtest <- function(mat,...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j],method=cor_method,exact=exact)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  }
  corrplot(cor(data,method=cor_method),method = "number",type = "lower",
           sig.level = 0.05,
           insig = "blank",
           tl.col = "black",
           tl.srt = 45,
           tl.cex = 0.8,
           number.cex = 0.6,
           diag = T,
           p.mat=cor.mtest(data))
  
}

corPlot(trait)

library(corrplot)
trimat<-as.matrix(tri[,-1])
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
M<-cor(trimat,method="spearman")
res1 <- cor.mtest(trimat, conf.level = .95)
res1
corrplot(M, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = res1$p, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

##########################################
#####Structural equation modeling(SEM)####
##########################################
library(lavaan)
growth<-'growth=~Shoot.number+Shoot.diameter+Greatest.height+Accumulative.height+
Rhizome.length+Rhizome.diameter+Stem.DM
Shoot.number ~~ Stem.DM
Shoot.number ~~  Rhizome.length
Shoot.number ~~ Greatest.height
Rhizome.length ~~          Stem.DM
Greatest.height ~~   Rhizome.length
Greatest.height ~~ Stem.DM
Shoot.number ~~ Accumulative.height
Shoot.diameter ~~ Accumulative.height 
Greatest.height ~~ Accumulative.height
Accumulative.height ~~ Stem.DM
growth~Lat+bio36+bio37+bio38+bio39+bio40
'
Growth.ci<-cfa(growth, data=both2,meanstructure=TRUE,fixed.x=F)
summary(Growth.ci,fit.measures=TRUE,rsq=TRUE,standardized=TRUE)
mi<-modindices(Growth.ci)
print(mi[mi$mi>3.84,])

growth<-'growth=~Shoot.number+Shoot.diameter+Greatest.height+
Rhizome.length+Rhizome.diameter+Stem.DM+Stem.allocation

Shoot.number ~~ Stem.DM
Shoot.number ~~  Rhizome.length
Shoot.number ~~ Greatest.height
Rhizome.length ~~          Stem.DM
Greatest.height ~~   Rhizome.length
Greatest.height ~~ Stem.DM
Greatest.height ~~ Stem.allocation
Rhizome.diameter ~~ Stem.allocation

growth~Lat+Long+bio36+0*bio37+bio38+bio39+bio40
bio36~Lat+Long
bio38~Lat+Long
bio39~Lat+Long
bio40~Long+Lat

bio36~~bio39
bio38~~bio39
Lat~~Long
Lat~~bio37
Long~~bio37
bio38~~bio37
'
Growth.ci<-cfa(growth, data=both2,meanstructure=TRUE,fixed.x=F)
summary(Growth.ci,fit.measures=TRUE,rsq=TRUE,standardized=TRUE)
mi<-modindices(Growth.ci)
print(mi[mi$mi>3.84,])

library(semPlot)
pdf("SEM.pdf", width=18.5, height=10)
par(mfrow=c(1,2))
semPaths(Growth.ci,"std")
dev.off()

##############################################
########repeated measured mixed model#########
##############################################
if(!require(psych)){install.packages("psych")}
if(!require(nlme)){install.packages("nlme")}
if(!require(car)){install.packages("car")}
if(!require(multcompView)){install.packages("multcompView")}
if(!require(lsmeans)){install.packages("lsmeans")}
if(!require(rcompanion)){install.packages("rcompanion")}

library(psych)
library(nlme)
library(car)
library(multcompView)
library(emmeans)
library(ggplot2)
library(rcompanion)

model.b = lm(highest ~ group*time, 
             random = ~1|genotype,
             data=rep,na.action = na.exclude)
ACF(model.b)

model = lme(highest ~ group*time, 
            random = ~1|genotype,
            correlation = corAR1(form = ~ time | genotype,
                                 value = 0.15),
            data=rep,na.action = na.exclude,
            method="REML")
Anova(model)

model.fixed = gls(highest~ group*time,
                  data=rep,na.action = na.exclude,
                  method="REML")

anova(model,
      model.fixed)

model.null = lme(highest ~ 1,
                 random = ~1|genotype,
                 data =rep,na.action = na.exclude)

nagelkerke(model, 
           model.null)

marginal = emmeans(model, 
                   ~time:group)

cld(marginal,
    alpha   = 0.05, 
    Letters = letters,     ### Use lower-case letters for .group
    adjust  = "tukey")

x = residuals(model)
plotNormalHistogram(x)

################################################
###Generalized linear mixed model (GLMM)########
################################################
library(lme4)
library(nlme)
library(MASS)
library(MuMIn)
library(RLRsim)
library(MuMIn)

###Shoot number
hist(both$Shoot.number)
plot(abs(both$Lat),both$Shoot.number)
shapiro.test(both$Shoot.number)
mshoot.number<-lme(Shoot.number~ Lat*Long,random=~1|group,method = "ML",
                   data2)
summary(mshoot.number)
plot(mshoot.number)
anova(mshoot.number)

mshoot.number1<-lmer(Shoot.number~ Lat*Long+(1|group),both,REML=FALSE)
mshoot.number0<-lm(Shoot.number~Lat*Long,both)
exactLRT(m=mshoot.number,m0=mshoot.number0)
anova(mshoot.number,mshoot.number0)

r.squaredGLMM(mshoot.number)
# library(lmerTest)
# mshoot.number2<-lmer(Shoot.number~ Lat*Long+(1|group),data2)
# a<-step(mshoot.number2,reduce.fixed = FALSE,reduce.random = FALSE)
# plot(a)

##################################
#####LMM+AIcc selection###########
##################################
plot(abs(both$Lat),both$Shoot.number)
shapiro.test(both$Shoot.number)
mShoot.number<-lm(Shoot.number~ Lat+Long+bio36+bio37+bio38+bio39+bio40,data=both)
summary(mShoot.number)
anova(mShoot.number)
options(na.action="na.fail")
ms2 <- dredge(mShoot.number, trace = TRUE, rank = "AICc", REML = FALSE)
fmList <- get.models(ms2, 1:4)
summary(model.avg(fmList))
mShoot.number2<-lm(Shoot.number~ Lat+Long+bio36,data=both)
plot(mShoot.number2)
summary(mShoot.number2)
anova(mShoot.number2)


###############################
#######random forest###########
###############################
library("randomForest")
library("plyr") # for the "arrange" function
library("rfUtilities") # to test model significance
library("caret") # to get leave-one-out cross-validation accuracies and also contains the nearZeroVar function 
library("vegan")#data standard normalize
library("magrittr")

tri<-read.csv("traitsnew3.csv")
both<-cbind(loc[,-c(2,6:7)],tri[,-1],env[,-1])
str(both)

set.seed(2018)
index = sample(1:nrow(both),nrow(both)*0.7)
train = both[index,]
test = both[-index,]

#########################
#grestest height
pgrid = expand.grid(2:10,seq(200,1000,100))
colnames(pgrid) = c("mtry","n.trees")

R2 = rep(NA,nrow(pgrid))
for(i in 1:nrow(pgrid))
{
  rf = randomForest(train$Greatest.height ~ ., data = train[,c(34:ncol(both)-5)],mtry = pgrid[i,1],n.trees = pgrid[i,2])
  pred = predict(rf,test)
  R2[i] = cor(test$Greatest.height,pred)^2 #Not bad - R2 of 0.74 in the test set
}

plot(pgrid,cex = R2*5)
i = which.max(R2)

set.seed(2018)
rf2 = randomForest(both$Greatest.height ~ ., data = both[,c(34:ncol(both)-5)],mtry = pgrid[i,1],n.trees = pgrid[i,2],keep.inbag = TRUE)
cor(test$Greatest.height,predict(rf2,test))^2 #R2 = 0.87 - not a huge amount of improvement

RF.mor.sig <- rf.significance( x=rf2,  xdata=both[,c(34:ncol(both)-5)], nperm=1000 , ntree=501 )  
RF.mor.sig 
fit_control <-trainControl( method = "LOOCV" )    
RF.mor.loocv <- train(both[,c(34:ncol(both)-5)] ,y= both$Greatest.height,
                      method="rf", ntree=501 , tuneGrid=data.frame( mtry=5 ) , trControl=fit_control )
RF.mor.loocv
varImpPlot(rf2)

library("dplyr")

RF.mor.imp <- as.data.frame( rf2$importance )
RF.mor.imp$features <- rownames( RF.mor.imp )
RF.mor.imp.sorted <- arrange( RF.mor.imp  , desc(`IncNodePurity`)  )
str(RF.mor.imp.sorted)
RF.mor.imp.sorted$features

library(ggplot2)
library(ggthemes)
RF.mor.imp.sorted10<-RF.mor.imp.sorted[1:10,]
RF.mor.imp.sorted10$features<-factor(RF.mor.imp.sorted10$features,levels=RF.mor.imp.sorted10$features[order(RF.mor.imp.sorted10$`IncNodePurity`,decreasing=TRUE)])
RF.mor.imp.sorted10$features

pdf("heightrf.pdf", width=8.5, height=5)
p<-ggplot(RF.mor.imp.sorted10,aes(x=features,y= `IncNodePurity`,fill=features))+geom_bar(stat="identity")
p+ylab("% Increase in Mean Squared Error\n (Variable Importance)")+xlab("")+
  theme_few()+scale_fill_brewer(type='div', palette=2,guide=FALSE)+
  theme(axis.text=element_text(face="bold", size=rel(1.5)),
        axis.title=element_text(face="bold", size=rel(1.5)))
dev.off()

library(forestFloor)
ff = forestFloor(rf2,both)
color = fcol(ff,cols=1)
pdf("grestestheight3.pdf", width=8.5, height=5)
plot(ff,plot_seq =1:9,orderByImportance=TRUE,col=color,plot_GOF = TRUE)
dev.off()


# ######################################################
# ################RDA+Forward selection#################
# ######################################################
# ###library required packages
# library(xlsx) # data import
# library(dplyr) # data manipulation
# library(vegan) # ordination methods
# library(packfor) # forward selection
# library(fossil) # PCNM analysis
# library(magrittr)###for %>%
# library(car)###test vif
# library(rda)###rda analysis
# library(CCA)###cca analysis
# library(fields)
# library(stats)
# 
# rda0<-data.frame(tri,envnew)
# rda1<-scale(rda0)
# rda1<-as.data.frame(rda1)
# str(rda1)
# rda2<-rda(rda1[,1:20]~.,data=envnew)
# summary(rda2)
# vif.cca(rda2)
# (R2a.all<-RsquareAdj(rda2)$adj.r.squared)
# library(vegan)
# step.forward<-ordistep(rda(rda1[,1:20]~1,data=envnew),scope=formula(rda2),
# direction="forward",pstep=1000)
# RsquareAdj(rda(rda1[,1:20]~bio23+bio05+bio14+bio22,data=envnew))$adj.r.squared
# 
# summary(step.forward)
# vif.cca(step.forward)
# anova.cca(step.forward)
# permutest(step.forward)

