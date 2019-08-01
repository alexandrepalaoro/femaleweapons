##################################################
###      SEXUAL DIMORPHISM ANALYSES (1/2)      ###
###    Allometry in male and female weapons    ###
##################################################

####
# The analyses used in this paper will be split
# in two parts. The first (this one) will contain
# all analyses related to sexual dimorphism in the
# allometry of male and female weapons. On the 
# second code, you will find the analyses regarding
# sexual dimorphism in contest variables. Cheers.
###

## Cleaning up any used space
rm(list=ls())
##

## Loading the csv containing the morphological variables
aegla<-read.csv("morpho-aegla1.csv",h=T)

head(aegla)

## Loading the packages we are going to need
library(scales)
library(nlme)
library(phia)
library(beeswarm)

## The following code is an extension for the phia package.
## It allows the package to work with gls() output as well.

model.matrix.gls <- function(object, ...) {
  model.matrix(terms(object), data = getData(object), ...)
}
model.frame.gls <- function(object, ...) {
  model.frame(formula(object), data = getData(object), ...)
}
terms.gls <- function(object, ...) {
  terms(model.frame(object), ...)
}


## Since we are going to perform our analyses only on the major claw (the left claw),
## we split the data in different objects. (Note. "esquerdo" = left, "direito" = right
## in Portuguese).

aegla.esq<-aegla[aegla$side=="esquerdo",]
aegla.dir<-aegla[aegla$side=="direito",]

## Separating both species

longi<-aegla[aegla$species=="longirostri",]
manu<-aegla[aegla$species=="manuinflata",]

#------ CENTERING-----
## Now, we are going to centralize the cephalothorax length for the sexes of both species
head(longi)

longi.fem<-longi[longi$sex=="female",]
fem.body.centr<-as.vector(scale(log(longi.fem$cc),scale=F))
longi.male<-longi[longi$sex=="male",]
male.body.centr<-as.vector(scale(log(longi.male$cc),scale=F))

longi$cc.centr<-c(fem.body.centr,male.body.centr)

manu.fem<-manu[manu$sex=="female",]
manu.fem.body.centr<-as.vector(scale(log(manu.fem$cc),scale=F))
manu.male<-manu[manu$sex=="male",]
manu.male.body.centr<-as.vector(scale(log(manu.male$cc),scale=F))

manu$cc.centr<-c(manu.fem.body.centr,manu.male.body.centr)

longi.esq<-longi[longi$side=="esquerdo",]
manu.esq<-manu[manu$side=="esquerdo",]

## As a reviewer noted, we have no biological reason to analyze two species within
## the same model (i.e., adding species as an explanatory variable). 
## Thus, we are going to analyze both species separately.

# First, just checking out the data without performing any tests.

par(mfrow=c(2,2))
plot(log(cp)~cc.centr,data=longi.esq,cex=1.5,las=1,bty='l',
     pch=21,bg=c(alpha("grey",0.75),alpha("black",0.75))[as.numeric(sex)],
     ylab="Propodus length (log)",
     xlab="Cephalothorax length (log-centered)")

plot(log(ap)~cc.centr,data=longi.esq,cex=1.5,las=1,bty='l',
     pch=21,bg=c(alpha("grey",0.75),alpha("black",0.75))[as.numeric(sex)],
     ylab="Propodus height (log)",
     xlab="Cephalothorax length (log-centered)")

plot(log(cp)~cc.centr,data=manu.esq,cex=1.5,las=1,bty='l',
     pch=21,bg=c(alpha("grey",0.75),alpha("black",0.75))[as.numeric(sex)],
     ylab="Propodus length (log)",
     xlab="Cephalothorax length (log-centered)")

plot(log(ap)~cc.centr,data=manu.esq,cex=1.5,las=1,bty='l',
     pch=21,bg=c(alpha("grey",0.75),alpha("black",0.75))[as.numeric(sex)],
     ylab="Propodus height (log)",
     xlab="Cephalothorax length (log-centered)")

par(mfrow=c(1,1))

#---- CLAW LENGTH ANALYSES----

# First, an OLS regression

m1<-lm(log(cp)~cc.centr*sex,data=longi.esq)
plot(m1,col=longi.esq$sex,which=1)
summary.aov(m1)
plot(hatvalues(m1),col=longi.esq$sex)
plot(cooks.distance(m1),col=longi.esq$species)
text(cooks.distance(m1),row.names(cooks.distance(m1)),pos=3)

m2<-lm(log(cp)~cc.centr*sex,data=manu.esq)
plot(m2,col=manu.esq$sex,which=1)
summary.aov(m2)
plot(hatvalues(m2),col=manu.esq$sex)
plot(cooks.distance(m2),col=manu.esq$sex)
text(cooks.distance(m2),row.names(cooks.distance(m2)),pos=3)

## For longirostri the residuals were a bit off. I will give a further look
## with a GLS. However, manuinflata seems fine and an OLS should work fine.

m1.g<-gls(log(cp)~cc.centr*sex,data=manu.esq,
          weights=varExp(form=~cc.centr))
plot(m1.g)
summary(m1.g)
confint(m1.g)

# It didn't change much regardless of the weight I used (i.e., varPower, varExp with cc, or Ident with sex).
# So, I will do the most parsimonious analyses and use the OLS.
# For more info, please check Pékar & Brabec (2016) in Ethology.

summary(m1)
summary.aov(m1) #Aegla longirostri model
summary(m2)
summary.aov(m2) #Aegla manuinflata model

coef.m1<-coef(m1)
ci.m1<-confint(m1)

coef.m2<-coef(m2)
ci.m2<-confint(m2)

## Since we saved the coefficients and CIs, we need to 
## calculate the values for each level of the factors.
## I wrote the following code to build a table containing
## all the values in that analyses. It is not scalable to 
## other scenarios, and it is highly redundant,
## but it does the job for this analyses.

estimates = data.frame(rep(NA, 8),
                       rep(NA, 8),
                       rep(NA, 8),
                       rep(NA, 8),
                       rep(NA, 8))
colnames(estimates) = c("lower.95","estimate","upper.95","measure","groups")
estimates$measure=sort(factor(rep(c("intercept","slope"),1)))
estimates$groups=factor(c("female.longi","female.longi",
                          "female.manu","female.manu","male.longi","male.longi",
                          "male.manu","male.manu"))
estimates$sex=factor(c("female","female",
                          "female","female","male","male",
                          "male","male"))
estimates$species=factor(c("longirostri","longirostri",
                          "manuinflata","manuinflata","longirostri","longirostri",
                          "manuinflata","manuinflata"))

#intercepts
estimates[1,2]=coef.m1[1]
estimates[5,2]=coef.m1[1]+coef.m1[3]
estimates[3,2]=coef.m2[1]
estimates[7,2]=coef.m2[1]+coef.m2[3]

#slopes
estimates[2,2]=coef.m1[2]
estimates[6,2]=coef.m1[2]+coef.m1[4]
estimates[4,2]=coef.m2[2]
estimates[8,2]=coef.m2[2]+coef.m2[4]

#CI - female.longi
estimates[1,1]=ci.m1[1,1]
estimates[1,3]=ci.m1[1,2]
estimates[2,1]=ci.m1[2,1]
estimates[2,3]=ci.m1[2,2]

#CI - female.manu
estimates[3,1]=ci.m2[1,1] #intercept
estimates[3,3]=ci.m2[1,2]
estimates[4,1]=ci.m2[2,1] #slope
estimates[4,3]=ci.m2[2,2]

#CI - male.longi
estimates[5,1]=ci.m1[1,1]+ci.m1[3,1] #intercept
estimates[5,3]=ci.m1[1,2]+ci.m1[3,2]
estimates[6,1]=ci.m1[2,1]+ci.m1[4,1] #slope 
estimates[6,3]=ci.m1[2,2]+ci.m1[4,2]

#CI - male.manu
estimates[7,1]=ci.m2[1,1]+ci.m2[3,1] #intercept
estimates[7,3]=ci.m2[1,2]+ci.m2[3,2]
estimates[8,1]=ci.m2[2,1]+ci.m2[4,1] #slope 
estimates[8,3]=ci.m2[2,2]+ci.m2[4,2]

inters=estimates[estimates$measure=="intercept",]
slopes=estimates[estimates$measure=="slope",]

inters<-inters[order(inters$species),]
slopes<-slopes[order(slopes$species),]

inter.longi<-inters[inters$species=='longirostri',]
inter.manu<-inters[inters$species=='manuinflata',]

slope.longi<-slopes[slopes$species=='longirostri',]
slope.manu<-slopes[slopes$species=='manuinflata',]

# This is what you will get in the end.
estimates
# As you can see, the uncertainty regarding the manuinflata
# slope is high. I have tried to work around this issue, but 
# I had no luck. Since the data for this species was retrieved 
# from the literature, there is not much we can do about it. 
# However, we also got high uncertainty even using data I collected.
# So, the slopes might indeed be more variable in males (which could
# indicate towards a condition-dependent trait).

#---- PLOT----
# This is the plot for the paper. I am going to plot a three part figure.
# First, intercepts and CIs. Second, slopes and CIs. Third, claw length 
# against body size and regression lines.

pdf("FIGURE2.pdf",width=12,height=8)
#png(file="allometry_claw-length.png", units='mm',width=300,height=220,res=600)
par(mfrow=c(2,3),las=1,bty='l')
plot(inter.longi$estimate,cex=4,las=1,bty='l',pch=19,col=c(alpha("grey",0.75),alpha("black",0.75))[inters$sex],
     ylab="Intercept",xlab="sex",ylim=c(2.3,2.9),xlim=c(0.5,2.5),xaxt='n',cex.axis=1.5,cex.lab=1.5)
axis(side=1, at=c(1,2), labels=c("Female","Male"),cex.axis=1.5)

segments(x0=1:2,y0=inter.longi$lower.95,y1=inter.longi$upper.95,
         col=c(alpha("grey",0.7),alpha("black",0.7))[inters$sex],
         lwd=5)

points(inter.longi$estimate,cex=4,pch=19,col=c(alpha("grey",0.95),alpha("black",0.95))[inter.longi$sex])
legend('topleft',"A",bty='n',cex=1.5)

plot(slope.longi$estimate,cex=4,xaxt='n',las=1,bty='l',pch=19,cex.axis=1.5,cex.lab=1.5,
     col=c(alpha("grey",0.85),alpha("black",0.85))[slope.longi$sex],
     ylab="Slope",xlab="Sex",ylim=c(0.3,2.2),xlim=c(0.5,2.5))
axis(side=1, at=c(1,2), labels=c("Female","Male"),cex.axis=1.5)
segments(x0=1:2,y0=slope.longi$lower.95,y1=slope.longi$upper.95,
         col=c(alpha("grey",0.7),alpha("black",0.7))[slope.longi$sex],lwd=5)
points(slope.longi$estimate,cex=4,pch=19,col=c(alpha("grey",0.9),
                                               alpha("black",0.9))[slope.longi$sex])
legend('topleft',"B",bty='n',cex=1.5)

plot(log(cp)~cc.centr,data=longi.esq,cex=3,las=1,bty='l',cex.axis=1.5,cex.lab=1.5,
     pch=21,bg=c(alpha("grey",0.75),alpha("black",0.75))[as.numeric(sex)],
     ylab="Claw length (log)",
     xlab="Cephalothorax length (log-centered)")
curve(inter.longi[1,2]+(slope.longi[1,2]*x),add=T,lwd=5,col=alpha('black',0.75),lty=2,from=-0.2,to=0.3)
curve(inter.longi[2,2]+(slope.longi[2,2]*x),add=T,lwd=5,col=alpha('grey',0.75),from=-0.5,to=0.3)
legend('topleft',"C",bty='n',cex=1.5)

####

plot(inter.manu$estimate,cex=4,las=1,bty='l',pch=19,col=c(alpha("grey",0.75),alpha("black",0.75))[inters$sex],
     ylab="Intercept",xlab="Sex",ylim=c(2,2.9),xlim=c(0.5,2.5),xaxt='n',cex.lab=1.5,cex.axis=1.5)
axis(side=1, at=c(1,2), labels=c("Female","Male"),cex.axis=1.5)

segments(x0=1:2,y0=inter.manu$lower.95,y1=inter.manu$upper.95,
         col=c(alpha("grey",0.7),alpha("black",0.7))[inters$sex],
         lwd=5)

points(inter.manu$estimate,cex=4,pch=19,col=c(alpha("grey",0.9),alpha("black",0.9))[inter.manu$sex])
legend('topleft',"D",bty='n',cex=1.5)

plot(slope.manu$estimate,cex=4,xaxt='n',las=1,bty='l',pch=19,cex.axis=1.5,cex.lab=1.5,
     col=c(alpha("grey",0.85),alpha("black",0.85))[slope.manu$sex],
     ylab="Slope",xlab="Sex",ylim=c(0.4,1.5),xlim=c(0.5,2.5))
axis(side=1, at=c(1,2), labels=c("Female","Male"),cex.axis=1.5)
segments(x0=1:2,y0=slope.manu$lower.95,y1=slope.manu$upper.95,
         col=c(alpha("grey",0.7),alpha("black",0.7))[slope.longi$sex],lwd=5)
points(slope.manu$estimate,cex=4,pch=19,col=c(alpha("grey",0.9),
                                               alpha("black",0.9))[slope.manu$sex])
legend('topleft',"E",bty='n',cex=1.5)

plot(log(cp)~cc.centr,data=manu.esq,cex=3,las=1,bty='l',cex.axis=1.5,cex.lab=1.5,
     pch=21,bg=c(alpha("grey",0.75),alpha("black",0.7))[as.numeric(sex)],
     ylab="Claw length (log)",
     xlab="Cephalothorax length (log-centered)")
curve(inter.manu[1,2]+(slope.manu[1,2]*x),add=T,lwd=5,col=alpha('black',0.75),lty=2,from=-0.2,to=0.3)
curve(inter.manu[2,2]+(slope.manu[2,2]*x),add=T,lwd=5,col=alpha('grey',0.75),from=-0.15,to=0.25)
legend('topleft',"F",bty='n',cex=1.5)

dev.off()


#---- CLAW HEIGHT ANALYSES----
# This is pretty much the same analyses you saw before, but now
# I am using claw height as the response variable, instead of 
# claw length. The procedures are exactly the same.

m3<-lm(log(ap)~cc.centr*sex,data=longi.esq)
plot(m3,which=1,col=longi.esq$sex)
summary(m3)
summary.aov(m3)

m3.g<-gls(log(ap)~cc.centr*sex,data=longi.esq,
          weights=varPower(form=~cc.centr))
plot(m3.g)
summary(m3.g)
anova(m3.g)
 
m4<-lm(log(ap)~cc.centr*sex,data=manu.esq)
plot(m4,which=1,col=manu.esq$sex)
summary(m4)
summary.aov(m4)

m4.g<-gls(log(ap)~cc.centr*sex,data=manu.esq,
          weights=varExp(form=~cc.centr))
plot(m4.g)


# Hopefully you noticed that using the body size as a weight in this analysis
# "fixed" the residuals better than in the previous analyses. Maybe claw height
# has a stronger correlation to body size than claw length? (Quite possibly, 
# but not our goal here). 


coef.m3<-coef(m3)
ci.m3<-confint(m3)
coef.m4<-coef(m4)
ci.m4<-confint(m4)

estimates2 = data.frame(rep(NA, 8),
                       rep(NA, 8),
                       rep(NA, 8),
                       rep(NA, 8),
                       rep(NA, 8))
colnames(estimates2) = c("lower.95","estimate","upper.95","measure","groups")
estimates2$measure=sort(factor(rep(c("intercept","slope"),1)))
estimates2$groups=factor(c("female.longi","female.longi",
                          "female.manu","female.manu","male.longi","male.longi",
                          "male.manu","male.manu"))
estimates2$sex=factor(c("female","female",
                       "female","female","male","male",
                       "male","male"))
estimates2$species=factor(c("longirostri","longirostri",
                           "manuinflata","manuinflata","longirostri","longirostri",
                           "manuinflata","manuinflata"))

#intercepts
estimates2[1,2]=coef.m3[1]
estimates2[5,2]=coef.m3[1]+coef.m3[3]
estimates2[3,2]=coef.m4[1]
estimates2[7,2]=coef.m4[1]+coef.m4[3]

#slopes
estimates2[2,2]=coef.m3[2]
estimates2[6,2]=coef.m3[2]+coef.m3[4]
estimates2[4,2]=coef.m4[2]
estimates2[8,2]=coef.m4[2]+coef.m4[4]

#CI - female.longi
estimates2[1,1]=ci.m3[1,1]
estimates2[1,3]=ci.m3[1,2]
estimates2[2,1]=ci.m3[2,1]
estimates2[2,3]=ci.m3[2,2]

#CI - female.manu
estimates2[3,1]=ci.m4[1,1] #intercept
estimates2[3,3]=ci.m4[1,2]
estimates2[4,1]=ci.m4[2,1] #slope
estimates2[4,3]=ci.m4[2,2]

#CI - male.longi
estimates2[5,1]=ci.m3[1,1]+ci.m3[3,1] #intercept
estimates2[5,3]=ci.m3[1,2]+ci.m3[3,2]
estimates2[6,1]=ci.m3[2,1]+ci.m3[4,1] #slope (coef.m1[1]+coef.m1[6])
estimates2[6,3]=ci.m3[2,2]+ci.m3[4,2]

#CI - male.manu
estimates2[7,1]=ci.m4[1,1]+ci.m4[3,1] #intercept
estimates2[7,3]=ci.m4[1,2]+ci.m4[3,2]
estimates2[8,1]=ci.m4[2,1]+ci.m4[4,1] #slope (coef.m1[1]+coef.m1[5]+coef.m1[6]+coef.m1[8])
estimates2[8,3]=ci.m4[2,2]+ci.m4[4,2]


inters2=estimates2[estimates2$measure=="intercept",]
slopes2=estimates2[estimates2$measure=="slope",]

inters2<-inters2[order(inters2$species),]
slopes2<-slopes2[order(slopes2$species),]

inter2.longi<-inters2[inters2$species=='longirostri',]
inter2.manu<-inters2[inters2$species=='manuinflata',]

slope2.longi<-slopes2[slopes2$species=='longirostri',]
slope2.manu<-slopes2[slopes2$species=='manuinflata',]


# Same plot as before

#---- PLOT 2----

pdf("FIGURE3.pdf",width=12,height=8)
#png(file="allometry_claw-height.png", units='mm',width=300,height=220,res=600)
par(mfrow=c(2,3),las=1,bty='l')

plot(inter2.longi$estimate,cex=4,las=1,bty='l',pch=19,col=c(alpha("grey",0.75),alpha("black",0.75))[inters2$sex],
     ylab="Intercept",xlab="sex",ylim=c(1.7,2.5),xlim=c(0.5,2.5),xaxt='n',cex.axis=1.5,cex.lab=1.5)
axis(side=1, at=c(1,2), labels=c("Female","Male"),cex.axis=1.5)

segments(x0=1:2,y0=inter2.longi$lower.95,y1=inter2.longi$upper.95,
         col=c(alpha("grey",0.7),alpha("black",0.7))[inters2$sex],
         lwd=5)

points(inter2.longi$estimate,cex=4,pch=19,col=c(alpha("grey",0.9),alpha("black",0.9))[inter2.longi$sex])
legend('topleft',"A",bty='n',cex=1.5)

plot(slope2.longi$estimate,cex=4,xaxt='n',las=1,bty='l',pch=19,cex.axis=1.5,cex.lab=1.5,
     col=c(alpha("grey",0.85),alpha("black",0.85))[slope2.longi$sex],
     ylab="Slope",xlab="Sex",ylim=c(0.3,2.6),xlim=c(0.5,2.5))
axis(side=1, at=c(1,2), labels=c("Female","Male"),cex.axis=1.5)
segments(x0=1:2,y0=slope2.longi$lower.95,y1=slope2.longi$upper.95,
         col=c(alpha("grey",0.7),alpha("black",0.7))[slope2.longi$sex],lwd=5)
points(slope2.longi$estimate,cex=4,pch=19,col=c(alpha("grey",0.9),
                                               alpha("black",0.9))[slope2.longi$sex])
legend('topleft',"B",bty='n',cex=1.5)

plot(log(ap)~cc.centr,data=longi.esq,cex=3,las=1,bty='l',cex.axis=1.5,cex.lab=1.5,
     pch=21,bg=c(alpha("grey",0.75),alpha("black",0.75))[as.numeric(sex)],
     ylab="Claw height (log)",
     xlab="Cephalothorax length (log-centered)")
curve(inter2.longi[1,2]+(slope2.longi[1,2]*x),add=T,lwd=5,col=alpha('black',0.75),lty=2,from=-0.2,to=0.3)
curve(inter2.longi[2,2]+(slope2.longi[2,2]*x),add=T,lwd=5,col=alpha('grey',0.75),from=-0.3,to=0.3)
legend('topleft',"C",bty='n',cex=1.5)

####

plot(inter2.manu$estimate,cex=4,las=1,bty='l',pch=19,col=c(alpha("grey",0.75),alpha("black",0.75))[inters2$sex],
     ylab="Intercept",xlab="Sex",ylim=c(1.4,2.5),xlim=c(0.5,2.5),xaxt='n',cex.lab=1.5,cex.axis=1.5)
axis(side=1, at=c(1,2), labels=c("Female","Male"),cex.axis=1.5)

segments(x0=1:2,y0=inter2.manu$lower.95,y1=inter2.manu$upper.95,
         col=c(alpha("grey",0.7),alpha("black",0.7))[inters2$sex],
         lwd=5)

points(inter2.manu$estimate,cex=4,pch=19,col=c(alpha("grey",0.9),alpha("black",0.9))[inter2.manu$sex])
legend('topleft',"D",bty='n',cex=1.5)

plot(slope2.manu$estimate,cex=4,xaxt='n',las=1,bty='l',pch=19,cex.axis=1.5,cex.lab=1.5,
     col=c(alpha("grey",0.85),alpha("black",0.85))[slope2.manu$sex],
     ylab="Slope",xlab="Sex",ylim=c(0.8,1.8),xlim=c(0.5,2.5))
axis(side=1, at=c(1,2), labels=c("Female","Male"),cex.axis=1.5)
segments(x0=1:2,y0=slope2.manu$lower.95,y1=slope2.manu$upper.95,
         col=c(alpha("grey",0.7),alpha("black",0.7))[slope2.longi$sex],lwd=5)
points(slope2.manu$estimate,cex=4,pch=19,col=c(alpha("grey",0.9),
                                              alpha("black",0.9))[slope2.manu$sex])
legend('topleft',"E",bty='n',cex=1.5)

plot(log(ap)~cc.centr,data=manu.esq,cex=3,las=1,bty='l',cex.axis=1.5,cex.lab=1.5,
     pch=21,bg=c(alpha("grey",0.75),alpha("black",0.75))[as.numeric(sex)],
     ylab="Claw height (log)",
     xlab="Cephalothorax length (log-centered)")
curve(inter2.manu[1,2]+(slope2.manu[1,2]*x),add=T,lwd=5,col=alpha('black',0.75),lty=2,from=-0.23,to=0.3)
curve(inter2.manu[2,2]+(slope2.manu[2,2]*x),add=T,lwd=5,col=alpha('grey',0.75),from=-0.15,to=0.25)
legend('topleft',"F",bty='n',cex=1.5)

dev.off()


##----HETEROCHELY ANALYSES----


###
# The following analyses is used specific to crusties.
# Crusties have two claws that can be adapted to different functions.
# Once claw might be used to tear and the other claw used to crush, for instance.
# Whenever this occurs, we called it heterochely. Heterochely
# might entail in differences in claw shape, strength, and even size.
# (Imagine a fiddler crab here). Aeglids are not that different.
# The males and females have a big left claw and a small right claw.
# The left one is used to crush/grasp stuff, and the small claw to tear.
# In a fight, individuals use mainly their left claws because they 
# grasp their opponents. Thus, a fair assumption would be that individuals
# might not necessarily invest in larger claws; individuals might invest
# less in one claw. Thus, the degree of heterochely might tells us something
# about how much a sex invests in crushing.
# The following analyses is about that. We built and index of heterochely
# (large claw minus smaller claw) and tested whether the sexes and species
# differed in their degree of heterochely.
###

manu.dir<-manu[manu$side=="direito",]
longi.dir<-longi[longi$side=="direito",]

# Building our "complex" index of heterochely

aegla.esq$dif<-((aegla.esq$cp-aegla.dir$cp))

longi.esq$dif<-((longi.esq$cp-longi.dir$cp))
manu.esq$dif<-((manu.esq$cp-manu.dir$cp))

# We used a two-way anova for this analyses.
# Despite the skewness of the residuals, we could
# not find a way to circumvent those issues.
# We tried several gls() models, but to no avail.
# If you have any idea of how to "fix" this issue,
# please let me know. Your feedback will be immensely
# appreciated!!

m5<-lm(dif~sex,data=longi.esq)
plot(m5,which=1)
summary(m5)
summary.aov(m5)

m6<-lm(dif~sex,data=manu.esq)
plot(m6,which=1)
summary(m6)
summary.aov(m6)

# Now, the cute plot.

#---- PLOT 3----

pdf("FIGURE4.pdf",width=10,height=5)
#png(file="magnitude-of-heterochely.png", units='mm',width=220,height=160,res=600)
par(mfrow=c(1,2),las=1,bty='l')
beeswarm(longi.esq$dif~longi.esq$sex,
         pch=21,cex=1.3,bg=alpha('grey',0.5),
         ylab="Difference between left and right claw",xaxt='n',
         xlab='Sex')
boxplot(longi.esq$dif~longi.esq$sex,outline=F,boxwex=0.3,
        lwd=2,add=T,col=c(alpha("grey",0.75),alpha("grey30",0.75)),xaxt='n')        
axis(side=1, at=c(1,2), 
     labels=c("Female","Male"))
legend("topleft","A",cex=1,bty='n')

beeswarm(manu.esq$dif~manu.esq$sex,
         pch=21,cex=1.3,bg=alpha('grey',0.5),
         ylab="Difference between left and right claw",xaxt='n',
         xlab='Sex')
boxplot(manu.esq$dif~manu.esq$sex,outline=F,boxwex=0.3,
        lwd=2,add=T,col=c(alpha("grey",0.75),alpha("grey30",0.75)),xaxt='n')        
axis(side=1, at=c(1,2), 
     labels=c("Female","Male"))
legend("topleft","B",cex=1,bty='n')
dev.off()


## Now I will plot histograms of the raw values of body size and claw size.
# Since we used centred data in our paper, it is hard to visualize the true
# extent of sexual dimorphism in some variables. Thus, this plot will illustrate
# the sex dimorphism in both species.

longi.male<-longi.esq[longi.esq$sex=="male",]
longi.fem<-longi.esq[longi.esq$sex=="female",]
manu.male<-manu.esq[manu.esq$sex=="male",]
manu.fem<-manu.esq[manu.esq$sex=="female",]

png("sex-dimorph_histograms.png",width=220,height=180,res=600,units='mm')
par(las=1,bty='l',mfrow=c(2,3))
plot(density(aegla.esq$cc),xlim=c(12,32),ylim=c(0,0.3),main="",
     xlab="Cephalothorax length (mm)",col=NA)
polygon(density(longi.male$cc),col=alpha("black",0.75),border=NA)
polygon(density(longi.fem$cc),col=alpha("grey",0.75),border=NA)
legend("topright",c("Female","Male"),cex=1.3,bty='n',
       pch=21,pt.bg= c(alpha("gray",0.75), alpha("black",0.75)))
mtext("A",4,las=1,padj=-8,adj=0)

plot(density(aegla.esq$cp),xlim=c(5,28),ylim=c(0,0.35),main="",
     xlab="Claw length (mm)",col=NA)
polygon(density(longi.male$cp),col=alpha("black",0.75),border=NA)
polygon(density(longi.fem$cp),col=alpha("grey",0.75),border=NA)
mtext("B",4,las=1,padj=-8,adj=0)

plot(density(aegla.esq$ap),xlim=c(2,18),ylim=c(0,0.5),main="",
     xlab="Claw height (mm)",col=NA)
polygon(density(longi.male$ap),col=alpha("black",0.75),border=NA)
polygon(density(longi.fem$ap),col=alpha("grey",0.75),border=NA)
mtext("C",4,las=1,padj=-8,adj=0)

###

plot(density(aegla.esq$cc),xlim=c(10,30),ylim=c(0,0.3),main="",
     xlab="Cephalothorax length (mm)",col=NA)
polygon(density(manu.male$cc),col=alpha("black",0.75),border=NA)
polygon(density(manu.fem$cc),col=alpha("grey",0.75),border=NA)
mtext("D",4,las=1,padj=-8,adj=0)

plot(density(aegla.esq$cp),xlim=c(0,28),ylim=c(0,0.35),main="",
     xlab="Claw length (mm)",col=NA)
polygon(density(manu.male$cp),col=alpha("black",0.75),border=NA)
polygon(density(manu.fem$cp),col=alpha("grey",0.75),border=NA)
mtext("E",4,las=1,padj=-8,adj=0)

plot(density(aegla.esq$ap),xlim=c(2,16),ylim=c(0,0.5),main="",
     xlab="Claw height (mm)",col=NA)
polygon(density(manu.male$ap),col=alpha("black",0.75),border=NA)
polygon(density(manu.fem$ap),col=alpha("grey",0.75),border=NA)
mtext("F",4,las=1,padj=-8,adj=0)
dev.off()

# DONE :D