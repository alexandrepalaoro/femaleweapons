##################################################
###      SEXUAL DIMORPHISM ANALYSES (2/2)      ###
###         Male and female contests           ###
##################################################

####
# The analyses used in this paper will be split
# in two parts. The first contained all analyses
# related to sexual dimorphism in the allometry
# of male and female weapons. On the second code
# (this one), you will find the analyses regarding
# sexual dimorphism in contest variables. Cheers.
###

## Cleaning up any used space
rm(list=ls())
##

## Loading the csv containing the contest variables
fights<-read.csv("aegla-fights2.csv",h=T)

## Loading the packages needed
library(scales)
library(beeswarm)
library(betareg)
library(lmtest)

## The first test is regarding the difference in contest duration.
## I first tried to use an OLS regression, but the residuals were 
## all over the place.
manu<-fights[fights$species=="manuinflata",]
longi<-fights[fights$species=="longirostri",]


#---- FIGHT DURATION----

m1<-lm(duration~sex,data=longi)
plot(m1,which=1)


## The plot below show that the longirostri species has a large 
## variance.

plot(duration~sex,data=longi)

## Instead of using a log-transformation, I used a Gamma distribution.
## The gamma is useful because the variance increases together with 
## mean value. Thus, as you can see below, it fitted the data well.
m2<-glm(duration+1~sex,family=Gamma,data=longi)
plot(m2,which=1)
summary(m2)
anova(m2,test="Chisq")

m3<-glm(duration+1~sex,family=Gamma,data=manu)
plot(m3,which=1)
summary(m3)
anova(m3,test="Chisq")

#---- PROPORTION OF HIGHLY AGGRESSIVE BEHAVIORS----

## For the proportion of highly aggressive acts, I summed +0.001
## to all observations. Why? Beta regressions do not allow zeroes nor ones.
## This, I simply added a little nudge to be able to perform the analyses
## without having to rely on a binomial distribution.

manu$prop.aggr<-(manu$highl.aggr+0.001)/(manu$duration+1)
longi$prop.aggr<-(longi$highl.aggr+0.001)/(longi$duration+1)


plot(prop.aggr~sex,data=longi)

# Since I did not state what to use as precision parameter,
# The warning message states that it is using one. No biggie.

m4<-betareg(prop.aggr~sex,data=longi)
plot(m4,which=1)
summary(m4)

m5<-betareg(prop.aggr~sex,data=manu)
plot(m5,which=1)
summary(m5)

# Unfortunately, the beta regression does not have an "aov table"
# that shows the test for the factor. Thus, I built mine using
# lrtest(). The differences between the models with interactions
# (or without them) are the same thing you will find in the classic
# aov table.

m4.1<-betareg(prop.aggr~sex,data=longi)
m4.2<-betareg(prop.aggr~1,data=longi)
m5.2<-betareg(prop.aggr~1,data=manu)

lrtest(m4,m4.2)
lrtest(m5,m5.2)

# Overall, these tests indicate that the species factor 
# has more information than models without species. But, 
# it also says that females are fighting similarly to males 
# - models that have the sex factor are not different from
# models with no factors (m3.4).

#------ LATENCY-----

# Lastly, we tested whether the latency to start a contest 
# differed between species and sex.
# I tried to use the Gamma distribution again, but the model
# did not fit very well. Thus, I turned to good-and-old
# log-transformation.

m6.1<-glm(latency~sex,family=Gamma,data=longi)
plot(m6.1,which=1)
summary(m6.1)
anova(m6.1,test="Chisq")

m6.2<-lm(log(latency)~sex,data=longi)
plot(m6.2,which=1)
summary.aov(m6.2)

m7.1<-glm(latency~sex,family=Gamma,data=manu)
plot(m7.1,which=1)
summary(m7.1)
anova(m7.1,test="Chisq")

m7.2<-lm(log(latency)~sex,data=manu)
plot(m7.2,which=1)
summary(m7.2)
summary.aov(m7.2)

# Once again, no difference.

#---- PLOT----

# The last part is the plot. It is a three-part figure 
# that follows the logic of the analyses above.

pdf("FIGURE5.pdf",width=10,height=8)
#png(file="fight-dynamics.png", units='mm',width=300,height=200,res=600)
par(las=1,bty='l',mfrow=c(2,3))

beeswarm(duration+1~sex,data=longi,pch=21,cex=1.5,bg=alpha('grey',0.5),
         ylab="Fight duration (s)",xaxt='n',cex.axis=1.3,cex.lab=1.3,
         xlab='Sex',las=1,bty='l')
boxplot(duration+1~sex,data=longi,outline=F,yaxt='n',boxwex=0.3,
        lwd=2,add=T,col=c(alpha("grey",0.5),alpha("black",0.5)),xaxt='n')   
axis(side=1, at=c(1,2),cex.axis=1.3, 
     labels=c("Female","Male"))
legend('topleft',"A",bty='n',cex=1.2)

beeswarm(prop.aggr~sex,data=longi,pch=21,cex=1.5,bg=alpha('grey',0.5),
         ylab="Proportion of time spent in highly aggressive behaviors",xaxt='n',
         xlab='Sex',ylim=c(0,1),cex.axis=1.3,cex.lab=1.3,las=1,bty='l')
boxplot(prop.aggr~sex,data=longi,outline=F,ylim=c(0,1),yaxt='n',boxwex=0.3,
        lwd=2,add=T,col=c(alpha("grey",0.5),alpha("black",0.5)),xaxt='n')   
axis(side=1, at=c(1,2), 
     labels=c("Female","Male"),cex.axis=1.3)
legend('topleft',"B",bty='n',cex=1.2)

beeswarm(log(latency)~sex,data=longi,pch=21,cex=1.5,bg=alpha('grey',0.5),
         ylab="Latency to start contest (log)",xaxt='n',las=1,bty='l',
         xlab='Sex',las=1,bty='l',cex.axis=1.3,cex.lab=1.3)
boxplot(log(latency)~sex,data=longi,outline=F,yaxt='n',boxwex=0.3,
        lwd=2,add=T,col=c(alpha("grey",0.5),alpha("black",0.5)),xaxt='n')
axis(side=1, at=c(1,2),cex.axis=1.3,
     labels=c("Female","Male"))
legend('topright',"C",bty='n',cex=1.2)

####

beeswarm(duration+1~sex,data=manu,pch=21,cex=1.5,bg=alpha('grey',0.5),
         ylab="Fight duration (s)",xaxt='n',cex.axis=1.3,cex.lab=1.3,
         xlab='Sex',las=1,bty='l')
boxplot(duration+1~sex,data=manu,outline=F,yaxt='n',boxwex=0.3,
        lwd=2,add=T,col=c(alpha("grey",0.5),alpha("black",0.5)),xaxt='n')   
axis(side=1, at=c(1,2),cex.axis=1.3,
     labels=c("Female","Male"))
legend('topleft',"D",bty='n',cex=1.2)

beeswarm(prop.aggr~sex,data=manu,pch=21,cex=1.5,bg=alpha('grey',0.5),
         ylab="Proportion of time spent in highly aggressive behaviors",xaxt='n',
         xlab='Sex',ylim=c(0,1),cex.axis=1.3,cex.lab=1.3,las=1,bty='l')
boxplot(prop.aggr~sex,data=manu,outline=F,ylim=c(0,1),yaxt='n',boxwex=0.3,
        lwd=2,add=T,col=c(alpha("grey",0.5),alpha("black",0.5)),xaxt='n')   
axis(side=1, at=c(1,2),cex.axis=1.3, 
     labels=c("Female","Male"))
legend('topleft',"E",bty='n',cex=1.2)

beeswarm(log(latency)~sex,data=manu,pch=21,cex=1.5,bg=alpha('grey',0.5),
         ylab="Latency to start contest (log)",xaxt='n',
         xlab='Sex',las=1,bty='l',cex.axis=1.3,cex.lab=1.3,las=1,bty='l')
boxplot(log(latency)~sex,data=manu,outline=F,yaxt='n',boxwex=0.3,
        lwd=2,add=T,col=c(alpha("grey",0.5),alpha("black",0.5)),xaxt='n')
axis(side=1, at=c(1,2), 
     labels=c("Female","Male"),cex.axis=1.3)
legend('topright',"F",bty='n',cex=1.2)

dev.off()

#DONE :D
