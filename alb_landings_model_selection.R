library(mgcv)
library(corrplot)
library(robustHD)
library(tidyverse)
library(pscl)
library(broom)
library(ROCR)
library(MuMIn)
library(car)
library(pammtools)
library(mgcViz)

#develop model of albacore landings in season 3 for major ports on US West Coast

#read the updated data frames for all vessel lengths
Sdat = read.csv("/home/desiree/Documents/COCA/Albacore/Sdatgam18.csv")
Mdat = read.csv("/home/desiree/Documents/COCA/Albacore/Mdatgam18.csv")
Ldat = read.csv("/home/desiree/Documents/COCA/Albacore/Ldat18gam.csv")

#read the cost data and calcualte the distnace to cog
Sdatc = read.csv("/home/desiree/Documents/COCA/Albacore/Sdat18bipred.csv")
Mdatc = read.csv("/home/desiree/Documents/COCA/Albacore/Mdat18bipred.csv")
Ldatc = read.csv("/home/desiree/Documents/COCA/Albacore/Ldat18bipred.csv")

#compute distance to max costs - opportunity costs
Sdatc$dist = sqrt((Sdatc$cogya-Sdatc$lat)^2+(Sdatc$cogxa2-Sdatc$lon)^2)
Mdatc$dist = sqrt((Mdatc$cogya-Mdatc$lat)^2+(Mdatc$cogxa2-Mdatc$lon)^2)
Ldatc$dist = sqrt((Ldatc$cogya-Ldatc$lat)^2+(Ldatc$cogxa2-Ldatc$lon)^2)

#compute distance weighted by fuel price
Sdatc$distw = Sdatc$dist*Sdatc$fuelp
Mdatc$distw = Mdatc$dist*Mdatc$fuelp
Ldatc$distw = Ldatc$dist*Ldatc$fuelp

#compute fuel price adjusted by inflation
Sdatc$fuelpa = Sdatc$fuelp/(Sdatc$cpi/100)
Mdatc$fuelpa = Mdatc$fuelp/(Mdatc$cpi/100)
Ldatc$fuelpa = Ldatc$fuelp/(Ldatc$cpi/100)

#compute distance weighted by adjusted fuel price
Sdatc$distw = Sdatc$dist*Sdatc$fuelpa
Mdatc$distw = Mdatc$dist*Mdatc$fuelpa
Ldatc$distw = Ldatc$dist*Ldatc$fuelpa#code to compute if landings are 0 or not based on distance weighted by fuel costs

#Put all in one table
Sdat$distw=Sdatc$distw
Sdat$NVes[is.na(Sdat$NVes)]=0
Mdat$distw=Mdatc$distw
Mdat$NVes[is.na(Mdat$NVes)]=0
Ldat$distw=Ldatc$distw
Ldat$NVes[is.na(Ldat$NVes)]=0

#generate catch per effort
Sdat$cpue = Sdat$LBS/Sdat$NVes
Mdat$cpue = Mdat$LBS/Mdat$NVes
Ldat$cpue = Ldat$LBS/Ldat$NVes
#the 0/0 turn to NaN so change to 0
Sdat$cpue[is.na(Sdat$cpue)]=0
Mdat$cpue[is.na(Mdat$cpue)]=0
Ldat$cpue[is.na(Ldat$cpue)]=0

#use lognormal model, select positive catches
Sdat12=Sdat[Sdat$cpue>0,]
Mdat12=Mdat[Mdat$cpue>0,]
Ldat12=Ldat[Ldat$cpue>0,]

#construct a price deviation by port metric
Pfun <- function(x){
  pm=gam(PPP~s(LBS,k=3),data=x)
  pdev=pm$residuals
  return(pdev)
}
#apply the function by port for each vessel type
pde=by(Sdat12, Sdat12$PORT_NAME,Pfun)
pde.m=by(Mdat12, Mdat12$Port,Pfun)
#for southern CA ports, the number of landings by port are too small have to combine, create a port 2 var
Ldat12$Port2=c(as.character(Ldat12$Port[1:123]),rep("CAPort",6))
pde.l=by(Ldat12, Ldat12$Port2,Pfun)

#add price deviation to the data frame
Sdat12$PPD=c(pde[[1]],pde[[2]],pde[[3]],
             pde[[4]],pde[[5]],pde[[9]],
             pde[[11]],pde[[6]],pde[[7]],
             pde[[8]],pde[[10]])
Mdat12$PPD=c(pde.m[[1]],pde.m[[2]],pde.m[[3]],
             pde.m[[5]],pde.m[[6]],pde.m[[7]],
             pde.m[[8]],pde.m[[4]])
Ldat12$PPD=c(pde.l[[1]],pde.l[[2]],pde.l[[4]],
             pde.l[[5]],pde.l[[6]],pde.l[[7]],pde.l[[3]])


# use crossvalidation and backward stepwise selection to select final model for each vessel size

#################Small Vessels###########################################################################################
#fit full model on 20 random datasets that contian 75% of the data 
mlist.sc.f <- vector(mode = "list", length = 20)
mlist.sc2 <- vector(mode = "list", length = 20)
mlist.sc3 <- vector(mode = "list", length = 20)
mlist.sc4 <- vector(mode = "list", length = 20)
mlist.sc5 <- vector(mode = "list", length = 20)
mlist.sc6 <- vector(mode = "list", length = 20)
mlist.sc7 <- vector(mode = "list", length = 20)
mlist.sc8 <- vector(mode = "list", length = 20)
dev.mat=matrix(0,nrow = 8, ncol=20)
aic.mat=matrix(0,nrow = 8, ncol=20)

for (j in 1:20){
  tinds= sample(1:182, 136, replace = FALSE)
  train.s = Sdat12[tinds,]
  #build full model for cpue
  cfits = gam(cpue ~s(vb602,k=3)+s(biomass,k=3)+s(chinook,k=3)+s(dunge,k=3)+s(net,k=3)+s(tz,k=3)+s(PPD,k=3)+PORT_NAME,family = Gamma(link="log") ,data=train.s)
  cfit2 = gam(cpue ~s(vb602,k=3)+s(biomass,k=3)+s(chinook,k=3)+s(dunge,k=3)+s(tz,k=3)+s(PPD,k=3)+PORT_NAME,family = Gamma(link="log") ,data=train.s)
  cfit3=gam(cpue ~s(biomass,k=3)+s(chinook,k=3)+s(dunge,k=3)+s(tz,k=3)+s(PPD,k=3)+PORT_NAME,family = Gamma(link="log"), data=train.s)
  cfit4=gam(cpue ~s(biomass)+s(chinook,k=3)+s(dunge,k=3)+s(PPD,k=3)+PORT_NAME,family = Gamma(link="log"), data=train.s)
  cfit5=gam(cpue ~s(biomass)+s(chinook,k=3)+s(PPD,k=3)+PORT_NAME,family = Gamma(link="log"), data=train.s)
  cfit6=gam(cpue ~s(chinook,k=3)+s(PPD,k=3)+PORT_NAME,family = Gamma(link="log"), data=train.s)
  cfit7=gam(cpue ~s(chinook,k=3)+PORT_NAME,family = Gamma(link="log"), data=train.s)
  cfit8=gam(cpue ~PORT_NAME,family = Gamma(link="log"), data=train.s)

  #save model object
  mlist.sc.f[[j]]=cfits
  mlist.sc2[[j]]=cfit2
  mlist.sc3[[j]]=cfit3
  mlist.sc4[[j]]=cfit4
  mlist.sc5[[j]]=cfit5
  mlist.sc6[[j]]=cfit6
  mlist.sc7[[j]]=cfit7
  mlist.sc8[[j]]=cfit8
  
  #save the deviance explained
  dev.mat[1,j]=as.numeric(summary(mlist.sc.f[[j]])[14])
  dev.mat[2,j]=as.numeric(summary(mlist.sc2[[j]])[14])
  dev.mat[3,j]=as.numeric(summary(mlist.sc3[[j]])[14])
  dev.mat[4,j]=as.numeric(summary(mlist.sc4[[j]])[14])
  dev.mat[5,j]=as.numeric(summary(mlist.sc5[[j]])[14])
  dev.mat[6,j]=as.numeric(summary(mlist.sc6[[j]])[14])
  dev.mat[7,j]=as.numeric(summary(mlist.sc7[[j]])[14])
  dev.mat[8,j]=as.numeric(summary(mlist.sc8[[j]])[14])
  
  #save the aic
  aic.mat[1,j]=mlist.sc.f[[j]]$aic
  aic.mat[2,j]=mlist.sc2[[j]]$aic
  aic.mat[3,j]=mlist.sc3[[j]]$aic
  aic.mat[4,j]=mlist.sc4[[j]]$aic
  aic.mat[5,j]=mlist.sc5[[j]]$aic
  aic.mat[6,j]=mlist.sc6[[j]]$aic
  aic.mat[7,j]=mlist.sc7[[j]]$aic
  aic.mat[8,j]=mlist.sc8[[j]]$aic
}

#check which model has an AIC that is comparable (less than 2 AIC units) or lower to the full model
rowMeans(aic.mat) #model 2 has an AIC within 2 units of full model, this will be the final model

#fit final model to full dataset
m1.s=gam(cpue ~s(vb602,k=3)+s(biomass,k=3)+s(chinook,k=3)+s(dunge,k=3)+s(tz,k=3)+s(PPD,k=3)+PORT_NAME,family = Gamma(link="log"), data=Sdat12)
summary(m1.s)
AICc(m1.s)

#acfs and checks are fine, now we check for collinearity using VIFs
tmps= lm(cpue ~vb602+biomass+ chinook+dunge+tz+PPD+PORT_NAME+year,data=Sdat12)
vif(tmps)
#All are less than 10, so we are ok. 

#do a loop to plot ACf
Sdat12$resid=resid(m1.s)
Sports=unique(Sdat12$PORT_NAME)
for (j in 1:length(Sports)){
  tmp=Sdat12 %>% filter(PORT_NAME==Sports[j])
  dev.new(width=7, height=6)
  tiff(paste("acfs",Sports[j],".tiff", sep=""), width = 6, height = 5, units = 'in', res = 300, compression = 'lzw')
  acf(tmp$resid)
  dev.off()
}

dev.new(width=7, height=6)
tiff("Chekcs", width = 6, height = 5, units = 'in', res = 300, compression = 'lzw')
par(mfrow=c(2,2))
gam.check(m1.s)
dev.off()

par(mfrow=c(3,2))
plot(m1.s)

Sdat3=Sdat12
Sdat3$cpue = (fitted(m1.s,type="response"))
Sdat3$model="fit"
Sdat12$model="obs"
Sdat4=rbind(Sdat12,Sdat3)

dev.new(width=4, height=6)
tiff("fits", width = 3, height = 5, units = 'in', res = 300, compression = 'lzw')
ggplot(Sdat4, aes(x=year, y=cpue, colour=model,group=model)) + 
  geom_line() + geom_point()+
  ylab("Landings (lbs) per vessel")+xlab("") +
  ggtitle("") + theme_bw()+facet_wrap(.~PORT_NAME)
dev.off

#Calculate total landings given participation
Sdat4$LBS2=Sdat4$LBS
Sdat4$LBS2[183:364]=Sdat4$cpue[183:364]*Sdat4$NVes[183:364]

ggplot(Sdat4, aes(x=year, y=LBS2, colour=model,group=model)) + 
  geom_line() + geom_point()+
  ylab("Landings (lbs)")+xlab("") +
  ggtitle("") + theme_bw()+facet_wrap(.~PORT_NAME)

#################Medium Vessels###########################################################################################
#preliminary anlayses showed that for medium vessels, there remained some autocorrelation so allowed an autorrelated error strucutre and used a gamm model
mlist <- vector(mode = "list", length = 20)
mlist2 <- vector(mode = "list", length = 20)
mlist3 <- vector(mode = "list", length = 20)
mlist4 <- vector(mode = "list", length = 20)
mlist5 <- vector(mode = "list", length = 20)
mlist6<- vector(mode = "list", length = 20)
mlist7<- vector(mode = "list", length = 20)
mdev.mat=matrix(0,nrow = 7, ncol=20)
maic.mat=matrix(0,nrow = 7, ncol=20)

for (j in 1:20){
  tindm= sample(1:142, 107, replace = FALSE)
  train.m = Mdat12[tindm,]
  #build full model for cpue
  cfitm =gam(cpue ~s(vb1506,k=3)+s(biomass,k=3)+s(tz,k=3)+s(chinook,k=3)+s(dunge,k=3)+s(ground,k=3)+Port,data=train.m,family = Gamma(link = "log"))
  cfitm2 =gam(cpue ~s(vb1506,k=3)+s(biomass,k=3)+s(tz,k=3)+s(chinook,k=3)+s(dunge,k=3)+Port,data=train.m,family = Gamma(link = "log"))
  cfitm3 =gam(cpue ~s(vb1506,k=3)+s(biomass,k=3)+s(tz,k=3)+s(dunge,k=3)+Port,data=train.m,family = Gamma(link = "log"))
  cfitm4 =gam(cpue ~s(vb1506,k=3)+s(biomass,k=3)+s(tz,k=3)+Port,data=train.m,family = Gamma(link = "log"))
  cfitm5 = gam(cpue ~s(biomass,k=3)+s(tz,k=3)+Port,data=train.m,family = Gamma(link = "log"))
  cfitm6 = gam(cpue ~s(tz,k=3)+Port,data=train.m,family = Gamma(link = "log"))
  cfitm7 = gam(cpue ~Port,data=train.m,family = Gamma(link = "log"))
  #save model object
  mlist[[j]]=cfitm
  mlist2[[j]]=cfitm2
  mlist3[[j]]=cfitm3
  mlist4[[j]]=cfitm4
  mlist5[[j]]=cfitm5
  mlist6[[j]]=cfitm6
  mlist7[[j]]=cfitm7
  
  #save the deviance explained
  mdev.mat[1,j]=as.numeric(summary(mlist[[j]])[14])
  mdev.mat[2,j]=as.numeric(summary(mlist2[[j]])[14])
  mdev.mat[3,j]=as.numeric(summary(mlist3[[j]])[14])
  mdev.mat[4,j]=as.numeric(summary(mlist4[[j]])[14])
  mdev.mat[5,j]=as.numeric(summary(mlist5[[j]])[14])
  mdev.mat[6,j]=as.numeric(summary(mlist6[[j]])[14])
  mdev.mat[7,j]=as.numeric(summary(mlist7[[j]])[14])

  #save the aic
  maic.mat[1,j]=mlist[[j]]$aic
  maic.mat[2,j]=mlist2[[j]]$aic
  maic.mat[3,j]=mlist3[[j]]$aic
  maic.mat[4,j]=mlist4[[j]]$aic
  maic.mat[5,j]=mlist5[[j]]$aic
  maic.mat[6,j]=mlist6[[j]]$aic
  maic.mat[7,j]=mlist7[[j]]$aic

}

#check which model has an AIC that is comparable (less than 2 AIC units) or lower to the full model
rowMeans(maic.mat) #model 2 has an AIC within 2 units of full model, this will be the final model

#there was some autocorrelation in the residuals for some port, so need to add an autocorrelated error structure to the best model
m1.ma = gamm(cpue ~s(vb1506,k=3)+s(biomass,k=3)+s(tz,k=3)+s(dunge,k=3)+
               Port,data=Mdat12,family = Gamma(link = "log"),
             correlation = corAR1(form=~year|Port))
summary(m1.ma$gam)

par(mfrow=c(2,2))
gam.check(m1.ma$gam)

dev.new(width=9, height=4)
tiff("effects_m", width = 8, height = 3, units = 'in', res = 300, compression = 'lzw')
par(mfrow=c(2,2))
plot(m1.ma$gam)
dev.off()

Mdat3=Mdat12
Mdat3$cpue = fitted(m1.ma$gam,type="response")
Mdat3$model="fit"
Mdat12$model="obs"
Mdat4=rbind(Mdat12,Mdat3)

ggplot(Mdat4, aes(x=year, y=cpue, colour=model,group=model)) + 
  geom_line() + geom_point()+
  ylab("Landings (lbs) per vessel")+xlab("") +
  ggtitle("") + theme_bw()+facet_wrap(.~Port)

#Calculate total landings
Mdat4$LBS2=Mdat4$LBS
Mdat4$LBS2[143:284]=Mdat4$cpue[143:284]*Mdat4$NVes[143:284]

ggplot(Mdat4, aes(x=year, y=LBS2, colour=model,group=model)) + 
  geom_line() + geom_point()+
  ylab("Landings (lbs)")+xlab("") +
  ggtitle("") + theme_bw()+facet_wrap(.~Port)

#################Large Vessels###########################################################################################
llist <- vector(mode = "list", length = 20)
llist2 <- vector(mode = "list", length = 20)
llist3 <- vector(mode = "list", length = 20)
llist4 <- vector(mode = "list", length = 20)
ldev.mat=matrix(0,nrow = 4, ncol=20)
laic.mat=matrix(0,nrow = 4, ncol=20)

for (j in 1:20){
  tindl= sample(1:129, 97, replace = FALSE)
  train.l = Ldat12[tindl,]
  #build full model for cpue
  cfitm =gam(cpue ~s(biomass,k=3)+s(tz,k=3)+s(distw,k=3)+Port,data=train.l,family = Gamma(link = "log"))
  cfitm2 =gam(cpue ~s(biomass,k=3)+s(distw,k=3)+Port,data=train.l,family = Gamma(link = "log"))
  cfitm3 =gam(cpue ~s(biomass,k=3)+Port,data=train.l,family = Gamma(link = "log"))
  cfitm4 =gam(cpue ~Port,data=train.l,family = Gamma(link = "log"))
  #save model object
  llist[[j]]=cfitm
  llist2[[j]]=cfitm2
  llist3[[j]]=cfitm3
  llist4[[j]]=cfitm4

  #save the deviance explained
  ldev.mat[1,j]=as.numeric(summary(llist[[j]])[14])
  ldev.mat[2,j]=as.numeric(summary(llist2[[j]])[14])
  ldev.mat[3,j]=as.numeric(summary(llist3[[j]])[14])
  ldev.mat[4,j]=as.numeric(summary(llist4[[j]])[14])

  #save the aic
  laic.mat[1,j]=llist[[j]]$aic
  laic.mat[2,j]=llist2[[j]]$aic
  laic.mat[3,j]=llist3[[j]]$aic
  laic.mat[4,j]=llist4[[j]]$aic

}

#check which model has an AIC that is comparable (less than 2 AIC units) or lower to the full model
rowMeans(laic.mat) 

#final model
m1.la = gamm(cpue ~s(biomass,k=3)+s(distw,k=3)+Port,data=Ldat12,family = Gamma(link = "log"),
             correlation = corAR1(form=~year|Port))
summary(m1.la$gam)

dev.new(width=7, height=6)
tiff("Chekcl", width = 6, height = 5, units = 'in', res = 300, compression = 'lzw')
par(mfrow=c(2,2))
gam.check(m1.la$gam)
dev.off()

#check for collinearity using VIFs
tmp= lm(cpue ~biomass+ distw+ Port+year,data=Ldat12)
vif(tmp)
#All are less than 5, so we are ok.

Ldat3=Ldat12
Ldat3$cpue = fitted(m1.la$gam,type="response")
Ldat3$model="fit"
Ldat12$model="obs"
Ldat4=rbind(Ldat12,Ldat3)

ggplot(Ldat4, aes(x=year, y=cpue, colour=model,group=model)) + 
  geom_line() + geom_point()+
  ylab("Landings (lbs) per vessel")+xlab("") +
  ggtitle("") + theme_bw()+facet_wrap(.~Port)

#Calculate total landings
Ldat4$LBS2=Ldat4$LBS
Ldat4$LBS2[130:258]=Ldat4$cpue[130:258]*Ldat4$NVes[130:258]

ggplot(Ldat4, aes(x=year, y=LBS2, colour=model,group=model)) + 
  geom_line() + geom_point()+
  ylab("Landings (lbs)")+xlab("") +
  ggtitle("") + theme_bw()+facet_wrap(.~Port)

#check fit for total landings
Ldat5 = as.data.frame(Ldat4 %>% group_by(year,model) %>% summarize(LBS=mean(LBS2), cpue=mean(cpue)))

dev.new(width=5, height=5)
png("l_hist_l.png", width = 5, height = 5, units = 'in', res = 300)
ggplot(Ldat5, aes(x=year, y=LBS, colour=model,group=model)) + 
  geom_line() + geom_point()+
  ylab("Landings (lbs)")+xlab("") +
  ggtitle("") + theme_bw()
dev.off()


dev.new(width=5, height=5)
png("lv_hist_l.png", width = 5, height = 5, units = 'in', res = 300)
ggplot(Ldat5, aes(x=year, y=cpue, colour=model,group=model)) + 
  geom_line() + geom_point()+
  ylab("Landings (lbs)/Large Vessel")+xlab("") +
  ggtitle("") + theme_bw()
dev.off()

