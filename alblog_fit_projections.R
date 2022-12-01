#code to compute if landings are 0 or not based on distance weighted by fuel costs
library(ROCR)
library(mgcv)
library(matrixStats)
library(tidyverse)
library(ggplot2)
library(ape)
#read the updated data frames for all vessel lengths
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

#logistic fit to see if port operational or not -Small vessels####################################################
#create a training and testing set
#remove NaNs
Sdatc = Sdatc[-( which(is.na(Sdatc$oper))),]
auc.s = 1:20
fit.s = matrix(NA,ncol=20,nrow = 259)
ci5.s = matrix(NA,ncol=20,nrow = 259)
ci95.s = matrix(NA,ncol=20,nrow = 259)
mlist.s = vector(mode = "list", length = 20)
moran.s = vector(mode = "list", length = 20)

for (j in 1:20){
  tinds= sample(1:259, 194, replace = FALSE)
  train.s = Sdatc[tinds,]
  test.s = Sdatc[-tinds,]
  
  #Generate  the inverse distance matrix to calculate Moran's I
  Sdat.dists = as.matrix(dist(cbind(train.s$lon, train.s$lat)))
  Sdat.dists.inv = 1/Sdat.dists
  diag(Sdat.dists.inv) = 0
  Sdat.dists.inv[Sdat.dists.inv==Inf]=0 #some ports have same longitude so set to 0 
  
  #build binomial model
  lfits = gam(oper~s(distw,k=4), data=train.s, family = binomial)
  #save model object
  mlist.s[[j]]=lfits
  # make prediction for test set
  pred.test=as.numeric(predict(lfits,test.s,type='response')) 
  
  #make a ROCR prediction object using the predicted values from
  #our model and the true values from the real data
  rp = prediction(pred.test,test.s$oper) 
  #now calculate AUC
  auc.s[j] <- performance( rp, "auc")@y.values[[1]]
  
  #compute fitted values for the whole data set
  #note we need the answer on the response scale, not the default link scale
  fit.s[,j]=as.numeric(predict(lfits,Sdatc,type='response')) 
  #as the model is not linear, we can't use the simple +1.96*SE to get the 95% CI on the response scale
  #have to do it on the link scale and map it back to the response by taking the inverse of the link function
  #extract family object
  fam.s = family(lfits)
  #extract the inverse of the link function from the family object
  ilink.s = fam.s$linkinv
  #generate fitted values and standard errors on the link scale 
  #using predict(...., type = 'link'), which is the default
  lfit.s = predict(lfits,Sdatc,se.fit=TRUE)
  fit.link.s=lfit.s$fit
  se.link.s=lfit.s$se.fit
  #compute the confidence interval using these fitted values and standard error 
  #backtransform them to the response scale using the inverse of the link function
  ci5.s[,j] = ilink.s(fit.link.s - (2 * se.link.s))
  ci95.s[,j] = ilink.s(fit.link.s + (2 * se.link.s))
  
  #save the Moran's I output
  moran.s[[j]]=Moran.I(resid(mlist.s[[j]]), Sdat.dists.inv)
}

#compute average AUC of test set
mean(auc.s) 

#compute ensemble mean of fitted values
mfit.s = rowMeans(fit.s)

#compute mean ci
mci5.s = rowMeans(ci5.s)
mci95.s = rowMeans(ci95.s)

#create data frame for plotting
#pfdats = data.frame(fit=mfit.s,q5 = qfit.s$`5%`, q95 = qfit.s$`95%`, distw=Sdatc$distw)
pfdats = data.frame(fit=mfit.s,q5 = mci5.s, q95 = mci95.s, distw=Sdatc$distw)
pfdats2 = pfdats[order(pfdats$distw),]

#create logical vector of port beign active or not for plotting
Sdatc = mutate(Sdatc, loper = as.logical(oper))
setwd("/home/desiree/Documents/COCA/Albacore/")
dev.new(width=7, height=6)
tiff("probop.s_pred.tiff", width = 6, height = 5, units = 'in', res = 300, compression = 'lzw')
ggplot(pfdats2, aes(x = distw, y = fit)) +
  geom_line() +
  geom_rug(aes(y = oper, colour = loper), data = Sdatc) +
  scale_colour_discrete(name = 'Observed \nPort Active') +
  labs(x = 'Cost Index', y = 'Probability of Port Being Active')+
  ggtitle("Estimated probability of port being active \n as a function of cost - Vessels <45 ft")+
  geom_ribbon(data = pfdats2, aes(ymin = q5, ymax = q95),alpha = 0.1)+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+geom_text(x=23, y=0.9, label="AUC = 0.93")
dev.off()

#compare mean observed and predicted probability
pfdats2$port = Sdatc$PORT_NAME
pfdats2$lat=Sdatc$lat
pfdats2$type=rep("Estimated",259)
pfdats3=rbind(pfdats2,pfdats2)
pfdats3$fit[260:518]=Sdatc$oper
pfdats3$type[260:518]=rep("Observed",259)

ggplot(pfdats3, aes(x=as.factor(lat), y=fit, fill=type)) + 
  geom_boxplot() + scale_fill_discrete(name="")+
  labs(x = '', y = 'Probability of Port Being Active')+
  theme_bw()+#scale_x_discrete(labels=portlab)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Historical (1995-2018) Observed and Estimated Probability of port being active - Vessels < 45 ft")


##########################################################################################################################
#medium vessels
#create a training and testing set
auc.m = 1:20
fit.m = matrix(NA,ncol=20,nrow = 192)
ci5.m = matrix(NA,ncol=20,nrow = 192)
ci95.m = matrix(NA,ncol=20,nrow = 192)
mlist.m <- vector(mode = "list", length = 20)
moran.m = vector(mode = "list", length = 20)
for (j in 1:20){
  #since there are few 0s have to make sure there are enough zeros
  Mdatc1 = Mdatc %>% filter(oper>0)
  Mdatc0 = Mdatc %>% filter(oper==0)
  tindm1= sample(1:157, 118, replace = FALSE)
  tindm0= sample(1:35, 25, replace = FALSE)
  train1.m = Mdatc1[tindm1,]
  test1.m = Mdatc1[-tindm1,]
  train0.m = Mdatc0[tindm0,]
  test0.m = Mdatc0[-tindm0,]
  train.m = rbind(train1.m,train0.m)
  test.m = rbind(test1.m,test0.m)
  
  #Generate  the inverse distance matrix to calculate Moran's I
  Mdat.dists = as.matrix(dist(cbind(train.m$lon, train.m$lat)))
  Mdat.dists.inv = 1/Mdat.dists
  diag(Mdat.dists.inv) = 0
  Mdat.dists.inv[Mdat.dists.inv==Inf]=0 #some ports have same longitude so set to 0 
  
  lfitm = gam(oper~s(distw,k=4), data=train.m, family = binomial)
  #save model object
  mlist.m[[j]]=lfitm
  
  pred.test=as.numeric(predict(lfitm,test.m,type='response'))                  
  rp = prediction(pred.test,test.m$oper) 
  auc.m[j] <- performance( rp, "auc")@y.values[[1]]
  fit.m[,j]=as.numeric(predict(lfitm,Mdatc,type='response')) 
  #as the model is not linear, we can't use the simple +1.96*SE to get the 95% CI on the response scale
  #have to do it on the link scale and map it back to the response by taking the inverse of the link function
  #extract family object
  fam.m = family(lfitm)
  #extract the inverse of the link function from the family object
  ilink.m = fam.m$linkinv
  #generate fitted values and standard errors on the link scale 
  #using predict(...., type = 'link'), which is the default
  lfit.m = predict(lfitm,Mdatc,se.fit=TRUE)
  fit.link.m=lfit.m$fit
  se.link.m=lfit.m$se.fit
  #compute the confidence interval using these fitted values and standard 
  #backtransform them to the response scale using the inverse of the link function
  ci5.m[,j] = ilink.m(fit.link.m - (2 * se.link.m))
  ci95.m[,j] = ilink.m(fit.link.m + (2 * se.link.m))
  
  #save the Moran's I output
  moran.m[[j]]=Moran.I(resid(mlist.m[[j]]), Mdat.dists.inv)
  
}

mean(auc.m) #0.98
mfit.m = rowMeans(fit.m)

#compute mean ci
mci5.m = rowMeans(ci5.m)
mci95.m = rowMeans(ci95.m)

#create data frame for plotting
#pfdats = data.frame(fit=mfit.s,q5 = qfit.s$`5%`, q95 = qfit.s$`95%`, distw=Sdatc$distw)
pfdatm = data.frame(fit=mfit.m,q5 = mci5.m, q95 = mci95.m, distw=Mdatc$distw)
pfdatm2 = pfdatm[order(pfdatm$distw),]

#create logical vector of port beign active or not for plotting
Mdatc = mutate(Mdatc, loper = as.logical(oper))
setwd("/home/desiree/Documents/COCA/Albacore/")
dev.new(width=7, height=6)
tiff("probop.m_pred.tiff", width = 6, height = 5, units = 'in', res = 300, compression = 'lzw')
ggplot(pfdatm2, aes(x = distw, y = fit)) +
  geom_line() +
  geom_rug(aes(y = oper, colour = loper), data = Mdatc) +
  scale_colour_discrete(name = 'Observed \nPort Active') +
  labs(x = 'Cost Index', y = 'Probability of Port Being Active')+
  ggtitle("Estimated probability of port being active \n as a function of cost - Vessels 45-60 ft")+
  geom_ribbon(data = pfdatm2, aes(ymin = q5, ymax = q95),alpha = 0.1)+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+geom_text(x=23, y=0.9, label="AUC = 0.94")
dev.off()

#################################################################################################################
#Large vessels
#create a training and testing set
auc.l = 1:20
fit.l = matrix(NA,ncol=20,nrow = 205)
ci5.l = matrix(NA,ncol=20,nrow = 205)
ci95.l = matrix(NA,ncol=20,nrow = 205)
mlist.l <- vector(mode = "list", length = 20)
moran.l = vector(mode = "list", length = 20)

for (j in 1:20){
  #since there are few 0s have to make sure there are enough zeros
  Ldatc1 = Ldatc %>% filter(oper>0)
  Ldatc0 = Ldatc %>% filter(oper==0)
  tindl1= sample(1:138, 103, replace = FALSE)
  tindl0= sample(1:67, 50, replace = FALSE)
  train1.l = Ldatc1[tindl1,]
  test1.l = Ldatc1[-tindl1,]
  train0.l = Ldatc0[tindl0,]
  test0.l = Ldatc0[-tindl0,]
  train.l = rbind(train1.l,train0.l)
  test.l = rbind(test1.l,test0.l)
  
  #Generate  the inverse distance matrix to calculate Moran's I
  Ldat.dists = as.matrix(dist(cbind(train.l$lon, train.l$lat)))
  Ldat.dists.inv = 1/Ldat.dists
  diag(Ldat.dists.inv) = 0
  Ldat.dists.inv[Ldat.dists.inv==Inf]=0 #some ports have same longitude so set to 0 
  
  
  lfitl = gam(oper~s(distw,k=4), data=train.l, family = binomial)
  #save model object
  mlist.l[[j]]=lfitl
  
  pred.test=as.numeric(predict(lfitl,test.l,type='response'))                  
  rp = prediction(pred.test,test.l$oper) 
  auc.l[j] <- performance( rp, "auc")@y.values[[1]]
  fit.l[,j]=as.numeric(predict(lfitl,Ldatc,type='response')) 
  #as the model is not linear, we can't use the simple +1.96*SE to get the 95% CI on the response scale
  #have to do it on the link scale and map it back to the response by taking the inverse of the link function
  #extract family object
  fam.l = family(lfitl)
  #extract the inverse of the link function from the family object
  ilink.l = fam.l$linkinv
  #generate fitted values and standard errors on the link scale 
  #using predict(...., type = 'link'), which is the default
  lfit.l = predict(lfitl,Ldatc,se.fit=TRUE)
  fit.link.l=lfit.l$fit
  se.link.l=lfit.l$se.fit
  #compute the confidence interval using these fitted values and standard 
  #backtransform them to the response scale using the inverse of the link function
  ci5.l[,j] = ilink.l(fit.link.l - (2 * se.link.l))
  ci95.l[,j] = ilink.l(fit.link.l + (2 * se.link.l))
  
  #save the Moran's I output
  moran.l[[j]]=Moran.I(resid(mlist.l[[j]]), Ldat.dists.inv)
  
}

mean(auc.l) #0.99
mfit.l = rowMeans(fit.l)

#compute mean ci
mci5.l = rowMeans(ci5.l)
mci95.l = rowMeans(ci95.l)

#create data frame for plotting
#pfdats = data.frame(fit=mfit.s,q5 = qfit.s$`5%`, q95 = qfit.s$`95%`, distw=Sdatc$distw)
pfdatl = data.frame(fit=mfit.l,q5 = mci5.l, q95 = mci95.l, distw=Ldatc$distw)
pfdatl2 = pfdatl[order(pfdatl$distw),]

#create logical vector of port beign active or not for plotting
Ldatc = mutate(Ldatc, loper = as.logical(oper))
setwd("/home/desiree/Documents/COCA/Albacore/")
dev.new(width=7, height=6)
tiff("probop.l_pred.tiff", width = 6, height = 5, units = 'in', res = 300, compression = 'lzw')
ggplot(pfdatl2, aes(x = distw, y = fit)) +
  geom_line() +
  geom_rug(aes(y = oper, colour = loper), data = Ldatc) +
  scale_colour_discrete(name = 'Observed \nPort Active') +
  labs(x = 'Cost Index', y = 'Probability of Port Being Active')+
  ggtitle("Estimated probability of port being active \n as a function of cost - Vessels >60 ft")+
  geom_ribbon(data = pfdatl2, aes(ymin = q5, ymax = q95),alpha = 0.1)+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+geom_text(x=23, y=0.9, label="AUC = 0.90")
dev.off()

#############load COGs from Barb's analysis#####################################################
setwd("/home/desiree/Documents/COCA/Albacore/")
cogpred1 = read.csv("projectedFutureCOGs_10iter_noSeed_1991_2005.csv")
cogpred2 = read.csv("projectedFutureCOGs_10iter_noSeed.csv")

cogpred = rbind(cogpred1,cogpred2)

##############Small Vessels####################################################
#need to generate a table of probabilities for each port/vessel length combination
sport=as.character(unique(Sdatc$PORT_NAME))
slat=unique(Sdatc$lat)
slon1=unique(Sdatc$lon)
#Eureka and Coos Bay ahave the same lon and Ilwaco and Newport, so have to create new vector
slon=c(slon1[1:2],slon1[2],slon1[3:4],slon1[4],slon1[5:9])
#small vessel
#create a matrix for 11 ports x 110 years x 10 ensemble members x 3 GCMs
scogp = do.call(rbind, replicate(11, cogpred, simplify=FALSE))
scogp$cogx2=scogp$COGx-360
scogp$port=c(rep(sport[1],3300),rep(sport[2],3300),rep(sport[3],3300),rep(sport[4],3300),rep(sport[5],3300),
             rep(sport[6],3300),rep(sport[7],3300),rep(sport[8],3300),rep(sport[9],3300),rep(sport[10],3300),
             rep(sport[11],3300))
scogp$lat=c(rep(slat[1],3300),rep(slat[2],3300),rep(slat[3],3300),rep(slat[4],3300),rep(slat[5],3300),
             rep(slat[6],3300),rep(slat[7],3300),rep(slat[8],3300),rep(slat[9],3300),rep(slat[10],3300),
             rep(slat[11],3300))
scogp$lon=c(rep(slon[1],3300),rep(slon[2],3300),rep(slon[3],3300),rep(slon[4],3300),rep(slon[5],3300),
             rep(slon[6],3300),rep(slon[7],3300),rep(slon[8],3300),rep(slon[9],3300),rep(slon[10],3300),
             rep(slon[11],3300))
#for fuel price use the mean of the last 10 years. Note it differs by state and is adjusted by inflation
mean(Sdatc$fuelpa[15:24])#OR-1.26
mean(Sdatc$fuelpa[87:96])#WA - 1.29
mean(Sdatc$fuelpa[63:72])#CA - 1.4
scogp$fuelpa=c(rep(1.26,3300),rep(1.26,3300),rep(1.4,3300),rep(1.29,3300),rep(1.29,3300),
              rep(1.26,3300),rep(1.29,3300),rep(1.4,3300),rep(1.4,3300),rep(1.4,3300),
              rep(1.4,3300))
scogp$dist=sqrt((scogp$COGy-scogp$lat)^2+(scogp$cogx2-scogp$lon)^2)
scogp$distw = scogp$dist*scogp$fuelpa

pred.s=matrix(NA,ncol=14,nrow = 3135000)
pred.s= scogp[FALSE,-c(1,2,4)]
pred.s$iter2=double()
pred.s$mod=double()
pred.s$prob=double()
ind=seq(1,3135000,100)
ct=1
#100 different models to account for model uncertianity and 100 for parameter uncertainity
for (j in 1:20){
  #generated predicted values for the 100 models
  #to capture uncertainity we sample 100 outcomes from that distribution
  #extract family object
  fam = family(mlist.s[[j]])
  #extract the inverse of the link function from the family object
  ilink = fam$linkinv
  #generate fitted values and standard errors on the link scale 
  #using predict(...., type = 'link'), which is the default
  lfit = predict(mlist.s[[j]],scogp,se.fit=TRUE)
  fit.link=lfit$fit
  se.link=lfit$se.fit
  for (p in 1:31350){
    dis.link=rnorm(100,fit.link[p],se.link[p])
    dis=ilink(dis.link)
    pred.s[(ind[ct]:(ind[ct]+99)),1:11]=do.call(rbind, replicate(100, scogp[p,-c(1,2,4)], simplify=FALSE))
    pred.s[ind[ct]:(ind[ct]+99),12]=1:100
    pred.s[ind[ct]:(ind[ct]+99),13]=rep(j,100)
    pred.s[ind[ct]:(ind[ct]+99),14]=dis
    ct=ct+1
  }
}

#Need to add the 30 year chunks
pred.s1=pred.s %>%filter(yr>2010)
pred.s1$time = "2011-2040"
for (j in 1:2969900){
  if (pred.s1$yr[j] %in% (2041:2070) ) {pred.s1$time[j] = "2041-2070"}
  if (pred.s1$yr[j] >2070 ) {pred.s1$time[j] = "2071-2100"}
}
#mean model output only
predm.s= scogp[FALSE,-c(1,2,4)]
predm.s$iter2=double()
predm.s$mod=double()
predm.s$prob=double()

#100 different models to account for model uncertianity and 100 for parameter uncertainity
for (j in 1:20){
  #generated predicted values for the 100 models
  #to capture uncertainity we sample 100 outcomes from that distribution
  #extract family object
  fam = family(mlist.s[[j]])
  #extract the inverse of the link function from the family object
  ilink = fam$linkinv
  #generate fitted values and standard errors on the link scale 
  #using predict(...., type = 'link'), which is the default
  lfit = predict(mlist.s[[j]],scogp,se.fit=TRUE)
  fit.link=lfit$fit
  se.link=lfit$se.fit
  for (p in 1:36300){
    dis=ilink(fit.link[p])
    predm.s[p,1:11]=scogp[p,-c(1,2,4)]
    predm.s[p,12]=1
    predm.s[p,13]=j
    predm.s[p,14]=dis
  }
}

write.csv(predm.s,"predmat_small.csv")

predm.s=read.csv("predmat_small.csv")

predm.s$time = "2021-2040"
predm.s$time[which(predm.s$yr %in% (2041:2070))] = "2041-2070"
predm.s$time[which(predm.s$yr >2070)] = "2071-2100"
predm.s$time[which(predm.s$yr <2021)] = "Historical"

#historical (1995-2018) probabilites by port
probo.s=as.data.frame(Sdatc%>%group_by(PORT_NAME)%>%summarize(probo=mean(oper)))
names(probo.s)[1]="port"
#merge the two datasets by port
predm.s2=merge(predm.s,probo.s)
#calcualte % change in probability from historical
predm.s3=predm.s2%>%mutate(probr=((prob-probo)/probo)*100)

#do the same calcualation relative to the mean historical
probh.s = predm.s %>% filter(time=="Historical")%>%group_by(port)%>%summarize(probh=mean(prob))
#merge the two datasets by port
predm.s4=merge(predm.s3,probh.s)
#calcualte % change in probability from historical
predm.s5=predm.s4%>%mutate(probrh=((prob-probh)/probh)*100)

predm.s1=predm.s5 %>%filter(yr>2040)

#port label
portlab=c("San Diego","Long Beach","Morro Bay","Moss Landing","Eureka","Coos Bay","Newport","Tillamook","Astoria","Ilwaco","Westport")

dev.new(width=9, height=8)
tiff("proj_probs.tiff", width = 7, height = 6, units = 'in', res = 300, compression = 'lzw')
ggplot(predm.s1, aes(x=as.factor(lat), y=prob, fill=time)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Projection Period")+
  facet_grid(esm~ .)+
  labs(x = '', y = 'Estimated Probability of Port Being Active')+
  theme_bw()+scale_x_discrete(labels=portlab)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Projected probability of port being active - Vessels < 45 ft")+
  geom_hline(yintercept=0, linetype="dashed",color="black")
dev.off()

dev.new(width=9, height=8)
tiff("proj_probs_esm", width = 7, height = 6, units = 'in', res = 300, compression = 'lzw')
ggplot(predm.s1, aes(x=as.factor(lat), y=prob, fill=esm)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Earth System Model")+
  facet_grid(time~ .)+
  labs(x = '', y = 'Probability of Port Being Active')+
  theme_bw()+scale_x_discrete(labels=portlab)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Projected probability of port being active - Vessels < 45 ft")
dev.off()

dev.new(width=9, height=8)
tiff("proj_probrs.tiff", width = 7, height = 6, units = 'in', res = 300, compression = 'lzw')
ggplot(predm.s1, aes(x=as.factor(lat), y=probrh, fill=time)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Projection Period")+
  facet_grid(esm~ .)+
  labs(x = '', y = '% Change in Probability of Port Being Active \n Relative to Historical (1991-2020)')+
  theme_bw()+scale_x_discrete(labels=portlab)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Projected change in probability of port being active - Vessels < 45 ft")+
  geom_hline(yintercept=0, linetype="dashed",color="black")
dev.off()

#compare observed with historical
pred9518=predm.s5%>%filter(yr>1994&yr<2019)
pred9518h=pred9518[,1:16]
pred9518o=pred9518[,c(1:14,17)]
pred9518o$time="Observed"
names(pred9518o)[15]="prob"
predho=rbind(pred9518o,pred9518h)
ggplot(predho, aes(x=as.factor(lat), y=prob, fill=time)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="")+
  facet_grid(esm~ .)+
  labs(x = '', y = 'Probability of Port Being Active')+
  theme_bw()+scale_x_discrete(labels=portlab)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Probability of port being active - Vessels < 45 ft")+
  geom_hline(yintercept=0, linetype="dashed",color="black")


##############################Medium sized vessels
#need to generate a table of probabilities for each port/vessel length combination
mport=as.character(unique(Mdatc$Port))
mlat=unique(Mdatc$lat)
mlon1=unique(Mdatc$lon)
#Eureka and Coos Bay ahave the same lon and Ilwaco and Newport, so have to create new vector
mlon=c(mlon1[1:4],mlon1[3],mlon1[5:7])
#create a matrix for 8 ports x 110 years x 10 ensemble members x 3 GCMs
mcogp = do.call(rbind, replicate(8, cogpred, simplify=FALSE))
mcogp$cogx2=mcogp$COGx-360
mcogp$port=c(rep(mport[1],3300),rep(mport[2],3300),rep(mport[3],3300),rep(mport[4],3300),rep(mport[5],3300),
             rep(mport[6],3300),rep(mport[7],3300),rep(mport[8],3300))
mcogp$lat=c(rep(mlat[1],3300),rep(mlat[2],3300),rep(mlat[3],3300),rep(mlat[4],3300),rep(mlat[5],3300),
            rep(mlat[6],3300),rep(mlat[7],3300),rep(mlat[8],3300))
mcogp$lon=c(rep(mlon[1],3300),rep(mlon[2],3300),rep(mlon[3],3300),rep(mlon[4],3300),rep(mlon[5],3300),
            rep(mlon[6],3300),rep(mlon[7],3300),rep(mlon[8],3300))
#for fuel price use the mean of the last 10 years. Note it differs by state and is adjusted by inflation
mcogp$fuelpa=c(rep(1.26,3300),rep(1.26,3300),rep(1.29,3300),rep(1.4,3300),rep(1.26,3300),
               rep(1.4,3300),rep(1.29,3300),rep(1.4,3300))
mcogp$dist=sqrt((mcogp$COGy-mcogp$lat)^2+(mcogp$cogx2-mcogp$lon)^2)
mcogp$distw = mcogp$dist*mcogp$fuelpa

#mean model output only
predm.m= mcogp[FALSE,-c(1,2,4)]
predm.m$iter2=double()
predm.m$mod=double()
predm.m$prob=double()

#100 different models to account for model uncertianity and 100 for parameter uncertainity
for (j in 1:20){
  #generated predicted values for the 100 models
  #to capture uncertainity we sample 100 outcomes from that distribution
  #extract family object
  fam = family(mlist.m[[j]])
  #extract the inverse of the link function from the family object
  ilink = fam$linkinv
  #generate fitted values and standard errors on the link scale 
  #using predict(...., type = 'link'), which is the default
  lfit = predict(mlist.m[[j]],mcogp,se.fit=TRUE)
  fit.link=lfit$fit
  se.link=lfit$se.fit
  for (p in 1:26400){
    dis=ilink(fit.link[p])
    predm.m[p,1:11]=mcogp[p,-c(1,2,4)]
    predm.m[p,12]=1
    predm.m[p,13]=j
    predm.m[p,14]=dis
  }
}

write.csv(predm.m,"predmat_med.csv")

##############################Medium sized vessels
#need to generate a table of probabilities for each port/vessel length combination
lport=as.character(unique(Ldatc$Port))
llat=unique(Ldatc$lat)
llon1=unique(Ldatc$lon)
#ahave the same lon and Ilwaco and Newport, so have to create new vector
llon=c(llon1[1:4],llon1[4],llon1[5:8])
#create a matrix for 9 ports x 110 years x 10 ensemble members x 3 GCMs
lcogp = do.call(rbind, replicate(9, cogpred, simplify=FALSE))
lcogp$cogx2=lcogp$COGx-360
lcogp$port=c(rep(lport[1],3300),rep(lport[2],3300),rep(lport[3],3300),rep(lport[4],3300),rep(lport[5],3300),
             rep(lport[6],3300),rep(lport[7],3300),rep(lport[8],3300),rep(lport[9],3300))
lcogp$lat=c(rep(llat[1],3300),rep(llat[2],3300),rep(llat[3],3300),rep(llat[4],3300),rep(llat[5],3300),
            rep(llat[6],3300),rep(llat[7],3300),rep(llat[8],3300),rep(llat[9],3300))
lcogp$lon=c(rep(llon[1],3300),rep(llon[2],3300),rep(llon[3],3300),rep(llon[4],3300),rep(llon[5],3300),
            rep(llon[6],3300),rep(llon[7],3300),rep(llon[8],3300),rep(llon[9],3300))
#for fuel price use the mean of the last 10 years. Note it differs by state and is adjusted by inflation
lcogp$fuelpa=c(rep(1.26,3300),rep(1.29,3300),rep(1.26,3300),rep(1.29,3300),rep(1.26,3300),
               rep(1.29,3300),rep(1.4,3300),rep(1.4,3300),rep(1.4,3300))
lcogp$dist=sqrt((lcogp$COGy-lcogp$lat)^2+(lcogp$cogx2-lcogp$lon)^2)
lcogp$distw = lcogp$dist*lcogp$fuelpa

#mean model output only
predm.l= lcogp[FALSE,-c(1,2,4)]
predm.l$iter2=double()
predm.l$mod=double()
predm.l$prob=double()

#100 different models to account for model uncertianity and 100 for parameter uncertainity
for (j in 1:20){
  #generated predicted values for the 100 models
  #to capture uncertainity we sample 100 outcomes from that distribution
  #extract family object
  fam = family(mlist.l[[j]])
  #extract the inverse of the link function from the family object
  ilink = fam$linkinv
  #generate fitted values and standard errors on the link scale 
  #using predict(...., type = 'link'), which is the default
  lfit = predict(mlist.l[[j]],lcogp,se.fit=TRUE)
  fit.link=lfit$fit
  se.link=lfit$se.fit
  for (p in 1:29700){
    dis=ilink(fit.link[p])
    predm.l[p,1:11]=lcogp[p,-c(1,2,4)]
    predm.l[p,12]=1
    predm.l[p,13]=j
    predm.l[p,14]=dis
  }
}

write.csv(predm.l,"predmat_lar.csv")


  