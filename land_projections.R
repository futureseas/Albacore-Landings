# load input data to project landings

#Small vessels
years=1995:2100
sdmat.m=read.csv("/home/desiree/Documents/COCA/Albacore/sdmat_mproj.csv")
sdmat.s=read.csv("/home/desiree/Documents/COCA/Albacore/sdmat_sproj.csv")
sdmat.s=sdmat.s[,2:5]
sdmat.m=sdmat.m[,2:5]

PPDm=as.data.frame(Sdat12 %>% group_by(PORT_NAME)%>%summarize(PPD=mean(PPD)))
names(sdmat.s)[1]="yr"
names(sdmat.s)[4]="esm"
names(sdmat.s)[3]="PORT_NAME"
sdmat.s$PORT_NAME[sdmat.s$PORT_NAME=="CHARLESTON"] = "CHARLESTON (COOS BAY)"
sdmat.s$PORT_NAME[sdmat.s$PORT_NAME=="GARIBALDI"] = "GARIBALDI (TILLAMOOK)"

datds=read.csv("/home/desiree/Documents/COCA/Albacore/predmat_small.csv")
names(datds)[7]="PORT_NAME"
datds$PORT_NAME[datds$PORT_NAME=="CHARLESTON"] = "CHARLESTON (COOS BAY)"
datds$PORT_NAME[datds$PORT_NAME=="GARIBALDI"] = "GARIBALDI (TILLAMOOK)"
tmp1=datds%>%filter(yr>1994)

smatp=merge(tmp1,PPDm)
smatp2=merge(smatp,sdmat.s)

smatp2$chinook=rep(2340,34980)
smatp2$dunge=rep(20896,34980)
smatp2$biomass=rep(60172,34980)

tzp=read.csv("/home/desiree/Documents/COCA/Albacore/tzIndex_Projections.csv")
names(tzp)=c("yr","gfdl","ipsl","hadl")
tzpl <- gather(tzp, esm, tz, gfdl:hadl)
#note that we used detrended tz in fitting
tzpl$tzd=tzpl$tz-(5.172474+0.013208*tzpl$yr)
#only starts in 2006 so add the historical
tzplh=data.frame(yr=c(rep(Sdat12$year[1:11],3)),esm=c(rep("gfdl",11),rep("ipsl",11),rep("hadl",11)),tz=c(rep(Sdat12$tz[1:11],3)),tzd=c(rep(Sdat12$tz[1:11],3)))
tzpl2=rbind(tzpl,tzplh)
tzpl3=tzpl2[,c(1,2,4)]
names(tzpl3)[3]="tz"
#TZ  moves out of observed historical range, cap it at 0.922
tzpl3$tz[which(tzpl3$tz>0.922)]=0.922
#tzpl3$tz[which(tzpl3$tz>0.7)]=0.7

smatp3=merge(smatp2,tzpl3)
names(smatp3)[1]="year"

smatp3$vb602[which(smatp3$vb602<20)]=20


#medium Vessels
datdm=read.csv("/home/desiree/Documents/COCA/Albacore/predmat_med.csv")
names(datdm)[7]="Port"
tmp2=datdm%>%filter(yr>1994)

names(sdmat.m)[1]="yr"
names(sdmat.m)[4]="esm"

mmatp=merge(tmp2,sdmat.m)
mmatp2=merge(mmatp,tzpl3)
mmatp2$biomass=rep(60172,25440)
mmatp2$dunge=rep(20896,25440)

#since there is uncertainty on the effects of biomass over the fishing grounds below 140, to prevent extrapolation, all those less than 140 made 140.
mmatp2$vb1506[which(mmatp2$vb1506<140)]=140

#Large Vessels
datdl=read.csv("/home/desiree/Documents/COCA/Albacore/predmat_lar.csv")
names(datdl)[7]="Port"
datdl$Port[datdl$Port=="CHARLESTON"] = "CHARLESTON (COOS BAY)"
tmp3=datdl%>%filter(yr>1994)

tmp3$biomass=rep(60172,28620)

#fit
# model comes from the alb_landings_model_selection.R model
##############################SMALL VESSELS###############################################################################
cpuep.s=as.numeric(predict(m1.s,smatp3,type='response'))
smatp3$cpuep=cpuep.s
#find average vessels per port
Nves.s=as.data.frame(Sdat12 %>% group_by(PORT_NAME)%>%summarize(NVes=mean(NVes)))
#max cpue assuming a 22.5 mt capacity
Nves.s$cpuemax=Nves.s$NVesm*22.5*2204.62*91
#combine the two data sets
proj.s=merge(smatp3,Nves.s)
#find landings in pounds
proj.s$LBS=proj.s$cpuep*proj.s$NVes*proj.s$prob

#calculate the average landings by port for the baseline period - 1995-2018
lbsh.s=as.data.frame(proj.s%>%filter(year<2019)%>% group_by(PORT_NAME)%>%summarize(lbshm=mean(LBS)))

#merge data sets to have column of mean hisotrical landings by port
proj.s2=merge(proj.s,lbsh.s)
proj.s2$lbsa=proj.s2$LBS-proj.s2$lbshm

#need to also be able to compare historical landings to the landings observed assuming the same Nves
names(Nves.s)[2]="NVesm"
obs.s=merge(Sdat,Nves.s)
obs.s$LBS2=obs.s$NVesm*obs.s$cpue
#compute mean observed lbs by port
lbs.opm=as.data.frame(obs.s%>% group_by(PORT_NAME)%>%summarize(lbsom=mean(LBS2)))

#create a historical comparison mat
#add observed mean as column in hmat.s
proj.s3=merge(lbs.opm,proj.s2)
hmat.s=as.data.frame(proj.s3%>%filter(year<2019))
hmat.s$LBS2=hmat.s$lbsa+hmat.s$lbsom
hmat.s=hmat.s[,-c(2,5:27)]

tmp=data.frame(PORT_NAME=obs.s$PORT_NAME,year=obs.s$year,esm=rep("obs",264),LBS2=obs.s$LBS2)
hmato.s=rbind(hmat.s,tmp)

#need to add lat by port for plotting
latm=obs.s%>%group_by(PORT_NAME)%>%summarize(lat=mean(lat))
hmato2.s=merge(hmato.s,latm)

portlab=c("San Diego","Long Beach","Morro Bay","Moss Landing","Eureka","Coos Bay","Newport","Tillamook","Astoria","Ilwaco","Westport")

setwd("/home/desiree/Documents/COCA/Albacore/")
dev.new(width=9, height=8)
png("hcomp_lbs_s.png", width = 7, height = 6, units = 'in', res = 300)
ggplot(hmato2.s, aes(x=as.factor(lat), y=LBS2, fill=esm)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Earth System Model")+
  labs(x = '', y = 'Landings (LBS)')+
  theme_bw()+scale_x_discrete(labels=portlab)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Historical landings by port - Vessels < 45 ft")
dev.off()

#compute the relative change in landings by port
proj.s3$lbsr=proj.s3$lbsa/proj.s3$lbshm
proj.s3$time = "2041-2070"
proj.s3$time[which(proj.s3$year >2070)] = "2071-2100"
proj.s3$time[which(proj.s3$year <2020)] = "historical"
proj.s3$time[which(proj.s3$year >2019&proj.s3$year<2041)] = "2020-2040"

#sselect years for plotting
proj.s4100=proj.s3 %>%filter(year>2040)

#compute a histogram of future landings per vessel to compare to historical
dev.new(width=9, height=8)
png("proj_lbs_vsl_s_esm.png", width = 7, height = 6, units = 'in', res = 300)
ggplot(proj.s4100, aes(x=cpuep, color=esm, fill=esm)) +
  geom_histogram() +
  theme_bw() +
  xlab("Projected landings (LBS) per vessel") +
  ylab("count") +
  ggtitle("Vessels < 45 ft")+
  facet_wrap(~esm)
dev.off()


dev.new(width=9, height=8)
png("proj_lbs_s_esm.png", width = 7, height = 6, units = 'in', res = 300)
ggplot(proj.s4100, aes(x=as.factor(lat), y=lbsr*100, fill=time)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Projection Period")+
  facet_grid(esm~ .)+
  labs(x = '', y = '% Change in Landings by Port \n Relative to 1995-2018')+
  theme_bw()+scale_x_discrete(labels=portlab)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Projected change in landings by port - Vessels < 45 ft")+
  geom_hline(yintercept=0, linetype="dashed",color="black")
dev.off()

dev.new(width=9, height=8)
png("proj_lbs_s.png", width = 7, height = 6, units = 'in', res = 300)
ggplot(proj.s4100, aes(x=as.factor(lat), y=lbsr*100, fill=esm)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Earth System Model")+
  facet_grid(time~ .)+
  labs(x = '', y = '% Change in Landings by Port \n Relative to 1995-2018')+
  theme_bw()+scale_x_discrete(labels=portlab)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Projected change in landings by port - Vessels < 45 ft")+
  geom_hline(yintercept=0, linetype="dashed",color="black")
dev.off()

dev.new(width=9, height=8)
png("proj_lbs_all_s.png", width = 7, height = 6, units = 'in', res = 300)
ggplot(proj.s4100, aes(x=as.factor(lat), y=lbsr*100,fill=time)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Earth System Model")+
  labs(x = '', y = '% Change in Landings by Port \n Relative to 1995-2018')+
  theme_bw()+scale_x_discrete(labels=portlab)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Projected change in landings by port - Vessels < 45 ft")+
  geom_hline(yintercept=0, linetype="dashed",color="black")
dev.off()

#Plot of the input variables
proj.s3p = proj.s3[!(proj.s3$time=="2020-2040"),]

dev.new(width=9, height=8)
png("proj_distw_s.png", width = 7, height = 6, units = 'in', res = 300)
ggplot(proj.s3p, aes(x=as.factor(lat), y=distw, fill=time)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Projection Period")+
  facet_grid(esm~ .)+
  labs(x = '', y = 'Distance Index')+
  theme_bw()+scale_x_discrete(labels=portlab)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Vessels < 45 ft")
dev.off()

dev.new(width=9, height=8)
png("proj_vb602_s.png", width = 7, height = 6, units = 'in', res = 300)
ggplot(proj.s3p, aes(x=as.factor(lat), y=vb602, fill=esm)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Projection Period")+
  facet_grid(time~ .)+
  labs(x = '', y = 'Availability over fishing grounds')+
  theme_bw()+scale_x_discrete(labels=portlab)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Vessels < 45 ft")
dev.off()

dev.new(width=9, height=8)
png("proj_tz_s.png", width = 7, height = 6, units = 'in', res = 300)
ggplot(proj.s3p, aes(x=as.factor(lat), y=tz, fill=time)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Projection Period")+
  facet_grid(esm~ .)+
  labs(x = '', y = 'Transition zone location')+
  theme_bw()+scale_x_discrete(labels=portlab)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Vessels < 45 ft")
dev.off()

proj.s.smry = as.data.frame(proj.s4100 %>% group_by(PORT_NAME,time) %>% summarize(LBS=round(mean(LBS)),lbsa=round(mean(lbsa)),lbsr=round(mean(lbsr),digits=2),Nves=round(mean(NVes))))
write.csv(proj.s.smry,"mean_proj_2041_2100_s.csv")


############################Medium Vessels
cpuep.m=as.numeric(predict(m1.ma,mmatp2,type='response'))
mmatp2$cpuep=cpuep.m
#find average vessels per port
Nves.m=as.data.frame(Mdat12 %>% group_by(Port)%>%summarize(NVes=mean(NVes)))
#combine the two data sets
proj.m=merge(mmatp2,Nves.m)
#find landings in pounds
proj.m$LBS=proj.m$cpuep*proj.m$NVes*proj.m$prob

#calculate the average landings,vb150,pc1, and distw by port for the baseline period - 1995-2018
lbsh.m=as.data.frame(proj.m%>%filter(yr<2019)%>% group_by(Port)%>%summarize(lbshm=mean(LBS),
                                                                            vb1506hm=mean(vb1506),
                                                                            tzhm=mean(tz),
                                                                            distwhm=mean(distw)))

#merge data sets to have column of mean hisotrical landings by port
proj.m2=merge(proj.m,lbsh.m)
proj.m2$lbsa=proj.m2$LBS-proj.m2$lbshm

#need to also be able to compare historical landings to the landings observed assuming the same Nves
names(Nves.m)[2]="NVesm"
obs.m=merge(Mdat,Nves.m)
obs.m$LBS2=obs.m$NVesm*obs.m$cpue
#compute mean observed lbs by port
lbs.opmm=as.data.frame(obs.m%>% group_by(Port)%>%summarize(lbsom=mean(LBS2)))

#create a historical comparison mat
#add observed mean as column in hmat.s
proj.m3=merge(lbs.opmm,proj.m2)
hmat.m=as.data.frame(proj.m3%>%filter(yr<2019))
hmat.m$LBS2=hmat.m$lbsa+hmat.m$lbsom
hmat.m=hmat.m[,c(1,3,4,29)]

tmp.m=data.frame(Port=obs.m$Port,yr=obs.m$year,esm=rep("obs",192),LBS2=obs.m$LBS2)
hmato.m=rbind(hmat.m,tmp.m)

#need to add lat by port for plotting
latmm=obs.m%>%group_by(Port)%>%summarize(lat=mean(lat))
hmato2.m=merge(hmato.m,latmm)

portlabm=c("San Diego","Long Beach","Moss Landing","Coos Bay","Newport","Astoria","Ilwaco","Westport")

dev.new(width=9, height=8)
png("hcomp_lbs_m.png", width = 7, height = 6, units = 'in', res = 300)
ggplot(hmato2.m, aes(x=as.factor(lat), y=LBS2, fill=esm)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Earth System Model")+
  labs(x = '', y = 'Landings (LBS)')+
  theme_bw()+scale_x_discrete(labels=portlabm)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Historical landings by port - Vessels 45-60 ft")
dev.off()

#compute the relative change in landings by port
proj.m3$lbsr=proj.m3$lbsa/proj.m3$lbshm
proj.m3$vb1506r=(proj.m3$vb1506-proj.m3$vb1506hm)/proj.m3$vb1506hm
proj.m3$tzr=(proj.m3$tz-proj.m3$tzhm)/proj.m3$tzhm
proj.m3$distwr=(proj.m3$distw-proj.m3$distwhm)/proj.m3$distwhm
proj.m3$time = "2041-2070"
proj.m3$time[which(proj.m3$yr >2070)] = "2071-2100"
proj.m3$time[which(proj.m3$yr <2020)] = "historical"
proj.m3$time[which(proj.m3$year >2019&proj.m3$year<2041)] = "2020-2040"

#sselect years for plotting
proj.m4100=proj.m3 %>%filter(yr>2040)

#compute a histogram of future landings per vessel to compare to historical
dev.new(width=9, height=8)
png("proj_lbs_vsl_m_esm.png", width = 7, height = 6, units = 'in', res = 300)
ggplot(proj.m4100, aes(x=cpuep, color=esm, fill=esm)) +
  geom_histogram() +
  theme_bw() +
  xlab("Projected landings (LBS) per vessel") +
  ylab("count") +
  ggtitle("Vessels 45-60 ft")+
  facet_wrap(~esm)
dev.off()

dev.new(width=9, height=8)
png("proj_lbs_m_esm.png", width = 7, height = 6, units = 'in', res = 300)
ggplot(proj.m4100, aes(x=as.factor(lat), y=lbsr*100, fill=time)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Projection Period")+
  facet_grid(esm~ .)+
  labs(x = '', y = '% Change in Landings by Port \n Relative to 1995-2018')+
  theme_bw()+scale_x_discrete(labels=portlabm)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Projected change in landings by port - Vessels 45-60 ft")+
  geom_hline(yintercept=0, linetype="dashed",color="black")
dev.off() #soutern ports, was low remains low

dev.new(width=9, height=8)
png("proj_lbs_m.png", width = 7, height = 6, units = 'in', res = 300)
ggplot(proj.m4100, aes(x=as.factor(lat), y=lbsr, fill=esm)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Earth System Model")+
  facet_grid(time~ .)+
  labs(x = '', y = 'Change in Landings by Port \n Relative to 1995-2018')+
  theme_bw()+scale_x_discrete(labels=portlabm)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Projected change in landings by port - Vessels 45-60 ft")+
  geom_hline(yintercept=0, linetype="dashed",color="black")
dev.off()

ggplot(proj.m4100, aes(x=as.factor(lat), y=lbsr, fill=esm)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Earth System Model")+
  facet_grid(time~ .)+
  labs(x = '', y = 'Change in Landings by Port \n Relative to 1995-2018')+
  theme_bw()+scale_x_discrete(labels=portlabm)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Projected change in landings by port - Vessels 45-60 ft")+
  geom_hline(yintercept=0, linetype="dashed",color="black")
ggsave("proj_lbs_m.eps")



dev.new(width=9, height=8)
png("proj_lbs_m.png", width = 7, height = 6, units = 'in', res = 300)
ggplot(proj.m4100, aes(x=as.factor(lat), y=lbsr, fill=esm)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Earth System Model")+
  facet_grid(time~ .)+
  labs(x = '', y = 'Change in Landings by Port \n Relative to 1995-2018')+
  theme_bw()+scale_x_discrete(labels=portlabm)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Projected change in landings by port - Vessels 45-60 ft")+
  geom_hline(yintercept=0, linetype="dashed",color="black")
dev.off()

dev.new(width=9, height=8)
png("proj_lbs_all_m.png", width = 7, height = 6, units = 'in', res = 300)
ggplot(proj.m4100, aes(x=as.factor(lat), y=lbsr*100, fill=time)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Earth System Model")+
  labs(x = '', y = '% Change in Landings by Port \n Relative to 1995-2018')+
  theme_bw()+scale_x_discrete(labels=portlabm)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Projected change in landings by port - Vessels 45-60 ft")+
  geom_hline(yintercept=0, linetype="dashed",color="black")
dev.off()

#Plot of the input variables
proj.m3p = proj.m3[!(proj.m3$time=="2020-2040"),]

dev.new(width=9, height=8)
png("proj_distw_m.png", width = 7, height = 6, units = 'in', res = 300)
ggplot(proj.m3p, aes(x=as.factor(lat), y=distw, fill=time)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Projection Period")+
  facet_grid(esm~ .)+
  labs(x = '', y = 'Distance Index')+
  theme_bw()+scale_x_discrete(labels=portlab)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Vessels 45-60 ft")
dev.off()

dev.new(width=9, height=8)
png("proj_vb1506_m.png", width = 7, height = 6, units = 'in', res = 300)
ggplot(proj.m3p, aes(x=as.factor(lat), y=vb1506, fill=time)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Projection Period")+
  facet_grid(esm~ .)+
  labs(x = '', y = 'Availability over fishing grounds')+
  theme_bw()+scale_x_discrete(labels=portlab)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Vessels 45-60 ft")
dev.off()

dev.new(width=9, height=8)
png("proj_tz_m.tiff", width = 7, height = 6, units = 'in', res = 300)
ggplot(proj.m3p, aes(x=as.factor(lat), y=tz, fill=time)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Projection Period")+
  facet_grid(esm~ .)+
  labs(x = '', y = 'Transition zone location')+
  theme_bw()+scale_x_discrete(labels=portlab)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Vessels 45-60 ft")
dev.off()


proj.m.smry = as.data.frame(proj.m4100 %>% group_by(Port,time) %>% summarize(LBS=round(mean(LBS)),lbsa=round(mean(lbsa)),lbsr=round(mean(lbsr),digits=2),Nves=round(mean(NVes))))
write.csv(proj.m.smry,"mean_proj_2041_2100_m.csv")

#################################################Large Vessels
cpuep.l=as.numeric(predict(m1.la$gam,tmp3,type='response'))
tmp3$cpuep=cpuep.l
#find average vessels per port
Nves.l=as.data.frame(Ldat12 %>% group_by(Port)%>%summarize(NVes=mean(NVes)))
#combine the two data sets
proj.l=merge(tmp3,Nves.l)
#find landings in pounds
proj.l$LBS=proj.l$cpuep*proj.l$NVes*proj.l$prob

#calculate the average landings,vb150,pc1, and distw by port for the baseline period - 1995-2018
lbsh.l=as.data.frame(proj.l%>%filter(yr<2019)%>% group_by(Port)%>%summarize(lbshm=mean(LBS),
                                                                            distwhm=mean(distw)))

#merge data sets to have column of mean hisotrical landings by port
proj.l2=merge(proj.l,lbsh.l)
proj.l2$lbsa=proj.l2$LBS-proj.l2$lbshm

#need to also be able to compare historical landings to the landings observed assuming the same Nves
names(Nves.l)[2]="NVesm"
obs.l=merge(Ldat,Nves.l)
obs.l$LBS2=obs.l$NVesm*obs.l$cpue
#compute mean observed lbs by port
lbs.opml=as.data.frame(obs.l%>% group_by(Port)%>%summarize(lbsom=mean(LBS2)))

#create a historical comparison mat
#add observed mean as column in hmat.s
proj.l3=merge(lbs.opml,proj.l2)
hmat.l=as.data.frame(proj.l3%>%filter(yr<2019))
hmat.l$LBS2=hmat.l$lbsa+hmat.l$lbsom
hmat.l=hmat.l[,c(1,5,7,24)]

tmp.l=data.frame(Port=obs.l$Port,yr=obs.l$year,esm=rep("obs",205),LBS2=obs.l$LBS2)
hmato.l=rbind(hmat.l,tmp.l)

#need to add lat by port for plotting
latml=obs.l%>%group_by(Port)%>%summarize(lat=mean(lat))
hmato2.l=merge(hmato.l,latml)

portlabl=c("San Diego","San Pedro","Long Beach","Coos Bay","Newport","Astoria","Ilwaco","Westport","Bellingham Bay")

dev.new(width=9, height=8)
png("hcomp_lbs_l.tiff", width = 7, height = 6, units = 'in', res = 300)
ggplot(hmato2.l, aes(x=as.factor(lat), y=LBS2, fill=esm)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Earth System Model")+
  labs(x = '', y = 'Landings (LBS)')+
  theme_bw()+scale_x_discrete(labels=portlabl)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Historical landings by port - Vessels >60 ft")
dev.off()

#compute the relative change in landings by port
proj.l3$lbsr=proj.l3$lbsa/proj.l3$lbshm
proj.l3$distwr=(proj.l3$distw-proj.l3$distwhm)/proj.l3$distwhm
#sselect years for plotting
proj.l4100=proj.l3 %>%filter(yr>2040)
proj.l4100$time = "2041-2070"
proj.l4100$time[which(proj.l4100$yr >2070)] = "2071-2100"

# combine San Pedro 
ind = which(proj.l4100$Port=="SAN PEDRO")
proj.l41002=proj.l4100[-ind,]

#round the number of vessels
proj.l41002$NVes=round(proj.l41002$NVes)
proj.m4100$NVes=round(proj.m4100$NVes)
proj.s4100$NVes=round(proj.s4100$NVes)
write.csv(proj.l41002,"proj_2041_2100_l.csv")
write.csv(proj.m4100,"proj_2041_2100_m.csv")
write.csv(proj.s4100,"proj_2041_2100_s.csv")

#compute a histogram of future landings per vessel to compare to historical
dev.new(width=9, height=8)
png("proj_lbs_vsl_l_esm.png", width = 7, height = 6, units = 'in', res = 300)
ggplot(proj.l41002, aes(x=cpuep, color=esm, fill=esm)) +
  geom_histogram() +
  theme_bw() +
  xlab("Projected landings (LBS) per vessel") +
  ylab("count") +
  ggtitle("Vessels > 60 ft")+
  facet_wrap(~esm)
dev.off()

dev.new(width=9, height=8)
png("proj_lbs_l_esm.png", width = 7, height = 6, units = 'in', res = 300)
ggplot(proj.l41002, aes(x=as.factor(lat), y=lbsr, fill=time)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Projection Period")+
  facet_grid(esm~ .)+
  labs(x = '', y = '% Change in Landings by Port \n Relative to 1995-2018')+
  theme_bw()+scale_x_discrete(labels=portlabl[-2])+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Projected change in landings by port - Vessels >60 ft")+
  geom_hline(yintercept=0, linetype="dashed",color="black")
dev.off() #soutern ports, was low remains low

dev.new(width=9, height=8)
png("proj_lbs_l.png", width = 7, height = 6, units = 'in', res = 300)
ggplot(proj.l41002, aes(x=as.factor(lat), y=lbsr*100, fill=esm)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Earth System Model")+
  facet_grid(time~ .)+
  labs(x = '', y = '% Change in Landings by Port \n Relative to 1995-2018')+
  theme_bw()+scale_x_discrete(labels=portlabl[-2])+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Projected change in landings by port - Vessels >60 ft")+
  geom_hline(yintercept=0, linetype="dashed",color="black")
dev.off()

ggplot(proj.l41002, aes(x=as.factor(lat), y=lbsr, fill=esm)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Earth System Model")+
  facet_grid(time~ .)+
  labs(x = '', y = 'Change in Landings by Port \n Relative to 1995-2018')+
  theme_bw()+scale_x_discrete(labels=portlabl[-2])+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Projected change in landings by port - Vessels >60 ft")+
  geom_hline(yintercept=0, linetype="dashed",color="black")
ggsave("proj_lbs_l.eps")

dev.new(width=9, height=8)
png("proj_lbs_all_l.png", width = 7, height = 6, units = 'in', res = 300)
ggplot(proj.l41002, aes(x=as.factor(lat), y=lbsr, fill=time)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Earth System Model")+
  labs(x = '', y = '% Change in Landings by Port \n Relative to 1995-2018')+
  theme_bw()+scale_x_discrete(labels=portlabl[-2])+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Projected change in landings by port - Vessels >60 ft")+
  geom_hline(yintercept=0, linetype="dashed",color="black")
dev.off()

proj.l.smry = as.data.frame(proj.l4100 %>% group_by(Port,time) %>% summarize(LBS=round(mean(LBS)),lbsa=round(mean(lbsa)),lbsr=round(mean(lbsr),digits=2),Nves=round(mean(NVes))))
write.csv(proj.l.smry,"mean_proj_2041_2100_l.csv")

#Combine landings from all vessels
#need to change port levels for medium vessels
proj.m3p=proj.m3
levels(proj.m3p$Port)[levels(proj.m3p$Port)=="Charleston"] = "CHARLESTON (COOS BAY)"
levels(proj.m3p$Port)[levels(proj.m3p$Port)=="Astoria"] = "ASTORIA"
levels(proj.m3p$Port)[levels(proj.m3p$Port)=="Ilwaco"] = "ILWACO"
levels(proj.m3p$Port)[levels(proj.m3p$Port)=="LongBeach"] = "LONG BEACH"
levels(proj.m3p$Port)[levels(proj.m3p$Port)=="MossLanding"] = "MOSS LANDING"
levels(proj.m3p$Port)[levels(proj.m3p$Port)=="Newport"] = "NEWPORT"
levels(proj.m3p$Port)[levels(proj.m3p$Port)=="SanDiego"] = "SAN DIEGO"
levels(proj.m3p$Port)[levels(proj.m3p$Port)=="Westport"] = "WESTPORT"

#make factors consistent
levels(proj.s3$PORT_NAME) <- c(levels(proj.s3$PORT_NAME), "BELLINGHAM BAY","SAN PEDRO") 
levels(proj.m3p$Port) <- c(levels(proj.m3p$Port), "BELLINGHAM BAY","MOSS LANDING","EUREKA","GARIBALDI (TILLAMOOK)","MORRO BAY") 
levels(proj.l3$Port) <- c(levels(proj.l3$Port), "BELLINGHAM BAY","MOSS LANDING","EUREKA","GARIBALDI (TILLAMOOK)","MORRO BAY")

#also have to check esm factors
proj.st=proj.s3[,c(1:4,7,9,26)]
names(proj.st)[1]="Port"
proj.st$vtype="small"
proj.mt=proj.m3p[,c(1:4,9,11,21)]
proj.mt$vtype="medium"
proj.lt=proj.l3[,c(1:2,5:7,9,21)]
proj.lt$vtype="large"

ptot=rbind(proj.st,proj.mt,proj.lt)

#summarize by port
lbstot=as.data.frame(ptot%>%group_by(Port,yr,iter,esm)%>%summarize(LBS=sum(LBS),lat=mean(lat)))

#find the historical total mean landings and calculate anomalies
htot=as.data.frame(lbstot%>%filter(yr<2019)%>%group_by(Port,esm)%>%summarize(LBSmh=mean(LBS)))
lbstoth=merge(lbstot,htot)
lbstoth$lbsta=lbstoth$LBS-lbstoth$LBSmh

#obs.l
obs.l3=obs.l
obs.s3=obs.s
obs.m3=obs.m
levels(obs.m3$Port)[levels(obs.m3$Port)=="Charleston"] = "CHARLESTON (COOS BAY)"
levels(obs.m3$Port)[levels(obs.m3$Port)=="Astoria"] = "ASTORIA"
levels(obs.m3$Port)[levels(obs.m3$Port)=="Ilwaco"] = "ILWACO"
levels(obs.m3$Port)[levels(obs.m3$Port)=="LongBeach"] = "LONG BEACH"
levels(obs.m3$Port)[levels(obs.m3$Port)=="MossLanding"] = "MOSS LANDING"
levels(obs.m3$Port)[levels(obs.m3$Port)=="Newport"] = "NEWPORT"
levels(obs.m3$Port)[levels(obs.m3$Port)=="SanDiego"] = "SAN DIEGO"
levels(obs.m3$Port)[levels(obs.m3$Port)=="Westport"] = "WESTPORT"

#make factors consistent
levels(obs.s3$PORT_NAME) <- c(levels(obs.s3$PORT_NAME), "BELLINGHAM BAY","SAN PEDRO") 
levels(obs.m3$Port) <- c(levels(obs.m3$Port), "BELLINGHAM BAY","MOSS LANDING","EUREKA","GARIBALDI (TILLAMOOK)","MORRO BAY") 
levels(obs.l3$Port) <- c(levels(obs.l3$Port), "BELLINGHAM BAY","MOSS LANDING","EUREKA","GARIBALDI (TILLAMOOK)","MORRO BAY")

#also have to check esm factors
obs.st=obs.s3[,c(1:2,7,32)]
names(obs.st)[1]="Port"
obs.st$vtype="small"
obs.mt=obs.m3[,c(1,3,8,33)]
obs.mt$vtype="medium"
obs.lt=obs.l3[,c(1:2,7,29)]
obs.lt$vtype="large"

otot=rbind(obs.st,obs.mt,obs.lt)

#summarize by port across vessel types
tmp7=as.data.frame(otot%>%group_by(Port,year)%>%summarize(LBS=sum(LBS2),lat=mean(lat)))
#find mean by port for total 
totp=as.data.frame(tmp7%>%group_by(Port)%>%summarize(LBSo=mean(LBS)))

tmp8=merge(lbstoth,totp)
tmp8$LBS2=tmp8$lbsta+tmp8$LBSo
tmp8$lbsr=tmp8$lbsta/tmp8$LBSmh

#sselect years for plotting
proj.4100=tmp8 %>%filter(yr>2040)
proj.4100$time = "2041-2070"
proj.4100$time[which(proj.4100$yr >2070)] = "2071-2100"

portlabt=c("San Diego","San Pedro","Long Beach","Morro Bay","Moss Landing","Eureka","Coos Bay","Newport","Tillamook","Astoria","Ilwaco","Westport","Bellingham Bay")

dev.new(width=9, height=8)
png("proj_lbs_tot_esm.png", width = 7, height = 6, units = 'in', res = 300)
ggplot(proj.4100, aes(x=as.factor(lat), y=lbsr, fill=time)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Projection Period")+
  facet_grid(esm~ .)+scale_x_discrete(labels=portlabt)+
  labs(x = '', y = '% Change in Landings by Port \n Relative to 1995-2018')+
  theme_bw()+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Projected change in landings by port")+
  geom_hline(yintercept=0, linetype="dashed",color="black")
dev.off() #soutern ports, was low remains low

dev.new(width=9, height=8)
png("proj_lbs_tot.png", width = 7, height = 6, units = 'in', res = 300)
ggplot(proj.4100, aes(x=as.factor(lat), y=lbsr, fill=esm)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_discrete(name="Earth System Model")+
  facet_grid(time~ .)+
  labs(x = '', y = '% Change in Landings by Port \n Relative to 1995-2018')+
  theme_bw()+scale_x_discrete(labels=portlabl)+
  theme(axis.text.x  = element_text(angle=90))+
  ggtitle("Projected change in landings by port")+
  geom_hline(yintercept=0, linetype="dashed",color="black")
dev.off()

lbsa_proj_4100=proj.4100[,c(1:4,8,11,12)]
write.csv(lbsa_proj_4100,"alb_la_proj_2041_2100.csv")
