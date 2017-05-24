set.seed(4951)

#figures
dir.create("figures")
libs <- c("ggplot2","colorspace","RColorBrewer","vegan","oce","xtable","gridExtra","plyr")
lapply(libs, require, character.only=TRUE)

#load data
chem <- read.csv("ChemData_DeepDOM.csv",header=TRUE)
maxes <- read.csv("FlaMaxRates_DeepDOM.csv",header=TRUE,row.names=1)
monmaxes <- read.csv("MonMaxRates_DeepDOM.csv",header=TRUE)
#format factor levels
maxes$stn <- factor(maxes$stn,levels=c("stn2","stn4","stn7","stn10","stn15","stn18"),labels=c("stn2","stn4","stn7","stn10","stn15","stn18"))
maxes$depthlabel <- factor(maxes$depthlabel,levels=c("SuW","DCM","meso","AAIW","NADW","bot"))
monmaxes$stn <- factor(monmaxes$stn,levels=c("stn2","stn4","stn7","stn10","stn15","stn18","stn21","stn22","stn23"))
monmaxes$depthlabel <- factor(monmaxes$depthlabel,levels=c("SuW","DCM","meso","AAIW","NADW","bot"))
#rename substrate factor to single character
maxes$substrate <- revalue(maxes$substrate, c("ara"="A","chon"="C","fuc"="F","lam"="L","pul"="P","xyl"="Y"))
monmaxes$substrate <- revalue(monmaxes$substrate, c("aglu"="G","Leu"="Le"))

#define color palette for substrates, stn and depth
#FlaColors <- brewer.pal(n=6,name="Set3")
FlaColors <- c("white","turquoise","green3","lightgoldenrod1","dodgerblue4","firebrick")
FlaColorsBold <- c("#8DD3C7","lightgoldenrod1","#BEBADA","#FB8072","#80B1D3","#FDB462") 
names(FlaColors) <- levels(maxes$substrate)
names(FlaColorsBold) <- levels(maxes$substrate)
MonColors <- c("maroon3","yellowgreen")
names(MonColors) <- levels(monmaxes$substrate)
#define how want color code stn and depth factors
cols.stn <- c("palegreen2","firebrick","plum4","salmon2","skyblue3","goldenrod3")
cols.dep <- c("light blue","cornflowerblue","midnightblue")

###################Figure 2 - FLA max hydrolysis panels by site###################
maxSub <- subset(maxes,(mean.kcRate.nM.hr)>0)
ann <- data.frame(stn="stn15",depthlabel="meso",substrate="Y",mean.kcRate.nM.hr=10)
#add value labels to low activity you can barely see in 
low_values <- data.frame(stn=c("stn2","stn2","stn4","stn10","stn10","stn15","stn18","stn18"),depthlabel=c("AAIW","NADW","AAIW","NADW","bot","NADW","AAIW","NADW"),substrate=rep("L",8),mean.kcRate.nM.hr=rep(4.3,8),label=c("0.4","0.1","0.4","0.2","0.1","0.2","0.6","0.6"))
errorbars <- geom_errorbar(aes(ymax=(mean.kcRate.nM.hr+sd.kcRate.nM.hr), ymin=(mean.kcRate.nM.hr-sd.kcRate.nM.hr)),width=0.5,color="grey25",alpha=0.6,size=0.4)
theme_barplot <- theme_bw() + theme(axis.text.y=element_text(size=9),strip.background=element_rect(fill="grey93"),axis.title.x=element_text(vjust=0.1),axis.title.y=element_text(vjust=0.5),legend.position="none",panel.grid.major=element_line(color="grey93"),panel.grid.minor=element_blank()) 

tiff("figures/Fig2_FLAbarplotPanels.tiff",width=5,height=4.2,units="in",res=1200)
fig2 <- ggplot(maxes,aes(x=substrate,y=mean.kcRate.nM.hr,fill=substrate)) + geom_bar(position="dodge",stat="identity",color="black",size=0.2) + facet_grid(depthlabel~stn) + errorbars + scale_fill_manual(name="substrate",values=FlaColors) + coord_cartesian(ylim=c(0,20)) + scale_y_continuous(breaks=c(0,5,10,15,20),labels=c("0","5","10","15","")) + labs(y=substitute(paste("Maximum Hydrolysis Rate (nM ",h^-1,")"))) + geom_text(data=ann, label="*") + geom_text(data=low_values,aes(label=label),size=2.2,family="Arial") + geom_point(data=maxSub,aes(x=substrate,y=kcRate.1.nM.hr),alpha=0.4,size=1.3) + geom_point(data=maxSub, aes(x=substrate,y=kcRate.2.nM.hr),alpha=0.4,size=1.3) + geom_point(data=maxSub, aes(x=substrate,y=kcRate.3.nM.hr),alpha=0.4,size=1.3) + theme_barplot
print(fig2)
dev.off()
#b/w version
tiff("figures/BWFig2_FLAbarplotPanels.tiff",width=5,height=4.2,units="in",res=1200)
bwfig2 <- ggplot(maxes,aes(x=substrate,y=mean.kcRate.nM.hr,fill=substrate)) + geom_bar(position="dodge",stat="identity",color="black",size=0.2) + facet_grid(depthlabel~stn) + errorbars + scale_fill_grey(start=1,end=0) + coord_cartesian(ylim=c(0,20)) + scale_y_continuous(breaks=c(0,5,10,15,20),labels=c("0","5","10","15","")) + labs(y=substitute(paste("Maximum Hydrolysis Rate (nM ",h^-1,")"))) + geom_text(data=ann, label="*") + geom_text(data=low_values,aes(label=label),size=2.2,family="Arial") + geom_point(data=maxSub,aes(x=substrate,y=kcRate.1.nM.hr),alpha=0.4,size=1.3) + geom_point(data=maxSub, aes(x=substrate,y=kcRate.2.nM.hr),alpha=0.4,size=1.3) + geom_point(data=maxSub, aes(x=substrate,y=kcRate.3.nM.hr),alpha=0.4,size=1.3) + theme_barplot
print(bwfig2)
dev.off()


###################Figure 3 - Shannon Diversity###################
shannon <- data.frame()
for (site in levels(maxes$site)) {
    #extract maxes rows in one Stn.Depth group
    stndep <- maxes[maxes$site==site,]
    if (sum(stndep$mean.kcRate.nM.hr,na.rm=TRUE)==0) {
        H <- NA
    } else {
        H <- diversity(x=stndep$mean.kcRate.nM.hr) 
    }
    stndep$H <- rep(H,nrow(stndep))
    shannon <- rbind(stndep[1,c("stn","depthlabel","site","H")],shannon)
    
    #shannon[site,] <- c(stndep[1,c("stn","depthlabel","site")],"H"=H)
}

shannon_theme <- theme_bw() + theme(strip.background=element_rect(fill="grey93"),strip.text=element_text(size=4.5),axis.title=element_text(size=4),axis.title.x=element_text(vjust=0.1),axis.text=element_text(size=3.5),axis.ticks=element_line(size=0.1),panel.grid.major=element_line(color="grey93"),panel.grid.minor=element_blank())
ann <- shannon[complete.cases(shannon)==FALSE,]
tiff("figures/Fig3_ShannonDivBySite.tiff",width=2.5,height=1.5,units="in",res=1200)
fig3 <- ggplot(shannon,aes(x=depthlabel,y=H)) + geom_bar(fill="firebrick",color="black",stat="identity",size=0.1) + facet_grid(~stn) + coord_flip() + scale_y_continuous(breaks=c(0,0.5,1,1.5),labels=c("0.0","","1.0","")) + scale_x_discrete(limits=rev(levels(shannon$depthlabel))) + geom_hline(aes(yintercept=1.79),linetype=2,size=0.1) + geom_text(data=ann,aes(x=depthlabel,y=0.2),label="*",size=2,vjust=0.8) + labs(x="Depth Level",y="Shannon Diversity (H)") + shannon_theme
print(fig3)
dev.off()
#b/w version
#tiff("figures/BWFig3_ShannonDivBySite.tiff",width=2.5,height=1.5,units="in",res=1200)
#fig3 <- ggplot(shannon,aes(x=depthlabel,y=H)) + geom_bar(fill="grey",color="black",stat="identity",size=0.1) + facet_grid(~stn) + coord_flip() + scale_y_continuous(breaks=c(0,0.5,1,1.5),labels=c("0.0","","1.0","")) + scale_x_discrete(limits=rev(levels(shannon$depthlabel))) + geom_hline(aes(yintercept=1.79),linetype=2,size=0.1) + geom_text(data=ann,aes(x=depthlabel,y=0.2),label="*",size=2,vjust=0.8) + labs(x="Depth Level",y="Shannon Diversity (H)") + shannon_theme
#print(fig3)
#dev.off()


#################Figure 4 - Percent Hydrolysis Above Pycnocline##################
CL <- list()
#ID only full surface-bottom casts sampled d5s from stns 2,4,7,10,15,18,21,22,23.
castsd5 <- chem[chem$depthid=="d5",]
#for each cast, read in .asc ctd file from chem-data/ folder
for (cast in as.character(castsd5$cast)) {
    CL[[cast]] <- read.table(paste("chem-data/kn210-04",sub(pattern="cast([0-9]+)",replacement="\\1",cast),".asc",sep=""),header=TRUE)
}
#calculate potential temp, pycnocline (potential density), and buoyancy frequency, append as new column for each cast. (uses oce pkg) 
for (k in 1:length(CL)) {
    CL[[k]]$theta <- swTheta(salinity=CL[[k]]$Sal11,temperature=CL[[k]]$T190C,pressure=CL[[k]]$PrDM)
    CL[[k]]$SigmaTheta <- swSigmaTheta(salinity=CL[[k]]$Sal11,temperature=CL[[k]]$T190C,pressure=CL[[k]]$PrDM)
    CL[[k]]$BuoyFreqN2 <- swN2(pressure=CL[[k]]$PrDM,sigmaTheta=CL[[k]]$SigmaTheta)
}
#get rows of max buoyancy frequency from each cast
maxN2 <- data.frame()
for (j in 1:length(CL)) {
    max <- CL[[j]][CL[[j]]["BuoyFreqN2"]==max(CL[[j]]["BuoyFreqN2"]),]
    max$cast <- names(CL[j])
    max2 <- data.frame(maxn2=max$BuoyFreqN2,cast=max$cast,chl=max$FlECO.AFL)
    #use second row, which is the upcast when samples were actually taken
    maxN2 <- rbind(max2[2,],maxN2)
}
#add second pycnocline max N2 at stn23 where first max is amazon plume
cast071n2 <- CL[["cast071"]][50:1000,]
pyc <- cast071n2[cast071n2[,"BuoyFreqN2"]==max(cast071n2$BuoyFreqN2),]
pyc$cast <- "cast071"
pyc2 <- data.frame(maxn2=pyc$BuoyFreqN2,cast=pyc$cast,chl=pyc$FlECO.AFL)
#pyc2$stn <- "stn23"
maxN2 <- rbind(maxN2,pyc2)
#add stn info to maxN2 from castsd5
maxN2 <- merge(maxN2,castsd5[,c("cast","stn","lat")],by="cast")

#calculate %hydrol in d0 and d1 vs total water column at each station
#total rates at each site
totals <- data.frame("stn"=NA,"shallowrangerate"=NA,"stntotalrate"=NA,"percentinshallow"=NA)
for (stn in levels(maxes$stn)) {
    onestn <- maxes[maxes$stn==stn,]
    #sum d0 and d1 depths
    shallow <- round(sum(onestn[onestn$depthid=="d0" | onestn$depthid== "d1","mean.kcRate.nM.hr"]), 2)
    wholestn <- round(sum(onestn$mean.kcRate.nM.hr,na.rm=TRUE), 2)
    percent <- round(shallow/wholestn, 3)
    insert <- c(stn,shallow,wholestn,percent)
    totals <- rbind(totals,insert)
}
totals <- totals[complete.cases(totals)==TRUE,]

#merge totals with maxN2 at that stn
totals <- merge(x=totals,y=maxN2,by="stn")
totals$percentinshallow <- as.numeric(totals$percentinshallow)
#lm max N2 vs percent hydrolysis occurring above the pycnocline. P<0.05, R2=0.66
n2lm <- lm(percentinshallow~maxn2,data=totals)
R2 <- round(summary(n2lm)$r.squared, 2)
Pvalue <- round(summary(n2lm)$coefficients[2,4], 3)

#plot
fig4_theme <- theme_bw() + theme(legend.position="none",axis.ticks=element_line(size=0.2),axis.title=element_text(size=5),axis.title.x=element_text(vjust=0.1),axis.text=element_text(size=4,color="black"),panel.grid.major=element_line(color="grey93"),panel.grid.minor=element_blank())
tiff("figures/Fig4_StratHydrol.tiff",width=2,height=2,units="in",res=1200)
fig4 <- ggplot(totals,aes(x=maxn2*1000,y=percentinshallow)) + geom_point(size=1.5,aes(color=stn),alpha=0.5) + geom_text(aes(label=stn),size=1.2,family="Arial") + scale_y_continuous(limits=c(0.5,1)) + geom_smooth(method="lm",aes(group=1),color="black",size=0.2,se=FALSE) + scale_color_manual(values=cols.stn) + labs(x=substitute(paste("Maximum Buoyance Frequency (x ",10^3,") ",s^-1)),y="Percent Hydrolysis Above the Pycnocline (%)") + fig4_theme
print(fig4) 
dev.off()
#b/w version
#tiff("figures/BWFig4_StratHydrol.tiff",width=2,height=2,units="in",res=1200)
#fig4 <- ggplot(totals,aes(x=maxn2*1000,y=percentinshallow)) + geom_point(size=1.5,color="grey60",alpha=0.5) + geom_text(aes(label=stn),size=1.2,family="Arial") + scale_y_continuous(limits=c(0.5,1)) + geom_smooth(method="lm",aes(group=1),color="black",size=0.2,se=FALSE) + scale_color_manual(values=cols.stn) + labs(x=substitute(paste("Maximum Buoyance Frequency (x ",10^3,") ",s^-1)),y="Percent Hydrolysis Above the Pycnocline (%)") + fig4_theme
#print(fig4) 
#dev.off()


################Figure 5 - BC distances by stn and depth################
#exclude d3, d4, d5 depth b/c can't calculate BC between two sites that both have zero activity
levs <- levels(maxes$site)[grep(pattern=".d0|.d1|.d2",levels(maxes$site))]
#make matrix of rates for each substrate at each site we're considering
comm <- matrix(nrow=length(levs),ncol=length(levels(maxes$substrate)),dimnames=list(levs,levels(maxes$substrate)))
for (site in levs) {
    #all rows from particular stn.depth site
    rows <- maxes[maxes$site==site,]
    #change Stn.Depth to class character so can insert into comm easily
    #sd$site <- as.character(sd$site)
    for (sub in rows$substrate) {
        #activity from row of sd of particular substrate, insert in comm
        comm[site,sub] <- rows[rows$substrate==sub,"mean.kcRate.nM.hr"]
    }
}
#define categories from rate table in data frame
comm <- na.omit(comm)
comm.group <- data.frame("stn"=sub(pattern="stn([0-9]+).d([0-9])",replacement="stn\\1",row.names(comm)),"depthid"=sub(pattern="stn([0-9]+).d([0-9])",replacement="d\\2",row.names(comm)))
comm.group$depthlabel <- factor(comm.group$depthid,levels=levels(comm.group$depthid),labels=c("SuW","DCM","meso"))

bc <- metaMDS(comm, "bray", 2)
stress <- bc$stress
#stress pretty good at 0.1
print(paste("Bray Curtis solution reached. Stress:",stress))
#are bc distances significantly different comparing across stations
#YES, P=0.004
permanova.stn <- adonis(comm~stn,data=comm.group)
stn.p <- permanova.stn$aov.tab[1,6]
print(paste("permanova distances by station:",stn.p))
#are bc distances significantly different comparing across depth levels
#NO, P=0.40
permanova.dep <- adonis(comm~depthlabel,data=comm.group)
dep.p <- permanova.dep$aov.tab[1,6]
print(paste("permanova distances by depth:",dep.p))

#5a
tiff("figures/Fig5a_NMDSbystn.tiff",width=5,height=5,units="in",res=1200)
#par(mfrow=c(2,1))
#plot by stn 
#make empty plot
plot(bc, type="n")
#add points 
points(bc, pch = c(16,18,0,15,17,3)[as.numeric(comm.group$stn)], col=cols.stn[comm.group$stn])
#add ellipses around single stn group colored by stn factor
for (i in 1:length(levels(comm.group$stn))) {
    ordiellipse(bc,groups=comm.group$stn,label=TRUE,col=cols.stn[i],show.groups=levels(comm.group$stn)[i],draw="polygon",alpha=40)
}
#mtext("a)",side=3,adj=0,line=2,font=2)
#add P-value
#text(c(-0.75,1),labels=paste("P=",permanova.stn$aov.tab$Pr[1]))
dev.off()

##Fig 5B
##plot by depth label
tiff("figures/Fig5b_NMDSbydep.tiff",width=5,height=5,units="in",res=1200)
plot(bc, type="n")
points(bc, pch = c(16,17,3)[as.numeric(comm.group$depthlabel)], col=cols.dep[comm.group$depthlabel])
for (i in 1:length(levels(comm.group$depthlabel))) {
    ordiellipse(bc,groups=comm.group$depthlabel,label=TRUE,col=cols.dep[i],show.groups=levels(comm.group$depthlabel)[i],draw="polygon",alpha=40)
}
dev.off()

#b/w version
#5a
#tiff("figures/BWFig5a_NMDSbystn.tiff",width=5,height=5,units="in",res=1200)
#plot(bc, type="n")
#points(bc, pch = c(16,18,0,15,17,3)[as.numeric(comm.group$stn)], col="black")
#add ellipses around single stn group colored by stn factor
#for (i in 1:length(levels(comm.group$stn))) {
#    ordiellipse(bc,groups=comm.group$stn,label=TRUE,col="grey93",show.groups=levels(comm.group$stn)[i],draw="polygon",alpha=50)
#}
#dev.off()
#5b
#tiff("figures/BWFig5b_NMDSbydep.tiff",width=5,height=5,units="in",res=1200)
#plot(bc, type="n")
#points(bc, pch = c(16,17,3)[as.numeric(comm.group$depthlabel)], col="black")
#for (i in 1:length(levels(comm.group$depthlabel))) {
#    ordiellipse(bc,groups=comm.group$depthlabel,label=TRUE,col="grey93",show.groups=levels(comm.group$depthlabel)[i],draw="polygon",alpha=50)
#}
#dev.off()


#################Figure 6 - Short Substrate Barplot####################
ann_mon <- monmaxes[monmaxes$stn=="stn7"&monmaxes$depthid=="d2"&monmaxes$substrate=="G",]
mon_lowvalues <- data.frame(stn=monmaxes_lowdeepvalues$stn,depthlabel=monmaxes_lowdeepvalues$depthlabel,substrate=monmaxes_lowdeepvalues$substrate,mean.kcRate.nM.hr=rep(6,nrow(monmaxes_lowdeepvalues)),label=c(0.02, 0.4, 0.1, 0.02, 0.3, 0.1, 0.2, 0.2, 0.2, 0.2, 0.1, 0.2, 0.2, 0.1, 0.02, 0.2, 0.1, 0.04, 0.4, 0.1, 0.2, 0.3, 0.2, 0.1, 0.3, 0.2, 0.3, 0.2))
mon_errorbars <- geom_errorbar(aes(ymin=mean.kcRate.nM.hr-sd.kcRate.nM.hr,ymax=mean.kcRate.nM.hr+sd.kcRate.nM.hr),width=0.5,color="grey25",alpha=0.6,size=0.3) 
theme_monbarplot <- theme_bw() + theme(axis.title=element_text(size=10),axis.text.x=element_text(size=6),axis.text.y=element_text(size=7.5),strip.background=element_rect(fill="grey93"),axis.title.x=element_text(vjust=0.5),axis.title.y=element_text(vjust=0.5),legend.position="none",panel.grid.major=element_line(color="grey93"),panel.grid.minor=element_blank()) 

tiff("figures/Fig6_ShortSubBarplot.tiff",width=5,height=4,units="in",res=1200)
fig6 <- ggplot(monmaxes,aes(x=substrate,y=mean.kcRate.nM.hr)) + geom_bar(stat="identity",aes(fill=substrate),color="black",size=0.1) + scale_y_continuous(breaks=c(0,5,10,15,20),labels=c("0","5","10","15","")) + scale_fill_manual(values=MonColors) + facet_grid(depthlabel~stn) + coord_flip(ylim=c(0,20)) + mon_errorbars + ylab(substitute(paste("Maximum Hydrolysis Rate (nM ",h^-1,")"))) + geom_text(data=ann_mon,y=15,label="43.0",size=2,vjust=-0.2) + geom_text(data=mon_lowvalues,aes(label=label),size=2,family="Arial") + theme_monbarplot
print(fig6)
dev.off()
#B/W version
tiff("figures/BWFig6_ShortSubBarplot.tiff",width=5,height=4,units="in",res=1200)
bwfig6 <- ggplot(monmaxes,aes(x=substrate,y=mean.kcRate.nM.hr)) + geom_bar(stat="identity",aes(fill=substrate),color="black",size=0.1) + scale_y_continuous(breaks=c(0,5,10,15,20),labels=c("0","5","10","15","")) + scale_fill_grey(start=0.7,end=0.3) + facet_grid(depthlabel~stn) + coord_flip(ylim=c(0,20)) + mon_errorbars + ylab(substitute(paste("Maximum Hydrolysis Rate (nM ",h^-1,")"))) + geom_text(data=ann_mon,y=15,label="43.0",size=2,vjust=-0.2) + geom_text(data=mon_lowvalues,aes(label=label),size=2,family="Arial") + theme_monbarplot
print(bwfig6)
dev.off()
