#FurtherProcessing_DeepDOM.R

#take final rates, add some factor labels and save just kill-corrected mean and sd rate as new table FlaRatesWithFactors_DeepDOM.csv
#find timepoint with max rate for each incubation set from FlaRatesWithFactors, save as FlaMaxRates_DeepDOM.csv
#take short substrate results master sheet, find timepoint with max rate and save as MonMaxRates_DeepDOM.csv

master <- read.csv("FLAMasterRates_DeepDOM.csv",row.names=1)
time <- read.csv("FLAElapsedTime_DeepDOM.csv",row.names=1)
mon <- read.csv("MonResultsMaster_DeepDOM.csv",header=TRUE)

#add stn, latitude, longitude, depth id, site, depth sampled (m), substrate, timepoint, elapsed time factor columns

factors <- data.frame(kcRate.1.nM.hr=master$kcRate.1.nM.hr, kcRate.2.nM.hr=master$kcRate.2.nM.hr, kcRate.3.nM.hr=master$kcRate.3.nM.hr, mean.kcRate.nM.hr=master$mean.kcRate.nM.hr, sd.kcRate.nM.hr=master$sd.kcRate.nM.hr,row.names=row.names(master))
#stn
factors$stn <- factor(sub(pattern="stn([0-9]+)-d([0-9])-([a-z]+)-t([0-9])",replacement="stn\\1",row.names(factors)),levels=c("stn2","stn4","stn7","stn10","stn15","stn18"))
#depthid
factors$depthid <- factor(sub(pattern="stn([0-9]+)-d([0-9])-([a-z]+)-t([0-9])",replacement="d\\2",row.names(factors)))
#depthlabel
factors$depthlabel <- factor(factors$depthid,levels=levels(factors$depthid),labels=c("SuW","DCM","meso","AAIW","NADW","bot"))
#site (stn.depth id)
factors$site <- factor(sub(pattern="stn([0-9]+)-d([0-9])-([a-z]+)-t([0-9])",replacement="stn\\1.d\\2",row.names(factors)))

#lat, long, and depth (m) retrieved from ship-board ctd log sheets:
coordinates <- data.frame(stn=rep(levels(factors$stn),each=6),depthid=rep(levels(factors$depthid),6),
                          lat=c(rep(-38,6),rep(-31.25,6),rep(-22.5,6),rep(-9.5,6),rep(-2.7,6),rep(3,6)),
                          long=c(rep(-45,6),rep(-41,6),rep(-32.7,6),rep(-26,6),rep(-28.5,6),rep(-39,6)),
                          depthsampled=c(5,80,220,750,2502,5109,9,121,253,850,2503,3781,6,131,251,753,2505,4513,5,126,250,849,2500,5225,5,63,251,750,2500,5009,6,65,251,760,2500,4511))
#stn2 was at 38S, 45W, sampled d0=5m,d1=80m,d2=220m,d3=750m,d4=2502m,d5=5109m
#stn4 was at 31.25S,41W, d0=9m,d1=121m,d2=253m,d3=850m,d4=2503m,d5=3781m
#stn7 was at 22.5S,32.7W, d0=6m,d1=131m,d2=251m,d3=753m,d4=2505m,d5=4513m
#stn10 was at 9.5S,26W, d0=5m,d1=126m, d2=250m,d3=849m,d4=2500m,d5=5225m
#stn15 was at 2.7S,28.5W, d0=5m,d1=63m,d2=251m,d3=750m,d4=2500m,d5=5009m
#stn18 was at 3N,39W, d0=6m,d1=65m,d2=251m,d3=760m,d4=2500m,d5=4511m

#latitude
for (stn in levels(factors$stn)) {
    factors[factors$stn==stn,"lat"] <- coordinates[coordinates$stn==stn,"lat"]
}
factors$lat <- as.factor(factors$lat)
#longitude
for (stn in levels(factors$stn)) {
    factors[factors$stn==stn,"long"] <- coordinates[coordinates$stn==stn,"long"]
}
factors$long <- factor(factors$long, levels=c("-45","-41","-32.7","-26","-28.5","-39"))
#depthsampled
for(stn in levels(factors$stn)) {
    for (dep in levels(factors$depthid)) {
        factors[factors$stn==stn&factors$depthid==dep,"depthsampled"] <- coordinates[coordinates$stn==stn&coordinates$depthid==dep,"depthsampled"]
    }
}
#substrate
factors$substrate <- factor(sub(pattern="stn([0-9]+)-d([0-9])-([a-z]+)-t([0-9])",replacement="\\3",row.names(factors)))
#timepoint
factors$timepoint <- factor(sub(pattern="stn([0-9]+)-d([0-9])-([a-z]+)-t([0-9])",replacement="t\\4",row.names(factors)))
#elapsed time
row <- gsub(pattern="stn([0-9]+)-d([0-9])-([a-z]+)-t([0-9])",replacement="stn\\1-d\\2-\\3-1-t\\4",x=row.names(master))
factors$elapsedtime <- time[row,"elapsed.time.hrs"]

#save finished factors data frame as FlaRatesWithFactors.csv
write.csv(factors,"FlaRatesWithFactors_DeepDOM.csv",row.names=TRUE)


#####FlaMaxRates_DeepDOM.csv
#take factors, extract only timepoint with maximum hydrolysis rate for each incubation set; save as FlaMaxRates_DeepDOM.csv
FlaMaxes <- matrix(ncol=ncol(factors),dimnames=list(NULL,colnames(factors)))

#loop through each stn, then depth, then substrate, find max activity timepoint for that substrate, and rbind to FLAmax matrix
#for one stn...
for (stn in levels(factors$stn)) {
    #and one depth...
    for (dep in levels(factors$depthid)) {
        #and one substrate...
        for (sub in levels(factors$substrate)) {
            subset <- factors[factors$stn==stn&factors$depthid==dep&factors$substrate==sub,]
            #...find row with max activity
            m <- max(subset$mean.kcRate.nM.hr)
            maxsub <- subset[subset$mean.kcRate.nM.hr==m,]
            #if >1 row is max (e.g. rates are 0 at every timepoint), use first row
            if(nrow(maxsub)>1) {
                maxsub <- maxsub[1,]
            }
            
            #Add max to FLAmax
            FlaMaxes <- rbind(FlaMaxes,maxsub)
        }
    }
}

#remove row with NA
FlaMaxes <- FlaMaxes[complete.cases(FlaMaxes)==TRUE,]
#save as FlaMaxes_DeepDOM.csv
write.csv(FlaMaxes,"FlaMaxRates_DeepDOM.csv",row.names=TRUE)

######short substrate results maxes
#for each stn/depth/substrate, find max, save in monmaxes df
monmaxes <- data.frame()
for (stn in levels(mon$stn)) {
    #and one depth...
    for (dep in levels(mon$depthid)) {
        #and one substrate...
        for (sub in levels(mon$substrate)) {
            subset <- mon[mon$stn==stn&mon$depthid==dep&mon$substrate==sub,]
            #...find row with max activity
            m <- max(subset$mean.kcRate.nM.hr)
            maxsub <- subset[subset$mean.kcRate.nM.hr==m,]
            #if >1 row is max (e.g. rates are 0 at every timepoint), use first row
            if(nrow(maxsub)>1) {
                maxsub <- maxsub[1,]
            }
            
            #Add max to FLAmax
            monmaxes <- rbind(maxsub,monmaxes)
        }
    }
}

write.csv(monmaxes,"MonMaxRates_DeepDOM.csv",row.names=FALSE)
