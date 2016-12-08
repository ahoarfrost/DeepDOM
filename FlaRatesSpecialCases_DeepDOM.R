##Calculates rates for special cases
#ignore FLA - where ignore FLA bin b/c too much degraded substrate
#change to zero - where looked at chroms and identified no hydrolysis, but due to baseline bumpiness might result in small non-zero rate
#output new, adjusted master rates as FlaMasterListFinal_DeepDOM

path <- getwd()

#name want output file to be
outputMaster <- "FLAMasterRates_DeepDOM.csv"
#read in master list of rates
master <- read.csv("FLAMasterRates_DeepDOM.csv",header=TRUE,row.names=1)

#name folder in your wd in which have slant corrected data want to calc rates with
CsvDir <- "csvs-for-rates" 
#define substrates used in this experiment (using abbreviations used in file names)
substrates <- c("ara","chon","fuc","lam","pul","xyl")
#Read in info about mw of substrates and the # of cuts needed to get to a certain std size bin
cuts <- read.csv("HydrolysisCutsInfo_DeepDOM.csv", skip=1, header=TRUE, row.names=1)
#Read in table of sampling times and the Elapsed Time since t0
FLAElapsed <- read.csv("FLAElapsedTime_DeepDOM.csv", header=TRUE, row.names=1)
#Read in info about which std bins to use for which samples
StdsForRates <- read.csv("StdsForRates_DeepDOM.csv",header=TRUE,row.names=1)

####change to zero correction#####
#subset only the rows where need to change rates to zero, change mean.kcRate.nM.hr and sd.kcRate.nM.hr values to zero
zero <- master[which(master$addend=="change to zero"),]
zero$mean.kcRate.nM.hr <- 0
zero$sd.kcRate.nM.hr <- 0
#insert new zero values into master data frame
master[row.names(zero),] <- zero
##################

####ignore FLA correction####
#have to recalculate these rates, but ignore FLA bin. So, percent of hydrolysis within a bin is as divided by sum of just 150kD, 10kD, 4kD, and monomer bins only
igFLA <- master[which(master$addend=="ignore FLA"),]
#define partialname of sets of incubations to recalculate rates for
partial <- unique(sub(pattern="stn([0-9]+)-d([0-9])-([a-z]+)-t([0-9])",replacement="stn\\1-d\\2-\\3",row.names(igFLA)))
#make a list of all the csv file names containing tagged with 'ignore FLA' in csvs-for-rates folder
pattern <- paste(partial,collapse="|")
CsvNameList <- list.files(path=CsvDir,pattern=pattern)

#Make parent lists that divide csvs into useful groups for later 
RatesMasterList <- list()
IncList <- list()
TimeList <- list()

setwd(paste(path,CsvDir,sep="/"))

for (i in 1:length(CsvNameList)) {
    CsvName <- CsvNameList[i]
    #takes each csv file name you read into CsvNameList, and define the PartialName, IncName, TimeName, and FullName for that file name
    PartialName <- sub("stn([0-9]+)-d([0-5])-([a-z]+)-([0-9a-z]{1,4})-t([0-3]).csv","stn\\1-d\\2-\\3",CsvName)
    IncName <- sub("stn([0-9]+)-d([0-5])-([a-z]+)-([0-9a-z]{1,4})-t([0-3]).csv","stn\\1-d\\2-\\3-\\4",CsvName)
    FullName <- sub("stn([0-9]+)-d([0-5])-([a-z]+)-([0-9a-z]{1,4})-t([0-3]).csv","stn\\1-d\\2-\\3-\\4-t\\5",CsvName)
    TimeName <- sub("stn([0-9]+)-d([0-5])-([a-z]+)-([0-9a-z]{1,4})-t([0-3]).csv","stn\\1-d\\2-\\3-t\\5",CsvName)
    #sort these references into useful lists for later
    #RatesMasterList, each object in list is a list named PartialName ('stn-depth-substrate'). Within each PartialName list, each object in the list is a data frame named FullName ('stn-depth-sub-incubation-timepoint') with two columns of time and slant corrected reads
    RatesMasterList[[PartialName]][[FullName]] <- read.csv(CsvName, header=TRUE)
    #IncList, goes PartialName -> Inc Name (grouped by incubation, e.g. stn2-d3-lam-x) -> FullName (each group contains FullName csv for all timepoints within that incubation)
    IncList[[PartialName]][[IncName]][[FullName]] <- read.csv(CsvName,header=TRUE)
    #TimeList, PartialName -> TimeName (grouped by timpoint, e.g. stn4-d0-ara-t0 all the t0s together, t1s together...tn) -> FullName
    TimeList[[PartialName]][[TimeName]][[FullName]] <- read.csv(CsvName,header=TRUE)
}

setwd(path)

#this will just print out a summary of your lists so you can doublecheck it looked all right
print("RatesMasterList:")
print(summary(RatesMasterList))
print("IncList:")
print(summary(IncList))
print("TimeList:")
print(summary(TimeList))

for (j in 1:length(RatesMasterList)) {
    print(paste("Processing",names(RatesMasterList[j]),"...",sep=" "))
    #Detect which std bins to use for this incubation set
    stdsID <- as.character(StdsForRates[names(RatesMasterList[[j]][1]),])
    #Read in info about cutoffs for std bins from csv with correct stds
    StdBins <- read.csv(paste("stdbins/",stdsID,".csv",sep=""),header=TRUE)
    #Define std bins from StdBins.csv
    kD150 <- StdBins$read.number[1] #Is the read number (e.g. row 145) that bin starts; to the right of this read everything in the '150kD' bin
    kD10 <- StdBins$read.number[2] #to the right of this in '10kD' bin
    kD4 <- StdBins$read.number[3] #to the right of this in '4kD' bin
    kDmon <- StdBins$read.number[4] #to the right of this in 'monomer' bin
    kDfla <- StdBins$read.number[5] #to the right of this in 'FLA' bin                                                2
    
    #Make matrix with 6 rows and correct number columns (depending on number timepoints) named FullName
    #IGNORING FLA BIN - don't include an FLA row or FLA in total
    BinSums <- matrix(nrow=5, ncol=(length(names(RatesMasterList[[j]]))+1), dimnames=list(c("150kD", "10kD", "4kD", "monomer", "total (to fla)"), c("BinReadStart", names(RatesMasterList[[j]]))))
  
    #Add bin starts from standards to "BinReadStart" column in matrix
    #IG FLA BIN - "total" is now at the point where FLA bin starts
    BinSums[,"BinReadStart"] <- c(kD150,kD10,kD4,kDmon,kDfla)
    
    #For each FullName dataframe in this list, sum up fluorescence reads within the bins and insert into BinSums matrix in corresponding row
    #IG FLA BIN - sum within each bin, "total" sum is just up to FLA bin read start
    for (k in 1:length(RatesMasterList[[j]])) {
        reads <- RatesMasterList[[j]][[k]]$SlantCorrectedReads
        #For each FullName dataframe in this list, sum up fluorescence reads within the bins and insert into BinSums matrix in corresponding row
        BinSums[,k+1] <- c(sum(reads[kD150:(kD10-1)]),sum(reads[kD10:(kD4-1)]),sum(reads[kD4:(kDmon-1)]),sum(reads[kDmon:(kDfla-1)]), sum(reads[kD150:(kDfla-1)]))
    }
    
    #Make new matrix, in which will calculate percent of total fluorescence is in each bin
    PercentTotal <- matrix(nrow=4,ncol=ncol(BinSums)-1, dimnames=list(c("%150kD", "%10kD", "%4kD", "%monomer"), c(colnames(BinSums[,2:ncol(BinSums)]))))
    
    #Calculate percent of total fluorescence from BinSums and insert into PercentTotal matrix in appropriate cell
    for (n in 1:ncol(PercentTotal)) {
        PercentTotal[,n] <- c(BinSums[1,n+1]/BinSums[5,n+1], BinSums[2,n+1]/BinSums[5,n+1], BinSums[3,n+1]/BinSums[5,n+1], BinSums[4,n+1]/BinSums[5,n+1])
    }
    
    #Make new matrix in which will calculate change from time zero for each sample
    ChangeZero <- matrix(nrow=4,ncol=ncol(PercentTotal), dimnames=list(c("150kD-0","10kD-0","4kD-0","monomer-0"), c(colnames(PercentTotal))))
    
    #remove columns from PercentTotal by IncName (e.g. all columns that start with stn10-d2-chon-1), calculate ChangeZero, insert into correct column in ChangeZero
    for (p in 1:length(IncList[[j]])) {
        #makes new matrix with one incubation set, each column is T0 T1, T2...Tn
        PT <- PercentTotal[,grep(paste(names(IncList[[j]][[p]]),collapse="|"),x=colnames(PercentTotal))]
        #calculate change zero for that incubation, subtract column 1 from all other columns
        CZ <- PT-PT[,1]
        #insert CZ into correct columns in ChangeZero matrix
        ChangeZero[,colnames(CZ)] <- CZ
    }
    
    #Correct any negative values to be zero
    for (t in 1:length(ChangeZero)) { 
        if (ChangeZero[t]<0) {
            ChangeZero[t]=0
        } 
    }
    
    #Calculate number of hydrolysis events in 10kD-FLA bins (excluding 150kD) for each sample (column), using correct hydrolysis cuts and nmol substrate info from cuts.csv read in above; record in new matrix hydrolEventPerBin
    for (aa in substrates) {
        if (grepl(paste("*",aa,sep=""),names(RatesMasterList[j]))==TRUE) {
            u <- cuts[4:6,aa]
            umat <- matrix(rep(u), ncol=ncol(ChangeZero), nrow=3, dimnames=list(c("# cuts 10kD","# cuts 4kD","# cuts monomer"), c(rep(aa, ncol(ChangeZero)))))
            nmolSub <- cuts["nmol substrate",aa]
            hydrolEventPerBin <- matrix(ChangeZero[2:4,]*umat*nmolSub, nrow=3, ncol=ncol(ChangeZero), dimnames=list(c("10kD events","4kD events","monomer events"), c(colnames(ChangeZero))))
        }
    }
    
    #Sum hydrolysis events from each bin (sum all rows in each column), new matrix with each column FullName (name of incubation and timepoint) and 1 row of total hydrolysis events 
    hydrolysisEvents <- matrix(data=colSums(hydrolEventPerBin), nrow=1, ncol=ncol(hydrolEventPerBin), dimnames=list(c("#cuts"),c(colnames(hydrolEventPerBin))))
    
    #make matrix for calculating rates
    Rates <- matrix(nrow=1,ncol=ncol(hydrolysisEvents),dimnames=list(c("rate #/nM.hr"),c(colnames(hydrolysisEvents))))
    
    #calculate rates from hydrolysis events - retrieve elapsed time from spreadsheet, divide hydrolysisEvents/(timeElapsed*volume incubation)
    for (v in 1:length(colnames(hydrolysisEvents))) {
        timeElapsed <- FLAElapsed[colnames(hydrolysisEvents)[v],"elapsed.time.hrs"]
        volumeL <- 0.0167
        Rates[,v] <- hydrolysisEvents[,v]/(timeElapsed*volumeL)
    }
    
    #if there are any NaN    or Inf values in Rates matrix (e.g. at t0 because dividing by 0), change value to 0
    for (w in 1:length(Rates)) { 
        if (is.nan(Rates[w])==TRUE) {
            Rates[w]=0
        } 
        if (is.infinite(Rates[w])==TRUE) {
            Rates[w]=0
        } 
    }
    
    #group Rates by timepoint using TimeList, then kill correct and record pre-kc and kc rates in FLAMaster
    for (x in 1:length(TimeList[[j]])) {
        #makes new matrix with one timepoint set, each column is rep1, rep2,...repn, x
        R2 <- Rates[,grep(paste(names(TimeList[[j]][[x]]),collapse="|"),x=colnames(Rates))]
        #if there are multiple x incubations, average rate, make one column for x. If there's only one x incubation, move that value to the first value to make subtracting kill control easier
        y <- R2[grep("-[x]-",x=names(R2),value=TRUE)]
        if (length(y)>1) {
            #mean of all the kills
            ymean <- mean(y)
            #identify position in R2 where the kills are
            yn <- which(names(R2) %in% names(y))
            #remove all the kill columns
            R2 <- R2[-yn]
            #Add the new kill column with the avg ymean in the first column
            R2 <- append(R2,values=c(ymean),after=0)
            #Name the appended column
            names(R2)[1] <- names(y)[1]
        } else {
            #Remove the one kill column, append it to front of vector. 
            yn <- which(names(R2) %in% names(y))
            R2 <- R2[-yn]
            R2 <- append(R2,values=c(y),after=0)
        }
        
        #calculate kill corrected rates by subtracting x rate from all other rates
        kcR2 <- R2 - R2[1]
        
        #If any kill corrected rates <0, turn to zero
        for (z in 1:length(kcR2)) { 
            if (kcR2[z]<0) {
                kcR2[z]=0
            } 
        }
        
        #Add mean and sd values to rates (R2) and kill corrected rates (kcR2)
        R2 <- append(R2,values=c("mean.Rate.nM.hr"=mean(R2[2:length(R2)]),"sd.Rate.nM.hr"=sd(R2[2:length(R2)])),after=length(R2))
        kcR2 <- append(kcR2,values=c("mean.kcRate.nM.hr"=mean(kcR2[2:length(kcR2)]),"sd.kcRate.nM.hr"=sd(kcR2[2:length(kcR2)])),after=length(kcR2))
        
        #rename columns in R2 and kcR2 to match what want in FLAMaster. takes just the replicate info, e.g. "stn4-d0-lam-2-t2" -> "2 rate nM.hr"
        names(R2) <- gsub(pattern="stn([0-9]+)-d([0-5])-([a-z]+)-([0-9a-z])-t([0-9])","Rate.\\4.nM.hr",names(R2))
        names(kcR2) <- gsub(pattern="stn([0-9]+)-d([0-5])-([a-z]+)-([0-9a-z])-t([0-9])","kcRate.\\4.nM.hr",names(kcR2))
        
        #append kcR2 to R2, so can insert into FLAMaster 
        finalR <- append(R2,kcR2,after=length(R2))
        
        #Write rates before kill corrected (R2) and kill corrected rates (kcR2 to correct row and columns in FLA Master Sheet)
        master[names(TimeList[[j]][x]),names(finalR)] <- finalR
        
        print(paste(names(TimeList[[j]][x]),"completed.",sep=" "))
        
    }
    
    #save the new FLAMaster with your calculated rates to a csv file in your wd named as defined by outputMaster
    write.csv(master,outputMaster,quote=FALSE)
    
}
