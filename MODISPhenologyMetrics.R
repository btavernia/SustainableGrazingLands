#Author: Brian G. Tavernia
#Date: 05.17.2015
#Purpose: This script is designed to accomplish the following:

#(1) to smoothed NDVI values using Savitzky-Golay filter approach outlined in Chen et al. (2004) and 
#(2) to calculate phenological metrics based on smoothed NDVI values to judge the effect of different ranch management practices
#on ecosystem productivity

#Resource: smoothing will be accomplished using the methods outlined in

#Chen et al. 2004. A simple method for reconstructing a high-quality NDVI time-series data set based on Savitzky-Golay
#filter. Remote Sensing of Environment. 91:332-344.  

#####
#Environment
#####
options(scipen=999) #setting precision of scientific notation

#####
#Libraries
#####
require(signal)
require(RODBC)

#####
#Input parameters
#####
indata <- "PtsFox_pModis.txt" #name of input MODIS data file
syr <- 2001 #first year
fyr <- 2015 #final year
nyr <- 15 #number of years
pref <- "StudyAreaThreshold2001to2015_20Percent_" #prefix for output pdf and csv files
ndvit <- 2510 #NDVI Threshold for calculating phenology metrics

#creating vector containing colors used to plot smoothed ndvi curves for 16 years
colStore = c("chartreuse","aquamarine3","chocolate1","blue1","darkgoldenrod",
             "brown","gold","darkseagreen1","darkviolet","greenyellow","lightcyan",
             "navy","orangered","tan","springgreen4","snow4")

#####
#Input Data
#####
#MODIS data for fox
rdata = read.csv(file=paste("W:\\TNC_Land_Protection\\2015_SustainableGrazinglands\\Data\\ModisPilots\\",indata,sep=""),stringsAsFactors=F) 

#creating vector containing unique ranch names; will loop through ranches
rnames = unique(rdata$RANCH_1)

#####
#Creating matrices to store metrics
#####
sgseas <- matrix(nrow=length(rnames),ncol=nyr) #start of growing season
colnames(sgseas) <- seq(syr,fyr,1)
rownames(sgseas) <- rnames

egseas <- matrix(nrow=length(rnames),ncol=nyr) #end of growing season
colnames(egseas) <- seq(syr,fyr,1)
rownames(egseas) <- rnames

lgseas <- matrix(nrow=length(rnames),ncol=nyr) #length of growing season
colnames(lgseas) <- seq(syr,fyr,1)
rownames(lgseas) <- rnames

auc <- matrix(nrow=length(rnames),ncol=nyr) #area under the curve for growing season
colnames(auc) <- seq(syr,fyr,1)
rownames(auc) <- rnames

mean.ndvi <- matrix(nrow=length(rnames),ncol=nyr) #mean ndvi during growing season
colnames(mean.ndvi) <- seq(syr,fyr,1)
rownames(mean.ndvi) <- rnames

median.ndvi <- matrix(nrow=length(rnames),ncol=nyr) #median ndvi during growing season
colnames(median.ndvi) <- seq(syr,fyr,1)
rownames(median.ndvi) <- rnames
 
min.ndvi <- matrix(nrow=length(rnames),ncol=nyr) #minimum ndvi during growing season
colnames(min.ndvi) <- seq(syr,fyr,1)
rownames(min.ndvi) <- rnames

min.ndvi.date <- matrix(nrow=length(rnames),ncol=nyr) #date of minimum ndvi during growing season
colnames(min.ndvi.date) <- seq(syr,fyr,1)
rownames(min.ndvi.date) <- rnames

max.ndvi <- matrix(nrow=length(rnames),ncol=nyr) #maximum ndvi during growing season
colnames(max.ndvi) <- seq(syr,fyr,1)
rownames(max.ndvi) <- rnames

max.ndvi.date <- matrix(nrow=length(rnames),ncol=nyr) #date of maximum ndvi during growing season
colnames(max.ndvi.date) <- seq(syr,fyr,1)
rownames(max.ndvi.date) <- rnames

#####
#Looping through ranches
#####

for (ranch in rnames){
  
  print(ranch)
  
  #keeping only records pertaining to a particular ranch
  rdata.sub = subset(rdata,RANCH_1==ranch)
  
  #####
  #Data Manipulation
  #####

  rdata.sub = rdata.sub[,-which(colnames(rdata.sub)=="RANCH_1")]
  
  #using for loop to convert character to numeric
  for (i in 1:(ncol(rdata.sub))){
    
    rdata.sub[,i]=as.numeric(gsub(",","",rdata.sub[,i]))
    
  }
  
  #creating a vector that contains column names
  cnames = colnames(rdata.sub)
  
  #creating vector containing years to be produced
  years = seq(syr,fyr,1)
  
  for (year in years){
    
    print(year)
    
    #####
    #Filtering NDVI values based on atmospheric conditions, soil, and NLCD cover classes
    #####
    
    #getting column numbers corresponding to ndvi for a particular year and for 5 extra days at the beginning and end of each year
    #extra days are to avoid edge effects when smoothing the curve
    cnumbers = grep(paste("^A",year,sep=""),cnames)
    cnumbers = c(seq((cnumbers[1]-5),(cnumbers[1]-1)),cnumbers,seq((cnumbers[length(cnumbers)]+1),(cnumbers[length(cnumbers+1)]+5)))
    
    #getting ndvi columns associated with specific year
    rdata.ndvi = rdata.sub[,cnames[cnumbers]]
    
    #getting column numbers corresponding to cloud data for a particular year and for 4 extra days at the beginning and end of each 
    #year extra days are to avoid edge effects when smoothing the curve
    cnumbers = grep(paste("^PR_A",year,sep=""),cnames)
    cnumbers = c(seq((cnumbers[1]-5),(cnumbers[1]-1)),cnumbers,seq((cnumbers[length(cnumbers)]+1),(cnumbers[length(cnumbers+1)]+5)))
    
    #getting cloud status columns associated with specific year
    rdata.cloud = rdata.sub[,cnames[cnumbers]]
    
    #binding ndvi and cloud data
    rdata.ndvi = cbind(rdata.ndvi,rdata.cloud)
    
    #filtering out pixels flagged as having snow or cloud cover (values>=2)
    for (i in 1:(ncol(rdata.ndvi)/2)){
      
      rdata.ndvi[,i] = ifelse(rdata.ndvi[,(i+ncol(rdata.ndvi)/2)]<2,rdata.ndvi[,i],NA)
      
    }
    
    #getting columns associated with specific year and extra days
    rdata.ndvi = rdata.ndvi[,colnames(rdata.ndvi)[grep("^A",colnames(rdata.ndvi))]]
    
    #getting columns associated with soils data
    soils = rdata.sub[,colnames(rdata.sub)[grep("^SGVALUE",colnames(rdata.sub))]]
    
    #determining dominant three soils based on percentages
    soilsn = names(sort((colSums(soils)/sum(colSums(soils)))*100,decreasing=T))[1:3] 
    
    #determining if dominant soils make up >50% of soil within a pixel
    check.soils = (rowSums(soils[,soilsn])/rowSums(soils))*100
    
    #replacing row values with NA for those pixels where dominant soils do not make >=50% of pixel
    for (i in 1:nrow(rdata.ndvi)){
      
      for (j in 1:ncol(rdata.ndvi)){
      
      rdata.ndvi[i,j] = ifelse(check.soils[i]>=80,rdata.ndvi[i,j],NA)
      
      }
      
    }
    
    #if open water or developed cover classes are present, will filter out pixels with either open water 
    #or developed cover present
    if (length(c(colnames(rdata.sub)[grep("^NLCDVALUE_2",colnames(rdata.sub))],colnames(rdata.sub)[grep("^NLCDVALUE_11",colnames(rdata.sub))]))>0){
      
      #getting columns associated with open water or developed areas
      nlcd = cbind(rdata.sub[,colnames(rdata.sub)[grep("^NLCDVALUE_2",colnames(rdata.sub))]],rdata.sub[,colnames(rdata.sub)[grep("^NLCDVALUE_11",colnames(rdata.sub))]])
      
      #summing across cover types
      nlcdrsums = rowSums(nlcd)
      
      #replacing pixel values with NA if open water or developed cover is present
      for (i in 1:nrow(rdata.ndvi)){
        
        for (j in 1:ncol(rdata.ndvi)){
        
        rdata.ndvi[i,j] = ifelse(nlcdrsums[i]==0,rdata.ndvi[i,j],NA)
        
        }
        
      }
      
      
    }
    
    #####
    #Creating data set containing mean daily NDVI values across pixels
    #####
    
    #determining mean NDVI value across pixels for each day
    ndvim = colMeans(rdata.ndvi,na.rm=T)
    
    #getting day numbers for NDVI plot
    dlab = c()
    
    for (i in 1:length(ndvim)){
      
      dlab = c(dlab,substr(names(ndvim)[i],6,8))
      
    }
    
    dlab = as.numeric(dlab)
    
    #creating data frame containing day and mean NDVI values
    ndvim = as.data.frame(cbind(dlab,ndvim))
    colnames(ndvim) = c("Day","No")
    
    #####
    #Chen et al. (2004) Step 1. Interpolating missing values in the NDVI data series
    #####
    
    #handling any missing means using average of adjacent points
    if (length(which(is.na(ndvim$No)))>0){
      
      nap = which(is.na(ndvim$No))
      nap = nap[which(nap!=1)]
      nap = nap[which(nap!=length(ndvim$No))]
      
      nnap = which(!is.na(ndvim$No))
      
      #using mean of nearest adjacent points to replace missing value
      for (v in 1:length(nap)){
        
        
        diff = nap[v]-nnap
        ndvim$No[nap[v]] = mean(c(ndvim$No[nnap[min(which(diff<0))]],ndvim$No[nnap[max(which(diff>0))]]),na.rm=T)
        
      }   
      
    }
    
    #removing first and last rows from data frame - do not need these dates any longer as they will be outside of
    #smoothing window
    ndvim = ndvim[-c(1,nrow(ndvim)),] 
    
    #####
    #using Savitzky-Golay filter to smooth NDVI data series
    #####
    
    #####
    #Chen et al. (2004) Step 2. Long-term change Trend Fitting
    #####
    
    #fitting long-term trends requires evaluating different paramters for the SG filter
    
    #creating a matrix to hold fit statistics for models with different parameters
    msel = matrix(nrow=12,ncol=3)
    colnames(msel) = c("span","degree","RSS")
    
    #creating index variable to populate the matrix
    m = 1
    
    #creating vectors holding a suite of span and polynomial degrees to be selected
    spanv = c(9,11,13,15)
    degreev = c(2,3,4)
    
    #using for loop to iteratively fit the loess models
    for (i in spanv){
      
      for (j in degreev){
        
        #using locally weighted scatterplot smoothing to produce a long-term change trend curve for NDVI values (Chen et al. 2004)
        sgpred = sgolayfilt(ndvim$No,p=j,n=i) #fitting Savitzky-Golay filter
              
        rss = sum((ndvim$No - sgpred)^2) #calculating residual sum of squares as a measure of model fit
        msel[m,1] = i
        msel[m,2] = j
        msel[m,3] = rss
        
        m = m + 1
        
      }
      
    }
    
    #Using best fit model to produce a long-term change trend curve for NDVI values (Chen et al. 2004)
    Ntr = sgolayfilt(ndvim$No,p = msel[,2][which(msel[,3] == min(msel[,3]))],n = msel[,1][which(msel[,3] == min(msel[,3]))]) #fitting modely using Savitzky-Golay filter
    ndvim$Ntr = Ntr
      
    #####
    #Chen et al. (2004) Step 3. Determination of weight for each point in NDVI time-series
    #####
    
    ndvim$wt = ifelse(ndvim$No>=ndvim$Ntr,1,(1-abs(ndvim$No-ndvim$Ntr)/max(abs(ndvim$No-ndvim$Ntr))))
    
    #####
    #Chen et al. (2004) Step 4. Generation of a new NDVI time-series
    #####
    
    ndvim$N1 = ifelse(ndvim$No>=ndvim$Ntr,ndvim$No,ndvim$Ntr)
    
    #####
    #Chen et al. (2004) Step 5/6. Fitting new NDVI time-series using the Savitzky-Golay filter/Calculation of a fitting-effect index
    #####
    
    NDVIStore = matrix(nrow=nrow(ndvim),ncol=100) #matrix to store iteratively produced smoothed NDVI time series, 100 iterations will be run after Chen et al. (2004)
    FkStore = c() #vector to store index to determine what NDVI time series to use
    
    Nk = ndvim$N1 #creating vector of time series to be smoothed by SG filter each time through the for loop
    
    for (i in 1:100){
      
      #using locally weighted scatterplot smoothing to produce a long-term change trend curve for NDVI values (Chen et al. 2004)
      Nk1 = sgolayfilt(Nk,p=6,n=9) #fitting Savitzky-Golay filter for smoothing, span and polynomial degress after Chen et al. (2004)
      Fk = sum(abs(Nk1-ndvim$No)*ndvim$wt) #calculating fit index to determine how well series fits new data
      FkStore = c(FkStore,Fk) #storing fit statistics
      NDVIStore[,i]=Nk1
      Nk = ifelse(ndvim$No>=Nk1,ndvim$No,Nk1)
      
    }
    
    #storing best smoothed NDVI values based on minimum value of fitting index
    ndvim$Nfinal = NDVIStore[,which(FkStore==min(FkStore))]
    
    #removing dates from other years used to avoid edge effects when smoothing
    ndvim = ndvim[rownames(ndvim)[grep(paste("^A",year,sep=""),rownames(ndvim))],]
    
    #####
    #linear interpolation used to determine daily NDVI values
    #####
    for (r in 1:(nrow(ndvim)-1)){
      
      lint = lm(Nfinal~Day,data=ndvim[r:(r+1),]) #creating linear model for two successive NDVI points
      days = seq(ndvim$Day[r],ndvim$Day[r+1],1) #creating sequence of dates between the two NDVI points
      days = as.data.frame(days[-c(1,length(days))]) #creating sequence of dates between the two NDVI points
      colnames(days)="Day" #renaming new dates so that model is able to recognize predictor
      ndvi.daily = predict.lm(lint,days) #predicting ndvi values for dates
      ndvi.daily = cbind(days,ndvi.daily)
      colnames(ndvi.daily)[2]=paste("y",year,sep="")
      
      if (r==1){
        
        daily=ndvi.daily
        
      }else{
        
        daily = rbind(daily,ndvi.daily)
        
      }
      
      if (r == (nrow(ndvim)-1)){ #adding back in days that were in the original data series
        
        mdata = ndvim[,c("Day","Nfinal")]
        colnames(mdata)[2] = paste("y",year,sep="")
        daily = rbind(mdata,daily)
        daily = daily[order(daily$Day),]
        
      }
      
    }
    
    #####
    #Storing smoothed data as well as interpolated daily data
    #####
    #storing output as data frame
    
    if (year==years[1]){
      
      yadd = ndvim[,c("Day","No","Nfinal")]
      
      colnames(yadd)[2:3] = c(paste("y",year,"orig",sep=""),paste("y",year,"smooth",sep="")) 
      
      yearStore = yadd
      
      dailyStore = daily
      
    } else {
      
      yadd = ndvim[,c("Day","No","Nfinal")]
      
      colnames(yadd)[2:3] = c(paste("y",year,"orig",sep=""),paste("y",year,"smooth",sep=""))
      
      yearStore = merge(yearStore,yadd,by="Day",all=T)
      
      dailyStore = merge(dailyStore,daily,by="Day",all=T)
    }
      
  }
  
  #creating plot containing smoothed NDVI curves for all years
  
  pdf(paste("W:\\TNC_Land_Protection\\2015_SustainableGrazinglands\\Output\\",pref,ranch,".pdf",sep=""))
  
  ymax = max(apply(dailyStore,2,max,na.rm=T))
  
  for (i in 1:nyr){
        
    if (i ==1){ #first year creates plots
      
      plot(dailyStore[,1],dailyStore[,i+1],ylim=c(0,ymax),ylab="NDVI",xlab="Day",col = colStore[i],type="l",lwd=2)
      
      
    } else { #all years after first year are added to existing plot
      
      lines(dailyStore[,1],dailyStore[,i+1],col = colStore[i],lwd=2)
      
    }
    
    if (i==nyr){ #after final year is added, a legend is created
      
      legend(50,1000,seq(syr,fyr,1),ncol=4,lwd=2,col=colStore)
      
    }
  }
  
  dailym = apply(dailyStore[,2:(nyr+1)],1,mean) #averaging annual average NDVI value for each day
  
  lines(dailyStore$Day,dailym,col="black",lwd=5) #adding annual average to plot of NDVI values
  
  abline(h=ndvit,lwd=2,lty=2)
  
  text(40,ndvit+1500,paste("Threshold =\n ",round(ndvit,1),sep=""))
  
  dev.off()
  
  #creating for loops to identify start and end of growing season, length of growing season, and integrated area under the curve
  for (i in 1:nyr){
    
    temp <- dailyStore[,(1+i)][which(dailyStore[,(1+i)]>ndvit)] #temp data set will be used to grab values associated with growing season
    
    if (length(temp)==0){ #if statement to assign values when no growing season occurred
      
      sgseas[ranch,i] = NA #first date above threshold is start of growing season
      egseas[ranch,i] = NA #last date above threshold is end of growing season
      lgseas[ranch,i] = 0 #length of growing season
      auc[ranch,i] = NA #integrated area under the curve for growing season
      mean.ndvi[ranch,i] = NA #Average daily NDVI
      median.ndvi[ranch,i] = NA #Median daily NDVI
      min.ndvi[ranch,i] = NA #Minimum daily NDVI
      min.ndvi.date[ranch,i] = NA #Date of Minimum daily NDVI
      max.ndvi[ranch,i] = NA #Maximum daily NDVI
      max.ndvi.date[ranch,i] = NA #Date of Maximum daily NDVI      
      
    } else {
      
      gseas <- dailyStore[,(1+i)][which(dailyStore[,(1+i)]==temp[1]):which(dailyStore[,(1+i)]==temp[length(temp)])] #grabbing values associated with growing season
      
      sgseas[ranch,i] = dailyStore[,1][which(dailyStore[,(1+i)]==gseas[1])]#first date above threshold is start of growing season
      egseas[ranch,i] = dailyStore[,1][which(dailyStore[,(1+i)]==gseas[length(gseas)])] #last date above threshold is end of growing season
      lgseas[ranch,i] = egseas[ranch,i]-sgseas[ranch,i] #length of growing season
      
      indvi = c() #creating vector to store area under the curve calculations
      
      for (f in 1:(length(gseas)-1)){ #using trapezoidal area under the curve approach, note that height is not shown because it is = 1 day for all intervals
        
        indvi[f] = (gseas[f]+gseas[f+1])/2
        
      }
      
      auc[ranch,i]= sum(indvi) #integrated area under the curve for growing season
      
      mean.ndvi[ranch,i] = mean(gseas) #mean daily NDVI within growing season
      median.ndvi[ranch,i] = median(gseas) #median daily NDVI within growing season
      min.ndvi[ranch,i] = min(gseas) #median daily NDVI within growing season
      min.ndvi.date[ranch,i] = dailyStore[,1][which(dailyStore[,(1+i)]==min(gseas))]#date of minimum NDVI value
      max.ndvi[ranch,i] = max(gseas) #median daily NDVI within growing season
      max.ndvi.date[ranch,i] = dailyStore[,1][which(dailyStore[,(1+i)]==max(gseas))] #date of maximum NDVI value
      
    }
  }
  
}

#####
#Creating output
#####

#list of data frames with output metrics
outfiles = list(sgseas, egseas, lgseas, auc, mean.ndvi, median.ndvi, min.ndvi, min.ndvi.date, max.ndvi, max.ndvi.date)

#list of names to assign to output files
outnames = c("StartOfGrowingSeason","EndOfGrowingSeason","LengthOfGrowingSeason","AreaUnderTheCurve","MeanNDVI","MedianNDVI","MinNDVI","MinNDVIDate",
             "MaxNDVI","MaxNDVIDate")

#as requested by project manager, creating separate output files for each metric
for (i in 1:10){
  
  out = outfiles[[i]]
  
  write.csv(out,paste("W:\\TNC_Land_Protection\\2015_SustainableGrazinglands\\Output\\",pref,outnames[i],".csv",sep=""))  
  
}
