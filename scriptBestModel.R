# Predict West Nile virus in mosquitos across the city of Chicago
# Kaggle competition https://www.kaggle.com/c/predict-west-nile-virus
# This script was ranked 85/1306 and achieved AUC 0.80562 in the 
# private part of the dataset
# Author: Andrey Bortsov
# May 24, 2015


setwd("~/Kaggle/westnile")

library(Metrics)
library(data.table)  
train <- fread("~/Kaggle/westnile/data/train.csv")
test <- fread("~/Kaggle/westnile/data/test.csv")
weather <- fread("~/Kaggle/westnile/data/weather.csv")


## prepare weather data
## assign closest weather stations to traps
distance <- function(longitude, latitude) { 
  dist1 <- sqrt((stations[1,]$Latitude-latitude)^2+(stations[1,]$Longitude-longitude)^2)
  dist2 <- sqrt((stations[2,]$Latitude-latitude)^2+(stations[2,]$Longitude-longitude)^2)  
  if(dist1<dist2){return(1)}
  return(2)}
stations<-data.frame(c(1,2),c(41.995,41.786),c(-87.933,-87.752))
names(stations)<-c("Station","Latitude","Longitude")

train$Station<-mapply(distance,train$Longitude,train$Latitude)
test$Station<-mapply(distance,test$Longitude,test$Latitude)

#calculate average daily temp in Celsius
weather$TavgC <- (as.numeric(weather$Tavg) - 32)*5/9 

weather[,dMonth:=as.factor(paste(substr(weather$Date,6,7)))]
weather[,dYear:=as.factor(paste(substr(weather$Date,1,4)))]
weather$Date1 = as.Date(weather$Date, format="%Y-%m-%d")
tsDate = as.Date(paste0(weather$dYear, "0101"), format="%Y%m%d")
weather$dWeek = as.numeric(paste(floor((weather$Date1 - tsDate + 1)/7)))
weather$dDay = as.numeric(paste(floor((weather$Date1 - tsDate + 1))))

weather$PrecipTotal[weather$PrecipTotal=='T'] <- '0' ##trace precip change to 0

# average temp and precip
tn <-      cbind(weather$TavgC[weather$dYear=='2007'], 
                 weather$TavgC[weather$dYear=='2008'],
                 weather$TavgC[weather$dYear=='2009'],
                 weather$TavgC[weather$dYear=='2010'],
                 weather$TavgC[weather$dYear=='2011'],
                 weather$TavgC[weather$dYear=='2012'],
                 weather$TavgC[weather$dYear=='2013'],
                 weather$TavgC[weather$dYear=='2014'])
mean1 <- function(x){mean(as.numeric(x), na.rm=TRUE)}
TmprNorm <- apply(tn, 1, mean1)

pn <-      cbind(weather$PrecipTotal[weather$dYear=='2007'], 
                 weather$PrecipTotal[weather$dYear=='2008'],
                 weather$PrecipTotal[weather$dYear=='2009'],
                 weather$PrecipTotal[weather$dYear=='2010'],
                 weather$PrecipTotal[weather$dYear=='2011'],
                 weather$PrecipTotal[weather$dYear=='2012'],
                 weather$PrecipTotal[weather$dYear=='2013'],
                 weather$PrecipTotal[weather$dYear=='2014'])
PrecipNorm <- apply(pn, 1, mean1)

time <- c(1:368)
PrecipNormSm <- loess(PrecipNorm ~ time)
PrecipNorm <- PrecipNormSm$fitted

TmprNormSm <- loess(TmprNorm ~ time)
TmprNorm <- TmprNormSm$fitted

weather$Tnorm <- rep(TmprNorm, 8)
weather$Pnorm <- rep(PrecipNorm, 8)


## calculate warm days (>22) and cool days (<22)
weather$DD <- (weather$TavgC - 22)*(weather$TavgC >= 22)
weather$CDD <- (weather$TavgC - 22)*(weather$TavgC < 22)
weather$DD[is.na(weather$DD)] <- 0
weather$CDD[is.na(weather$CDD)] <- 0

## precipitation
weather$PrecipTotal <- as.numeric(weather$PrecipTotal) 
weather$PrecipTotal[is.na(weather$PrecipTotal)] <- 0
weather$flood <- ifelse(weather$PrecipTotal > 0.75, 1, 0)
weather$flood[is.na(weather$flood)] <- 0

##wind
weather$AvgSpeed  <- as.numeric(weather$AvgSpeed) 
weather$AvgSpeed[is.na(weather$AvgSpeed)] <- 7.5

# create function to calculate lagged exposure
library(Hmisc)
ma <- function(x,days,dlag){as.numeric(Lag(filter(x,rep(1/days,days), method="convolution", sides=1), dlag))}

## calculate lagged and cumulative precip and warm days by weather station
weatherAll <- NULL
for (i in 1:2){  
  wst <- weather[weather$Station==i,]
  wst$DDcum <- unlist(tapply(wst$DD, wst$dYear, cumsum))
  wst$CDDcum <- unlist(tapply(wst$CDD, wst$dYear, cumsum))
  wst$PrecipCum <- unlist(tapply(wst$PrecipTotal, wst$dYear, cumsum))
  
  wst$Precip7lag0  <- ma(wst$PrecipTotal, 7, 0)
  wst$Precip14lag0 <- ma(wst$PrecipTotal, 14, 0)
  wst$Precip21lag0 <- ma(wst$PrecipTotal, 21, 0)
  wst$Precip28lag0 <- ma(wst$PrecipTotal, 28, 0)
  
  wst$Precip7lag7  <- ma(wst$PrecipTotal, 7, 7)
  wst$Precip14lag7 <- ma(wst$PrecipTotal, 14, 7)
  wst$Precip21lag7 <- ma(wst$PrecipTotal, 21, 7)
  wst$Precip28lag7 <- ma(wst$PrecipTotal, 28, 7)
  
  wst$Precip7lag14  <- ma(wst$PrecipTotal, 7, 14)
  wst$Precip14lag14 <- ma(wst$PrecipTotal, 14, 14)
  wst$Precip21lag14 <- ma(wst$PrecipTotal, 21, 14)
  wst$Precip28lag14 <- ma(wst$PrecipTotal, 28, 14)
  
  wst$Precip7lag21  <- ma(wst$PrecipTotal, 7, 21)
  wst$Precip14lag21 <- ma(wst$PrecipTotal, 14, 21)
  wst$Precip21lag21 <- ma(wst$PrecipTotal, 21, 21)
  wst$Precip28lag21 <- ma(wst$PrecipTotal, 28, 21)
  
  wst$DD7lag0  <- ma(wst$DD, 7, 0)
  wst$DD14lag0 <- ma(wst$DD, 14, 0)
  wst$DD21lag0 <- ma(wst$DD, 21, 0)
  wst$DD28lag0 <- ma(wst$DD, 28, 0)
  
  wst$DD7lag7  <- ma(wst$DD, 7, 7)
  wst$DD14lag7 <- ma(wst$DD, 14, 7)
  wst$DD21lag7 <- ma(wst$DD, 21, 7)
  wst$DD28lag7 <- ma(wst$DD, 28, 7)
  
  wst$DD7lag14  <- ma(wst$DD, 7, 14)
  wst$DD14lag14 <- ma(wst$DD, 14, 14)
  wst$DD21lag14 <- ma(wst$DD, 21, 14)
  wst$DD28lag14 <- ma(wst$DD, 28, 14)
  
  wst$DD7lag21  <- ma(wst$DD, 7, 21)
  wst$DD14lag21 <- ma(wst$DD, 14, 21)
  wst$DD21lag21 <- ma(wst$DD, 21, 21)
  wst$DD28lag21 <- ma(wst$DD, 28, 21)
  
  wst$CDD7lag0  <- ma(wst$CDD, 7, 0)
  wst$CDD14lag0 <- ma(wst$CDD, 14, 0)
  wst$CDD21lag0 <- ma(wst$CDD, 21, 0)
  wst$CDD28lag0 <- ma(wst$CDD, 28, 0)
  
  wst$CDD7lag7  <- ma(wst$CDD, 7, 7)
  wst$CDD14lag7 <- ma(wst$CDD, 14, 7)
  wst$CDD21lag7 <- ma(wst$CDD, 21, 7)
  wst$CDD28lag7 <- ma(wst$CDD, 28, 7)
  
  wst$CDD7lag14  <- ma(wst$CDD, 7, 14)
  wst$CDD14lag14 <- ma(wst$CDD, 14, 14)
  wst$CDD21lag14 <- ma(wst$CDD, 21, 14)
  wst$CDD28lag14 <- ma(wst$CDD, 28, 14)
  
  wst$CDD7lag21  <- ma(wst$CDD, 7, 21)
  wst$CDD14lag21 <- ma(wst$CDD, 14, 21)
  wst$CDD21lag21 <- ma(wst$CDD, 21, 21)
  wst$CDD28lag21 <- ma(wst$CDD, 28, 21)
  
  wst$Wind7lag0  <- ma(wst$AvgSpeed, 7, 0)
  wst$Wind14lag0 <- ma(wst$AvgSpeed, 14, 0)
  wst$Wind21lag0 <- ma(wst$AvgSpeed, 21, 0)
  wst$Wind28lag0 <- ma(wst$AvgSpeed, 28, 0)
  
  wst$Wind7lag7  <- ma(wst$AvgSpeed, 7, 7)
  wst$Wind14lag7 <- ma(wst$AvgSpeed, 14, 7)
  wst$Wind21lag7 <- ma(wst$AvgSpeed, 21, 7)
  wst$Wind28lag7 <- ma(wst$AvgSpeed, 28, 7)
  
  wst$Wind7lag14  <- ma(wst$AvgSpeed, 7, 14)
  wst$Wind14lag14 <- ma(wst$AvgSpeed, 14, 14)
  wst$Wind21lag14 <- ma(wst$AvgSpeed, 21, 14)
  wst$Wind28lag14 <- ma(wst$AvgSpeed, 28, 14)
  
  wst$Wind7lag21  <- ma(wst$AvgSpeed, 7, 21)
  wst$Wind14lag21 <- ma(wst$AvgSpeed, 14, 21)
  wst$Wind21lag21 <- ma(wst$AvgSpeed, 21, 21)
  wst$Wind28lag21 <- ma(wst$AvgSpeed, 28, 21)
  
  weatherAll <- rbind(weatherAll, wst)
}

## cumulative spring (up to 153rd day) temp and precip 
DDcumMay <- weatherAll[weatherAll$dDay==153 & weatherAll$Station==1,DDcum]
PrecipcumMay <- weatherAll[weatherAll$dDay==153 & weatherAll$Station==1,PrecipCum]
MayData <- data.frame(cbind(DDcumMay, PrecipcumMay, c(2007:2014)))
colnames(MayData)[3] <- 'dYear'
MayData$dYear <- as.factor(MayData$dYear)
weatherAll <- merge(weatherAll, MayData, by="dYear", all.x=TRUE)

## impute missing weather data
weatherAllnum <- apply(weatherAll,2,as.numeric)
weatherAllnum <- weatherAllnum[,-c(3,14,15,16,26)]

library(mice)
weatherAllnumi <- complete(mice(weatherAllnum,m=1))
weatherAllnumi$Date <- weatherAll$Date


#merge weather with train and test data
datechar <- as.character(train[,Date])
train$Date <- datechar
train <- merge(train, weatherAllnumi,by=c("Date","Station"), all.x=TRUE)

datechar1 <- as.character(test[,Date])
test$Date <- datechar1
test <- merge(test, weatherAllnumi,by=c("Date","Station"), all.x=TRUE)


## prep the species column
vSpecies<-c(as.character(train$Species),as.character(test$Species))
vSpecies[which(vSpecies == "CULEX ERRATICUS" |
                 vSpecies == "CULEX SALINARIUS" |
                 vSpecies == "CULEX TARSALIS" |
                 vSpecies == "CULEX TERRITANS" |
                 vSpecies == "UNSPECIFIED CULEX")] = "CULEX OTHER"
vSpecies<-factor(vSpecies,levels=unique(vSpecies))

vSpeciesDum <- model.matrix(~vSpecies)
vSpeciesDum3 <-vSpeciesDum[,2:4]
colnames(vSpeciesDum3) <- c("RESTUANS","PIPIENS","CulexOTHER")

train <-cbind(train, vSpeciesDum3[1:nrow(train),] )
test <-cbind(test, vSpeciesDum3[(nrow(train)+1):length(vSpecies),] )

#train[,Species2:=factor(vSpecies[1:nrow(train)],levels=unique(vSpecies))]
#test[,Species2:=factor(vSpecies[(nrow(train)+1):length(vSpecies)],levels=unique(vSpecies))]

train[,dMonth:=as.numeric(paste(substr(train$Date,6,7)))]
train[,dYear:=as.factor(paste(substr(train$Date,1,4)))]
train$Date = as.Date(train$Date, format="%Y-%m-%d")
xsDate = as.Date(paste0(train$dYear, "0101"), format="%Y%m%d")
train$dWeek = as.numeric(paste(floor((train$Date - xsDate + 1)/7)))
train$dDay = as.numeric(paste(floor((train$Date - xsDate + 1))))

test[,dMonth:=as.numeric(paste(substr(test$Date,6,7)))]
test[,dYear:=as.factor(paste(substr(test$Date,1,4)))]
test$Date = as.Date(test$Date, format="%Y-%m-%d")
tsDate = as.Date(paste0(test$dYear, "0101"), format="%Y%m%d")
test$dWeek = as.numeric(paste(floor((test$Date - tsDate + 1)/7)))
test$dDay = as.numeric(paste(floor((test$Date - tsDate + 1))))

# split predictors and outcome in train data
my.x = subset(train, select=-c(2:8,11:14))
my.y = train$WnvPresent


library(SuperLearner)

## modify SL wrappers
SL.Mygam <- function (Y, X, newX, family, deg.gam = 2, cts.num = 4, ...) 
{ require(gam)
  fit.gam <- gam::gam(
    
    ## edit this GAM model as needed
    Y ~  lo(Latitude, Longitude, span=0.7) + RESTUANS + PIPIENS + CulexOTHER + Precip14lag0 +
      DD7lag0 + AvgSpeed + DDcum*PrecipCum,
    
    data = X, family = family, control = gam::gam.control(maxit = 50, bf.maxit = 50))
  pred <- gam::predict.gam(fit.gam, newdata = newX, type = "response")
  fit <- list(object = fit.gam)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.gam")
  return(out)
}


screen.randomForest.5 <- function(..., nVar=5){screen.randomForest(..., nVar=nVar)}
screen.randomForest.20 <- function(..., nVar=20){screen.randomForest(..., nVar=nVar)}
screen.randomForest.30 <- function(..., nVar=30){screen.randomForest(..., nVar=nVar)}

create.SL.glmnet <- function(alpha = c(0.25, 0.50, 0.75)) {
  for(mm in seq(length(alpha))){
    eval(parse(text = paste('SL.glmnet.', alpha[mm], '<- function(..., alpha = ', alpha[mm], ') SL.glmnet(..., alpha = alpha)', sep = '')), envir = .GlobalEnv)
  }
  invisible(TRUE)
}
create.SL.glmnet(alpha=c(0.25, 0.50, 0.75))

## put all wrappers in SL library
SL.library <- list(c("SL.randomForest","screen.randomForest"),
                   c("SL.randomForest","screen.randomForest.5"),
                   c("SL.randomForest","screen.randomForest.20"),
                   c("SL.glmnet","screen.randomForest.20"),
                   c("SL.glmnet.0.25","screen.randomForest.20"),
                   c("SL.glmnet.0.5","screen.randomForest.20"),
                   c("SL.glmnet.0.75","screen.randomForest.20"),
                   c("SL.gam","screen.glmnet"),
                     "SL.Mygam")

## fit SuperLearner model
fitSL <- SuperLearner (Y = my.y,
                       X = my.x,
                       SL.library=SL.library,
                       family = binomial(),
                       method = "method.NNLS",
                       verbose = TRUE,
                       cvControl = list(V = 10))
summary(fitSL)

#selected <- fitSL$varNames[fitSL$whichScreen]

## cross-validate SuperLearner model
fitSL.CV <- CV.SuperLearner (Y = my.y,
                             X = my.x,
                             SL.library=SL.library,
                             V = 10, family = binomial(),
                             method = "method.NNLS",
                             cvControl = list(stratifyCV=TRUE), parallel = 'multicore',
                             verbose=TRUE)
summary(fitSL.CV)
predictions2 <- fitSL.CV$SL.predict

## check cross-validated AUC
library(pROC)
ROC1 <- roc(my.y=='1' ~ predictions2)
plot(ROC1) 

## obtain predictions for test data using the fitted model
cols <- colnames(my.x)
test2 <- subset(test, select=c(cols))
predictions <- predict(fitSL, test2)

## prepare submission file
submissionFile<-cbind(test$Id,predictions[[1]])
colnames(submissionFile)<-c("Id","WnvPresent")
options("scipen"=100, "digits"=8)
write.csv(submissionFile,"submitSLv6-4.csv",row.names=FALSE,quote=FALSE)

