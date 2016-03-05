#Initialisation


# prior distributions
prior=list()
prior$K$distribution = "uniform"
prior$K$model = c(K="proportional")
prior$K$p = data.frame( a=c(min=0.001,max=0.5))
prior$R$distribution = "fixed"
prior$R$model = c(K="constant")
prior$R$p = data.frame(a=20)
prior$mutation_rate$model = "stepwise"
prior$mutation_rate$distribution = "loguniform"
prior$mutation_rate$p =data.frame(a=c(min=1E-6,max=1E-2))
prior$dispersion$model="contiguous"
prior$dispersion$distribution="uniform"
prior$dispersion$p=data.frame(a=c(min=0.001,max=0.5))

# fixed uniform normal loguniform lognormal



wd="~/Documents/M1bi 2015 2016/stage/busseola/"
setwd(wd)
source("Laurianne.R")
source("generic.R");source("method.R")
library(raster)
library(rgdal) # necessary to use function raster()
environmentalData <- raster("busseola.tif")
genetData <- read.table("WBf16genelandcoord.txt")
genetData <- cbind(genetData,read.table("WBf16genelandgeno.txt"))
genetSP <- SpatialPoints(genetData[,c("x","y")])
bbox(genetSP)
environmentalData <- crop(environmentalData, bbox(genetSP)+res(environmentalData)*c(-2,-2,2,2))
genetData$Cell_numbers <- cellFromXY(environmentalData,genetData)
environmentalData[environmentalData==0] <- NA
ncellA(environmentalData)
any(is.na(extract(environmentalData,genetData[,c("x","y")])))
plot(environmentalData)
plot(genetSP,add=TRUE)
min(extract(environmentalData,genetData[,c("x","y")]))
genetData[which(extract(environmentalData,genetData[,c("x","y")])<0.01),]
genetData=cleanerData(genetData)
transitionMatrixBackward(environmentalData,prior)

