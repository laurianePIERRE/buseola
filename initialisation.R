#Initialisation


# prior distributions
prior=list()
prior$K$busseola$a$distribution = "uniform"
prior$K$busseola$model = data.frame(busseola="proportional")
prior$K$busseola$a$p = data.frame(busseola=c(min=0.001,max=0.5))
prior$R$busseola$a$distribution = "fixed"
prior$R$busseola$model =data.frame(busseola="constant")
prior$R$busseola$a$p = 20
prior$mutation_rate$busseola$model =data.frame(busseola= "stepwise")
prior$mutation_rate$busseola$a$distribution = "loguniform"
prior$mutation_rate$busseola$a$p =data.frame(busseola=c(min=1E-6,max=1E-2))
prior$dispersion$busseola$model=data.frame(busseola="contiguous")
prior$R$busseola$model = c(K="constant")

prior$K$busseola$model =c("proportional")
prior$K$busseola$a$p = data.frame(busseola=c(min=0.001,max=0.5))
prior$R$busseola$a$distribution = "fixed"
prior$R$busseola$model =c("constant")
prior$R$busseola$a$p = 20
prior$mutation_rate$busseola$model =data.frame(busseola= "stepwise")
prior$mutation_rate$busseola$a$distribution = "loguniform"
prior$mutation_rate$busseola$a$p =data.frame(busseola=c(min=1E-6,max=1E-2))

prior$dispersion$busseola$model="contiguous"
prior$dispersion$busseola$model=data.frame(busseola="contiguous")
prior$R$busseola$model = c("constant")
prior$R$busseola$a$p = 20
prior$mutation_rate$busseola$model =c("stepwise")
prior$mutation_rate$busseola$a$distribution = "loguniform"
prior$mutation_rate$busseola$a$p =data.frame(busseola=c(min=1E-6,max=1E-2))
prior$dispersion$busseola$model=c("contiguous")
prior$dispersion$busseola$a$distribution="uniform"
prior$dispersion$busseola$a$p=data.frame(busseola=c(min=0.001,max=0.5))

# fixed uniform normal loguniform lognormal



wd="~/Documents/M1bi 2015 2016/stage/busseola/"

wd="~/Documents/Lauriane/busseola/"

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
matriceMigration=migrationMatrixA(environmentalData,listePrior$dispersion$busseola$model,listePrior$dispersion$busseola$p)
matriceTrans=transitionMatrixA(environmentalData,prior)
listePrior=sampleP(prior) 
plot(raster(matriceMigration))
plot(raster(matriceTrans))
AbsorbingTransition(matriceTrans[1:50,1:50],valuesA(environmentalData))
