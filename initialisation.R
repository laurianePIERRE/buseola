#Initialisation


# prior distributions
prior=list()
prior$K$distribution = "uniform"
prior$K$model = c(K="proportional")
prior$K$p = c(min=0.001,max=0.5)
prior$R$distribution = "fixed"
prior$R$model = c(K="constant")
prior$R$p = 20
prior$mutation_rate$model = "stepwise"
prior$mutation_rate$distribution = "loguniform"
prior$mutation_rate$p = c(min=1E-6,max=1E-2)
prior$dispersion$model="contiguous"
prior$dispersion$distribution="uniform"
prior$dispersion$p=c(min=0.001,max=0.5)


rm(list=ls())
wd="/home/dupas/Bureau/" # fixe
setwd(wd)
source.file=c("Laurianne.R")

library(raster)

environmentalData_ <- raster("AfrBfGAMActA.asc");environmentalData_[environmentalData_>200]=200
environmentalData_2 <- trim(environmentalData_, padding=(8*dim(environmentalData_)%/%8)[1:2]) 
environmentalData_2 <- aggregate(environmentalData_2,8)
genetData <- read.table("WBf16genelandcoord.txt")
genetData <- cbind(genetData,read.table("WBf16genelandgeno.txt"))
genetSP <- SpatialPoints(genetData[,c("x","y")])
bbox(genetSP)
environmentalData <- crop(environmentalData_2, bbox(genetSP)+res(environmentalData_2)*c(-2,-2,2,2))
genetData$Cell_numbers <- cellFromXY(environmentalData,genetData)
environmentalData[environmentalData==0] <- NA
ncellA(environmentalData)
any(is.na(extract(environmentalData,genetData[,c("x","y")])))
plot(environmentalData)
plot(genetSP,add=TRUE)
min(extract(environmentalData,genetData[,c("x","y")]))
genetData[which(extract(environmentalData,genetData[,c("x","y")])<0.01),]

genetData[is.na(genetData)] <- as.integer(1E9)
genetData=genetData[rowSums(genetData[,grep("ocus",colnames(genetData),value=TRUE)])<7E9,]
genetData[genetData==1E9]=NA
genetData <- TwoCols2OneCol(genetData)
