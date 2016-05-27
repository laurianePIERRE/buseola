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



wd="/M1bi 2015 2016/stage/busseola/"

wd="~/Documents/Lauriane/busseola/"

setwd(wd)
library(raster)
library(rgdal) # necessary to use function raster()
source("Laurianne.R")
source("class.R")
source("generic.R");source("method.R")
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
listePrior=sampleP(prior) 
matriceMigration=migrationMatrixA(environmentalData,listePrior$dispersion$busseola$model,listePrior$dispersion$busseola$p)
matriceTrans=transitionMatrixA(environmentalData,prior)
listePrior=sampleP(prior) 
plot(raster(matriceMigration))
plot(raster(matriceTrans))
#diminution des pixels
environmentalDataSimple=Aggregate_and_adjust_raster_to_data(environmentalData,xy=genetData[,c("x","y")],extend_band_size=1,aggregate_index=4)
plot(environmentalDataSimple)
matriceTrans=transitionMatrixA(environmentalDataSimple,prior)
plot(raster(matriceTrans))
matriceAbsorbante=absorbingTransitionA(matriceTrans,valuesA(environmentalDataSimple))
plot(raster(matriceAbsorbante[[1]]))

# avec l utilisation de la class

Matrice_de_transition=new("Transition_Matrix",matriceTrans,populationsize=environmentalDataSimple)
absorbante=absorbingTransitionA(Matrice_de_transition)
# pour simplifier 

ma=raster(matrix(c(20,30,40,10,NA,NA,20,30,40),nrow=3,ncol=3))
extent(ma)=c(0,3,0,3)
names(ma)="busseola"
values(ma)=c(20,30,40,10,NA,NA,20,30,40)
plot(ma)
dim(ma)
matriceAbsMa=migrationMatrixA(ma,listePrior$dispersion$busseola$model,listePrior$dispersion$busseola$p)

plot(raster(matriceAbsMa))
matriceTransMa=transitionMatrixA(ma,prior)
plot(raster(matriceTransMa))
matriceAbsorbanteMa=absorbingTransitionA(matriceTransMa,valuesA(ma))
plot(raster(matriceAbsorbanteMa[[1]]))
pdf("exemplesimple.pdf")
par(mfrow=c(2,2))
plot(ma,main="raster de depart")
plot(raster(matriceAbsMa),main="matrice de migration")
plot(raster(matriceTransMa),main="matrice de transition")
dev.off()

pdf("comparaison.pdf")
par(mfrow=c(4,2))
plot(ma,main="raster de depart simple")
plot(environmentalDataSimple,main="raster de depart busseola")
plot(raster(matriceAbsMa),main="matrice de migration simple")
plot(raster(matriceMigration),main="matrice de migration busseola")
plot(raster(matriceTransMa),main="matrice de transition simple")
plot(raster(matriceTrans),main="matrice de transition busseola")
plot(raster(matriceAbsorbanteMa[[1]]),main="matrice absorbante simple")
plot(raster(matriceAbsorbante[[1]]),main="matrice absorbante busseola")
dev.off()


#genealogie 
load("genealogy.rda")

# donnes genetique pour le mini jeu de donnees genotypes
load("genotypes.rda")
gen <- new("genetic",data.frame(Locus1=genotypes[,"Locus1"]),ploidy=as.integer(1), ploidyByrow=FALSE)
spgen <- new("spatialGenetic",gen,x=genotypes[,"x"],y=genotypes[,"y"],Cell_numbers=genotypes[,"Cell_numbers"])

genealogy[[2]]
library(ape)
library(stringr)
plot(coalescent_2_phylog(genealogy))
coalescent=new("Genealogy",genealogy,genotypes[,2])
plotgenealogy(coalescent,tipcols = NA)
geneasimple=new("LandGenealogy",ma,genealogy=coalescent)
plotLandG(geneasimple,rasK = NULL)
geneaandgenet=new("LandGenetGenalogy",geneasimple,Genotype=spgen)


########## modify gealogy to colorate this

genealogymod=change_genealogy(genealogy)


# try to use phtools for colorate but no my data
# our data = ma (raster simplify dim (3,3,1))

## set seed for reproducibility
set.seed(100)
## load phytools
library(phytools)
## simulate a tree & some data
tree <- pbtree(n = 26, scale = 1)
tree$tip.label <- LETTERS[26:1]
x <- fastBM(tree)

Q <- matrix(c(-2, 1, 1, 1, -2, 1, 1, 1, -2), 3, 3)
colnames(Q) <- rownames(Q) <- letters[1:3]
x <- sim.history(tree, Q, anc = "a")$states
print(x)
# estimate ancestral states under a ER model
fitER <- rerootingMethod(tree, x, model = "ER")
print(fitER)
plotTree(tree, setEnv = TRUE, offset = 0.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "red", "yellow"), cex = 0.6)
tiplabels(pie = to.matrix(x, sort(unique(x))), piecol = c("blue", "red", "yellow"), 
          cex = 0.3)
