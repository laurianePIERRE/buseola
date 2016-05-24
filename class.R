prior <- setClass("prior",
                  contains="list") # add validity (contains K, R and mutation_rate)

parameters <- setClass("parameters",
                       contains="list") # add validity (contains K, R and mutation_rate)


transitionModel <- setClass("transitionModel",
                           slots = c(demicTransition="matrix",
                                     allelicTransition="matrix",
                                     Ne="numeric", 
                                     demesNames="character", 
                                     allelesNames="character", 
                                     demicStatusOfStartingIndividuals="character", 
                                     allelicStatusOfStartingIndividuals="character"),
                           validity = function(object){
                                if(length(object@demicStatusOfStartingIndividuals)!=length(object@allelicStatusOfStartingIndividuals))  {stop("length of demicStatusOfStartingIndividuals differs from length of allelicStatusOfStartingIndividuals")}
                                if(length(object@Ne)!=length(object@demesNames))  {stop("length of Ne differs from length of demesnames")}
                                if(length(object@Ne)!=nrow(object@demicTransition))  {stop("length of Ne differs from number of rows and columns in demic transition matrix")}
                                if(ncol(object@demicTransition)!=nrow(object@demicTransition))  {stop("demic transition matrix is not square")}
                                if(length(object@allelesNames)!=nrow(object@demicTransition))  {stop("length of allele names differs from nrow and ncols of demicTransition")}
                              }
) 
# nrow replace dim() because object is a matrix 


setClassUnion("numericOrNULL", c("numeric", "NULL"))

branchTransition <- setClass("branchTransition",
                             slots=c(tipAge="numeric",ancestorAge="numeric",
                                     statusDemes="integer",agesDemes="numeric",
                                     statusAlleles="integer",agesAlleles="numeric"))


Node <- setClass("Node",
                      contains="branchTransition",
                      slots = c(nodeNo="integer",descendant="integer")
                      )

listOfNodes <- setClass("listOfNodes",
                      contains ="list",
                      validity = function(object){
                        if (all(lapply(object,"class")=="Node")) TRUE else FALSE
                      }
                      )

spatialListOfNodes <- setClass("spatialListOfNodes",
                        contains ="listOfNodes",
                        slots=c(populationsSizes="RasterLayer"),
                        validity = function(object){
                          if (all(state(object,Inf,"Demes")%in%1:ncellA(object@populationsSizes))) TRUE else FALSE
                        }
)

# start to remove
LandGenealogy <- setClass("LandGenealogy",
                           slots=c(genealogy="Node"),
                           contains="RasterLayer")
#stop to remove

genetic <- setClass("genetic",
                    slots = c(ploidy="integer",ploidyByrow="logical"),
                    contains="data.frame",
                    prototype = prototype(data.frame(locus1.1=c(200,202),locus1.2=c(204,200)),ploidy=as.integer(2), ploidyByrow=FALSE),
                    validity = function(object){
                      if (all(grepl("ocus",names(object)))) TRUE else stop("col names of genetic data.frame do not contain 'ocus'")
                      if ((object@ploidy==2)&(object@ploidyByrow==FALSE)) {#tip.color=tipcols
                        if (length(grep("\\.1",names(object)))==0|length(grep("\\.2",names(object)))==0) {
                          if ((grep("\\.1",names(object))%%2!=1)|(grep("\\.2",names(object))%%2!=0)){
                            stop("Columns of diploid by row FALSE data frame have to be named as follows: 'c('.1','.2','.1','.2')'")
                          }
                        }
                      }
                    }
)


spatialGenetic <- setClass("spatialGenetic",
                           slots = c(x="numeric", y="numeric",Cell_numbers="integer"),
                           contains = "genetic",
                           prototype = prototype(genetic(),x=c(1,2),y=c(1,1),Cell_numbers=as.integer(c(1,2))),
                           validity = function(object){
                             if (length(object@x)!=length(object@y)) stop("slots x and y do not have the same length")
                             if (length(object@x)!=nrow(object)) stop("slots x and genetic do not have the same number of individuals")
                           }
)

GenealPopGenet <- setClass("GenealPopGenet",
                   contains="listOfNodes",
                   slots=c(Genet="spatialGenetic",Pop="RasterBrick"),

                   )

# to remove start
LandGenetGenealogy <- setClass( "LandGenetGenealogy",
                                contains = "LandGenealogy",
                                slot= c(Genotype="spatialGenetic"))
# stop to remove

genetic <- function(df=data.frame(locus1.1=c(200,202),locus1.2=c(204,200)),ploidy=NA, ploidyByrow=NA){
  if (is.na(ploidy)) { 
    if (any(grep("\\.2", colnames(df)))) 
    {
      P <- (unlist(strsplit(grep("\\.",colnames(df),value=TRUE),"\\.")))
      ploidy <- max(suppressWarnings(as.integer(P))[!is.na(suppressWarnings(as.integer(P)))])
    }
    if (any(grep("\\.2", rownames(df)))) 
    {
      P <- (unlist(strsplit(grep("\\.",rownames(df),value=TRUE),"\\.")))
      ploidy <- max(suppressWarnings(as.integer(P))[!is.na(suppressWarnings(as.integer(P)))])
    }
  }
  if (is.na(ploidyByrow)) ploidyByrow = !(any(grep(paste("\\.",ploidy,sep=""), colnames(df))))
  new("genetic",df,ploidy=ploidy,ploidyByrow=ploidyByrow)
}

spatialGenetic <- function(df=NA,x=NA,y=NA,Cell_numbers=NA)
{
  if(is.na(df)) df=cbind(data.frame(genetic()),x=c(1,2),y=c(1,1),Cell_numbers=c(1,2))
  if (!(all(c("x","y","Cell_numbers")%in%colnames(df)))){
    if (is.na(Cell_numbers)){
      if ("Cell_numbers"%in%colnames(df)) {
        Cell_numbers=df$Cell_numbers
        df <- df[,-which(colnames(df)=="Cell_numbers")]
      }
    }
    if (is.na(x)){
      if ("x"%in%colnames(df)) {
        x=df$x
        df <- df[,-which(colnames(df)=="x")]
      }
    }
    if (is.na(y)){
      if ("x"%in%colnames(df)) {
        y=df$y
        df <- df[,-which(colnames(df)=="y")]
      }
    }
  }
  new("spatialGenetic",genetic(df[,-which(colnames(df)%in%c("x","y","Cell_numbers"))]),x=df$x,y=df$y,Cell_numbers=df$Cell_numbers)
}

