prior <- setClass("prior",
                  contains="list") # add validity (contains K, R and mutation_rate)

parameters <- setClass("parameters",
                       contains="list") # add validity (contains K, R and mutation_rate)


Transition_Matrix <- setClass("Transition_Matrix",
                              slots = c(populationsize="RasterLayer"),
                              contains="matrix",
                  #  prototype = prototype(matrix(c(0.1,0.3,0.9,0.7),nrow=2,ncol=2),populationsize=raster(matrix(c(20,10)))),
                              validity = function(object){
                      if(nrow(object)!=ncellA(object@populationsize))  {stop("error in dimentions")}
                    }
) 
# nrow replace dim() because object is a matrix 

branchTransition <- setClass("branchTransition",
                             slots=c(statusDemes="integer",agesDemes="numeric",
                                     statusAlleles="integer",agesAlleles="numeric"))


Genealogy <- setClass("Genealogy",
                      contains="branchTransition",
                      slots = c(age="numeric",nodeNo="integer",descendantList="list")
                      )

listOfGenealogies <- setClass("listOfGenealogies",
                      contains ="list",
                      validity = function(object){
                        if (all(lapply(object,"class")=="Genealogy")) TRUE else FALSE
                      }
                      )


# start to remove
LandGenealogy <- setClass("LandGenealogy",
                           slots=c(genealogy="Genealogy"),
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
                           slots = c(x="numeric", y="numeric",Cell_numbers="numeric"),
                           contains = "genetic",
                           prototype = prototype(genetic(),x=c(1,2),y=c(1,1),Cell_numbers=c(1,2)),
                           validity = function(object){
                             if (length(object@x)!=length(object@y)) stop("slots x and y do not have the same length")
                             if (length(object@x)!=nrow(object)) stop("slots x and genetic do not have the same number of individuals")
                           }
)

GenealPopGenet <- setClass("GenealPopGenet",
                   contains="listOfGenealogies",
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

