
Transition_Matrix <- setClass("Transition_Matrix",
                    slots = c(transition="matrix", populationsize="raster"),
                    contains="matrix",
                    prototype = prototype(data.frame(locus1.1=c(200,202),locus1.2=c(204,200)),ploidy=as.integer(2), ploidyByrow=FALSE),
                    validity = function(object){
                      if (all(grepl("ocus",names(object)))) TRUE else stop("col names of genetic data.frame do not contain 'ocus'")
                      if ((object@ploidy==2)&(object@ploidyByrow==FALSE)) {
                        if (length(grep("\\.1",names(object)))==0|length(grep("\\.2",names(object)))==0) {
                          if ((grep("\\.1",names(object))%%2!=1)|(grep("\\.2",names(object))%%2!=0)){
                            stop("Columns of diploid by row FALSE data frame have to be named as follows: 'c('.1','.2','.1','.2')'")
                          }
                        }
                      }
                    }
)

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