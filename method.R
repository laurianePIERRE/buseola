setMethod(
  f = "ncellA",
  signature = "RasterLayer",
  definition = function(object){
    length(na.omit(values(object)))
  }
)
setMethod(
  f = "ncellA",
  signature = "RasterStack",
  definition = function(object){
    ncellA(object[[1]])
  }
)
setMethod(
  f = "ncellA",
  signature = "RasterBrick",
  definition = function(object){
    ncellA(object[[1]])
  }
)

setMethod(
  f = "valuesA",
  signature = "RasterLayer",
  definition = function(object){
    x=data.frame(variable=na.omit(values(object)))
    colnames(x)=names(object)
  x
  }
  )

setMethod(
  f = "valuesA",
  signature = "RasterBrick",
  definition = function(object){
    x=na.omit(values(object))
    colnames(x)=names(object)
    x
    }
)


    
    
setMethod(
      f = "valuesA",
      signature = "RasterStack",
      definition = function(object){
        x=na.omit(values(object))
        colnames(x)=names(object)
        x
      }
)
        