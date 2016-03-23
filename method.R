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
  f = "xyFromCellA",
  signature = "RasterLayer",
  definition = function(object){
    xyFromCell(object,cellNumA(object))
  }
)
setMethod(
  f = "xyFromCellA",
  signature = "RasterStack",
  definition = function(object){
    xyFromCell(object,cellNumA(object))
  }
)
setMethod(
  f = "xyFromCellA",
  signature = "RasterBrick",
  definition = function(object){
    xyFromCell(object,cellNumA(object))
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

# copier la fonction distanceMatrix  dans la method en  l'adaptant --> idem pour transitionMatrix        


setMethod(
  f="distanceMatrixA",
  signature="RasterLayer",
  definition = function(object) {
      coords = xyFromCellA(objet, 1:length(values(object[[1]])), spatial=FALSE)
      distance = as.matrix(dist(coords)) # distance matrix of coordinates
      distance
  }
)

setMethod(
  f="migrationMatrixA",
  signature="RasterLayer",
  definition = function(object,shape,p)
)