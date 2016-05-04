setGeneric(
  name = "ncellA",
  def = function(object) { return(standardGeneric("ncellA"))})

setGeneric(
  name = "valuesA",
  def = function(object) { return(standardGeneric("valuesA"))})

setGeneric(
  name="nodes",
  def= function(object) { return(standardGeneric("nodes"))})

setGeneric(
  name="leaves",
  def= function(object) { return(standardGeneric("leaves"))})

setGeneric(
  name="statesOfLeaves",
  def= function(object,type) { return(standardGeneric("statesOfLeaves"))})

setGeneric(
  name="currentState",
  def= function(object,type,nodes) { return(standardGeneric("currentState"))})

setGeneric(
  name="varnames",
  def=function(object){return(standardgeneric("varnames"))})

setGeneric(
  name = "cellNumA",
  def = function(object) { return(standardGeneric("cellNumA"))})


setGeneric(
  name = "xyFromCellA",
  def = function(object,cellNum) { return(standardGeneric("xyFromCellA"))})

setGeneric(
  name="distA",
  def=function(object)
  { return(standardGeneric("distA"))
  }
)

setGeneric( 
  name="distanceMatrixA",
  def = function(object) { return(standardGeneric("distanceMatrixA"))}
)

setGeneric(
  name="migrationMatrixA",
  def= function(object,shapeDisp,pDisp) { return(standardGeneric("migrationMatrixA"))}
)

setGeneric(
  name="transitionMatrixA",
  def= function(object1,object2)  { return(standardGeneric("transitionMatrixA"))}
)

setGeneric(
  name="absorbingTransitionA",
  def= function(object,...)  { return(standardGeneric("absorbingTransitionA"))}
)

setGeneric(
  name="absorbingTransitionA",
  def= function(object,...)  { return(standardGeneric("absorbingTransitionA"))}
)

setGeneric(
  name="plotgenealogy",
  def=function(object,tipcols) {return(standardGeneric("plotgenealogy"))}
)

setGeneric(
  name = "plotLandG",
  def=function(object,rasK) {return(standardGeneric("plotLandG"))}
)
