
Transition_Matrix <- setClass("Transition_Matrix",
                              slots = c(populationsize="RasterLayer"),
                              contains="matrix",
                  #  prototype = prototype(matrix(c(0.1,0.3,0.9,0.7),nrow=2,ncol=2),populationsize=raster(matrix(c(20,10)))),
                              validity = function(object){
                      if(nrow(object)!=ncellA(object@populationsize))  {stop("error in dimentions")}
                    }
) 
# nrow replace dim() because object is a matrix 
