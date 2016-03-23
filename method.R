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
  f = "cellNumA",
  signature = "RasterLayer",
  definition = function(object){
    which(!is.na(values(object)))
  }
)

setMethod(
  f = "cellNumA",
  signature = "RasterStack",
  definition = function(object){
    cellNumA(object[[1]])
  }
)

setMethod(
  f = "cellNumA",
  signature = "RasterBrick",
  definition = function(object){
    cellNumA(object[[1]])
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
      distance = as.matrix(dist(coords)) 
      distance
  }
)

setMethod(
  f="migrationMatrixA",
  signature="RasterLayer",
  definition = function(object,shapeDisp,pDisp){
      distMat<-distanceMatrixA(object)
      Ndim = 1+all(ncellA(object)!=dim(o)[1:2])
      migration = apply(distMat, c(1,2), 
                        function(x)(switch(shapeDisp,
                                           # 1: alphaDisp   2: betaDisp ; note: esperance = 1/alphaDisp
                                           fat_tail1 = 1/(1+x^pDisp[,2]/pDisp[,1]), # Molainen et al 2004
                                           # 1: sigmaDisp
                                           gaussian = (dnorm(x, mean = 0, sd = pDisp[,1], log = FALSE)),
                                           # 1: sigma Disp
                                           exponential = (dexp(x, rate = 1/pDisp[,1], log = FALSE)),
                                           # 1: sigmaDisp
                                           contiguous = (x==0)*(1-pDisp[,1])+((x>0)-(x>1.4*res(object)[1]))*(pDisp[,1]/(2*Ndim)),
                                           # 1: sigmaDisp
                                           contiguous8 = (x==0)*(1-pDisp[,1])+((x>0)-(x>2*res(object)[1]))*(pDisp[,1]/(4*Ndim)),
                                           #island = (x==0)*(1-pDisp[1])+(x>0)*(pDisp[1]/(nCell-1)),
                                           island = (x==0)*(1-pDisp[,1])+(x>0)*(pDisp[,1]),
                                           #1: sigmaDisp    2: gammaDisp
                                           fat_tail2 = x^pDisp[,2]*exp(-2*x/(pDisp[,1]^0.5)),
                                           #1: sigmaDisp, 2: 
                                           contiguous_long_dist_mixt = pDisp["plongdist"]/ncellA(object)+(x==0)*(1-pDisp["pcontiguous"]-pDisp["plongdist"])+((x>0)-(x>1.4*res(rasterStack)[1]))*(pDisp["pcontiguous"]/2),
                                           gaussian_long_dist_mixt = pDisp[,2]/ncellA(object) + (dnorm(x, mean = 0, sd = pDisp[,1], log = FALSE))
                        )))
      migration
    }
    
    

)

setMethod(
  f="transitionMatrixA",
  signature="RasterLayer",
  definition = function(object,prior){
      listeSample=sampleP(prior)
      X=valuesA(object)
      K = ReactNorm(X,listeSample$K$busseola$p,listeSample$K$busseola$model)[,"Y"]
      r = ReactNorm(X,listeSample$R$busseola$p,listeSample$R$busseola$model)[,"Y"] 
      migration <- migrationMatrixA(object,listeSample$dispersion$busseola$model, listeSample$dispersion$busseola$p)
      if ((length(r)==1)&(length(K)==1)){transition = r * K * t(migration)}
      if ((length(r)>1)&(length(K)==1)){transition = t(matrix(r,nrow=length(r),ncol=length(r))) * K * t(migration)}
      if ((length(r)==1)&(length(K)>1)){transition = r * t(matrix(K,nrow=length(K),ncol=length(K))) * t(migration)}
      if ((length(r)>1)&(length(K)==1)){transition = t(matrix(r,nrow=length(r),ncol=length(r))) * K * t(migration)}
      if ((length(r)>1)&(length(K)>1)) {transition = t(matrix(r,nrow=length(r),ncol=length(r))) * t(matrix(K,nrow=length(K),ncol=length(K))) * t(migration)}
      transition = transition / rowSums(transition)  # Standardisation
      transition = transition * as.numeric(transition>1e-6) # removal of values below 1e-6
      transition = transition / rowSums(transition)  # Standardisation again
      transition
    }

  )

