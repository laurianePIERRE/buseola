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
  signature = "Transition_Matrix",  
  definition = function(object){
    x=na.omit(values(object@populationsize))
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
      coords = xyFromCellA(object)
      distance = as.matrix(dist(coords)) 
      distance
  }
)

setMethod(
  f="transitionMatrixA",
  signature=c("RasterLayer","prior"),
  definition = function(object1,object2){
    object2 <- sampleP(object2)
    transitionMatrixA(object1,object2)
  })

setMethod(
  f="transitionMatrixA",
  signature=c("RasterLayer","parameters"),
  definition = function(object1,object2){
  X=valuesA(object1)
  K = ReactNorm(X,object2$K[[names(object1)]]$p,object2$K[[names(object1)]]$model)[,"Y"]
  r = ReactNorm(X,object2$R[[names(object1)]]$p,object2$R[[names(object1)]]$model)[,"Y"] 
  migration <- migrationMatrixA(object1,object2$dispersion[[names(object1)]]$model, object2$dispersion[[names(object1)]]$p)
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

setMethod(
  f="varnames",
  signature="parameters",
  definition=function(object){
    names(object$K)
  })

setMethod(
  f="varnames",
  signature="prior",
  definition=function(object){
    names(object$K)
  })

setMethod(
  f="nodes",
  signature="listOfGenealogies",
  definition=function(object){
    unlist(lapply(object,function(sub) sub@nodeNo))
  }
)

setMethod(
  f="leaves",
  signature="listOfGenealogies",
  definition=function(object){
    which(lapply((lapply(object,function(x) x@descendantList)),"length")==0)  }
)

setMethod(
  f="leavesDemes",
  signature="listOfGenealogies",
  definition=function(object){
    unlist(lapply(object,function(object) object@States$demes))[leaves(object)]}
)

setMethod(
  f="migrationMatrixA",
  signature="RasterLayer",
  definition = function(object,shapeDisp,pDisp){
      distMat<-distanceMatrixA(object)
      Ndim = 1+all(ncellA(object)!=dim(object)[1:2])
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
  f="absorbingTransitionA",
  signature="matrix",
  definition = function(object,N){
    N[N<1]=1
    Ndeme <- dim(object)[1]
    # States of pairs of genes
    # N Homodemic not coalesced states
    # Heterodemic states (ij)!=(kl)
    Nhetero <- Ndeme*(Ndeme-1)/2
    # Ndeme Coalesced demic state, can be merged to 1 if we don't mind the deme
    #
    # First: 
    # Calculate the transient states (not coalesced) transition matrix
    # for transition among Q among N*(N+1)/2 not coalesced states 
    # It is composed of a submatrixes of transition between not coalesced heterodemic and
    # homodemic states
    #
    #    /          \
    #    |HeHe  HeHo|
    # Q= |HoHe  HoHo|
    #    \          /
    # where He are heterodemic states, and Ho are homodemic states
    # In submatrix HeHe:
    # lines ij are ordered for (i in 2:Ndeme) for (j in 1:(i-1))
    # columns kl are ordered for (k in 2:Ndeme) for (l in 1:(k-1))
    # In submatrix HoHe the lines are from 1 to Ndeme
    Qline <- matrix(NA,Ndeme,Ndeme) # this matrix of the same dimension as the transition matrix
    # contains the line number or column number 
    # in transient Q matrix for each pair of deme {i,j} 
    # presented as in the transition matrix
    # 
    QheteroHetero <- matrix(NA,Nhetero,Nhetero)
    ij=0;kl=0
    Check=TRUE
    for (i in 2:Ndeme){
      for(j in 1:(i-1)) {
        ij=ij+1
        Qline[i,j] <- Qline[j,i] <- ij # this matrix aims to find where the pair of demes {i,j}
        # is in Q matrix lines or columns
        #      i_j_Q_lines[ij,]<-c(i,j)
        kl=0
        for (k in 2:Ndeme){
          for (l in 1:(k-1)){
            kl=kl+1
            QheteroHetero[ij,kl] <- object[i,k]*object[j,l]
            #          Check <- Check*(P_i_k_j_l[P_ij_kl[ij,kl,"i"],P_ij_kl[ij,kl,"j"],P_ij_kl[ij,kl,"k"],P_ij_kl[ij,kl,"l"]]==Q[ij,kl])
          }
        }
      }
    }
    
    # then transient matrix 
    QhomoHetero <- matrix(NA,Ndeme,Nhetero)
    kl=0
    Check=TRUE
    for (i in 1:Ndeme){ # only homodemic sources are considered (i=j)
      Qline[i,i] <- Nhetero+i # the homodemic states are after the  heterodemic states
      # in the lines of Q matrix
      for (k in 2:Ndeme){
        for (l in 1:(k-1)){ # only heterodemic targets
          kl=kl+1
          QhomoHetero[i,kl] <- object[i,k]*object[i,l]    # i=j (homodemic sources)
          #          Check <- Check*(P_i_k_j_l[P_ij_kl[ij,kl,"i"],P_ij_kl[ij,kl,"j"],P_ij_kl[ij,kl,"k"],P_ij_kl[ij,kl,"l"]]==Q[ij,kl])
        }
      }
      kl=0
    }
    QheteroHomo <- matrix(NA,Nhetero,Ndeme)
    ij=0
    for (i in 2:Ndeme){
      for (j in 1:(i-1)){ # only heterodemic sources
        ij=ij+1
        for(k in 1:Ndeme){ # only homodemic targets
          QheteroHomo[ij,k] <- object[i,k]*object[j,k]*(1-1/(2*N[k,]))
          # homodemic targets that have not coalesced
        }
      }
    }
    QhomoHomo <- object*object*matrix(1-1/(2*N[,1]),nrow=Ndeme,ncol=Ndeme,byrow=TRUE)
    Q <- cbind(rbind(QheteroHetero,QhomoHetero),rbind(QheteroHomo,QhomoHomo))
    result=list(Q=Q,Qline=Qline)
    result
  } 
  
  )


setMethod(
  f="absorbingTransitionA",
  signature="Transition_Matrix",
  definition = function(object){
    N=valuesA(object)
    N[N<1]=1
    Ndeme <- dim(object)[1]
    # States of pairs of genes
    # N Homodemic not coalesced states
    # Heterodemic states (ij)!=(kl)
    Nhetero <- Ndeme*(Ndeme-1)/2
    # Ndeme Coalesced demic state, can be merged to 1 if we don't mind the deme
    #
    # First: 
    # Calculate the transient states (not coalesced) transition matrix
    # for transition among Q among N*(N+1)/2 not coalesced states 
    # It is composed of a submatrixes of transition between not coalesced heterodemic and
    # homodemic states
    #
    #    /          \
    #    |HeHe  HeHo|
    # Q= |HoHe  HoHo|
    #    \          /
    # where He are heterodemic states, and Ho are homodemic states
    # In submatrix HeHe:
    # lines ij are ordered for (i in 2:Ndeme) for (j in 1:(i-1))
    # columns kl are ordered for (k in 2:Ndeme) for (l in 1:(k-1))
    # In submatrix HoHe the lines are from 1 to Ndeme
    Qline <- matrix(NA,Ndeme,Ndeme) # this matrix of the same dimension as the transition matrix
    # contains the line number or column number 
    # in transient Q matrix for each pair of deme {i,j} 
    # presented as in the transition matrix
    # 
    QheteroHetero <- matrix(NA,Nhetero,Nhetero)
    ij=0;kl=0
    Check=TRUE
    for (i in 2:Ndeme){
      for(j in 1:(i-1)) {
        ij=ij+1
        Qline[i,j] <- Qline[j,i] <- ij # this matrix aims to find where the pair of demes {i,j}
        # is in Q matrix lines or columns
        #      i_j_Q_lines[ij,]<-c(i,j)
        kl=0
        for (k in 2:Ndeme){
          for (l in 1:(k-1)){
            kl=kl+1
            QheteroHetero[ij,kl] <- object[i,k]*object[j,l]
            #          Check <- Check*(P_i_k_j_l[P_ij_kl[ij,kl,"i"],P_ij_kl[ij,kl,"j"],P_ij_kl[ij,kl,"k"],P_ij_kl[ij,kl,"l"]]==Q[ij,kl])
          }
        }
      }
    }
    
    # then transient matrix 
    QhomoHetero <- matrix(NA,Ndeme,Nhetero)
    kl=0
    Check=TRUE
    for (i in 1:Ndeme){ # only homodemic sources are considered (i=j)
      Qline[i,i] <- Nhetero+i # the homodemic states are after the  heterodemic states
      # in the lines of Q matrix
      for (k in 2:Ndeme){
        for (l in 1:(k-1)){ # only heterodemic targets
          kl=kl+1
          QhomoHetero[i,kl] <- object[i,k]*object[i,l]    # i=j (homodemic sources)
          #          Check <- Check*(P_i_k_j_l[P_ij_kl[ij,kl,"i"],P_ij_kl[ij,kl,"j"],P_ij_kl[ij,kl,"k"],P_ij_kl[ij,kl,"l"]]==Q[ij,kl])
        }
      }
      kl=0
    }
    QheteroHomo <- matrix(NA,Nhetero,Ndeme)
    ij=0
    for (i in 2:Ndeme){
      for (j in 1:(i-1)){ # only heterodemic sources
        ij=ij+1
        for(k in 1:Ndeme){ # only homodemic targets
          QheteroHomo[ij,k] <- object[i,k]*object[j,k]*(1-1/(2*N[k,]))
          # homodemic targets that have not coalesced
        }
      }
    }
    QhomoHomo <- object*object*matrix(1-1/(2*N[,1]),nrow=Ndeme,ncol=Ndeme,byrow=TRUE)
    Q <- cbind(rbind(QheteroHetero,QhomoHetero),rbind(QheteroHomo,QhomoHomo))
    result=list(Q=Q,Qline=Qline)
    result
  } 
  
)

setMethod(
  f="plotgenealogy",
  signature="Genealogy",
  definition = function(object,tipcols=NA)    {
    if(is.na(tipcols)) { tipcols=1}
    plot(coalescent_2_phylog(object),direction="downward",tip.color=tipcols)
    }
)

setMethod(
  f="plotLandG",
  signature="LandGenetGenealogy",
  definition=function(object,rasK=NULL) {
    par(mfrow=c(1,2))
    plot(object)
    tipcells <- genotypes[,2][as.numeric(coalescent_2_phylog(object,direction="downward",tip.color=tip.colors)$tip.label)]
    tipcols = rainbow(ncell(rasK))[tipcells]
    plotgenealogy(object@genealogy,tip.color=tipcols)
    legend("topright", title="demes", cex=0.75, pch=16, col=tipcols[!duplicated(tipcols)], legend=tipcells[!duplicated(tipcols)], ncol=2, inset=c(legend_right_move,0))
    
    
  }
)


