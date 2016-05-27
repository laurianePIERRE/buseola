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
    #x=data.frame(variable=na.omit(values(object)))
    select <- !is.na(values(object))
    x=values(object)[select]
    names(x) <- which(select)
    #colnames(x)=names(object)
  x
  }
  )

setMethod(
  f = "valuesA",
  signature = "RasterBrick",  
  definition = function(object){
    x=na.omit(values(object))
    colnames(x)=names(object)
    rownames(x)=cellNumA(object)
    x
    }
)

setMethod(
  f = "valuesA",
  signature = "RasterStack",
  definition = function(object){
    x=na.omit(values(object))
    colnames(x)=names(x)
    rownames(x) <- cellNumA(object)
    x
  }
)

setMethod(
  f = "xyFromCellA",
  signature = "RasterLayer",
  definition = function(object){
    df=xyFromCell(object,cellNumA(object))
    rownames(df) <- cellNumA(object)
  }
)

setMethod(
  f = "xyFromCellA",
  signature = "RasterStack",
  definition = function(object){
    df =xyFromCell(object,cellNumA(object))
    rownames(df) <- cellNumA(object)
    df
  }
)

setMethod(
  f = "xyFromCellA",
  signature = "RasterBrick",
  definition = function(object){
    df=xyFromCell(object,cellNumA(object))
    rownames(df) <- cellNumA(object)
    df
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

    
    
# copier la fonction distanceMatrix  dans la method en  l'adaptant --> idem pour transitionMatrix        


setMethod(
  f="distanceMatrixA",
  signature="RasterLayer",
  definition = function(object) {
      coords = xyFromCellA(object)
      distance = as.matrix(dist(coords)) 
      dimnames(distance) <- list(which(!is.na(values(object))),which(!is.na(values(object))))
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
  f="sampleP",
  signature="prior",
  definition=function(prior) {
  Result=list()
  for(parametreBio in names(prior)) {
    for (variableEnvironnemental in names(prior[[parametreBio]])) {
      Result[[parametreBio]] <- list() 
      Result[[parametreBio]][[variableEnvironnemental]]$p <- switch(prior[[parametreBio]][[variableEnvironnemental]]$a$distribution,
                                                                    uniform=data.frame(variableEnvironnemental=runif(1,min=prior[[parametreBio]][[variableEnvironnemental]]$a$p[1,1],max=prior[[parametreBio]][[variableEnvironnemental]]$a$p[2,1])),
                                                                    fixed =data.frame(variableEnvironnemental=prior[[parametreBio]][[variableEnvironnemental]]$a$p),
                                                                    normal=data.frame(variableEnvironnemental=rnorm(1,mean=prior[[parametreBio]][[variableEnvironnemental]]$a$p[1,1],sd=prior[[parametreBio]][[variableEnvironnemental]]$a$p[2,1])),
                                                                    loguniform=data.frame(variableEnvironnemental=log(runif(1,min=prior[[parametreBio]][[variableEnvironnemental]]$a$p[1,1],max=prior[[parametreBio]][[variableEnvironnemental]]$a$p[2,1]))),
                                                                    uniform=data.frame(variableEnvironnemental= runif(1,min=prior[[parametreBio]][[variableEnvironnemental]]$a$p[1,1],max=prior[[parametreBio]][[variableEnvironnemental]]$a$p[2,1])),
                                                                    fixed =data.frame(variableEnvironnemental= prior[[parametreBio]][[variableEnvironnemental]]$a$p),
                                                                    normal=data.frame(variableEnvironnemental =rnorm(1,mean=prior[[parametreBio]][[variableEnvironnemental]]$a$p[1,1],sd=prior[[parametreBio]][[variableEnvironnemental]]$a$p[2,1])),
                                                                    loguniform=data.frame(variableEnvironnemental= log(runif(1,min=prior[[parametreBio]][[variableEnvironnemental]]$a$p[1,1],max=prior[[parametreBio]][[variableEnvironnemental]]$a$p[2,1]))))
      
      Result[[parametreBio]][[variableEnvironnemental]]$model <-  prior[[parametreBio]][[variableEnvironnemental]]$model
      colnames(Result[[parametreBio]][[variableEnvironnemental]]$p)=variableEnvironnemental
      rownames(Result[[parametreBio]][[variableEnvironnemental]]$p)=c("a")
      names(Result[[parametreBio]][[variableEnvironnemental]]$model)=variableEnvironnemental
      
    }
  }
  #  names(Result) <- names(prior)
  return(new("parameters",Result))
}
)

setMethod(
  f="nodes",
  signature="listOfNodes",
  definition=function(object){
    unlist(lapply(object,function(sub) sub@nodeNo))
  }
)

setMethod(
  f="currentNodes",
  signature=c("listOfNodes","numeric"),
  definition=function(object,age){
    which((sapply(object,function(object) object@ancestorAge)>=(age+1E-15))&(sapply(object,function(object) object@tipAge)<=age))
  }
)

setMethod(
  f="nodesByStates",
  signature=c("listOfNodes","numeric","character"),
  definition=function(object,age,Which){
    switch(Which,
           allByDemeAndAllele = {
             states=state(object,age,"DemesAndAlleles")
             NBS=split(states,paste(states$statusDemes,states$statusAlleles,sep="."))
             lapply(NBS,function(x) as.numeric(row.names(x)))},
           allByDeme = {
             states=state(object,age,"Demes")
             lapply(split(states,states),function(x) as.integer(names(x)))},
           allAllele = {
             states=state(object,age,"Alleles")
             lapply(split(states,states),function(x) as.integer(names(x)))},
           notAloneByDemeAndAllele = {
             states=state(object,age,"DemesAndAlleles")
             NBS=split(states,paste(states$statusDemes,states$statusAlleles,sep="."))
             whichNBS<- which(lapply(NBS,nrow)>1)
             lapply(whichNBS,function(x) as.integer(row.names(NBS[[x]])))},
           notAloneAndDemeAndAllele = {
             states=state(object,age,"DemesAndAlleles")
             NBS=split(states,paste(states$statusDemes,states$statusAlleles,sep="."))
             whichNBS<- which(lapply(NBS,nrow)>1)
             lapply(whichNBS,function(x) NBS[[x]])},
           notAloneByDeme = {
             states=state(object,age,"Demes")
             NBS=split(states,states)
             whichNBS<- which(lapply(NBS,length)>1)
             lapply(whichNBS,function(x) as.integer(names(NBS[[x]])))},
           notAloneByAllele = {
             states=state(object,age,"Alleles")
             NBS=split(states,states)
             whichNBS<- which(lapply(NBS,length)>1)
             lapply(whichNBS,function(x) as.integer(names(NBS[[x]])))},
           statesWhereCoalesce ={
             states=state(object,age,"DemesAndAlleles")
             NBS=split(states,paste(states$statusDemes,states$statusAlleles,sep="."))
             whichNBS<- which(lapply(NBS,nrow)>1)
             lapply(whichNBS,function(x) as.numeric(NBS[[x]][1,]))})
  }
  )    


setMethod(
  f="state",
  signature=c("listOfNodes","numeric","character"),
  definition=function(object,age,type){
    Nodes = currentNodes(object,age)
    switch(type,
           DemesAndAlleles = data.frame(
             statusDemes = unlist(lapply(Nodes,function(i) {
               demes <- object[[i]]@statusDemes[(object[[i]]@agesDemes<=age)]
               demes[length(demes)]})),
             statusAlleles = unlist(lapply(Nodes,function(i) {
               alleles <- object[[i]]@statusAlleles[(object[[i]]@agesAlleles<=age)]
               alleles[length(alleles)]}))),
           Demes = unlist(lapply(Nodes,function(i) {
             demes <- object[[i]]@statusDemes[(object[[i]]@agesDemes<=age)]
             demes[length(demes)]})),
           Alleles = unlist(lapply(Nodes,function(i) {
             alleles <- object[[i]]@statusAlleles[(object[[i]]@agesAlleles<=age)]
             alleles[length(alleles)]})),
           aparitionDemes = unlist(lapply(Nodes,function(i){
             apparition <- object[[i]]@agesDemes[(object[[i]]@agesDemes<=age)]
             apparition[length(apparition)]})),
           aparitionAlleles = unlist(lapply(Nodes,function(i){
             apparition <- object[[i]]@agesAlleles[(object[[i]]@agesAlleles<=age)]
             apparition[length(apparition)]})),
           tipAge =  unlist(lapply(Nodes,function(i){
             apparition <- object[[i]]@tipAge})),
           ancestorAge =  unlist(lapply(Nodes,function(i){
             apparition <- object[[i]]@ancestorAge})))
  }
)


setMethod("setState",
          signature=c("listOfNodes","integer","character","list"),
          definition=function(object,Nodes,attribut,newValues){
            switch(attribut,
                   ancestorAge={for (i in 1:length(Nodes)){
                     object[[Nodes[i]]]@ancestorAge = newValues[[i]]}},
                   tipAge={for (i in 1:length(Nodes)){
                     object[[Nodes[i]]]@tipAge = newValues[[i]]}},
                   statusDemes={for (i in 1:length(Nodes)){
                     object[[Nodes[i]]]@statusDemes = append(object[[Nodes[i]]]@statusDemes,newValues[[i]])}},
                   ageDemes={for (i in 1:length(Nodes)){
                     object[[Nodes[i]]]@ageDemes = append(object[[Nodes[i]]]@ageDemes,newValues[[i]])}},
                   statusAlleles={for (i in 1:length(Nodes)){
                     object[[Nodes[i]]]@statusAlleles = append(object[[Nodes[i]]]@statusAlleles,newValues[[i]])}},
                   ageAlleles={for (i in 1:length(Nodes)){
                     object[[Nodes[i]]]@ageAlleles = append(object[[Nodes[i]]]@ageAlleles,newValues[[i]])}},
                   nodeNo={for (i in 1:length(Nodes)){
                     object[[Nodes[i]]]@nodeNo = newValues[[i]]
                     names(object)[Nodes[i]] <- newValues[[i]]}},
                   descendants={for (i in 1:length(Nodes)){
                     object[[Nodes[i]]]@descendants = newValues[[i]]}})
          object
          }
)
            


setMethod(
  f="currentState",
  signature=c("listOfNodes","character","integer"),
  definition=function(object,type,nodes){
    switch(type,
           allStatus = data.frame(statusDemes = unlist(lapply(nodes,function(i) coalescent[[i]]@statusDemes[length(coalescent[[i]]@statusDemes)])),
                                  statusAlleles = unlist(lapply(nodes,function(i) coalescent[[i]]@statusAlleles[length(coalescent[[i]]@statusAlleles)]))),
           statusDemes = unlist(lapply(nodes,function(i) coalescent[[i]]@statusDemes[length(coalescent[[i]]@statusDemes)])),
           statusAlleles = unlist(lapply(nodes,function(i) coalescent[[i]]@statusAlleles[length(coalescent[[i]]@statusAlleles)])),
           agesDemes = unlist(lapply(nodes,function(i) coalescent[[i]]@agesDemes[length(coalescent[[i]]@agesDemes)])),
           agesAlleles = unlist(lapply(nodes,function(i) coalescent[[i]]@agesAlleles[length(coalescent[[i]]@agesAlleles)])))
  }
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

setAs("demeTransition", "matrix",
       function(from , to){
       })


setMethod(f="plot", 
          signature="demeTransition",
          definition = function(x, y , ...){
            
})

setMethod(
  f="absorbingTransitionA",
  signature="demeTransition",
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
  signature="listOfNodes",
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


setMethod(
  f="simul_coalescent",
  signature="transitionModel",
  definition = function(transitionMod)#transitionList,geneticData)
    {
      # transitionList =  list of transition matrix 
      #                   sublist demes contains list of demic transitions
      #                   sublist alleles contains list of allelic transitions for each allele
      # Ne = a data.frame with number of individuals in each deme
      # demes = all the possible deme status (attibutted cells in the raster lanscape of population sizes)
      # alleles
      # demeStatus= deme status of the nodes
      # alleleStatus = allele status of the nodes
      Ne <- round(transitionMod@Ne);Ne[Ne==0]<-1 
      coalescent <- list()
      for (i in 1:length(transitionMod@demicStatusOfStartingIndividuals)){
        coalescent[[i]] <- new("Node",nodeNo=i,descendant=integer(),new("branchTransition",tipAge=0,ancestorAge=Inf,
                                                                        statusDemes=transitionMod@demicStatusOfStartingIndividuals[i],
                                                                        agesDemes=0,
                                                                        statusAlleles=transitionMod@allelicStatusOfStartingIndividuals[i],
                                                                        agesAlleles=0))
        names(coalescent)[i]=i
      }
      numberOfNodes <- length(coalescent)
      coalescent= new("listOfNodes",coalescent)
      Age=0
      notCoalesced <- which(state(object = coalescent,age = 0,type = "ancestorAge")==Inf)
      while(length(notCoalesced)>1){
        Age=Age+1
        nodesThatCanCoalesce <- nodesByStates(object = coalescent,age = Age, Which = "notAloneByDemeAndAllele")
        nodesAndStatesThatCanCoalesce <- nodesByStates(object = coalescent,age = Age, Which = "notAloneAndDemeAndAllele")
        for (States in names(nodesAndStatesThatCanCoalesce)){
          # get get the Deme of the nodes that can coalesce in the list
          currentDeme = strsplit() nodesAndStatesThatCanCoalesce[[States]]$statusDemes[1]
          currentAllele = nodesAndStatesThatCanCoalesce[[States]]$statusAlleles[1]
          if (runif(1,0,1) < 1/(2*Ne[currentDeme])){
            numberOfNodes <- numberOfNodes+as.integer(1)
            coalescent[[numberOfNodes]] <- new("Node",tipAge=Age,ancestorAge=Inf,
                                               statusDemes=currentDeme,agesDemes=Age,
                                               statusAlleles=currentAllele,agesAlleles=Age,
                                               nodeNo=numberOfNodes,
                                               descendant=nodesThatCanCoalesce[[States]])
            names(coalescent)[numberOfNodes] <- numberOfNodes
            #lapply(coalescent ,function(x) modifyList(x,x[[x]]@ancestorAge=Age)
          }
        }
        for (node in nodes)#node = nodes[1];node = nodes[2];node = nodes[3]# parent_cell_number_of_nodes
        {
          # migrations
          parent_deme_status_of_nodes[node] = sample(demes,size=1,prob=c(transitionList$demes[as.character(deme_status_of_nodes[node]),]))
          # mutations
          parent_allele_status_of_nodes[node] = sample(alleles,size=1,prob=c(transitionList$alleles[as.character(allele_status_of_nodes[node]),]))
        }
        
        
        Node <- setClass("Node",
                         contains="branchTransition",
                         slots = c(nodeNo="integer",descendantList="list")
        )
        
      }
      }
  )
