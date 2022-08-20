
SMCbruit <- function(Dat, Param, Prior, ncores) {
  #avec le bruit qui diminue
  npart=Param$npart
  npas=Param$npas
  npasfin=Param$npasfin
  Nmin=Param$Nmin
  batches = Param$batches
  nfeuilles = ncol(Dat[[1]])
  state = inipopbruit(npart, Dat, Param, Prior,ncores)
  #for (i in 2:(nfeuilles - 1)) {
  #  state = lapply(state, function(x) {
  #    forward.topo(x, Dat, Prior, Param)
  #  })}
  pds = sapply(state, function(x) {
    x$weight
  })
  pds = pds - max(pds)
  pds = pds - log(sum(exp(pds)))
  pdshist = pds
  pds[is.na(pds)] = -Inf
  histgeneal = matrix(ncol = 1, nrow = npart)
  if (1 / (exp(max(pds)) ^ 2 * sum(exp(pds - max(pds)) ^ 2)) < Nmin) {
    state = lapply(sample(
      1:npart,
      npart,
      replace = T,
      prob = exp(pds - max(pds))
    ), function(x) {
      u = state[[x]]
      u$weight = -log(npart)
      return(u)
    })
    pds = rep(-log(npart), npart)
  }
  state = pasparalbruit(state, npas, Dat, Prior, Param, ncores,  nrow(Dat[[1]]),F)
  if (!(all(sapply(state,length)==length(state[[1]])))){
    return(state)
  }
  for (i in 2:length(Prior$bruit)) {
    state = lapply(state, function(x) {
      forward.bruit(x, Dat, Prior, Param,i)
    })
    pds = sapply(state, function(x) {
      x$weight
    })
    pds = pds - max(pds)
    pds = pds - log(sum(exp(pds)))
    pdshist = rbind(pdshist, pds)
    if (1 / (exp(max(pds)) ^ 2 * sum(exp(pds - max(pds)) ^ 2)) < Nmin) {
      quisample = sample(1:npart,
                         npart,
                         replace = T,
                         prob = exp(pds - max(pds)))
      histgeneal = cbind(histgeneal, quisample)
      state = lapply(quisample, function(x) {
        u = state[[x]]
        u$weight = -log(npart)
        return(u)
      })
      pds = rep(-log(npart), npart)
    }
    state = pasparalbruit(state, npas, Dat, Prior, Param, ncores, nrow(Dat[[1]]),F)
    if (!(all(sapply(state,length)==length(state[[1]])))){
      return(state)
    }
  }
  state = pasparalbruit(state, npasfin, Dat, Prior, Param, ncores, nrow(Dat[[1]]),T)
  if (!(all(sapply(state,length)==length(state[[1]])))){
    return(state)
  }
  #state=lapply(state,function(x){y=x;y$Tr$root=x$Tr$root[[2]];return(y)})
  statefinresamp = state
  renvoi = lapply(statefinresamp, function(x) {
    Z=sapply(Param$ReconstructClade,function(zz){nearestcommonancestor(x$Tr,zz)})
    list(Tr=x$Tr, La=x$La, P=x$P, Rho=x$Rho, X=x$X, Tps=x$Tps, beta=x$taux, noise=x$bruit,L=x$L, NL=x$NL, weight=x$weight,
         RootState=lapply(1:nch,function(z){x$Lin[[z]][[x$Tr$root[[2]]]]}),
         OtherReconstruct=lapply(Z,function(zz){lapply(1:nch,function(z){x$Lin[[z]][[zz]]})}))
  })
  renvoi=lapply(renvoi,function(x){y=x;y$Tr$root=y$Tr$root[[2]];return(y)})
  return(list(renvoi, pdshist, histgeneal))
}
