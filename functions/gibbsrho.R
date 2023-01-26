#étape pour la màj de rho sachant tout le bazar, ie ici seulement L
#voir calculs ailleurs

gibbsrhobis <- function(rho, NL, tr, Prirho) {
  N = length(tr$edge.length)
  rhot = rnorm(1, rho, .00005)
  if (rhot < 0) {
    return(rho)
  } else {
    lkld1 = Prirho[[1]](rhot) - Prirho[[1]](rho)
    Lbis = sapply(1:ncol(NL), function(x) {
      sum(dpois(NL[, x], tr$edge.length[x] * rhot, log = T))
    })
    Lter = sapply(1:ncol(NL), function(x) {
      sum(dpois(NL[, x], tr$edge.length[x] * rho, log = T))
    })
    lkld2 = sum(Lbis) - sum(Lter)
    u = log(runif(1, 0, 1))
    if (u < lkld1 + lkld2) {
      return(rhot)
    } else {
      return(rho)
    }
  }
}


gibbsrho <- function(rho, NL, tr, Prirho) {
  N = length(tr$edge.length)
  rhot = rnorm(1, rho, .00005)
  if (rhot < 0) {
    return(rho)
  } else {
    lkld1 = dgamma(rhot, Prirho[1], Prirho[2], log = T) - dgamma(rho, Prirho[1], Prirho[2], log = T)
    Lbis = sapply(1:ncol(NL), function(x) {
      sum(dpois(NL[, x], tr$edge.length[x] * rhot, log = T))
    })
    Lter = sapply(1:ncol(NL), function(x) {
      sum(dpois(NL[, x], tr$edge.length[x] * rho, log = T))
    })
    lkld2 = sum(Lbis) - sum(Lter)
    u = log(runif(1, 0, 1))
    if (u < lkld1 + lkld2) {
      return(rhot)
    } else {
      return(rho)
    }
  }
}

#calcul de la vraisemblance pour rho
rhoaux <- function(x, l, rho) {
  y = rep(NA, length(x))
  qui = (x == 0)
  y[qui] = -l * rho
  if (sum(!qui) > 0) {
    y[!qui] = dexp(1 - x[!qui], rho * l, log = T) - pnorm(1, rho * l, log =
                                                            T)
  }
  return(y)
}

#mise à jour de L
gibbsL <-
  function(tr,
           M,
           X,
           Dat,
           Trposs,
           nph,
           nnodes,
           nch,
           root,
           loiini,
           bruit,
           lkldinter,
           L,
           NL,
           rho,
           ncogn,
           loiapplistetot,
           Tps) {
    Lt = L
    NLt = NL
    listc = lkldinter
    i = sample(1:length(tr$edge.length), 1)
    ncogntot = nrow(Dat[[1]])
    treli = tr$edge.length[i]
    Lte = rpois(ncogntot, rho * treli)
    Lt[Lte == 0, i] = 0
    Lt[Lte != 0, i] = rbeta(sum(Lte != 0), Lte[Lte != 0], 1)
    path = cheminrac(tr, tr$edge[i, 1], root)
    lkld = 0
    listc = list()
    loiapplistetotc = loiapplistetot
    for (x in 1:nch) {
      loiapplistetotc[[x]][[i]] = loiapp(X[[x]][[i]], Lt[Lt[, i] > 0, i], loiini[[x]], Trposs[[x]], bruit, treli, Tps[[x]][[i]])
      listc[[x]] = pruningeco(
        tr,
        X[[x]],
        M[[x]],
        Dat[[x]],
        nph[x],
        path,
        lkldinter[[x]],
        Trposs[[x]],
        Lt,
        bruit,
        loiini[[x]],
        loiapplistetotc[[x]]
      )
      #listc[[x]]=pruninglkldlistini(tr,X[[x]],M[[x]],Lt,Dat[[x]],nph[x],nnodes,root,Trposs[[x]],loiini[[x]],bruit)
      if (ncogn>0){
              lkld = lkld + rowSums(sapply(root,function(xx){(log(loiini[[x]] %*% listc[[x]][[xx]][, 1:ncogn])) -
        log(loiini[[x]] %*% lkldinter[[x]][[xx]][, 1:ncogn])}))
      } else {
        lkld=numeric()
      }
    }
    #lkld2=-dpois(Lte,rho*tr$edge.length[i],log=T)+dpois(,log=T)
    AA = as.vector((log(runif(ncogntot, 0, 1)) < c(lkld, rep(0, ncogntot -
                                                               ncogn))))
    Lres = L
    Lres[AA, i] = Lt[AA, i]
    NLt[AA, i] = Lte[AA]
    listres = lkldinter
    loiapplistetotfin = loiapplistetot
    for (x in 1:nch) {
      loiapplistetotfin[[x]][[i]] = loiapp(X[[x]][[i]], Lres[Lres[, i] > 0, i], loiini[[x]], Trposs[[x]], bruit, tr$edge.length[i], Tps[[x]][[i]])
      listres[[x]] = pruningeco(
        tr,
        X[[x]],
        M[[x]],
        Dat[[x]],
        nph[x],
        path,
        lkldinter[[x]],
        Trposs[[x]],
        Lres,
        bruit,
        loiini[[x]],
        loiapplistetotfin[[x]]
      )
      #listres[[x]]=pruninglkldlistini(tr,X[[x]],M[[x]],Lt,Dat[[x]],nph[x],nnodes,root,Trposs[[x]],loiini[[x]],bruit)
    }
    return(list(Lres, listres, NLt, loiapplistetotfin,sum(AA)/(2*ncogntot)))
  }
