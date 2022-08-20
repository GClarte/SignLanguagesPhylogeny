      
rdirichlet <- function (n, alpha) 
{
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x/as.vector(sm))
}

#mise à jour du paramètre lambda
gibbslambda <- function(X, tr, alpha, beta,tronc)
{
  n = sum(sapply(X, length))
  l = sum(tr$edge.length)
  U=runif(1,0,1)
  rr=qgamma(pgamma(tronc[1], alpha + n, rate = beta + l) + (pgamma(tronc[2], alpha + n, rate = beta + l)-pgamma(tronc[1], alpha + n, rate = beta + l))*U, alpha + n, rate = beta + l)
  return(list(rr, 4))
}

gibbslambdaold <- function(X, tr, alpha, beta)
{
  n = sum(sapply(X, length))
  l = sum(tr$edge.length)
  return(list(rgamma(1, alpha + n, rate = beta + l), 4))
}


gibbsprobas <- function(X, inithyper, nbtransf) {
  q = unlist(X)
  return(rdirichlet(1, inithyper + table(factor(
    q, levels = 1:nbtransf
  ))))
}


gibbstaux <-
  function(tr,
           Da,
           la,
           Xt,
           M,
           Trposs,
           Trpossliste,
           nph,
           root,
           loiini,
           bruit,
           Lkldinter,
           L,
           NL,
           rho,
           ncogn,
           taux,
           nch,
           hyperbeta,
           loiappliste,
           Tps) {
    tauxt = rbeta(nch, hyperbeta[1,], hyperbeta[2,])
    tauxfin = taux
    Trpossfin = Trposs
    Lkldinterfin = Lkldinter
    loiapplistefin = loiappliste
    Trposst = lapply(1:nch, function(y) {
      lapply(Trpossliste[[y]], function(x) {
        A = diag(rep(1, nph[y]))
        A[x[1], x[2]] = tauxt[y]
        A[x[1], x[1]] = 1 - tauxt[y]
        return(A)
      })
    })
    loiapplistec = lapply(1:nch, function(x) {
      lapply(1:length(tr$edge.length), function(y) {
        loiapp(Xt[[x]][[y]], L[L[, y] > 0, y], loiini[[x]], Trposst[[x]], bruit, tr$edge.length[y], Tps[[x]][[y]])
      })
    })
    Mt = lapply(1:nch, function(z) {
      lapply(1:length(tr$edge.length), function(y) {
        transitionmatrix(Xt[[z]][[y]], Trposst[[z]], nph[z], bruit, tr$edge.length[y])
      })
    })
    lkldintert = lapply(1:length(Lkldinter), function(x) {
      pruninglkldlistini(
        tr,
        Xt[[x]],
        Mt[[x]],
        L,
        Da[[x]],
        nph[x],
        length(Lkldinter),
        root[[2]],
        Trposst[[x]],
        loiini[[x]],
        bruit,
        loiapplistec[[x]]
      )
    })
    Mfin=M
    bet = sapply(1:nch,function(z){lkldcogn(list(lkldintert[[z]]),1,root,ncogn,list(loiini[[z]]))-lkldcogn(list(Lkldinter[[z]]),1,root,ncogn,list(loiini[[z]]))})
    for (i in 1:nch) {
      if (log(runif(1)) < bet[i]) {
        tauxfin[i] = tauxt[i]
        Trpossfin[[i]] = Trposst[[i]]
        Lkldinterfin[[i]] = lkldintert[[i]]
        loiapplistefin[[i]] = loiapplistec[[i]]
        Mfin[[i]]=Mt[[i]]
      }
    }
    return(list(tauxfin, Trpossfin, Lkldinterfin, loiapplistefin,Mfin))
  }

gibbsbruit <-
  function(tr,
           Da,
           la,
           X,
           M,
           Trposs,
           nph,
           root,
           loiini,
           bruit,
           Lkldinter,
           L,
           NL,
           rho,
           ncogn,
           taux,
           nch,
           hyperbruit,
           loiappliste,
           Tps) {
    Mt=M
    bruitt = hyperbruit[[2]][[1]](bruit)
      gam = hyperbruit[[2]][[2]](bruit,bruitt) - hyperbruit[[2]][[2]](bruitt,bruit) + hyperbruit[[1]][[2]](bruitt) - hyperbruit[[1]][[2]](bruit)
    if (gam>-Inf){
      loiapplistec = loiappliste
      loiapplistec = lapply(1:nch, function(x) {
        lapply(1:length(tr$edge.length), function(y) {
          loiapp(X[[x]][[y]], L[L[, y] > 0, y], loiini[[x]], Trposs[[x]], bruitt, tr$edge.length[y], Tps[[x]][[y]])
        })
      })
      Mt = lapply(1:nch, function(z) {
        lapply(1:length(tr$edge.length), function(y) {
          transitionmatrix(X[[z]][[y]], Trposs[[z]], nph[z], bruitt, tr$edge.length[y])
        })
      })
      lkldintert = lapply(1:nch, function(x) {
        pruninglkldlistini(
          tr,
          X[[x]],
          Mt[[x]],
          L,
          Da[[x]],
          nph[x],
          length(Lkldinter),
          root[[2]],
          Trposs[[x]],
          loiini[[x]],
          bruitt,
          loiapplistec[[x]]
        )
      })
      bet = lkldcogn(lkldintert,nch,root,ncogn,loiini)- lkldcogn(Lkldinter,nch,root,ncogn,loiini)
      if (!is.numeric(bet)){
        return(list(bruit, Lkldinter, loiappliste,M,-12))
      } else  {
        if (log(runif(1)) < bet+gam) {
          return(list(bruitt, lkldintert, loiapplistec,Mt,12))
        } else {
          return(list(bruit, Lkldinter, loiappliste,M,-12))
        }
      }
    } else {
      return(list(bruit, Lkldinter, loiappliste,M,-12))
    }
  }

gibbsbruitfauxpost <-
  function(tr,
           Da,
           la,
           X,
           M,
           Trposs,
           nph,
           root,
           loiini,
           bruit,
           Lkldinter,
           L,
           NL,
           rho,
           ncogn,
           taux,
           nch,
           hyperbruit,
           loiappliste,
           Tps) {
    Mt=M
    bruitt = sample(Prior$pribruit[1,],1,prob=Prior$pribruit[2,])
    loiapplistec = loiappliste
    loiapplistec = lapply(1:nch, function(x) {
      lapply(1:length(tr$edge.length), function(y) {
        loiapp(X[[x]][[y]], L[L[, y] > 0, y], loiini[[x]], Trposs[[x]], bruitt, tr$edge.length[y], Tps[[x]][[y]])
      })
    })
    Mt = lapply(1:nch, function(z) {
      lapply(1:length(tr$edge.length), function(y) {
        transitionmatrix(X[[z]][[y]], Trposs[[z]], nph[z], bruitt, tr$edge.length[y])
      })
    })
    lkldintert = lapply(1:nch, function(x) {
      pruninglkldlistini(
        tr,
        X[[x]],
        Mt[[x]],
        L,
        Da[[x]],
        nph[x],
        length(Lkldinter),
        root[[2]],
        Trposs[[x]],
        loiini[[x]],
        bruitt,
        loiapplistec[[x]]
      )
    })
    if (log(runif(1,0,1))<lkldcogn(lkldintert,nch,root,ncogn,loiini)- lkldcogn(Lkldinter,nch,root,ncogn,loiini)){
      return(list(bruitt, lkldintert, loiapplistec,Mt))
    } else {
      return(list(bruit, Lkldinter, loiappliste,M))
    }
  }

