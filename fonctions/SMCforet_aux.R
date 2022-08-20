

pasparalbruit <- function(etats, pas, Dat, Prior, Param, ncores, ncogn,fin) {
  VV = mclapply(1:length(etats), function(x) {
    gibbsbr2(Dat, ncogn, c(1, pas), Prior, Param, etats[[x]],fin)[[3]]
  }, mc.cores = ncores)
  return(VV)
}


initreeold <- function(n,tiplabel,tipprior,cladeage){
  if (length(tipprior)==0){
    return(rcoal(n,tiplabel))
  } else {
    corr=0
    cpt=n+2
    act=1:n
    edge=matrix(ncol=2,nrow=0)
    edge.length=numeric(0)
    haut=numeric(2*n-1)
    for (i in 1:length(tipprior)){
      act=act[-which(act%in%tipprior[[i]])]
    }
    for (i in 1:length(tipprior)){
      actt=tipprior[[i]]
      for (j in 1:(length(tipprior[[i]])-1)){
        qui=sample(actt,2)
        corr=corr-log(length(actt))
        actt=actt[-which(actt%in%qui)]
        x=rexp(1)
        edge=rbind(edge,c(cpt,qui[1]))
        edge=rbind(edge,c(cpt,qui[2]))
        edge.length=c(edge.length,x-haut[qui]+max(haut[qui]))
        haut[cpt]=max(haut[qui])+x
        actt=c(actt,cpt)
        cpt=cpt+1
        corr=corr+dexp(x,log=T)
      }
      act=c(act,actt)
    }
    while(length(act)>2){
      qui=sample(act,2)
      corr=corr-log(length(act))
      act=act[-which(act%in%qui)]
      x=rexp(1)
      edge=rbind(edge,c(cpt,qui[1]))
      edge=rbind(edge,c(cpt,qui[2]))
      edge.length=c(edge.length,x-haut[qui]+max(haut[qui]))
      haut[cpt]=max(haut[qui])+x
      act=c(act,cpt)
      cpt=cpt+1
      corr=corr+dexp(x,log=T)
    }
    x=rexp(1)
    edge=rbind(edge,c(n+1,act[1]))
    edge=rbind(edge,c(n+1,act[2]))
    edge.length=c(edge.length,x-haut[act]+max(haut[act]))
    phy <- list(edge = edge, edge.length = edge.length)
    phy$tip.label <- tiplabel
    phy$Nnode <- n - 1
    class(phy) <- "phylo"
    phy$root=list(edge[edge[,1]==n+1,2],n+1)
    return(list(phy,corr))
  }
}


initree <- function(n,tiplabel,tipprior,cladeage){
  if (length(tipprior)==0){
    return(rcoal(n,tiplabel))
  } else {
    corr=0
    cpt=n+2
    act=1:n
    edge=matrix(ncol=2,nrow=0)
    edge.length=numeric(0)
    haut=numeric(2*n-1)
    ancetresclades=c()
    for (i in 1:length(tipprior)){
      if (length(act)>1){
        act=act[-which(act%in%tipprior[[i]])]
      } else {
        if (act%in%tipprior[[i]]){
          act=numeric(0)
        }
      }
    }
    for (i in 1:length(tipprior)){
      actt=tipprior[[i]]
      if(i>1){
        for (k in 1:(i-1)){
          if(all(tipprior[[k]]%in%tipprior[[i]])){
            actt=actt[-which(actt%in%tipprior[[k]])]
            actt=c(actt,ancetresclades[k])
            act=act[-which(act==ancetresclades[k])]
          } else if (any(tipprior[[k]]%in%tipprior[[i]])){
            stop("your constraints are not compatible")
          }
        }
        
      }
      for (j in 1:(length(actt)-1)){
        qui=sample(actt,2)
        corr=corr-log(length(actt))
        actt=actt[-which(actt%in%qui)]
        x=rexp(1)
        edge=rbind(edge,c(cpt,qui[1]))
        edge=rbind(edge,c(cpt,qui[2]))
        edge.length=c(edge.length,x-haut[qui]+max(haut[qui]))
        haut[cpt]=max(haut[qui])+x
        actt=c(actt,cpt)
        cpt=cpt+1
        corr=corr+dexp(x,log=T)
      }
      act=c(act,actt)
      ancetresclades=c(ancetresclades,actt)
    }
    while(length(act)>2){
      qui=sample(act,2)
      corr=corr-log(length(act))
      act=act[-which(act%in%qui)]
      x=rexp(1)
      edge=rbind(edge,c(cpt,qui[1]))
      edge=rbind(edge,c(cpt,qui[2]))
      edge.length=c(edge.length,x-haut[qui]+max(haut[qui]))
      haut[cpt]=max(haut[qui])+x
      act=c(act,cpt)
      cpt=cpt+1
      corr=corr+dexp(x,log=T)
    }
    x=rexp(1)
    edge=rbind(edge,c(n+1,act[1]))
    edge=rbind(edge,c(n+1,act[2]))
    edge.length=c(edge.length,x-haut[act]+max(haut[act]))
    phy <- list(edge = edge, edge.length = edge.length)
    phy$tip.label <- tiplabel
    phy$Nnode <- n - 1
    class(phy) <- "phylo"
    phy$root=list(edge[edge[,1]==n+1,2],n+1)
    return(list(phy,corr))
  }
}

rcoal <- function (n, tip.label = NULL, br = "coalescent", ...) 
{
  #copy paste of the rcoal function of phytools
  n <- as.integer(n)
  nbr <- 2 * n - 2
  edge <- matrix(NA, nbr, 2)
  x <- 2 * rexp(n - 1)/(as.double(n:2) * as.double((n - 1):1))
  corr=sum(dexp(x,log=T))
  if (n == 2) {
    edge[] <- c(3L, 3L, 1:2)
    edge.length <- rep(x, 2)
  } else if (n == 3) {
    edge[] <- c(4L, 5L, 5L, 4L, 5L, 1:3)
    edge.length <- c(x[c(2, 1, 1)], sum(x))
  } else {
    edge.length <- numeric(nbr)
    h <- numeric(2 * n - 1)
    node.height <- cumsum(x)
    pool <- 1:n
    nextnode <- 2L * n - 1L
    for (i in 1:(n - 1)) {
      y <- sample(pool, size = 2)
      corr=corr+log(1/length(pool))
      ind <- (i - 1) * 2 + 1:2
      edge[ind, 2] <- y
      edge[ind, 1] <- nextnode
      edge.length[ind] <- node.height[i] - h[y]
      h[nextnode] <- node.height[i]
      pool <- c(pool[!pool %in% y], nextnode)
      nextnode <- nextnode - 1L
    }
  }
  phy <- list(edge = edge, edge.length = edge.length)
  if (is.null(tip.label)) 
    tip.label <- paste("t", 1:n, sep = "")
  phy$tip.label <- tip.label
  phy$Nnode <- n - 1
  class(phy) <- "phylo"
  phy$root=list(edge[edge[,1]==n+1,2],n+1)
  return(list(phy,corr))
}


iniforet <- function(nfeuilles, agemaxt, tiplabel,tipprior,cladeage) {
  #initialises a forest of trees
  u = sample(1:nfeuilles, 2)
  VV=initree(nfeuilles,tiplabel,tipprior)
  tr=VV[[1]]
  corr=VV[[2]]
  if (length(agemaxt)==2){
    agemax=rgamma(1,agemaxt[1],agemaxt[2])
    corr=dgamma(agemax,agemaxt[1],agemaxt[2],log=T)
  } else {
    agemax=agemaxt
  }
  tr$edge.length=tr$edge.length/hauteur(tr,tr$root[[2]])*agemax
  while (!all(sapply(tipprior,function(tip){is.monophyletic.perso(tr,tip)})) || 
         !all(sapply(cladeage,function(clade){ageconstraints(tr,clade)})) ){
    VV=initree(nfeuilles,tiplabel,tipprior)
    tr=VV[[1]]
    corr=VV[[2]]
    #print(i)
    tr$edge.length=tr$edge.length/hauteur(tr,tr$root[[2]])*agemax
    #print("zut")
  }
  tr$edge.length=tr$edge.length/hauteur(tr,tr$root[[2]])*agemax
  return(list(tr,corr))
}


iniedges <- function(nfeuil){
  #initialise
  N=nfeuil+2
  actifs=1:nfeuil
  br=matrix(ncol=2,nrow=0)
  while(length(actifs)>2){
    t=sample(1:length(actifs),2)
    br=rbind(br,c(N,actifs[t[1]]),c(N,actifs[t[2]]))
    actifs=actifs[-t]
    actifs=c(actifs,N)
    N=N+1
  }
  br=rbind(br,c(nfeuil+1,actifs[1]),c(nfeuil+1,actifs[2]))
  return(br)
}
  

forward.bruit <- function(state, Dat, Prior, Param,j) {
  statet = state
  nph = Param$nph
  agemax = Param$agemax
  tipprior=Prior$tipprior
  bruit=Prior$bruittemp[j]
  statet$bruit=bruit
  taux=state$taux
  Prirho = Prior$Prirho
  Pri = Prior$Prila
  hyperpbini = Prior$hyperpbini
  Trposs = lapply(1:Param$nch, function(y) {
    lapply(Param$Trpossliste[[y]], function(x) {
      A = diag(rep(1, nph[y]))
      A[x[1], x[2]] = taux[y]
      A[x[1], x[1]] = 1 - taux[y]
      return(A)
    })
  })
  loiapplistetot = lapply(1:nch, function(x) {
    lapply(1:nrow(statet$Tr$edge), function(y) {
      loiapp(statet$X[[x]][[y]],
             statet$L[statet$L[, y] > 0, y],
             Param$loiini[[x]],
             Trposs[[x]],
             bruit,
             statet$Tr$edge.length[y],
             statet$Tps[[x]][[y]])
    })
  })
  statet$loiapptot = loiapplistetot
  M = lapply(1:nch, function(z) {
    lapply(1:length(statet$Tr$edge.length), function(y) {
      transitionmatrix(statet$X[[z]][[y]], Trposs[[z]], nph[z], bruit, statet$Tr$edge.length[y])
    })
  })
  statet$Lin = lapply(1:nch, function(x) {
    pruninglkldlistini(
      statet$Tr,
      statet$X[[x]],
      M[[x]],
      statet$L,
      Dat[[x]],
      nph[x],
      max(statet$Tr$edge),
      statet$Tr$root[[2]],
      Trposs[[x]],
      loiini[[x]],
      bruit,
      loiapplistetot[[x]]
    )
  })
  statet$Ltp = lkldtopotout(statet$Tr, statet$X, statet$La, statet$NL, statet$rho, statet$L,statet$P)
  weight = Prior$priortree(statet$Tr, agemax) - Prior$priortree(state$Tr, agemax) + #les lkld des variables latentes se simplifient entre elles ou avec le prior
    lkldcogn(statet$Lin,
             Param$nch,
             statet$Tr$root,
             nrow(Dat[[1]]),
             Param$loiini) -
    lkldcogn(state$Lin, Param$nch, state$Tr$root, nrow(Dat[[1]]), Param$loiini)
  statet$weight = weight + state$weight
  return(statet)
}

initialisation1partbruit <- function(Dat, Param, Prior) {
  #itiniailisation pour une seule particule
  nfeuilles = ncol(Dat[[1]])
  #ZZ = iniforet(nfeuilles, Param$agemax, Param$tiplabel,Prior$tipprior,Param$Cladeage)
  ZZ=iniaveccontraintes(Param,Prior)
  trt=ZZ[[1]]
  k = nrow(Dat[[1]])
  nph = Param$nph
  nch = Param$nch
  bruit=Prior$bruittemp[1]
  if (is.matrix(Prior$pribeta)) {
    #si Pribeta est une matrice ce sont les paramètres de la loi beta
    taux = rbeta(nch, Prior$pribeta[1,], Prior$pribeta[2,])
  } else {
    #sinon ce sont les paramètres fixés.
    taux = Prior$pribeta
  }
  Prirho = Prior$Prirho
  Pri = Prior$Prila
  hyperpbini = Prior$hyperpbini
  Trposs = lapply(1:Param$nch, function(y) {
    lapply(Param$Trpossliste[[y]], function(x) {
      A = diag(rep(1, nph[y]))
      A[x[1], x[2]] = taux[y]
      A[x[1], x[1]] = 1 - taux[y]
      return(A)
    })
  })
  nnodes = max(trt$edge)
  loiini = Param$loiini
  root = trt$root[[2]]
  Probtransf = lapply(hyperpbini, function(z) {
    rdirichlet(1, z)
  })
  la = sapply(1:nch, function(x) {
    rr=rgamma(1, Pri[[x]][1], Pri[[x]][2])
    while(rr>Pri[[x]][4] || rr<Pri[[x]][3]){
      rr=rgamma(1, Pri[[x]][1], Pri[[x]][2])
    }
    return(rr)
  })
  X = lapply(1:nch, function(x) {
    transf(trt, Trposs[[x]], la[x], Probtransf[[x]])
  })
  if(is.numeric(Prirho)){
    rho = rgamma(1,Prirho[1],Prirho[2])
    corrrho=1
  } else {
    rho = Prirho[[2]][[1]]()
    corrrho=Prirho[[1]](rho)-Prirho[[2]][[2]](rho)
  }
  NL = matrix(rpois(nrow(trt$edge) * k, rho * trt$edge.length),
              ncol = nrow(trt$edge),
              byrow = T)
  L = NL
  Tps = lapply(X, function(x) {
    lapply(x, function(y) {
      Q = length(y)
      if (Q > 0) {
        Y = cumsum(rexp(Q + 1, 1))
        return(Y[1:Q] / Y[Q + 1])
      } else{
        return(numeric(0))
      }
    })
  })
  L[NL != 0] = rbeta(length(which(NL != 0)), L[NL != 0], 1)
  loiapplistetot = lapply(1:nch, function(x) {
    lapply(1:length(trt$edge.length), function(y) {
      loiapp(X[[x]][[y]], L[L[, y] > 0, y], loiini[[x]], Trposs[[x]], bruit, trt$edge.length[y], Tps[[x]][[y]])
    })
  })
  M = lapply(1:nch, function(z) {
    lapply(1:length(trt$edge.length), function(y) {
      transitionmatrix(X[[z]][[y]], Trposs[[z]], nph[z], bruit, trt$edge.length[y])
    })
  })
  lkldinter = lapply(1:nch, function(x) {
    pruninglkldlistini(trt,
                       X[[x]],
                       M[[x]],
                       L,
                       Dat[[x]],
                       nph[x],
                       nnodes,
                       root,
                       Trposs[[x]],
                       loiini[[x]],
                       bruit,
                       loiapplistetot[[x]])
  })
  weight = lkldcogn(lkldinter, Param$nch, trt$root, Param$batches[[1]][1], Param$loiini) +
    log(nfeuilles * (nfeuilles - 1))+Prior$priortree(trt, Param$agemax) + corrrho +ZZ[[2]]
  return(
    list(
      X = X,
      L = L,
      La = la,
      Lin=lkldinter,
      Tr = trt,
      Ltp = lkldtopotout(trt, X, la, NL, rho, L,Probtransf),
      P = Probtransf,
      Rho = rho,
      NL = NL,
      taux = taux,
      bruit = bruit,
      Tps = Tps,
      weight = weight
    )
  )
}

inipopbruit <- function(npart, Dat, Param, Prior,ncores) {
  etats = mclapply(1:npart, function(x) {
    initialisation1partbruit(Dat, Param, Prior)
  }, mc.cores = ncores)
}


nearestcommonancestor = function(tr,tip){
  #if tip has length 1, it return the only ancestor, can be used to add contraints of the form 
  #"we know that at this date this language existed"
  if (length(tip)>1){
    k=lapply(tip,function(x){ancetres(tr,x)})
    rootmax=sapply(k,function(x){x[1]})
    W=unique(rootmax)
    ll=sapply(k,length)
    i=0
    while ((i)<min(ll) & all(sapply(k,function(x){x[i+1]})==k[[1]][i+1])){
      i=i+1
    }
    rr=k[[1]][i]
    return(rr)
  } else {
    return(tr$edge[tr$edge[,2]==tip,1])
  }
}

autreini <- function(tr,cladage,corr,agemaxt){
  trt=tr
  corrb=corr
  if (length(cladage)==0){
    agemax=rgamma(1,agemaxt[1],agemaxt[2])
    tr$edge.length=tr$edge.length/hauteur(tr,tr$root[[2]])*agemax
    return(list(tr,corr))
  } else {
      ZUT=T
      if (length(agemaxt)==2){
        QQQQ=sapply(cladage,function(x){nearestcommonancestor(tr,x[[1]])})
        qqqq=which(QQQQ==tr$root[[2]])
        cstrroot= sapply(qqqq,function(zzz){c(cladage[[zzz]][[2]][2],cladage[[zzz]][[2]][1])})
        if (length(cstrroot)>0){
          rootmin=max(cstrroot[2,])
          rootmax=min(cstrroot[1,])
        } else {
          rootmin=0
          rootmax=Inf
        }
        if (any(QQQQ==tr$root[[2]])){
          agemax=rgamma(1,agemaxt[1],agemaxt[2])
          if (rootmax>rootmin){
            while(!(agemax<rootmax & agemax>rootmin)){
              Fmin=pgamma(rootmin,agemaxt[1],agemaxt[2])
              Fmax=pgamma(rootmax,agemaxt[1],agemaxt[2])
              uuu=runif(1,min=Fmin,max=Fmax)
              agemax=qgamma(uuu,agemaxt[1],agemaxt[2])
            }
            corrb=corrb+dgamma(agemax,agemaxt[1],agemaxt[2],log=T)
            
          } else {
            return(list(tr,NaN))
          }
        } else {
          agemax=rgamma(1,agemaxt[1],agemaxt[2])
          while(agemax<min(sapply(cladage,function(x){x[[2]][1]}))){
            agemax=rgamma(1,agemaxt[1],agemaxt[2])
          }
          corrb=corrb+dgamma(agemax,agemaxt[1],agemaxt[2],log=T)
        }

      } else {
        agemax=agemaxt
      }
      hauttous=rep(NA,max(tr$edge))
      hauttous[1:length(tr$tip.label)]=0
      hauttous[tr$root[[2]]]=agemax
      parcours=parcoursdepth(tr,tr$root[[2]])
      quib=sapply(cladage,function(x){y=nearestcommonancestor(tr,x[[1]]);return(c(y,x[[2]][1],min(x[[2]][2],agemax)))})
      XX=unique(quib[1,])
      qui=matrix(nrow=3,ncol=length(XX))
      for (jj in 1:length(XX)){
        uu=which(quib[1,]==XX[jj])
        qui[1,jj]=XX[jj]
        qui[2,jj]=max(quib[2,uu])
        qui[3,jj]=min(quib[3,uu])
      }
        for (ii in intersect(parcours,qui[1,])){
          i=which(qui[1,]==ii)
          if (qui[1,i]!=tr$root[[2]]){
            hmin=max(qui[2,i],hauttous[intersect(tr$edge[sousarbre(tr,qui[1,i]),2],qui[1,])])
            part=i
            parts=numeric()
            while(part!=tr$root[[2]]){
              part=tr$edge[tr$edge[,2]==part,1]
              parts=c(parts,part)
            }
            hmax=min(qui[3,i],qui[3,qui[1,]%in%parts])
            if (hmin>hmax){
              corrb=NaN
              return(list(tr,corrb))
            } else {
              hauttous[qui[1,i]]=runif(1,hmin,hmax)
              corrb=corrb+log(1/(hmax-hmin))
            }
          }
        }
      j=1
      while (j<=length(parcours) & ZUT){
        u=parcours[j]
        if (is.na(hauttous[u])){
          part=tr$edge[tr$edge[,2]==u,1]
          cpt=1
          while (is.na(hauttous[part])){
            part=tr$edge[tr$edge[,2]==part,1]
            cpt=cpt+1
          }
          hhh=rbeta(1,1,cpt)
          hauttous[u]=max(hauttous[tr$edge[tr$edge[,1]==u,2]])+(hauttous[part] -max(hauttous[tr$edge[tr$edge[,1]==u,2]]))*hhh
          ZUT=ZUT&hauttous[u]>0
          corrb=corrb+dbeta(hhh,1,cpt,log=T)
        }
        j=j+1
      }
      if (ZUT){
        for (i in 1:length(trt$edge.length)){
          trt$edge.length[i]=hauttous[trt$edge[i,1]]-hauttous[trt$edge[i,2]]
        }
      }
    return(list(trt,corrb))
  }
}

parcoursdepth <- function(tr,j){
  if (length(which(tr$edge[,1]==j))==0){
    return(j)
  } else {
    u=which(tr$edge[,1]==j)
    return(c(parcoursdepth(tr,tr$edge[u[1],2]),parcoursdepth(tr,tr$edge[u[2],2]),j))
  }
}


iniaveccontraintes <- function(Param,Prior){
  ZUT=F
  while (!ZUT){
    V=initree(length(Param$tiplabel),Param$tiplabel,Prior$tipprior,Param$Cladeage)
    VV=autreini(V[[1]],Param$Cladeage,V[[2]],Param$agemax)
    ZUT=(all(VV[[1]]$edge.length>0))&!is.na(VV[[2]]) & 
      all(sapply(Prior$tipprior,function(tip){is.monophyletic.perso(VV[[1]],tip)})) &
      all(sapply(Param$Cladeage,function(clade){ageconstraints(VV[[1]],clade)}))
  }
  return(VV)
}

