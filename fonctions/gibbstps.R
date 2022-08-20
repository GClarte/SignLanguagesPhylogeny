

gibbstemps <-
  function(tr,
           Da,
           La,
           X,
           M,
           Trposs,
           nph,
           root,
           loiini,
           bruit,
           Lkldinter,
           L,
           pr,
           ncogn,
           NL,
           rho,
           priortree,
           loiappliste,
           Tps,
           agemax){
    #modifie l'âge des branches
    #recommence entièrement le calcul de lkld.
    a=nrow(tr$edge)
    if(length(agemax)==1){
      qq=sample(2:(2*a-1),1)
    } else {
      qq=sample(0:(2*a-1),1)
    }
    if (qq==1){
      VV=modifrootage(tr,agemax)
      trc = VV[[1]]
      corr=VV[[2]]
    } else if (qq==0) {
      VV=modifresc(tr,agemax)
      trc = VV[[1]]
      corr=VV[[2]]
    } else if (qq <a+1){
      trc = modif(tr,agemax)
      corr=0
    } else {
      VV=modifst(tr)
      trc=VV[[1]]
      corr=VV[[2]]
    }
    loiapplistec = lapply(1:nch, function(x) {
      lapply(1:length(trc$edge.length), function(y) {
        loiapp(X[[x]][[y]], L[L[, y] > 0, y], loiini[[x]], Trposs[[x]], bruit, trc$edge.length[y], Tps[[x]][[y]])
      })
    })
    lkldintert = lapply(1:length(Lkldinter), function(x) {
      pruninglkldlistini(
        trc,
        X[[x]],
        M[[x]],
        L,
        Da[[x]],
        nph[x],
        length(Lkldinter),
        root,
        Trposs[[x]],
        loiini[[x]],
        bruit,
        loiapplistec[[x]]
      )
    })
    bet = lkldcogn(lkldintert,nch,trc$root,ncogn,loiini)- lkldcogn(Lkldinter,nch,tr$root,ncogn,loiini)
    u = log(runif(1, 0, 1))
    if (u < lkldtopotout(trc, X, La, NL, rho, L,pr) - lkldtopotout(tr, X, La, NL, rho, L,pr) +
        priortree(trc,agemax) - priortree(tr,agemax) + bet + corr)
    {
      return(list(trc, lkldintert, 5, loiapplistec))
    } else {
      return(list(tr, Lkldinter, -5, loiappliste))
    }
  }


lkldtopo <- function(tr, X, lambda)
{
  #calcule la vrais poissonienne associée à une longueur de branches
  #Non utilisée
  S = sapply(X, length)
  return(sum(dpois(S, lambda * tr$edge.length, log = T)))
}

#fonction pour la muse à jour de tout, on n'a pas besoin de p puisqu'il n'intervient que dans la màj de X et pas avec la longuguer des branches
lkldtopotout <- function(tr, X, lambda, NL, rho, L,P)
{
  #calcul de la vraisemblance des transf et apparitions. Notons que p n'apparait pas, puisqu'il se simplifie toujours
  r = 0
  for (i in 1:length(X)) {
    S = sapply(X[[i]], length)
    QQ=sum(sapply(X[[i]],function(x){sum(log(P[[i]][x]))}))
    r = r + sum(dpois(S, lambda[i] * tr$edge.length, log = T)) + QQ
  }
  r = r + sum(sapply(1:length(tr$edge.length), function(x) {
    sum(dpois(NL[, x], rho * tr$edge.length[x], log = T))
  }))
  r=r + sum(dbeta(L[L>0],NL[L>0],1,log=T))
  return(r)
}

lkldtopotransf <- function(tr, X, lambda,P)
{
  #calcule la vaisemblance des transformations, non utilisée
  r = 0
  QQ=sum(sapply(X[[i]],function(x){sum(log(P[x]))}))
  for (i in 1:length(X)) {
    S = sapply(X[[i]], length)
    r = r + sum(dpois(S, lambda[i] * tr$edge.length[i], log = T)) + QQ
  }
  return(r)
}



modif <- function(tra,agemax)
{
  #moves an internal node
  tr = tra
  root=tra$root[[2]]
  nf=length(tra$tip.label)
  Q=(nf+1):max(tr$edge)
  if (length(tra$root[[2]])>1){
    j=Q[sample(length(Q),1)]
  } else {
    Q=Q[-which(Q==root)]
    j=Q[sample(1:length(Q),1)]
  }
  A = which(tr$edge[, 1] == j)#branches partant des changés
  if (j %in% root){
    B=numeric()
    v=agemax-hauteur(tra,j)
  } else{
    B = which(tr$edge[, 2] == j)#branches allant aux changés
    v = min(tr$edge.length[B])
  }
  u = min(tr$edge.length[A])
  w = runif(1, 0, u + v)
  tr$edge.length[A] = tr$edge.length[A] - u + w
  tr$edge.length[B] = tr$edge.length[B] - v + u + v - w
  if (length(tr$root[[2]]>1) || (length(tr$root)==1 & hauteur(tr,tr$root[[2]])==agemax)){
    return(tr)
  } else {
    return(tra)
  }
}

modifst <-function(tra){
  #rescales a subtree starting from j
  tr = tra
  root=tra$root[[2]]
  nf=length(tra$tip.label)
  Q=(nf+1):max(tr$edge)
  if (length(tra$root[[2]])>1){
    j=Q[sample(length(Q),1)]
  } else {
    Q=Q[-which(Q==root)]
    j=Q[sample(1:length(Q),1)]
  }
  edge_enfants=which(tr$edge[,1]==j)
  edge_parent=which(tr$edge[,2]==j)
  hmin=min(hauteur(tr,tr$edge[edge_enfants[1],2]),hauteur(tr,tr$edge[edge_enfants[2],2]))
  hmax=hauteur(tr,tr$edge[edge_parent,1])
  new_h=runif(1,hmin,hmax)
  st=sousarbre(tr,j)
  tr$edge.length[st]=tr$edge.length[st]*new_h/hauteur(tr,j)
  tr$edge.length[edge_parent]=hmax-new_h
  return(list(tr,(length(st)-1)*new_h/hauteur(tr,j)))
}


modifrootage <- function(tra,agemax)
{
  #modifies the root age
  tr = tra
  root=tra$root[[2]]
  m=which(tr$edge[,1]==root)
  rootage=hauteur(tr,tr$root[[2]])
  enfants=tr$edge[tr$edge[,1]==root,2]
  ageenfantsroot=c(hauteur(tr,enfants[1]),hauteur(tr,enfants[2]))
  newrootage=rgamma(1,agemax[1],agemax[2])
  while (newrootage<ageenfantsroot[1] || newrootage<ageenfantsroot[2]){
    newrootage=rgamma(1,agemax[1],agemax[2])
  }
  corr=dgamma(rootage,agemax[1],agemax[2],log=T)-dgamma(newrootage,agemax[1],agemax[2],log=T)
  tr$edge.length[m]=newrootage-ageenfantsroot
  if (!all(tra$edge.length>0)){
    return(list(tra,0))
  } else {
    return(list(tr,corr))
  }
}

modifresc <- function(tra,agemax){
  #rescales a tree according to rootage
  tr=tra
  root=tra$root[[2]]
  m=which(tr$edge[,1]==root)
  rootage=hauteur(tr,tr$root[[2]])
  newrootage=rgamma(1,agemax[1],agemax[2])
  tr$edge.length=tr$edge.length*newrootage/rootage
  corr=length(tr$edge.length-1)*(log(newrootage)-log(rootage) ) + dgamma(rootage,agemax[1],agemax[2],log=T) - dgamma(newrootage,agemax[1],agemax[2],log=T)
  if (all(tr$edge.length>0)){
    return(list(tr,corr))
  } else {
    return(list(tra,0))
  }
}
