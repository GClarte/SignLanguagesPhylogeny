
gibbstopopart2prior <-  function(tr,
                            roott,
                            priortree,
                            agemax,nnodes,feuilles) {
  #change la topologie (et les ages) de l'arbre, pour une partie des cognats
  if (runif(1)>1/2){VV = modificationtopo2prior(tr,
                         nnodes,
                         agemax,
                         feuilles)
  } else {
    VV = modificationtopofrereprior(tr,
                                nnodes,
                                agemax,
                                feuilles)
    
  }
  trt = VV[[1]]
  gam = VV[[2]]
  root = trt$root
  alph =  priortree(trt, agemax) - priortree(tr, agemax)
  u = log(runif(1, 0, 1))
  if (u < min(0, alph + gam)) {
    return(list(
      trt))
  } else {
    return(list(tr))
  }
}

modificationtopo2prior <-
  function(tr,
           nnodes,
           agemax,
           feuilles) {
    #on prend un noeud e non racine et on le déplace pour qu'il ait a comme frère
    root = tr$root
    trt = tr
    corrb = 0
    e = sample((1:nnodes)[-root[[2]]], 1)#noeud qui bouge
    i = which(tr$edge[, 2] == e)
    d = tr$edge[i, 1]#parent de celui ci
    haut=hauteur(tr,d)
    hauts=hauteurtous(tr,feuilles,agemax)
    st = sousarbre(tr, e)#sous arbre qui en part
    ifr = which(tr$edge[, 1] == d & !tr$edge[, 2] == e)#branche frère
    f = tr$edge[ifr, 2]#neoud frère
    if (!d %in% root[[2]]) {
      k = which(tr$edge[, 2] == d) #si ce n'est pas une racine alors on recolle directement
      dd = tr$edge[k, 1]
      trt$edge[k, ] = c(dd, f)
      trt$edge.length[k] = trt$edge.length[k] + trt$edge.length[ifr]
      } else {
      #si on rend le frère de $e$ racine, on compte que dans le sens réciproque il faudrait tirer du prior pour avoir les valeurs attendues
      trt$root[[2]] = c(trt$root[[2]][-which(trt$root[[2]] == d)], f)
      w = agemax - hauteur(tr, f)
      dd = numeric()
      }
    if (d %in% root[[1]]) {
      trt$root[[1]] = c(trt$root[[1]][-which(trt$root[[1]] == d)], f)
    }
    quiposs=which(hauts[[1]]<haut & hauts[[2]]>haut & (1:max(tr$edge)) != d)
    if (length(quiposs) > 0) {
      a=quiposs[sample(1:length(quiposs),1)]
    } else {
      return(list(tr, 0))
    }
    u=-(+haut - hauts[[1]][a])+(hauts[[2]][a]-hauts[[1]][a])
    if (!a %in% root[[2]]) {
      j = which(tr$edge[, 2] == a)
      ww=tr$edge.length[j]
      y = tr$edge[j, 1]
      trt$edge[j, ] = c(y, d)
    } else {
      ww=agemax-hauteur(tr,a)
      trt$root[[2]] = c(trt$root[[2]][-which(trt$root[[2]] == a)], d)
    }
    if (a %in% root[[1]]) {
      trt$root[[1]] = c(trt$root[[1]][-which(trt$root[[1]] == a)], d)
    }
    trt$edge[ifr, ] = c(d, a)
    trt$edge.length[ifr] = ww-u
    if (!a %in% root[[2]]) {
      trt$edge.length[j] = u
    }
    if (testarbre3(trt, trt$root[[2]], nnodes, feuilles)) {
      return(list(trt, corrb))
    } else {
      return(list(tr, 0))
    }
  }


modificationtopofrereprior <-
  function(tr,
           nnodes,
           agemax,
           feuilles) {
    #on prend un noeud e non racine et on le déplace pour qu'il ait a comme frère
    trt=tr
    root = tr$root
    corrb = 0
    etot = (1:nnodes)[-root[[2]]]#noeud qui bouge
    tousparents=sapply(etot,function(z){tr$edge[which(tr$edge[,2]==z),1]})
    tousgrandsparents=sapply(tousparents,function(z){tr$edge[which(tr$edge[,2]==z),1]})
    quiposs=sapply(tousgrandsparents,length)>0
    if (sum(quiposs)>0){
      ij=sample(1:length(quiposs),1,prob=as.numeric(quiposs))
      e=etot[ij]
      a=tr$edge[which(tr$edge[,1]==tousgrandsparents[[ij]] & tr$edge[,2]!=tousparents[ij]),2]
      i = which(tr$edge[, 2] == e)
      d = tr$edge[i, 1]#parent de celui ci
      haut=hauteur(tr,d)
      st = sousarbre(tr, e)#sous arbre qui en part
      ifr = which(tr$edge[, 1] == d & !tr$edge[, 2] == e)#branche frère
      f = tr$edge[ifr, 2]#neoud frère
      k = which(tr$edge[, 2] == d) #si ce n'est pas une racine alors on recolle directement
      dd = tr$edge[k, 1]
      trt$edge[k, ] = c(dd, f)
      trt$edge.length[k] = tr$edge.length[k] + tr$edge.length[ifr]
      w = trt$edge.length[k]
      j = which(tr$edge[, 2] == a)
      u=runif(1,0,tr$edge.length[j])#endroit où on recolle
      ww=tr$edge.length[j]
      y = tr$edge[j, 1]
      trt$edge[j, ] = c(y, d)
      trt$edge[ifr, ] = c(d, a)
      trt$edge.length[j]=u
      trt$edge.length[ifr] = ww-u
      trt$edge.length[c(st,i)]=tr$edge.length[c(st,i)]*(hauteur(tr,a)+ww-u)/(hauteur(tr,e)+tr$edge.length[i])
      
      corrb = corrb  +log(1/(tr$edge.length[k] + tr$edge.length[ifr]))- log(1/tr$edge.length[j])
      if (d %in% root[[1]]) {
        trt$root[[1]] = c(trt$root[[1]][-which(trt$root[[1]] == d)], f)
      }
      if (a %in% root[[1]]) {
        trt$root[[1]] = c(trt$root[[1]][-which(trt$root[[1]] == a)], d)
      }
      if (testarbre3(trt, trt$root[[2]], nnodes, feuilles)) {
        return(list(trt, corrb))
      } else {
        return(list(tr, 0))
      }} else {
        return(list(tr, 0))
      }
  }
