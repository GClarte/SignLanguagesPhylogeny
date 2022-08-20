

#for the jacobian when rescaling
jacrescale=-1

gibbstopopart2 <-  function(tr,
                            Da,
                            la,
                            X,
                            M,
                            Trposs,
                            nph,
                            roott,
                            loiini,
                            bruit,
                            lkldinter,
                            L,
                            NL,
                            rho,
                            quelcogn,
                            priortree,
                            loiappliste,
                            Tps,
                            agemax,
                            P,
                            tipprior,
                            Prior,
                            fin) {
  #updated the topology according to 5 possible moves
  #qqqq=3
  qqqq=sample(1:2,1)
  if( qqqq==1){
    if (runif(1,0,1)>1/2 || length(tr$root[[2]])>1){
      VV = modificationtopo2(tr,
                             X ,
                             NL,
                             L ,
                             Tps,
                             P,
                             rho,
                             la,
                             length(lkldinter[[1]]),
                             agemax,
                             1:ncol(Da[[1]]),
                             Da,
                             tipprior,
                             M,
                             loiappliste,
                             Trposs,
                             bruit,
                             nph)
      commentbis=1
    } else {
      VV = modificationtopofrere(tr,
                                 X ,
                                 NL,
                                 L ,
                                 Tps,
                                 P,
                                 rho,
                                 la,
                                 length(lkldinter[[1]]),
                                 agemax,
                                 1:ncol(Da[[1]]),
                                 Da,
                                 tipprior,
                                 M,
                                 loiappliste,
                                 Trposs,
                                 bruit,
                                 nph)
      commentbis=2
    }
  } else if (qqqq==2) {
    if (runif(1,0,1)>1/2 || length(tr$root[[2]])>1){
      VV = modificationtopo3(tr,
                             X ,
                             NL,
                             L ,
                             Tps,
                             P,
                             rho,
                             la,
                             length(lkldinter[[1]]),
                             agemax,
                             1:ncol(Da[[1]]),
                             Da,
                             tipprior,
                             M,
                             loiappliste,
                             Trposs,
                             bruit,
                             nph)
      commentbis=3
    } else {
      VV = modificationtopofrere2(tr,
                                  X ,
                                  NL,
                                  L ,
                                  Tps,
                                  P,
                                  rho,
                                  la,
                                  length(lkldinter[[1]]),
                                  agemax,
                                  1:ncol(Da[[1]]),
                                  Da,
                                  tipprior,
                                  M,
                                  loiappliste,
                                  Trposs,
                                  bruit,
                                  nph)
      commentbis=4
    }
  } else {
    if (runif(1)-1>1/2){
      VV = modificationtopooutgroup(tr,
                                    X ,
                                    NL,
                                    L ,
                                    Tps,
                                    P,
                                    rho,
                                    la,
                                    length(lkldinter[[1]]),
                                    agemax,
                                    1:ncol(Da[[1]]),
                                    Da,
                                    tipprior,
                                    M,
                                    loiappliste,
                                    Trposs,
                                    bruit,
                                    nph,
                                    Prior,
                                    fin)
      commentbis=5
    } else {
      VV = modificationtopooutgroupsansresc(tr,
                                            X ,
                                            NL,
                                            L ,
                                            Tps,
                                            P,
                                            rho,
                                            la,
                                            length(lkldinter[[1]]),
                                            agemax,
                                            1:ncol(Da[[1]]),
                                            Da,
                                            tipprior,
                                            M,
                                            loiappliste,
                                            Trposs,
                                            bruit,
                                            nph,
                                            Prior,
                                            fin)
      commentbis=6
    }
  }
  trt = VV[[1]]
  # iii=1
  # while (!all(sapply(tipprior,function(tip){is.monophyletic.perso(trt,tip)})) & iii<1000){
  #   VV = modificationtopo2(tr,
  #                          X ,
  #                          NL,
  #                          L ,
  #                          Tps,
  #                          P,
  #                          rho,
  #                          la,
  #                          length(lkldinter[[1]]),
  #                          agemax,
  #                          1:ncol(Da[[1]]),
  #                          Da)
  #   trt = VV[[1]]
  #   iii=iii+1
  # }
  Xt = VV[[2]]
  NLt = VV[[3]]
  Lt = VV[[4]]
  Tpst = VV[[5]]
  gam = VV[[6]]
  Mt=VV[[7]]
  loiapplistec=VV[[8]]
  alphabis =  priortree(trt, agemax) - priortree(tr, agemax)
  if (alphabis==-Inf || VV[[9]]==-1){
    return(list(tr, lkldinter, c(qqqq,0, 0, gam,-1),-6*(gam!=0), loiappliste, X, NL, L, Tps, M,rep(-2,nch),la,rho,bruit))
  } else if (VV[[9]]==0){
    lkldintert = lapply(1:nch, function(x) {
      pruninglkldlistini(trt,Xt[[x]],Mt[[x]],Lt,Da[[x]],nph[[x]],length(lkldinter[[1]]),
                         trt$root[[2]],
                         Trposs[[x]],
                         loiini[[x]],
                         VV[[14]],
                         loiapplistec[[x]])
    })
    alph = lkldtopotout(trt, Xt, la, NLt, rho, Lt) - lkldtopotout(tr, X, la, NL, rho, L) + alphabis
    bet = lkldcogn(lkldintert, nch, trt$root, quelcogn, loiini) - lkldcogn(lkldinter, nch, tr$root, quelcogn, loiini)
    #bet=0
    u = log(runif(1, 0, 1))
    alphaaccept=min(0, alph + bet + gam)
    if (is.nan(alphaaccept)){
      alphaaccept=-Inf
    }
    if (is.na(alphaaccept)){
      alphaaccept=-Inf
    }
    if (u < alphaaccept) {
      return(list(
        trt,
        lkldintert,
        c(qqqq,alph, bet, gam,1),
        commentbis,
        loiapplistec,
        Xt,
        NLt,
        Lt,
        Tpst,
        Mt,
        VV[[11]],
        VV[[12]],
        VV[[13]],
        VV[[14]]
      ))
    } else {
      return(list(tr, lkldinter, c(alph, bet, gam),-commentbis, loiappliste, X, NL, L, Tps,M,rep(-.2,nch),la,rho,bruit))
    }
  } else {
    path1=cheminrac(trt,VV[[9]],trt$root[[2]])
    path2=cheminrac(trt,VV[[10]],trt$root[[2]])
    root = trt$root
    nch = length(Xt)
    
    # Mt = lapply(1:nch, function(z) {
    #   lapply(1:length(tr$edge.length), function(y) {
    #     transitionmatrix(Xt[[z]][[y]], Trposs[[z]], nph[z], bruit, trt$edge.length[y])
    #   })
    # })
    # loiapplistec = lapply(1:nch, function(x) {
    #   lapply(1:length(trt$edge.length), function(y) {
    #     loiapp(Xt[[x]][[y]],
    #            Lt[Lt[, y] > 0, y],
    #            loiini[[x]],
    #            Trposs[[x]],
    #            bruit,
    #            trt$edge.length[y],
    #            Tpst[[x]][[y]])
    #   })
    # })
    lkldintert = lapply(1:nch, function(x) {
      pruningeco(trt,Xt[[x]],Mt[[x]],Da[[x]],nph[[x]],path1,
                 lkldinter[[x]],
                 Trposs[[x]],
                 Lt,
                 bruit,
                 loiini[[x]],
                 loiapplistec[[x]])
    })
    alph = lkldtopotout(trt, Xt, la, NLt, rho, Lt) - lkldtopotout(tr, X, la, NL, rho, L) + alphabis
    bet = lkldcogn(lkldintert, nch, trt$root, quelcogn, loiini) - lkldcogn(lkldinter, nch, tr$root, quelcogn, loiini)
    #bet=0
    u = log(runif(1, 0, 1))
    alphaaccept=min(0, alph + bet + gam)
    if (is.nan(alphaaccept)){
      alphaaccept=-Inf
    }
    if (is.na(alphaaccept)){
      alphaaccept=-Inf
    }
    if (u < alphaaccept) {
      return(list(
        trt,
        lkldintert,
        c(qqqq,alph, bet, gam,1),
        commentbis,
        loiapplistec,
        Xt,
        NLt,
        Lt,
        Tpst,
        Mt,
        VV[[11]],
        VV[[12]],
        VV[[13]],
        VV[[14]]
      ))
    } else {
      return(list(tr, lkldinter, c(qqqq,alph, bet, gam,-1),-commentbis, loiappliste, X, NL, L, Tps,M,rep(-2,nch),la,rho,bruit))
    }
    
  }
}

modificationtopo2 <-
  function(tr,
           X,
           NL,
           L,
           Tps,
           P,
           Rho,
           La,
           nnodes,
           agemax,
           feuilles,
           Dat,
           tipprior,
           M,
           loiappliste,
           Trposs,
           bruit,
           nph) {
    #on prend un noeud e non racine et on le déplace pour qu'il ait a comme frère
    root = tr$root
    Xt = X
    NLt = NL
    Lt = L
    Tpst = Tps
    trt = tr
    nch = length(X)
    corrb = 0
    Mt=M
    loiapplistet=loiappliste
    etot = sample((1:nnodes)[-root[[2]]])#noeud qui bouge
    iii=0
    ouiounon=T
    while(ouiounon & iii<length(etot)){
      iii=iii+1
      e=etot[iii]
      i = which(tr$edge[, 2] == e)
      d = tr$edge[i, 1]#parent de celui ci
      haut=hauteur(tr,d)
      hauts=hauteurtous(tr,feuilles,agemax)
      st = sousarbre(tr, e)#sous arbre qui en part
      ifr = which(tr$edge[, 1] == d & !tr$edge[, 2] == e)#branche frère
      f = tr$edge[ifr, 2]#neoud frère
      
      quiposs=which(hauts[[1]]<haut & hauts[[2]]>haut & (1:max(tr$edge)) != d)
      if (length(quiposs) > 0) {
        a=quiposs[sample(1:length(quiposs))]
        test=sapply(a,function(z){WW=c(feuillessousarbre(tr,z),feuillessousarbre(tr,e));return(all(sapply(tipprior,function(x){length(intersect(WW,x)) %in% c(0, length(WW), length(x))})) & z!=f)})
        if (any(test)){
          a=a[which(test)[1]]
          ouiounon=F
        } else {
          ouiounon=T
        }
      } else {
        ouiounon=T
      }
    }
    if(iii==length(etot) & ouiounon){
      return(list(tr, X, NL, L, Tps, 0,M ,loiappliste,-1,-1,rep(-1,nch),La,Rho,bruit))
    }
    
    if (!d %in% root[[2]]) {
      k = which(tr$edge[, 2] == d) #si ce n'est pas une racine alors on recolle directement
      dd = tr$edge[k, 1]
      trt$edge[k, ] = c(dd, f)
      trt$edge.length[k] = trt$edge.length[k] + trt$edge.length[ifr]
      for (ii in 1:nch) {
        Xt[[ii]][[k]] = c(X[[ii]][[k]], X[[ii]][[ifr]])
        Tpst[[ii]][[k]] = c(Tps[[ii]][[k]] * tr$edge.length[k], tr$edge.length[k]+ Tps[[ii]][[ifr]] *
                              tr$edge.length[ifr]) / (tr$edge.length[k] + tr$edge.length[ifr])
      }
      Lt[, k] = (
        (L[, ifr] * tr$edge.length[ifr] + tr$edge.length[k]) * (NL[, ifr] > 0) + (NL[, ifr] == 0 &
                                                                                    NL[, k] > 0) * L[, k] * tr$edge.length[k]
      ) / (tr$edge.length[k] + tr$edge.length[ifr])
      NLt[, k] = NL[, k] + NL[, ifr]
      w = trt$edge.length[k]
      #on rajoute le laplacien réciproque qui est lié à l'action de découpe de la branche en 2 cf ce qu'il y a après
      rrr=which(NLt[,k]>1 & Lt[,k]>(tr$edge.length[k])/trt$edge.length[k])#ceux pour lesquels il y a au moins une transformation qui pourrait être mise sur la branche en question
      for (ii in 1:nch){
        Mt[[ii]][[k]]=transitionmatrix(Xt[[ii]][[k]], Trposs[[ii]], nph[ii], bruit, trt$edge.length[k])
        loiapplistet[[ii]][[k]]=loiapp(Xt[[ii]][[k]],
                                       Lt[Lt[, k] > 0, k],
                                       loiini[[ii]],
                                       Trposs[[ii]],
                                       bruit,
                                       trt$edge.length[k],
                                       Tpst[[ii]][[k]])
      }
      corrb = corrb + sum(dbinom(NL[rrr,k],
                                 size=NLt[rrr,k]-1,
                                 prob=tr$edge.length[k]/(tr$edge.length[k] + L[rrr,ifr]*tr$edge.length[ifr])))+
        sum(dbeta(L[L[, k]>0 & L[,ifr] > 0, k], NL[NL[, k] > 0 & NL[,ifr]>0, k], 1, log = T)) 
    } else {
      #si on rend le frère de $e$ racine, on compte que dans le sens réciproque il faudrait tirer du prior pour avoir les valeurs attendues
      trt$root[[2]] = c(trt$root[[2]][-which(trt$root[[2]] == d)], f)
      w = agemax - hauteur(tr, f)
      dd = numeric()
      for (ii in 1:nch) {
        corrb = corrb + log(prod(P[[ii]][X[[ii]][[ifr]]])) + dpois(length(X[[ii]][[ifr]]), La[ii] *
                                                                     tr$edge.length[ifr], log = T)
      }
      corrb = corrb + sum(dbeta(L[L[, ifr] != 0, ifr], NL[L[, ifr] != 0, ifr], 1, log =
                                  T)) + sum(dpois(NL[, ifr], tr$edge.length[ifr] * Rho, log = T))
    }
    if (d %in% root[[1]]) {
      trt$root[[1]] = c(trt$root[[1]][-which(trt$root[[1]] == d)], f)
    }
    
    u=hauts[[2]][a]-haut
    if (!a %in% root[[2]]) {
      j = which(tr$edge[, 2] == a)
      ww=tr$edge.length[j]
      y = tr$edge[j, 1]
      trt$edge[j, ] = c(y, d)
      trt$edge[ifr, ] = c(d, a)
      trt$edge.length[ifr] = ww-u
      for (ii in 1:nch) {
        cbtr = Tps[[ii]][[j]] < (u / ww)
        Xt[[ii]][[j]] = X[[ii]][[j]][cbtr]
        Xt[[ii]][[ifr]] = X[[ii]][[j]][!cbtr]#puisque e a toujours un frère.
        Tpst[[ii]][[j]] = ww * (Tps[[ii]][[j]][cbtr]) / u
        Tpst[[ii]][[ifr]] = (Tps[[ii]][[j]][!cbtr]*ww - u)/(ww-u)
      }
      NLt[,j]=0
      Lt[,j]=0
      NLt[,ifr]=0
      Lt[,ifr]=0
      #quand il n'y a qu'une app ok
      sss1=which(NL[,j]==1 & L[,j]<(u/ww))
      sss2=which(NL[,j]==1 & L[,j]>(u/ww))
      NLt[sss1,j]=1
      NLt[sss2,ifr]=1
      Lt[sss1,j]=L[sss1,j]*ww/u
      Lt[sss2,ifr]=(L[sss2,j]*ww-u)/(ww-u)
      #maintenant les trucs problématiques, lorsqu'il y a plusieurs app
      sss1=which(NL[,j]>1 & L[,j]<u/ww)
      sss2=which(NL[,j]>1 & L[,j]>u/ww)
      NLt[sss1,j]=NL[sss1,j]
      NLt[sss2,j]=rbinom(length(sss2),size=NL[sss2,j]-1,prob=u/(L[sss2,j]*ww))
      NLt[sss2,ifr]=NL[sss2,j]-NLt[sss2,j]
      Lt[sss1,j]=L[sss1,j]*ww/u
      Lt[sss2,ifr]=(L[sss2,j]*ww-u)/(ww-u)
      sss3=which(NLt[sss2,j]>0)
      Lt[sss2[sss3],j]=rbeta(length(sss3),NLt[sss2[sss3],j],1)
      trt$edge.length[j]=u
      for (ii in 1:nch){
        Mt[[ii]][[j]]=transitionmatrix(Xt[[ii]][[j]], Trposs[[ii]], nph[ii], bruit, trt$edge.length[j])
        loiapplistet[[ii]][[j]]=loiapp(Xt[[ii]][[j]],
                                       Lt[Lt[, j] > 0, j],
                                       loiini[[ii]],
                                       Trposs[[ii]],
                                       bruit,
                                       trt$edge.length[j],
                                       Tpst[[ii]][[j]])
        Mt[[ii]][[ifr]]=transitionmatrix(Xt[[ii]][[ifr]], Trposs[[ii]], nph[ii], bruit, trt$edge.length[ifr])
        loiapplistet[[ii]][[ifr]]=loiapp(Xt[[ii]][[ifr]],
                                         Lt[Lt[, ifr] > 0, ifr],
                                         loiini[[ii]],
                                         Trposs[[ii]],
                                         bruit,
                                         trt$edge.length[ifr],
                                         Tpst[[ii]][[ifr]])
      }
      
      corrb=corrb-sum(dbinom(
        NLt[sss2,j],size=NL[sss2,j]-1,prob=u/(L[sss2,j]*ww),
        log = T
      ))
      corrb = corrb - sum(dbeta(Lt[sss2[sss3],j], NLt[sss2[sss3],j], 1, log = T))
    } else {
      ww=agemax-hauteur(tr,a)
      trt$root[[2]] = c(trt$root[[2]][-which(trt$root[[2]] == a)], d)
      trt$edge[ifr, ] = c(d, a)
      trt$edge.length[ifr] = ww-u
      for (ii in 1:nch) {
        Xt[[ii]][[ifr]] = sample(1:length(P[[ii]]),
                                 rpois(1, La[ii] * u),
                                 prob = P[[ii]],
                                 replace = T)
        Tpst[[ii]][[ifr]] = sort(runif(length(Xt[[ii]][[ifr]])))
        corrb = corrb - log(prod(P[[ii]][Xt[[ii]][[ifr]]])) - dpois(length(Xt[[ii]][[ifr]]), La[ii] *
                                                                      u, log = T)
      }
      NLt[, ifr] = rpois(nrow(Dat[[1]]), u * Rho)
      Lt[, ifr] = rbeta(nrow(Dat[[1]]), NLt[, ifr], 1)
      for (ii in 1:nch){
        Mt[[ii]][[ifr]]=transitionmatrix(Xt[[ii]][[ifr]], Trposs[[ii]], nph[ii], bruit, trt$edge.length[ifr])
        loiapplistet[[ii]][[ifr]]=loiapp(Xt[[ii]][[ifr]],
                                         Lt[Lt[, ifr] > 0, ifr],
                                         loiini[[ii]],
                                         Trposs[[ii]],
                                         bruit,
                                         trt$edge.length[ifr],
                                         Tpst[[ii]][[ifr]])
      }
      corrb = corrb - sum(dpois(NLt[, ifr], u * Rho, log = T)) - sum(dbeta(Lt[Lt[, ifr] !=
                                                                                0, ifr], NLt[Lt[, ifr] != 0, ifr], 1, log = T))
    }
    if (a %in% root[[1]]) {
      trt$root[[1]] = c(trt$root[[1]][-which(trt$root[[1]] == a)], d)
    }
    for (ii in 1:nch) {
      Xt[[ii]][[i]] = sample(1:length(P[[ii]]),
                             rpois(1, La[ii] * tr$edge.length[i]),
                             prob = P[[ii]],
                             replace = T)
      Tpst[[ii]][[i]] = sort(runif(length(Xt[[ii]][[i]])))
      corrb = corrb - sum(log(P[[ii]][Xt[[ii]][[i]]])) - dpois(length(Xt[[ii]][[i]]), La[ii] *
                                                                 trt$edge.length[i], log = T)
      corrb = corrb + sum(log(P[[ii]][X[[ii]][[i]]])) + dpois(length(X[[ii]][[i]]), La[ii] *
                                                                tr$edge.length[i], log = T)
    }
    NLt[, i] = rpois(nrow(Dat[[1]]), trt$edge.length[i] * Rho)
    Lt[, i] = rbeta(nrow(Dat[[1]]), NLt[, i], 1)
    
    for (ii in 1:nch){
      Mt[[ii]][[i]]=transitionmatrix(Xt[[ii]][[i]], Trposs[[ii]], nph[ii], bruit, trt$edge.length[i])
      loiapplistet[[ii]][[i]]=loiapp(Xt[[ii]][[i]],
                                     Lt[Lt[, i] > 0, i],
                                     loiini[[ii]],
                                     Trposs[[ii]],
                                     bruit,
                                     trt$edge.length[i],
                                     Tpst[[ii]][[i]])
    }
    
    corrb = corrb - sum(dpois(NLt[, i], trt$edge.length[i] * Rho, log =
                                T)) - sum(dbeta(Lt[Lt[, i] != 0, i], NLt[Lt[, i] != 0, i], 1, log = T))
    corrb = corrb + sum(dpois(NL[, i], tr$edge.length[i] * Rho, log =
                                T)) + sum(dbeta(L[L[, i] != 0, i], NL[L[, i] != 0, i], 1, log = T))
    if (a %in% root[[2]]){
      quoisur=rep(0,nch)
    } else {
      quoisur=sapply(1:nch,function(x){length(X[[x]][[j]])})+.5*(a<=ncol(Dat[[1]])&y %in%root[[2]])
    }
    if (testarbre3(trt, trt$root[[2]], nnodes, feuilles)) {
      return(list(trt, Xt, NLt, Lt, Tpst, corrb, Mt,loiapplistet,d,f,quoisur,La,Rho,bruit))
    } else {
      return(list(tr, X, NL, L, Tps, 0,M,loiappliste,-1,-1,rep(-1,nch),La,Rho,bruit))
    }
  }


modificationtopo3 <-
  function(tr,
           X,
           NL,
           L,
           Tps,
           P,
           Rho,
           La,
           nnodes,
           agemax,
           feuilles,
           Dat,
           tipprior,
           M,
           loiappliste,
           Trposs,
           bruit,
           nph) {
    #on prend un noeud e non racine et on le déplace pour qu'il ait a comme frère
    root = tr$root
    Xt = X
    NLt = NL
    Lt = L
    Tpst = Tps
    trt = tr
    nch = length(X)
    corrb = 0
    Mt=M
    loiapplistet=loiappliste
    etot = sample((1:nnodes)[-root[[2]]])#noeud qui bouge
    iii=0
    ouiounon=T
    while(ouiounon & iii<length(etot)){
      iii=iii+1
      e=etot[iii]
      i = which(tr$edge[, 2] == e)
      d = tr$edge[i, 1]#parent de celui ci
      haut=hauteur(tr,d)
      hauts=hauteurtous(tr,feuilles,agemax)
      st = sousarbre(tr, e)#sous arbre qui en part
      ifr = which(tr$edge[, 1] == d & !tr$edge[, 2] == e)#branche frère
      f = tr$edge[ifr, 2]#neoud frère
      
      quiposs=which(hauts[[1]]<haut & hauts[[2]]>haut & (1:max(tr$edge)) != d)
      if (length(quiposs) > 0) {
        a=quiposs[sample(1:length(quiposs))]
        test=sapply(a,function(z){WW=c(feuillessousarbre(tr,z),feuillessousarbre(tr,e));return(all(sapply(tipprior,function(x){length(intersect(WW,x)) %in% c(0, length(WW), length(x))})) & z!=f)})
        if (any(test)){
          a=a[which(test)[1]]
          ouiounon=F
        } else {
          ouiounon=T
        }
      } else {
        ouiounon=T
      }
    }
    if(iii==length(etot) & ouiounon){
      return(list(tr, X, NL, L, Tps, 0,M ,loiappliste,-1,-1,rep(-1,nch),La,Rho,bruit))
    }
    
    if (!d %in% root[[2]]) {
      k = which(tr$edge[, 2] == d) #si ce n'est pas une racine alors on recolle directement
      dd = tr$edge[k, 1]
      trt$edge[k, ] = c(dd, f)
      trt$edge.length[k] = trt$edge.length[k] + trt$edge.length[ifr]
      for (ii in 1:nch) {
        Xt[[ii]][[k]] = c(X[[ii]][[k]], X[[ii]][[ifr]])
        Tpst[[ii]][[k]] = c(Tps[[ii]][[k]] * tr$edge.length[k], tr$edge.length[k]+ Tps[[ii]][[ifr]] *
                              tr$edge.length[ifr]) / (tr$edge.length[k] + tr$edge.length[ifr])
      }
      Lt[, k] = (
        (L[, ifr] * tr$edge.length[ifr]+ tr$edge.length[k]) * (NL[, ifr] > 0) + (NL[, ifr] == 0 &
                                                                                   NL[, k] > 0) * L[, k] * tr$edge.length[k]
      ) / (tr$edge.length[k] + tr$edge.length[ifr])
      NLt[, k] = NL[, k] + NL[, ifr]
      w = trt$edge.length[k]
      #on rajoute le laplacien réciproque qui est lié à l'action de découpe de la branche en 2 cf ce qu'il y a après
      rrr=which(NLt[,k]>1 & Lt[,k]>(tr$edge.length[k])/trt$edge.length[k])#ceux pour lesquels il y a au moins une transformation qui pourrait être mise sur la branche en question
      for (ii in 1:nch){
        Mt[[ii]][[k]]=transitionmatrix(Xt[[ii]][[k]], Trposs[[ii]], nph[ii], bruit, trt$edge.length[k])
        loiapplistet[[ii]][[k]]=loiapp(Xt[[ii]][[k]],
                                       Lt[Lt[, k] > 0, k],
                                       loiini[[ii]],
                                       Trposs[[ii]],
                                       bruit,
                                       trt$edge.length[k],
                                       Tpst[[ii]][[k]])
      }
      corrb = corrb + sum(dbinom(NL[rrr,k],
                                 size=NLt[rrr,k]-1,
                                 prob=tr$edge.length[k]/(tr$edge.length[k] + L[rrr,ifr]*tr$edge.length[ifr])))+
        sum(dbeta(L[L[, k] > 0 & L[,ifr]>0, k], NL[NL[, k] > 0 & NL[,ifr]>0, k], 1, log = T))#en suis-je sûr ?
    } else {
      #si on rend le frère de $e$ racine, on compte que dans le sens réciproque il faudrait tirer du prior pour avoir les valeurs attendues
      trt$root[[2]] = c(trt$root[[2]][-which(trt$root[[2]] == d)], f)
      w = agemax - hauteur(tr, f)
      dd = numeric()
      for (ii in 1:nch) {
        corrb = corrb + log(prod(P[[ii]][X[[ii]][[ifr]]])) + dpois(length(X[[ii]][[ifr]]), La[ii] *
                                                                     tr$edge.length[ifr], log = T)
      }
      corrb = corrb + sum(dbeta(L[L[, ifr] != 0, ifr], NL[L[, ifr] != 0, ifr], 1, log =
                                  T)) + sum(dpois(NL[, ifr], tr$edge.length[ifr] * Rho, log = T))
    }
    if (d %in% root[[1]]) {
      trt$root[[1]] = c(trt$root[[1]][-which(trt$root[[1]] == d)], f)
    }
    
    u=-(+haut - hauts[[1]][a])+(hauts[[2]][a]-hauts[[1]][a])
    if (!a %in% root[[2]]) {
      j = which(tr$edge[, 2] == a)
      ww=tr$edge.length[j]
      y = tr$edge[j, 1]
      trt$edge[j, ] = c(y, d)
      trt$edge[ifr, ] = c(d, a)
      trt$edge.length[ifr] = ww-u
      trt$edge.length[j] = u
      for (ii in 1:nch) {
        cbtr = Tps[[ii]][[j]] < u / ww
        Xt[[ii]][[j]] = X[[ii]][[j]][cbtr]
        Xt[[ii]][[ifr]] = X[[ii]][[j]][!cbtr]#puisque e a toujours un frère.
        Tpst[[ii]][[j]] = ww * (Tps[[ii]][[j]][cbtr]) / u
        Tpst[[ii]][[ifr]] = (Tps[[ii]][[j]][!cbtr]*ww - u)/(ww-u)
      }
      NLt[,j]=0
      Lt[,j]=0
      NLt[,ifr]=0
      Lt[,ifr]=0
      #quand il n'y a qu'une app ok
      sss1=which(NL[,j]==1 & L[,j]<u/ww)
      sss2=which(NL[,j]==1 & L[,j]>u/ww)
      NLt[sss1,j]=1
      NLt[sss2,ifr]=1
      Lt[sss1,j]=L[sss1,j]*ww/u
      Lt[sss2,ifr]=(L[sss2,j]*ww-u)/(ww-u)
      #maintenant les trucs problématiques, lorsqu'il y a plusieurs app
      sss1=which(NL[,j]>1 & L[,j]<u/ww)
      sss2=which(NL[,j]>1 & L[,j]>u/ww)
      NLt[sss1,j]=NL[sss1,j]
      NLt[sss2,j]=rbinom(length(sss2),size=NL[sss2,j]-1,prob=u/(L[sss2,j]*ww))
      NLt[sss2,ifr]=NL[sss2,j]-NLt[sss2,j]
      Lt[sss1,j]=L[sss1,j]*ww/u
      Lt[sss2,ifr]=(L[sss2,j]*ww-u)/(ww-u)
      sss3=which(NLt[sss2,j]>0)
      Lt[sss2[sss3],j]=rbeta(length(sss3),NLt[sss2[sss3],j],1)
      
      for (ii in 1:nch){
        Mt[[ii]][[j]]=transitionmatrix(Xt[[ii]][[j]], Trposs[[ii]], nph[ii], bruit, trt$edge.length[j])
        loiapplistet[[ii]][[j]]=loiapp(Xt[[ii]][[j]],
                                       Lt[Lt[, j] > 0, j],
                                       loiini[[ii]],
                                       Trposs[[ii]],
                                       bruit,
                                       trt$edge.length[j],
                                       Tpst[[ii]][[j]])
        Mt[[ii]][[ifr]]=transitionmatrix(Xt[[ii]][[ifr]], Trposs[[ii]], nph[ii], bruit, trt$edge.length[ifr])
        loiapplistet[[ii]][[ifr]]=loiapp(Xt[[ii]][[ifr]],
                                         Lt[Lt[, ifr] > 0, ifr],
                                         loiini[[ii]],
                                         Trposs[[ii]],
                                         bruit,
                                         trt$edge.length[ifr],
                                         Tpst[[ii]][[ifr]])
      }
      
      corrb=corrb-sum(dbinom(
        NLt[sss2,j],size=NL[sss2,j]-1,prob=u/(L[sss2,j]*ww),
        log = T
      ))
      corrb = corrb - sum(dbeta(Lt[sss2[sss3],j], NLt[sss2[sss3],j], 1, log = T))
    } else {
      ww=agemax-hauteur(tr,a)
      trt$root[[2]] = c(trt$root[[2]][-which(trt$root[[2]] == a)], d)
      for (ii in 1:nch) {
        Xt[[ii]][[ifr]] = sample(1:length(P[[ii]]),
                                 rpois(1, La[ii] * u),
                                 prob = P[[ii]],
                                 replace = T)
        Tpst[[ii]][[ifr]] = sort(runif(length(Xt[[ii]][[ifr]])))
        corrb = corrb - log(prod(P[[ii]][Xt[[ii]][[ifr]]])) - dpois(length(Xt[[ii]][[ifr]]), La[ii] *
                                                                      u, log = T)
      }
      NLt[, ifr] = rpois(nrow(Dat[[1]]), u * Rho)
      Lt[, ifr] = rbeta(nrow(Dat[[1]]), NLt[, ifr], 1)
      for (ii in 1:nch){
        Mt[[ii]][[ifr]]=transitionmatrix(Xt[[ii]][[ifr]], Trposs[[ii]], nph[ii], bruit, trt$edge.length[ifr])
        loiapplistet[[ii]][[ifr]]=loiapp(Xt[[ii]][[ifr]],
                                         Lt[Lt[, ifr] > 0, ifr],
                                         loiini[[ii]],
                                         Trposs[[ii]],
                                         bruit,
                                         trt$edge.length[ifr],
                                         Tpst[[ii]][[ifr]])
      }
      corrb = corrb - sum(dpois(NLt[, ifr], u * Rho, log = T)) - sum(dbeta(Lt[Lt[, ifr] !=
                                                                                0, ifr], NLt[Lt[, ifr] != 0, ifr], 1, log = T))
    }
    if (a %in% root[[1]]) {
      trt$root[[1]] = c(trt$root[[1]][-which(trt$root[[1]] == a)], d)
    }
    if (a %in% root[[2]]){
      quoisur=rep(0,nch)
    } else {
      quoisur=sapply(1:nch,function(x){length(X[[x]][[j]])})+.5*(a<=ncol(Dat[[1]])&y %in%root[[2]])
    }
    if (testarbre3(trt, trt$root[[2]], nnodes, feuilles)) {
      return(list(trt, Xt, NLt, Lt, Tpst, corrb, Mt,loiapplistet,d,f,quoisur,La,Rho,bruit))
    } else {
      return(list(tr, X, NL, L, Tps, 0,M,loiappliste,-1,-1,rep(-1,nch),La,Rho,bruit))
    }
  }



modificationtopofrere <-
  function(tr,
           X,
           NL,
           L,
           Tps,
           P,
           Rho,
           La,
           nnodes,
           agemax,
           feuilles,
           Dat,
           tipprior,
           M,
           loiappliste,
           Trposs,
           bruit,
           nph) {
    #on prend un noeud e non racine et on le déplace pour qu'il ait a comme frère
    root = tr$root
    Xt = X
    NLt = NL
    Lt = L
    Tpst = Tps
    trt = tr
    nch = length(X)
    corrb = 0
    Mt=M
    loiapplistet=loiappliste
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
      #corrective term for the rescaling, it simplifies as the rescaling in one sense is the opposite of the rescaling in the other sense
      resc=jacrescale*(length(unique(tr$edge[st,1])))*log((hauteur(tr,a)+ww-u)/(hauteur(tr,e)+tr$edge.length[i]))
      for (ii in 1:nch) {
        Xt[[ii]][[k]] = c(X[[ii]][[k]], X[[ii]][[ifr]])
        Tpst[[ii]][[k]] = c(Tps[[ii]][[k]] * tr$edge.length[k], tr$edge.length[k]+ Tps[[ii]][[ifr]] *
                              tr$edge.length[ifr]) / (tr$edge.length[k] + tr$edge.length[ifr])
      }
      Lt[, k] = (
        (L[, ifr] * tr$edge.length[ifr]+ tr$edge.length[k]) * (NL[, ifr] > 0) + (NL[, ifr] == 0 &
                                                                                   NL[, k] > 0) * L[, k] * tr$edge.length[k]
      ) / (tr$edge.length[k] + tr$edge.length[ifr])
      NLt[, k] = NL[, k] + NL[, ifr]
      #on rajoute le laplacien réciproque qui est lié à l'action de découpe de la branche en 2 cf ce qu'il y a après
      rrr=which(NLt[,k]>1 & Lt[,k]>(tr$edge.length[k])/trt$edge.length[k])#ceux pour lesquels il y a au moins une transformation qui pourrait être mise sur la branche en question
      for (ii in 1:nch){
        Mt[[ii]][[k]]=transitionmatrix(Xt[[ii]][[k]], Trposs[[ii]], nph[ii], bruit, trt$edge.length[k])
        loiapplistet[[ii]][[k]]=loiapp(Xt[[ii]][[k]],
                                       Lt[Lt[, k] > 0, k],
                                       loiini[[ii]],
                                       Trposs[[ii]],
                                       bruit,
                                       trt$edge.length[k],
                                       Tpst[[ii]][[k]])
      }
      corrb = corrb  +log(1/(tr$edge.length[k] + tr$edge.length[ifr]))- log(1/tr$edge.length[j]) + sum(dbinom(NL[rrr,k],
                                                                                                              size=NLt[rrr,k]-1,
                                                                                                              prob=tr$edge.length[k]/(tr$edge.length[k] + L[rrr,ifr]*tr$edge.length[ifr])))+
        sum(dbeta(L[L[, k] > 0 & L[,ifr]>0, k], NL[NL[, k] >0 & NL[,ifr] > 0, k], 1, log = T))#en suis-je sûr ?
      if (d %in% root[[1]]) {
        trt$root[[1]] = c(trt$root[[1]][-which(trt$root[[1]] == d)], f)
      }
      for (ii in 1:nch) {
        cbtr = Tps[[ii]][[j]] < u / ww
        Xt[[ii]][[j]] = X[[ii]][[j]][cbtr]
        Xt[[ii]][[ifr]] = X[[ii]][[j]][!cbtr]#puisque e a toujours un frère.
        Tpst[[ii]][[j]] = ww * (Tps[[ii]][[j]][cbtr]) / u
        Tpst[[ii]][[ifr]] = (Tps[[ii]][[j]][!cbtr]*ww - u)/(ww-u)
      }
      NLt[,j]=0
      Lt[,j]=0
      NLt[,ifr]=0
      Lt[,ifr]=0
      #quand il n'y a qu'une app ok
      sss1=which(NL[,j]==1 & L[,j]<u/ww)
      sss2=which(NL[,j]==1 & L[,j]>u/ww)
      NLt[sss1,j]=1
      NLt[sss2,ifr]=1
      Lt[sss1,j]=L[sss1,j]*ww/u
      Lt[sss2,ifr]=(L[sss2,j]*ww-u)/(ww-u)
      #maintenant les trucs problématiques, lorsqu'il y a plusieurs app
      sss1=which(NL[,j]>1 & L[,j]<u/ww)
      sss2=which(NL[,j]>1 & L[,j]>u/ww)
      NLt[sss1,j]=NL[sss1,j]
      NLt[sss2,j]=rbinom(length(sss2),size=NL[sss2,j]-1,prob=u/(L[sss2,j]*ww))
      NLt[sss2,ifr]=NL[sss2,j]-NLt[sss2,j]
      Lt[sss1,j]=L[sss1,j]*ww/u
      Lt[sss2,ifr]=(L[sss2,j]*ww-u)/(ww-u)
      sss3=which(NLt[sss2,j]>0)
      Lt[sss2[sss3],j]=rbeta(length(sss3),NLt[sss2[sss3],j],1)
      
      for (ii in 1:nch){
        Mt[[ii]][[j]]=transitionmatrix(Xt[[ii]][[j]], Trposs[[ii]], nph[ii], bruit, trt$edge.length[j])
        loiapplistet[[ii]][[j]]=loiapp(Xt[[ii]][[j]],
                                       Lt[Lt[, j] > 0, j],
                                       loiini[[ii]],
                                       Trposs[[ii]],
                                       bruit,
                                       trt$edge.length[j],
                                       Tpst[[ii]][[j]])
        Mt[[ii]][[ifr]]=transitionmatrix(Xt[[ii]][[ifr]], Trposs[[ii]], nph[ii], bruit, trt$edge.length[ifr])
        loiapplistet[[ii]][[ifr]]=loiapp(Xt[[ii]][[ifr]],
                                         Lt[Lt[, ifr] > 0, ifr],
                                         loiini[[ii]],
                                         Trposs[[ii]],
                                         bruit,
                                         trt$edge.length[ifr],
                                         Tpst[[ii]][[ifr]])
      }
      
      corrb=corrb-sum(dbinom(
        NLt[sss2,j],size=NL[sss2,j]-1,prob=u/(L[sss2,j]*ww),
        log = T
      ))
      corrb = corrb - sum(dbeta(Lt[sss2[sss3],j], NLt[sss2[sss3],j], 1, log = T))
      if (a %in% root[[1]]) {
        trt$root[[1]] = c(trt$root[[1]][-which(trt$root[[1]] == a)], d)
      }
      for (ii in 1:nch) {
        Xt[[ii]][[i]] = sample(1:length(P[[ii]]),
                               rpois(1, La[ii] * tr$edge.length[i]),
                               prob = P[[ii]],
                               replace = T)
        Tpst[[ii]][[i]] = sort(runif(length(Xt[[ii]][[i]])))
        corrb = corrb - log(prod(P[[ii]][Xt[[ii]][[i]]])) - dpois(length(Xt[[ii]][[i]]), La[ii] *
                                                                    trt$edge.length[i], log = T)
        corrb = corrb + log(prod(P[[ii]][X[[ii]][[i]]])) + dpois(length(X[[ii]][[i]]), La[ii] *
                                                                   tr$edge.length[i], log = T)
      }
      NLt[, i] = rpois(nrow(Dat[[1]]), trt$edge.length[i] * Rho)
      Lt[, i] = rbeta(nrow(Dat[[1]]), NLt[, i], 1)
      
      for (ii in 1:nch){
        Mt[[ii]][[i]]=transitionmatrix(Xt[[ii]][[i]], Trposs[[ii]], nph[ii], bruit, trt$edge.length[i])
        loiapplistet[[ii]][[i]]=loiapp(Xt[[ii]][[i]],
                                       Lt[Lt[, i] > 0, i],
                                       loiini[[ii]],
                                       Trposs[[ii]],
                                       bruit,
                                       trt$edge.length[i],
                                       Tpst[[ii]][[i]])
      }
      
      corrb = corrb - sum(dpois(NLt[, i], trt$edge.length[i] * Rho, log =
                                  T)) - sum(dbeta(Lt[Lt[, i] != 0, i], NLt[Lt[, i] != 0, i], 1, log = T))
      corrb = corrb + sum(dpois(NL[, i], tr$edge.length[i] * Rho, log =
                                  T)) + sum(dbeta(L[L[, i] != 0, i], NL[L[, i] != 0, i], 1, log = T))
      
      if (a %in% root[[2]]){
        quoisur=rep(0,nch)
      } else {
        quoisur=sapply(1:nch,function(x){length(X[[x]][[j]])})+.5*(a<=ncol(Dat[[1]])&y %in%root[[2]])
      }
      if (testarbre3(trt, trt$root[[2]], nnodes, feuilles)) {
        return(list(trt, Xt, NLt, Lt, Tpst, corrb+resc, Mt,loiapplistet,d,f,quoisur,La,Rho,bruit))
      } else {
        return(list(tr, X, NL, L, Tps, 0,M,loiappliste,-1,-1,rep(-1,nch),La,Rho,bruit))
      }} else {
        return(list(tr, X, NL, L, Tps, 0,M,loiappliste,-1,-1,rep(-1,nch),La,Rho,bruit))
      }
  }


modificationtopofrere2 <-
  function(tr,
           X,
           NL,
           L,
           Tps,
           P,
           Rho,
           La,
           nnodes,
           agemax,
           feuilles,
           Dat,
           tipprior,
           M,
           loiappliste,
           Trposs,
           bruit,
           nph) {
    #on prend un noeud e non racine et on le déplace pour qu'il ait son cousin comme frère
    #version sans changement des variables latentes
    root = tr$root
    Xt = X
    NLt = NL
    Lt = L
    Tpst = Tps
    trt = tr
    nch = length(X)
    corrb = 0
    Mt=M
    loiapplistet=loiappliste
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
      resc=jacrescale*(length(unique(tr$edge[st,1])))*log((hauteur(tr,a)+ww-u)/(hauteur(tr,e)+tr$edge.length[i]))
      for (ii in 1:nch) {
        Xt[[ii]][[k]] = c(X[[ii]][[k]], X[[ii]][[ifr]])
        Tpst[[ii]][[k]] = c(Tps[[ii]][[k]] * tr$edge.length[k], tr$edge.length[k]+ Tps[[ii]][[ifr]] *
                              tr$edge.length[ifr]) / (tr$edge.length[k] + tr$edge.length[ifr])
      }
      Lt[, k] = (
        (L[, ifr] * tr$edge.length[ifr]+ tr$edge.length[k]) * (NL[, ifr] > 0) + (NL[, ifr] == 0 &
                                                                                   NL[, k] > 0) * L[, k] * tr$edge.length[k]
      ) / (tr$edge.length[k] + tr$edge.length[ifr])
      NLt[, k] = NL[, k] + NL[, ifr]
      #on rajoute le laplacien réciproque qui est lié à l'action de découpe de la branche en 2 cf ce qu'il y a après
      rrr=which(NLt[,k]>1 & Lt[,k]>(tr$edge.length[k])/trt$edge.length[k])#ceux pour lesquels il y a au moins une transformation qui pourrait être mise sur la branche en question
      for (ii in 1:nch){
        Mt[[ii]][[k]]=transitionmatrix(Xt[[ii]][[k]], Trposs[[ii]], nph[ii], bruit, trt$edge.length[k])
        loiapplistet[[ii]][[k]]=loiapp(Xt[[ii]][[k]],
                                       Lt[Lt[, k] > 0, k],
                                       loiini[[ii]],
                                       Trposs[[ii]],
                                       bruit,
                                       trt$edge.length[k],
                                       Tpst[[ii]][[k]])
      }
      corrb = corrb +log(1/(tr$edge.length[k] + tr$edge.length[ifr])) - log(1/tr$edge.length[j]) + sum(dbinom(NL[rrr,k],
                                                                                                              size=NLt[rrr,k]-1,
                                                                                                              prob=tr$edge.length[k]/(tr$edge.length[k] + L[rrr,ifr]*tr$edge.length[ifr])))+
        sum(dbeta(L[L[, k] > 0 & L[,ifr]>0, k], NL[NL[, k] > 0 & NL[,ifr]>0, k], 1, log = T))#en suis-je sûr ?
      if (d %in% root[[1]]) {
        trt$root[[1]] = c(trt$root[[1]][-which(trt$root[[1]] == d)], f)
      }
      for (ii in 1:nch) {
        cbtr = Tps[[ii]][[j]] < u / ww
        Xt[[ii]][[j]] = X[[ii]][[j]][cbtr]
        Xt[[ii]][[ifr]] = X[[ii]][[j]][!cbtr]#puisque e a toujours un frère.
        Tpst[[ii]][[j]] = ww * (Tps[[ii]][[j]][cbtr]) / u
        Tpst[[ii]][[ifr]] = (Tps[[ii]][[j]][!cbtr]*ww - u)/(ww-u)
      }
      NLt[,j]=0
      Lt[,j]=0
      NLt[,ifr]=0
      Lt[,ifr]=0
      #quand il n'y a qu'une app ok
      sss1=which(NL[,j]==1 & L[,j]<u/ww)
      sss2=which(NL[,j]==1 & L[,j]>u/ww)
      NLt[sss1,j]=1
      NLt[sss2,ifr]=1
      Lt[sss1,j]=L[sss1,j]*ww/u
      Lt[sss2,ifr]=(L[sss2,j]*ww-u)/(ww-u)
      #maintenant les trucs problématiques, lorsqu'il y a plusieurs app
      sss1=which(NL[,j]>1 & L[,j]<u/ww)
      sss2=which(NL[,j]>1 & L[,j]>u/ww)
      NLt[sss1,j]=NL[sss1,j]
      NLt[sss2,j]=rbinom(length(sss2),size=NL[sss2,j]-1,prob=u/(L[sss2,j]*ww))
      NLt[sss2,ifr]=NL[sss2,j]-NLt[sss2,j]
      Lt[sss1,j]=L[sss1,j]*ww/u
      Lt[sss2,ifr]=(L[sss2,j]*ww-u)/(ww-u)
      sss3=which(NLt[sss2,j]>0)
      Lt[sss2[sss3],j]=rbeta(length(sss3),NLt[sss2[sss3],j],1)
      
      for (ii in 1:nch){
        Mt[[ii]][[j]]=transitionmatrix(Xt[[ii]][[j]], Trposs[[ii]], nph[ii], bruit, trt$edge.length[j])
        loiapplistet[[ii]][[j]]=loiapp(Xt[[ii]][[j]],
                                       Lt[Lt[, j] > 0, j],
                                       loiini[[ii]],
                                       Trposs[[ii]],
                                       bruit,
                                       trt$edge.length[j],
                                       Tpst[[ii]][[j]])
        Mt[[ii]][[ifr]]=transitionmatrix(Xt[[ii]][[ifr]], Trposs[[ii]], nph[ii], bruit, trt$edge.length[ifr])
        loiapplistet[[ii]][[ifr]]=loiapp(Xt[[ii]][[ifr]],
                                         Lt[Lt[, ifr] > 0, ifr],
                                         loiini[[ii]],
                                         Trposs[[ii]],
                                         bruit,
                                         trt$edge.length[ifr],
                                         Tpst[[ii]][[ifr]])
      }
      
      corrb=corrb-sum(dbinom(
        NLt[sss2,j],size=NL[sss2,j]-1,prob=u/(L[sss2,j]*ww),
        log = T
      ))
      corrb = corrb - sum(dbeta(Lt[sss2[sss3],j], NLt[sss2[sss3],j], 1, log = T))
      if (a %in% root[[2]]){
        quoisur=rep(0,nch)
      } else {
        quoisur=sapply(1:nch,function(x){length(X[[x]][[j]])})+.5*(a<=ncol(Dat[[1]])&y %in%root[[2]])
      }
      if (testarbre3(trt, trt$root[[2]], nnodes, feuilles)) {
        return(list(trt, Xt, NLt, Lt, Tpst, corrb, Mt,loiapplistet,d,f,quoisur,La,Rho,bruit))
      } else {
        return(list(tr, X, NL, L, Tps, 0,M,loiappliste,-1,-1,rep(-1,nch),La,Rho,bruit))
      }} else {
        return(list(tr, X, NL, L, Tps, 0,M,loiappliste,-1,-1,rep(-1,nch),La,Rho,bruit))
      }
  }

testarbre <- function(tr, root, nnodes) {
  #voir la durée de la fonction et l'enlever si non necessaire.
  #teste l'arbre
  visites = numeric()
  mvt = root
  i = 1
  while (length(mvt) > 0 & i < 2 * nnodes) {
    u = sample(1:length(mvt), 1)
    w = tr$edge[tr$edge[, 1] == mvt[u], 2]
    visites = unique(c(visites, mvt[u]))
    mvt = c(mvt[-u], w[w != mvt[u]])
    i = i + 1
  }
  return(length(visites) == nnodes)
}

testarbre2 = function(tr, root, nnodes, feuilles) {
  # tentative de test plus rapide
  ok = all(dim(tr$edge) == c(nnodes - length(root), 2))
  ok = ok & all(1:nnodes == sort(c(tr$edge[, 2], root)))
  #ok = ok & (root == (nnodes + 3) / 2)
  ok = ok &
    all(sort(tr$edge[, 1]) == rep((1:nnodes)[-unique(c(feuilles, root))], each = 2))
  return(ok)
}

testarbre3 = function(tr, root, nnodes, feuilles) {
  qui = root
  actif = root
  for (i in 1:floor(nnodes / 2)) {
    actif = unlist(sapply(actif, function(x) {
      tr$edge[which(tr$edge[, 1] == x), 2]
    }))
    qui = c(qui, actif)
  }
  if (length(qui) == nnodes) {
    return(all(sort(qui) == 1:nnodes) & all(tr$edge.length>0))
  }
  else {
    return(FALSE)
  }
}

# testarbre3 = function(tr,root,nnodes,feuilles){
#   qui=root
#   actif=root
#   for (i in 1:floor(nnodes/2)){
#     actif=unlist(sapply(actif, function(x){tr$edge[which(tr$edge[,1]==x),2]}))
#     qui=c(qui,actif)
#   }
#     return(all(sort(qui)==1:nnodes))
# }

testhauteur <- function(tr, root, bon) {
  #teste que toutes les hauteurs sont bonnes
  VV = sapply(1:(root - 1), function(x) {
    remonteehauteur(tr, x, 0, root)
  })
  all(VV == bon)
}

remonteehauteur <- function(tr, i, h, root) {
  #fonction rec auxiliaire
  if (i == root) {
    return(h)
  } else {
    q = which(tr$edge[, 2] == i)
    return(remonteehauteur(tr, tr$edge[q, 1], h + tr$edge.length[q], root))
  }
}

nouvelletopo <- function(tr, X, root, nnodes, nch) {
  #génère un arbre aléatoire totalement, non utilisée
  TT = rcoal(root - 1, rooted = T)
  l = hauteur(TT, root)
  TT$edge.length = TT$edge.length / l * hauteur(tr, root)
  TT$tip.label = 1:(root - 1)
  return(list(TT, 0))
}

sousarbre <- function(tr, i) {
  #renvoie les indices des branches dans le sous arbre partant de i
  if (sum(tr$edge[, 1] == i) > 0) {
    w = which(tr$edge[, 1] == i)
    return(c(w[1], w[2], sousarbre(tr, tr$edge[w[1], 2]), sousarbre(tr, tr$edge[w[2], 2])))
  } else {
    return(numeric())
  }
}

feuillessousarbre <- function(tr,i){
  if (sum(tr$edge[, 1] == i) > 0) {
    w = which(tr$edge[, 1] == i)
    return(c(feuillessousarbre(tr, tr$edge[w[1], 2]), feuillessousarbre(tr, tr$edge[w[2], 2])))
  } else {
    return(i)
  }
}


hauteur <- function(tr, i) {
  #renvoie la hauteur du sous arbre partant de i
  if (sum(tr$edge[, 1] == i) > 0) {
    w = which(tr$edge[, 1] == i)
    return(tr$edge.length[w[1]] + hauteur(tr, tr$edge[w[1], 2]))
  } else {
    return(0)
  }
}

randomtree <- function(n) {
  #crée un arbre aléatoire.
  randomtreerec(1:n, n, n + 2, c())
}

randomtreerecpart <- function(noeuds, n, k, branches, nrac) {
  #on part d'une liste de noeuds et on en recolle deux uniformément (cf coalescence), avec n nombre de noeuds
  noeudst = noeuds
  branchest = branches
  if (length(noeuds) == nrac + 1) {
    branchest = rbind(branches, c(n + 1, noeuds[1]), c(n + 1, noeuds[2]))
    return(branchest)
  } else {
    r = sample(1:length(noeuds), 2)
    noeudst = noeuds[-r]
    noeudst = c(noeudst, k)
    branchest = rbind(branchest, c(k, noeuds[r[1]]), c(k, noeuds[r[2]]))
    return(randomtreerec(noeudst, n, k + 1, branchest))
  }
}

hauteurtous <- function(tr,feuilles,agemax) {
  hauts=numeric(max(tr$edge))
  hauts2=numeric(max(tr$edge))
  for (i in 1:length(feuilles)){
    act=feuilles[i]
    while(length(which(tr$edge[,2]==act))>0){
      j=which(tr$edge[,2]==act)
      uu=tr$edge[j,1]
      hauts[uu]=tr$edge.length[j]+hauts[act]
      hauts2[act]=hauts[uu]
      act=uu
    }
  }
  hauts2[act]=agemax
  return(list(hauts,hauts2))
}


modificationtopooutgroup <-
  function(tr,
           X,
           NL,
           L,
           Tps,
           P,
           Rho,
           La,
           nnodes,
           agemax,
           feuilles,
           Dat,
           tipprior,
           M,
           loiappliste,
           Trposs,
           bruit,
           nph,
           Prior,
           fin) {
    #on prend un noeud e non racine et on le déplace pour qu'il ait a comme frère
    root = tr$root
    Xt = X
    NLt = NL
    Lt = L
    Tpst = Tps
    trt = tr
    nch = length(X)
    corrb = 0
    Mt=M
    loiapplistet=loiappliste
    if (runif(1)<1/2){ #on transforme en outgroup
      #we do not change the latent variables on the subtree starting from e
      e=sample((1:nnodes)[-c(root[[2]],tr$edge[tr$edge[,1]==root[[2]],2])], 1)
      i = which(tr$edge[, 2] == e)
      d = tr$edge[i, 1]#parent de celui ci
      st = sousarbre(tr, e)#sous arbre qui en part
      ifr = which(tr$edge[, 1] == d & !tr$edge[, 2] == e)#branche frère
      f = tr$edge[ifr, 2]#noeud frère
      g=tr$edge[tr$edge[,2]==d,1]#noeud grand père
      j=which(tr$edge[,2]==d) #branche du grand père 
      trt$edge[j,]=c(d,root[[2]])
      trt$edge[ifr,]=c(g,f)
      trt$edge.length[ifr]=tr$edge.length[j]+tr$edge.length[ifr]
      trt$edge.length[j]=rexp(1,1/agemax)
      trt$edge.length[i]=trt$edge.length[j]+agemax-hauteur(tr,e)
      scale=agemax/(trt$edge.length[j]+agemax)
      corrb=corrb-dexp(trt$edge.length[j],1/agemax,log=T)
      trt$edge.length=trt$edge.length*scale
      trt$root[[2]]=d
      corrb=corrb+log(1/(hauteur(tr,g)-max(hauteur(tr,e),hauteur(tr,f))))
      haut=hauteur(trt,e)
      hauts=hauteurtous(trt,feuilles,agemax)
      quiposs=which(hauts[[2]]>haut & (1:max(tr$edge)) != d & (1:max(tr$edge)) != f)
      corrb = corrb - log(nrow(tr$edge)-2) + log(length(quiposs))
      resc=jacrescale*((nnodes-length(feuilles)-1)*log(scale))
      for (ii in 1:nch) {
        Xt[[ii]][[ifr]] = c(X[[ii]][[j]], X[[ii]][[ifr]])
        Tpst[[ii]][[ifr]] = c(Tps[[ii]][[j]] * tr$edge.length[j], tr$edge.length[j]+ Tps[[ii]][[ifr]] *
                                tr$edge.length[ifr]) / (tr$edge.length[j] + tr$edge.length[ifr])
      }
      Lt[, ifr] = ((tr$edge.length[j] + L[, ifr] * tr$edge.length[ifr]) * as.integer(NL[, ifr] > 0) + 
                     (NL[, ifr] == 0 &NL[, j] > 0) * L[, j] * tr$edge.length[j]
      ) / (tr$edge.length[j] + tr$edge.length[ifr])
      NLt[, ifr] = NL[, j] + NL[, ifr]
      rrr=which(NLt[,ifr]>1 & Lt[,ifr]>(tr$edge.length[j])/(tr$edge.length[ifr]+tr$edge.length[j]))#ceux pour lesquels il y a au moins une transformation qui pourrait être mise sur la branche en question
      for (ii in 1:nch){
        Mt[[ii]][[ifr]]=transitionmatrix(Xt[[ii]][[ifr]], Trposs[[ii]], nph[ii], bruit, trt$edge.length[ifr])
        loiapplistet[[ii]][[ifr]]=loiapp(Xt[[ii]][[ifr]],
                                         Lt[Lt[, ifr] > 0, ifr],
                                         loiini[[ii]],
                                         Trposs[[ii]],
                                         bruit,
                                         trt$edge.length[ifr],
                                         Tpst[[ii]][[ifr]])
      }
      corrb = corrb + sum(dbinom(NL[rrr,ifr],
                                 size=NLt[rrr,ifr]-1,
                                 prob=tr$edge.length[j]/(tr$edge.length[j] + L[rrr,ifr]*tr$edge.length[ifr])))+
        sum(dbeta(L[L[rrr, j] > 0, j], NL[NL[rrr, j] > 0, j], 1, log = T))
      for (ii in 1:nch) {
        for (jj in c(i)){
          Xt[[ii]][[jj]] = sample(1:length(P[[ii]]),
                                  rpois(1, La[ii] * trt$edge.length[jj]),
                                  prob = P[[ii]],
                                  replace = T)
          Tpst[[ii]][[jj]] = sort(runif(length(Xt[[ii]][[jj]])))
          corrb = corrb - log(prod(P[[ii]][Xt[[ii]][[jj]]])) - dpois(length(Xt[[ii]][[jj]]), La[ii] *
                                                                       trt$edge.length[jj], log = T)
          corrb = corrb + log(prod(P[[ii]][X[[ii]][[jj]]])) + dpois(length(X[[ii]][[jj]]), La[ii] *
                                                                      tr$edge.length[jj], log = T)
          
          
        }
        Xt[[ii]][[j]] = sample(1:length(P[[ii]]),
                               rpois(1, La[ii] * trt$edge.length[j]),
                               prob = P[[ii]],
                               replace = T)
        Tpst[[ii]][[j]] = sort(runif(length(Xt[[ii]][[j]])))
        corrb = corrb - log(prod(P[[ii]][Xt[[ii]][[j]]])) - dpois(length(Xt[[ii]][[j]]), La[ii] *
                                                                    trt$edge.length[j], log = T)
        
      }
      for (jj in c(i)){
        NLt[, jj] = rpois(nrow(Dat[[1]]), trt$edge.length[jj] * Rho)
        Lt[, jj] = rbeta(nrow(Dat[[1]]), NLt[, jj], 1)
        corrb=corrb-sum(dpois(NLt[,jj], trt$edge.length[jj] * Rho,log=T)) + sum(dpois(NL[,jj], tr$edge.length[jj] * Rho,log=T))
        corrb=corrb-sum(dbeta(Lt[Lt[,jj]>0,jj],  NLt[NLt[,jj]>0, jj], 1,log=T)) + sum(dbeta(L[L[,jj]>0,jj],  NL[NL[,jj]>0, jj], 1,log=T))
      }
      NLt[, j] = rpois(nrow(Dat[[1]]), trt$edge.length[j] * Rho)
      Lt[, j] = rbeta(nrow(Dat[[1]]), NLt[, j], 1)
      corrb=corrb-sum(dpois(NLt[,j], trt$edge.length[j] * Rho,log=T))
      corrb=corrb-sum(dbeta(Lt[Lt[,j]>0,j],  NLt[NLt[,j]>0, j], 1,log=T))
      for (ii in 1:nch){
        for (jj in 1:length(M[[1]])){
          Mt[[ii]][[jj]]=transitionmatrix(Xt[[ii]][[jj]], Trposs[[ii]], nph[ii], bruit, trt$edge.length[jj])
          loiapplistet[[ii]][[jj]]=loiapp(Xt[[ii]][[jj]],
                                          Lt[Lt[, jj] > 0, jj],
                                          loiini[[ii]],
                                          Trposs[[ii]],
                                          bruit,
                                          trt$edge.length[jj],
                                          Tpst[[ii]][[jj]])
        }
      }
    } else { #on bouge l'outgroup
      aaa=which(tr$edge[,1]==root[[2]])#the node that moves is the youngest son of the root
      e=tr$edge[aaa[which.max(tr$edge.length[aaa])],2]#the youngest is the one with longest branch at the root
      f=tr$edge[aaa[which.min(tr$edge.length[aaa])],2]
      st = sousarbre(tr, e)#sous arbre qui en part
      haut_e=hauteur(tr,e)
      hauts=hauteurtous(tr,feuilles,agemax)
      quiposs=which(hauts[[2]]>haut_e & (1:max(tr$edge)) != tr$root[[2]] & (1:max(tr$edge))!=f)
      if (length(tr$edge[,1]==f)==0 || length(quiposs)<1){
        return(list(tr, X, NL, L, Tps, 0,M,loiappliste,-1,-1,rep(-1,nch),La,Rho,bruit))
      } else {
        a=quiposs[sample(1:length(quiposs),1)]
        z=tr$edge[tr$edge[,2]==a,1]
        hk=which(tr$edge[,2]==a)#branche sur laquelle on va recoller le sous arbre e
        corrb = corrb - log(length(quiposs)) + log(nrow(tr$edge) - 2)
        i = which(tr$edge[, 2] == e)
        d = tr$edge[i, 1]#parent de celui ci
        ifr = which(tr$edge[, 1] == d & !tr$edge[, 2] == e)#branche frère
        f = tr$edge[ifr, 2]#neoud frère
        j=which(tr$edge[,2]==e)
        scale=agemax/hauteur(tr,f)
        trt$edge[ifr,]=c(z,root[[2]])
        trt$edge[hk,]=c(root[[2]],a)
        haut=runif(1,0,hauteur(tr,z)-max(hauteur(tr,e),hauteur(tr,a)))
        trt$edge.length[ifr]=haut
        trt$edge.length[hk]=tr$edge.length[hk]-haut
        trt$edge.length[i]=hauteur(tr,z)-haut-hauteur(tr,e)
        trt$edge.length=trt$edge.length*scale
        corrb=corrb+dexp(tr$edge.length[ifr]*scale,1/agemax,log=T)#la longeur de la branche qui n'existe plus à la racine avant transformation, on ne passe jamais par cet état donc il faut corriger
        trt$root[[2]]=f
        corrb=corrb-log(1/(min(trt$edge.length[ifr]+trt$edge.length[i],trt$edge.length[hk]+trt$edge.length[ifr])))
        resc=jacrescale*((nnodes-length(feuilles)-1)*log(scale))
        for (ii in 1:nch) {
          cbtr = Tps[[ii]][[hk]] < haut / tr$edge.length[hk]
          Xt[[ii]][[ifr]] = X[[ii]][[hk]][cbtr]
          Xt[[ii]][[hk]] = X[[ii]][[hk]][!cbtr]#puisque e a toujours un frère.
          Tpst[[ii]][[ifr]] = tr$edge.length[hk] * (Tps[[ii]][[hk]][cbtr])/haut
          Tpst[[ii]][[hk]] = (Tps[[ii]][[hk]][!cbtr]*tr$edge.length[hk] - haut)/(tr$edge.length[hk]-haut)
        }
        for (ii in 1:nch) {  
          corrb = corrb + log(prod(P[[ii]][X[[ii]][[ifr]]])) + dpois(length(X[[ii]][[ifr]]), La[ii] *
                                                                       tr$edge.length[ifr], log = T)
          
        }
        NLt[,hk]=0
        Lt[,hk]=0
        NLt[,ifr]=0
        Lt[,ifr]=0
        #quand il n'y a qu'une app ok
        sss1=which(NL[,hk]==1 & L[,hk]<haut / tr$edge.length[hk])
        sss2=which(NL[,hk]==1 & L[,hk]>haut / tr$edge.length[hk])
        NLt[sss1,ifr]=1
        NLt[sss2,hk]=1
        Lt[sss1,ifr]=L[sss1,hk]*tr$edge.length[hk]/haut
        Lt[sss2,hk]=(L[sss2,hk]*tr$edge.length[hk]-haut)/(tr$edge.length[hk]-haut)
        #maintenant les trucs problématiques, lorsqu'il y a plusieurs app
        sss1=which(NL[,hk]>1 & L[,hk]<haut / tr$edge.length[hk])
        sss2=which(NL[,hk]>1 & L[,hk]>haut / tr$edge.length[hk])
        NLt[sss1,ifr]=NL[sss1,hk]
        NLt[sss2,hk]=rbinom(length(sss2),size=NL[sss2,hk]-1,prob=(L[sss2,hk]*tr$edge.length[hk] - haut)/(L[sss2,hk]*tr$edge.length[hk]))+1
        NLt[sss2,ifr]=NL[sss2,hk]-NLt[sss2,hk]
        Lt[sss1,ifr]=L[sss1,hk]*haut/tr$edge.length[hk]
        Lt[sss2,hk]=(L[sss2,hk]*tr$edge.length[hk]-haut)/(tr$edge.length[hk]-haut)
        sss3=which(NLt[sss2,ifr]>0)
        Lt[sss2[sss3],ifr]=rbeta(length(sss3),NLt[sss2[sss3],ifr],1)
        
        corrb=corrb-sum(dbinom(
          NLt[sss2,hk]-1,size=NL[sss2,hk]-1,
          prob=(L[sss2,hk]*tr$edge.length[hk] - haut)/(L[sss2,hk]*tr$edge.length[hk]),
          log = T
        ))
        corrb = corrb - sum(dbeta(Lt[sss2[sss3],ifr], NLt[sss2[sss3],ifr], 1, log = T))
        
        corrb=corrb+ sum(dpois(NL[,ifr], tr$edge.length[ifr] * Rho,log=T))
        corrb=corrb+ sum(dbeta(L[L[,ifr]>0,ifr],  NL[NL[,ifr]>0, ifr], 1,log=T))
        for (ii in 1:nch) {
          for (jj in c(i)){
            Xt[[ii]][[jj]] = sample(1:length(P[[ii]]),
                                    rpois(1, La[ii] * trt$edge.length[jj]),
                                    prob = P[[ii]],
                                    replace = T)
            Tpst[[ii]][[jj]] = sort(runif(length(Xt[[ii]][[jj]])))
            corrb = corrb - log(prod(P[[ii]][Xt[[ii]][[jj]]])) - dpois(length(Xt[[ii]][[jj]]), La[ii] *
                                                                         trt$edge.length[jj], log = T)
            corrb = corrb + log(prod(P[[ii]][X[[ii]][[jj]]])) + dpois(length(X[[ii]][[jj]]), La[ii] *
                                                                        tr$edge.length[jj], log = T)
            
            
          }
        }
        for (jj in c(i)){
          NLt[, jj] = rpois(nrow(Dat[[1]]), trt$edge.length[jj] * Rho)
          Lt[, jj] = rbeta(nrow(Dat[[1]]), NLt[, jj], 1)
          corrb=corrb-sum(dpois(NLt[,jj], trt$edge.length[jj] * Rho,log=T)) + sum(dpois(NL[,jj], tr$edge.length[jj] * Rho,log=T))
          corrb=corrb-sum(dbeta(Lt[Lt[,jj]>0,jj],  NLt[NLt[,jj]>0, jj], 1,log=T)) + sum(dbeta(L[L[,jj]>0,jj],  NL[NL[,jj]>0, jj], 1,log=T))
        }
        for (ii in 1:nch){
          for (jj in 1:length(M[[1]])){
            Mt[[ii]][[jj]]=transitionmatrix(Xt[[ii]][[jj]], Trposs[[ii]], nph[ii], bruit, trt$edge.length[jj])
            loiapplistet[[ii]][[jj]]=loiapp(Xt[[ii]][[jj]],
                                            Lt[Lt[, jj] > 0, jj],
                                            loiini[[ii]],
                                            Trposs[[ii]],
                                            bruit,
                                            trt$edge.length[jj],
                                            Tpst[[ii]][[jj]])
          }
        }
        
      }
      
    }
    if (testarbre3(trt, trt$root[[2]], nnodes, feuilles) & !is.na(corrb)) {
      return(list(trt, Xt, NLt, Lt, Tpst, corrb+resc, Mt,loiapplistet,0,0,rep(.2,nch),La,Rho,bruit))
    } else {
      return(list(tr, X, NL, L, Tps, 0,M,loiappliste,-1,-1,rep(-1,nch),La,Rho,bruit))
    }
  }


modificationtopooutgroupsansresc <-
  function(tr,
           X,
           NL,
           L,
           Tps,
           P,
           Rho,
           La,
           nnodes,
           agemax,
           feuilles,
           Dat,
           tipprior,
           M,
           loiappliste,
           Trposs,
           bruit,
           nph,
           Prior,
           fin) {
    #on prend un noeud e non racine et on le déplace pour qu'il ait a comme frère
    root = tr$root
    Xt = X
    NLt = NL
    Lt = L
    Tpst = Tps
    trt = tr
    nch = length(X)
    corrb = 0
    Mt=M
    resc=0
    loiapplistet=loiappliste
    if (runif(1)<1/2){ #on transforme en outgroup
      #we do not change the latent variables on the subtree starting from e
      e=sample((1:nnodes)[-c(root[[2]],tr$edge[tr$edge[,1]==root[[2]],2])])[1]
      i = which(tr$edge[, 2] == e)
      d = tr$edge[i, 1]#parent de celui ci
      st = sousarbre(tr, e)#sous arbre qui en part
      ifr = which(tr$edge[, 1] == d & !tr$edge[, 2] == e)#branche frère
      f = tr$edge[ifr, 2]#neoud frère
      g=tr$edge[tr$edge[,2]==d,1]#noeud grand père
      j=which(tr$edge[,2]==d) #branche du grand père 
      trt$edge[j,]=c(d,root[[2]])
      trt$edge[ifr,]=c(g,f)
      trt$edge.length[ifr]=tr$edge.length[j]+tr$edge.length[ifr]
      bbb=which(trt$edge[,1]==root[[2]])
      hautmax=min(trt$edge.length[bbb])#minimum length at the root, that we split in two
      trt$edge.length[j]=runif(1,0,hautmax)
      corrb=corrb-log(1/hautmax)
      trt$edge.length[bbb]=trt$edge.length[bbb]-trt$edge.length[j]
      trt$edge.length[i]=agemax-hauteur(tr,e)
      trt$root[[2]]=d
      #corrb=corrb+log(1/(hauteur(trt,g)-max(hauteur(trt,e),hauteur(trt,f))))#legth of the possible branh inside
      haut=hauteur(trt,e)
      hauts=hauteurtous(trt,feuilles,agemax)
      quiposs=which(hauts[[2]]>haut & (1:max(tr$edge)) != d & (1:max(tr$edge))!=tr$root[[2]]) 
      corrb=corrb+log((nrow(tr$edge)-2))-log(length(quiposs))#on the other direction, correction for the new position of e
      
      for (ii in 1:nch) {
        Xt[[ii]][[ifr]] = c(X[[ii]][[j]], X[[ii]][[ifr]])
        Tpst[[ii]][[ifr]] = c(Tps[[ii]][[j]] * tr$edge.length[j], tr$edge.length[j]+ Tps[[ii]][[ifr]] *
                                tr$edge.length[ifr]) / (tr$edge.length[j] + tr$edge.length[ifr])
      }
      Lt[, ifr] = ((tr$edge.length[j] + L[, ifr] * tr$edge.length[ifr]) * as.integer(NL[, ifr] > 0) + 
                     (NL[, ifr] == 0 &NL[, j] > 0) * L[, j] * tr$edge.length[j]
      ) / (tr$edge.length[j] + tr$edge.length[ifr])
      NLt[, ifr] = NL[, j] + NL[, ifr]
      rrr=which(NLt[,ifr]>1 & Lt[,ifr]>(tr$edge.length[j])/(tr$edge.length[ifr]+tr$edge.length[j]))#ceux pour lesquels il y a au moins une transformation qui pourrait être mise sur la branche en question
      corrb = corrb + sum(dbinom(NL[rrr,j],
                                 size=NLt[rrr,ifr]-1,
                                 prob=tr$edge.length[j]/(tr$edge.length[j] + L[rrr,ifr]*tr$edge.length[ifr])))+
        sum(dbeta(L[L[, j] > 0 & L[,ifr]>0, j], NL[NL[, j] > 0 & NL[,ifr]>0, j], 1, log = T))
      for (ii in 1:nch) {
        Xt[[ii]][[i]] = sample(1:length(P[[ii]]),
                               rpois(1, La[ii] * trt$edge.length[i]),
                               prob = P[[ii]],
                               replace = T)
        Tpst[[ii]][[i]] = sort(runif(length(Xt[[ii]][[i]])))
        corrb = corrb - sum(log(P[[ii]][Xt[[ii]][[i]]])) - dpois(length(Xt[[ii]][[i]]), La[ii] *
                                                                   trt$edge.length[i], log = T)
        corrb = corrb + sum(log(P[[ii]][X[[ii]][[i]]])) + dpois(length(X[[ii]][[i]]), La[ii] *
                                                                  tr$edge.length[i], log = T)
        Xt[[ii]][[j]] = sample(1:length(P[[ii]]),
                               rpois(1, La[ii] * trt$edge.length[j]),
                               prob = P[[ii]],
                               replace = T)
        Tpst[[ii]][[j]] = sort(runif(length(Xt[[ii]][[j]])))
        corrb = corrb - sum(log(P[[ii]][Xt[[ii]][[j]]])) - dpois(length(Xt[[ii]][[j]]), La[ii] *
                                                                   trt$edge.length[j], log = T)
        
      }
      NLt[, i] = rpois(nrow(Dat[[1]]), trt$edge.length[i] * Rho)
      Lt[, i] = rbeta(nrow(Dat[[1]]), NLt[, i], 1)
      corrb=corrb-sum(dpois(NLt[,i], trt$edge.length[i] * Rho,log=T)) + sum(dpois(NL[,i], tr$edge.length[i] * Rho,log=T))
      corrb=corrb-sum(dbeta(Lt[Lt[,i]>0,i],  NLt[NLt[,i]>0, i], 1,log=T)) + sum(dbeta(L[L[,i]>0,i],  NL[NL[,i]>0, i], 1,log=T))
      NLt[, j] = rpois(nrow(Dat[[1]]), trt$edge.length[j] * Rho)
      Lt[, j] = rbeta(nrow(Dat[[1]]), NLt[, j], 1)
      corrb=corrb-sum(dpois(NLt[,j], trt$edge.length[j] * Rho,log=T))
      corrb=corrb-sum(dbeta(Lt[Lt[,j]>0,j],  NLt[NLt[,j]>0, j], 1,log=T))
      for (ii in 1:nch){
        for (jj in 1:length(M[[1]])){
          Mt[[ii]][[jj]]=transitionmatrix(Xt[[ii]][[jj]], Trposs[[ii]], nph[ii], bruit, trt$edge.length[jj])
          loiapplistet[[ii]][[jj]]=loiapp(Xt[[ii]][[jj]],
                                          Lt[Lt[, jj] > 0, jj],
                                          loiini[[ii]],
                                          Trposs[[ii]],
                                          bruit,
                                          trt$edge.length[jj],
                                          Tpst[[ii]][[jj]])
        }
      }
      #dans l'autre sens on fait sur trc les modif qu'on aurait faites, seul e ne change pas
      at=f
      zt=tr$edge[trt$edge[,2]==at,1]#nœud parent de a
      hkt=which(trt$edge[,2]==at)#branche sur laquelle on va recoller le sous arbre e
      it = which(trt$edge[, 2] == e)
      dt = trt$edge[i, 1]#parent de celui ci, c'est en fait la racine
      ifrt = which(trt$edge[, 1] == dt & !trt$edge[, 2] == e)#branche frère
      ft = tr$edge[ifrf, 2]#neoud frère
      trc$edge[ifrt,]=c(zt,dt)
      trc$edge[hkt,]=c(dt,at)
      hautt=runif(1,0,hauteur(trt,zt)-max(hauteur(trt,e),hauteur(trt,at)))
      trc$edge.length[ifrt]=hautt
      trc$edge.length[hkt]=tr$edge.length[hkt]-hautt
      trc$edge.length[it]=hauteur(trt,zt)-hautt-hauteur(trt,e) 
      ccc=which(trc$edge[,1]==ft)
      trc$edge.length[ccc]=trc$edge.length[ccc]+agemax-hauteur(trc,ft)    
      corrb=corrb+log(1/(hauteur(trt,zt)-max(hauteur(trt,e),hauteur(trt,at))))
      
      
    } else { #on bouge l'outgroup
      aaa=which(tr$edge[,1]==root[[2]])#the node that moves is the youngest son of the root
      e=tr$edge[aaa[which.max(tr$edge.length[aaa])],2]#the youngest is the one with longest branch at the root
      f=tr$edge[aaa[which.min(tr$edge.length[aaa])],2]
      st = sousarbre(tr, e)#sous arbre qui en part
      haut_e=hauteur(tr,e)
      hauts=hauteurtous(tr,feuilles,agemax)
      quiposs=which(hauts[[2]]>haut_e & (1:max(tr$edge)) != tr$root[[2]] & (1:max(tr$edge))!=f) 
      if (length(tr$edge[,1]==f)==0 || length(quiposs)<1){
        return(list(tr, X, NL, L, Tps, 0,M,loiappliste,-1,-1,rep(-1,nch),La,Rho,bruit))
      } else {
        a=quiposs[sample(1:length(quiposs),1)]
        z=tr$edge[tr$edge[,2]==a,1]#nœud parent de a
        hk=which(tr$edge[,2]==a)#branche sur laquelle on va recoller le sous arbre e
        corrb=corrb -log((nrow(tr$edge)-2)) +log(length(quiposs)) #number of nodes to choose
        i = which(tr$edge[, 2] == e)
        d = tr$edge[i, 1]#parent de celui ci, c'est en fait la racine
        st = sousarbre(tr, e)#sous arbre qui en part
        ifr = which(tr$edge[, 1] == d & !tr$edge[, 2] == e)#branche frère
        f = tr$edge[ifr, 2]#neoud frère
        trt$edge[ifr,]=c(z,root[[2]])
        trt$edge[hk,]=c(root[[2]],a)
        haut=runif(1,0,hauteur(tr,z)-max(hauteur(tr,e),hauteur(tr,a)))
        trt$edge.length[ifr]=haut
        trt$edge.length[hk]=tr$edge.length[hk]-haut
        trt$edge.length[i]=hauteur(tr,z)-haut-hauteur(tr,e) 
        ccc=which(trt$edge[,1]==f)
        trt$edge.length[ccc]=trt$edge.length[ccc]+agemax-hauteur(trt,f)
        #bbb=which(tr$edge[,1]==f)
        corrb=corrb-log(1/(hauteur(tr,z)-max(hauteur(tr,e),hauteur(tr,a))))
        #+log(1/(min(tr$edge.length[bbb])+tr$edge.length[ifr]))
        trt$root[[2]]=f
        trc=trt
        a=
          for (ii in 1:nch) {
            cbtr = Tps[[ii]][[hk]] < haut / tr$edge.length[hk]
            Xt[[ii]][[ifr]] = X[[ii]][[hk]][cbtr]
            Xt[[ii]][[hk]] = X[[ii]][[hk]][!cbtr]
            Tpst[[ii]][[ifr]] = tr$edge.length[hk] * (Tps[[ii]][[hk]][cbtr])/haut
            Tpst[[ii]][[hk]] = (Tps[[ii]][[hk]][!cbtr]*tr$edge.length[hk] - haut)/(tr$edge.length[hk]-haut)
          }
        for (ii in 1:nch) {  
          corrb = corrb + log(prod(P[[ii]][X[[ii]][[ifr]]])) + dpois(length(X[[ii]][[ifr]]), La[ii] *
                                                                       tr$edge.length[ifr], log = T)
          
        }
        NLt[,hk]=0
        Lt[,hk]=0
        NLt[,ifr]=0
        Lt[,ifr]=0
        #quand il n'y a qu'une app ok
        sss1=which(NL[,hk]==1 & L[,hk]<haut / tr$edge.length[hk])
        sss2=which(NL[,hk]==1 & L[,hk]>haut / tr$edge.length[hk])
        NLt[sss1,ifr]=1
        NLt[sss2,hk]=1
        Lt[sss1,ifr]=L[sss1,hk]*tr$edge.length[hk]/haut
        Lt[sss2,hk]=(L[sss2,hk]*tr$edge.length[hk]-haut)/(tr$edge.length[hk]-haut)
        #maintenant les trucs problématiques, lorsqu'il y a plusieurs app
        sss1=which(NL[,hk]>1 & L[,hk]<haut / tr$edge.length[hk])
        sss2=which(NL[,hk]>1 & L[,hk]>haut / tr$edge.length[hk])
        NLt[sss1,ifr]=NL[sss1,hk]
        NLt[sss2,hk]=rbinom(length(sss2),size=NL[sss2,hk]-1,prob=(L[sss2,hk]*tr$edge.length[hk] - haut)/(L[sss2,hk]*tr$edge.length[hk]))+1
        NLt[sss2,ifr]=NL[sss2,hk]-NLt[sss2,hk]
        Lt[sss1,ifr]=L[sss1,hk]*tr$edge.length[hk]/haut
        Lt[sss2,hk]=(L[sss2,hk]*tr$edge.length[hk]-haut)/(tr$edge.length[hk]-haut)
        sss3=which(NLt[sss2,ifr]>0)
        Lt[sss2[sss3],ifr]=rbeta(length(sss3),NLt[sss2[sss3],ifr],1)
        corrb=corrb-sum(dbinom(
          NLt[sss2,hk]-1,size=NL[sss2,hk]-1,
          prob=(L[sss2,hk]*tr$edge.length[hk] - haut)/(L[sss2,hk]*tr$edge.length[hk]),
          log = T
        ))
        corrb = corrb - sum(dbeta(Lt[sss2[sss3],ifr], NLt[sss2[sss3],ifr], 1, log = T))
        corrb=corrb+ sum(dpois(NL[,ifr], tr$edge.length[ifr] * Rho,log=T))
        corrb=corrb+ sum(dbeta(L[L[,ifr]>0,ifr],  NL[NL[,ifr]>0, ifr], 1,log=T))
        
        for (ii in 1:nch) {
          Xt[[ii]][[i]] = sample(1:length(P[[ii]]),
                                 rpois(1, La[ii] * trt$edge.length[i]),
                                 prob = P[[ii]],
                                 replace = T)
          Tpst[[ii]][[i]] = sort(runif(length(Xt[[ii]][[i]])))
          corrb = corrb - log(prod(P[[ii]][Xt[[ii]][[i]]])) - dpois(length(Xt[[ii]][[i]]), La[ii] *
                                                                      trt$edge.length[i], log = T)
          corrb = corrb + log(prod(P[[ii]][X[[ii]][[i]]])) + dpois(length(X[[ii]][[i]]), La[ii] *
                                                                     tr$edge.length[i], log = T)
          
          
        }
        NLt[, i] = rpois(nrow(Dat[[1]]), trt$edge.length[i] * Rho)
        Lt[, i] = rbeta(nrow(Dat[[1]]), NLt[, i], 1)
        corrb=corrb-sum(dpois(NLt[,i], trt$edge.length[i] * Rho,log=T)) + sum(dpois(NL[,i], tr$edge.length[i] * Rho,log=T))
        corrb=corrb-sum(dbeta(Lt[Lt[,i]>0,i],  NLt[NLt[,i]>0, i], 1,log=T)) + sum(dbeta(L[L[,i]>0,i],  NL[NL[,i]>0, i], 1,log=T))
        for (ii in 1:nch){
          for (jj in 1:length(M[[1]])){
            Mt[[ii]][[jj]]=transitionmatrix(Xt[[ii]][[jj]], Trposs[[ii]], nph[ii], bruit, trt$edge.length[jj])
            loiapplistet[[ii]][[jj]]=loiapp(Xt[[ii]][[jj]],
                                            Lt[Lt[, jj] > 0, jj],
                                            loiini[[ii]],
                                            Trposs[[ii]],
                                            bruit,
                                            trt$edge.length[jj],
                                            Tpst[[ii]][[jj]])
          }
        }
        trc=trt
        it = which(trt$edge[, 2] == e)
        dt = trt$edge[it, 1]#parent de celui ci
        ifrt = which(trt$edge[, 1] == dt & !trt$edge[, 2] == e)#branche frère
        ft = trt$edge[ifrt, 2]#neoud frère
        gt=trt$edge[trt$edge[,2]==dt,1]#noeud grand père
        jt=which(trt$edge[,2]==dt) #branche du grand père 
        trc$edge[j,]=c(dt,f)
        trc$edge[ifr,]=c(gt,ft)
        trc$edge.length[ifrt]=tr$edge.length[j]+tr$edge.length[ifr]
        bbb=which(trt$edge[,1]==root[[2]])
        hautmax=min(trt$edge.length[bbb])#minimum length at the root, that we split in two
        trt$edge.length[j]=runif(1,0,hautmax)
        corrb=corrb-log(1/hautmax)
      }
      
    }
    #as we have updated the parameters we recompute their values rwt the prior
    if (testarbre3(trt, trt$root[[2]], nnodes, feuilles) & !is.na(corrb) 
        & (all(1:nnodes==sort(c(trt$edge[,2],trt$root[[2]]))))
        & (all((length(feuilles)+1):nnodes==sort(unique(trt$edge[,1]))))) {
      return(list(trt, Xt, NLt, Lt, Tpst, corrb, Mt,loiapplistet,0,0,rep(.2,nch),La,Rho,bruit))
    } else {
      return(list(tr, X, NL, L, Tps, 0,M,loiappliste,-1,-1,rep(-1,nch),La,Rho,bruit))
    }
  }
