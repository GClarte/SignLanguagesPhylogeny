
modificationtopooutgroup2 <-
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
           Prior) {
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
      e=sample((1:nnodes)[-c(root[[2]],tr$edge[tr$edge[,1]==root[[2]],2])])[1]
      i = which(tr$edge[, 2] == e)
      d = tr$edge[i, 1]#parent de celui ci
      haut=hauteur(tr,d)
      #hauts=hauteurtous(tr,feuilles,agemax)
      st = sousarbre(tr, e)#sous arbre qui en part
      corr=-log((nnodes - length(st)-2 )/(nnodes-3))
      ifr = which(tr$edge[, 1] == d & !tr$edge[, 2] == e)#branche frère
      f = tr$edge[ifr, 2]#neoud frère
      g=tr$edge[tr$edge[,2]==d,1]#noeud grand père
      j=which(tr$edge[,2]==d) #branche du grand père 
      scale=rbeta(1,10,1)
      corr=corr-dbeta(scale,10,1)
      Lat=1/scale*La
      Rhot=1/scale*Rho
      bruitt=bruit/scale
      trt$edge[j,]=c(d,root[[2]])
      trt$edge[ifr,]=c(g,f)
      trt$edge.length[ifr]=tr$edge.length[j]+tr$edge.length[ifr]
      corr=corr+log(1/(tr$edge.length[j]+tr$edge.length[ifr]))
      trt$edge.length[j]=agemax*(1-scale)
      trt$edge.length[-c(st,j,i)]=trt$edge.length[-c(st,j,i)]*scale
      scale2=runif(1,0,1)
      trt$edge.length[st]=(1-scale2)*tr$edge.length[st]/hauteur(tr,e)*agemax
      trt$edge.length[i]=agemax-hauteur(trt,e)
      trt$root[[2]]=d
      corr=corr+log(1/(tr$edge.length[j]+tr$edge.length[ifr]))
      for (ii in 1:nch) {
        Xt[[ii]][[ifr]] = c(X[[ii]][[j]], X[[ii]][[ifr]])
        Tpst[[ii]][[ifr]] = c(Tps[[ii]][[j]] * tr$edge.length[j], tr$edge.length[j]+ Tps[[ii]][[ifr]] *
                                tr$edge.length[ifr]) / (tr$edge.length[j] + tr$edge.length[ifr])
      }
      Lt[, ifr] = ((tr$edge.length[j] + L[, ifr] * tr$edge.length[ifr]) * as.integer(NL[, ifr] > 0) + 
                     (NL[, ifr] == 0 &NL[, j] > 0) * L[, j] * tr$edge.length[j]
      ) / (tr$edge.length[j] + tr$edge.length[ifr])
      NLt[, ifr] = NL[, j] + NL[, ifr]
      rrr=which(NLt[,ifr]>1 & Lt[,ifr]>(tr$edge.length[j])/trt$edge.length[ifr])#ceux pour lesquels il y a au moins une transformation qui pourrait être mise sur la branche en question
      for (ii in 1:nch){
        Mt[[ii]][[ifr]]=transitionmatrix(Xt[[ii]][[ifr]], Trposs[[ii]], nph[ii], bruitt, trt$edge.length[ifr])
        loiapplistet[[ii]][[ifr]]=loiapp(Xt[[ii]][[ifr]],
                                         Lt[Lt[, ifr] > 0, ifr],
                                         loiini[[ii]],
                                         Trposs[[ii]],
                                         bruitt,
                                         trt$edge.length[ifr],
                                         Tpst[[ii]][[ifr]])
      }
      corrb = corrb + sum(dbinom(NL[rrr,ifr],
                                 size=NLt[rrr,ifr]-1,
                                 prob=tr$edge.length[j]/(tr$edge.length[j] + L[rrr,ifr]*tr$edge.length[ifr])))+
        sum(dbeta(L[L[, j] > 0, j], NL[NL[, j] > 0, j], 1, log = T))
      for (ii in 1:nch) {
        Xt[[ii]][[i]] = sample(1:length(P[[ii]]),
                               rpois(1, Lat[ii] * trt$edge.length[i]),
                               prob = P[[ii]],
                               replace = T)
        Tpst[[ii]][[i]] = sort(runif(length(Xt[[ii]][[i]])))
        Xt[[ii]][[j]] = sample(1:length(P[[ii]]),
                               rpois(1, Lat[ii] * trt$edge.length[j]),
                               prob = P[[ii]],
                               replace = T)
        Tpst[[ii]][[j]] = sort(runif(length(Xt[[ii]][[j]])))
        corrb = corrb - log(prod(P[[ii]][Xt[[ii]][[i]]])) - dpois(length(Xt[[ii]][[i]]), Lat[ii] *
                                                                    trt$edge.length[i], log = T)
        corrb = corrb + log(prod(P[[ii]][X[[ii]][[i]]])) + dpois(length(X[[ii]][[i]]), La[ii] *
                                                                   tr$edge.length[i], log = T)
        corrb = corrb - log(prod(P[[ii]][Xt[[ii]][[j]]])) - dpois(length(Xt[[ii]][[j]]), Lat[ii] *
                                                                    trt$edge.length[j], log = T)
        
      }
      
      NLt[, j] = rpois(nrow(Dat[[1]]), trt$edge.length[j] * Rhot)
      Lt[, j] = rbeta(nrow(Dat[[1]]), NLt[, j], 1)
      corr=corr-sum(dpois(NLt[,j], trt$edge.length[j] * Rhot,log=T)) 
      corr=corr-sum(dbeta(Lt[Lt[,j]>0,j],  NLt[NLt[,j]>0, j], 1,log=T)) 
      NLt[, i] = rpois(nrow(Dat[[1]]), trt$edge.length[i] * Rhot)
      Lt[, i] = rbeta(nrow(Dat[[1]]), NLt[, i], 1)
      corr=corr-sum(dpois(NLt[,i], trt$edge.length[i] * Rhot,log=T)) + sum(dpois(NL[,i], tr$edge.length[i] * Rho,log=T))
      corr=corr-sum(dbeta(Lt[Lt[,i]>0,i],  NLt[NLt[,i]>0, i], 1,log=T)) + sum(dbeta(L[L[,i]>0,i],  NL[NL[,i]>0, i], 1,log=T))
      
      for (ii in 1:nch){
        for (jj in c(i,j)){
          Mt[[ii]][[jj]]=transitionmatrix(Xt[[ii]][[jj]], Trposs[[ii]], nph[ii], bruitt, trt$edge.length[jj])
          loiapplistet[[ii]][[jj]]=loiapp(Xt[[ii]][[jj]],
                                          Lt[Lt[, jj] > 0, jj],
                                          loiini[[ii]],
                                          Trposs[[ii]],
                                          bruitt,
                                          trt$edge.length[jj],
                                          Tpst[[ii]][[jj]])
        }
      }
    } else { #on bouge l'outgroup
      aaa=sample(tr$edge[tr$edge[,1]==root[[2]],2]) #on choisit l'outgroup à bouger
      e=aaa[1]
      f=aaa[2]
      st = sousarbre(tr, e)#sous arbre qui en part
      pfff=(1:nnodes)[-c(tr$edge[st,2],e,root[[2]],f)]#noeuds où on peut mettre le sous arbre^
      if (length(tr$edge[,1]==f)==0 || length(pfff)<1){
        return(list(tr, X, NL, L, Tps, 0,M,loiappliste,-1,-1,rep(-1,nch),La,Rho,bruit))
      } else {a=sample(pfff,1)#direction, nouveau frère qui n'est pas dans st ni dans racine et noeud à bouger
      z=tr$edge[tr$edge[,2]==a,1]
      hk=which(tr$edge[,2]==a)#branche sur laquelle on va recoller le sous arbre e
      corr= -log((nnodes-3)/length(pfff)) #pour la correction
      i = which(tr$edge[, 2] == e)
      d = tr$edge[i, 1]#parent de celui ci
      st = sousarbre(tr, e)#sous arbre qui en part
      ifr = which(tr$edge[, 1] == d & !tr$edge[, 2] == e)#branche frère
      f = tr$edge[ifr, 2]#neoud frère
      j=which(tr$edge[,2]==e)
      trt$root[[2]]=f
      h1=hauteur(tr,f)
      h2=hauteur(tr,e)
      corr=corr+dbeta(h1/agemax,10,1)
      scale=h1/agemax
      corr=corr-dbeta(scale,10,1)
      Lat=scale*La
      Rhot=scale*Rho
      bruitt=bruit*scale
      trt$edge[hk,]=c(z,root[[2]])
      trt$edge[ifr,]=c(root[[2]],a)
      haut=runif(1,0,tr$edge.length[hk])
      corr=corr-length(1/tr$edge.length[hk])
      trt$edge.length[hk]=haut
      trt$edge.length[ifr]=tr$edge.length[hk]-haut
      trt$edge.length[-c(st,i)]=trt$edge.length[-c(st,i)]/scale
      trt$edge.length[c(st,i)]=(hauteur(tr,z)-haut)*trt$edge.length[c(st,i)]/agemax
      trt$edge.length[st]=(1-scale)*tr$edge.length[st]
      trt$edge.length[i]=scale*agemax
      trt$root[[2]]=f
      for (ii in 1:nch) {
        cbtr = Tps[[ii]][[hk]] < haut / tr$edge.length[hk]
        Xt[[ii]][[hk]] = X[[ii]][[hk]][cbtr]
        Xt[[ii]][[ifr]] = X[[ii]][[hk]][!cbtr]#puisque e a toujours un frère.
        Tpst[[ii]][[hk]] = haut / tr$edge.length[hk] * (Tps[[ii]][[hk]][cbtr])
        Tpst[[ii]][[ifr]] = (1 - haut / tr$edge.length[hk]) * (Tps[[ii]][[hk]][!cbtr])
      }
      NLt[,hk]=0
      Lt[,hk]=0
      NLt[,ifr]=0
      Lt[,ifr]=0
      #quand il n'y a qu'une app ok
      sss1=which(NL[,hk]==1 & L[,hk]<haut / tr$edge.length[hk])
      sss2=which(NL[,hk]==1 & L[,hk]>haut / tr$edge.length[hk])
      NLt[sss1,hk]=1
      NLt[sss2,ifr]=1
      Lt[sss1,hk]=L[sss1,hk]*haut / tr$edge.length[hk]
      Lt[sss2,ifr]=(L[sss2,hk]*tr$edge.length[hk]-haut)/(tr$edge.length[hk]-haut)
      #maintenant les trucs problématiques, lorsqu'il y a plusieurs app
      sss1=which(NL[,hk]>1 & L[,hk]<haut / tr$edge.length[hk])
      sss2=which(NL[,hk]>1 & L[,hk]>haut / tr$edge.length[hk])
      NLt[sss1,hk]=NL[sss1,hk]
      NLt[sss2,hk]=rbinom(length(sss2),size=NL[sss2,hk]-1,prob=haut/(L[sss2,hk]*tr$edge.length[hk]))
      NLt[sss2,ifr]=NL[sss2,hk]-NLt[sss2,hk]
      Lt[sss1,hk]=L[sss1,hk]*haut/tr$edge.length[hk]
      Lt[sss2,ifr]=(L[sss2,hk]*tr$edge.length[hk]-haut)/(tr$edge.length[hk]-haut)
      sss3=which(NLt[,hk]>0)
      Lt[sss3,hk]=rbeta(length(sss3),NLt[sss3,hk],1)
      
      #on recompte le fait que dans l'autre sens on doit rajouter les transf sur ifr
      corr=corr+ sum(dpois(NL[,ifr], tr$edge.length[ifr] * Rho,log=T))
      corr=corr+ sum(dbeta(L[L[,ifr]>0,ifr],  NL[NL[,ifr]>0, ifr], 1,log=T))
      for (ii in 1:nch) {  
        corrb = corrb + log(prod(P[[ii]][X[[ii]][[ifr]]])) + dpois(length(X[[ii]][[ifr]]), La[ii] *
                                                                     tr$edge.length[ifr], log = T)
        
      }
      
      for (ii in 1:nch){
        Mt[[ii]][[hk]]=transitionmatrix(Xt[[ii]][[hk]], Trposs[[ii]], nph[ii], bruit, trt$edge.length[hk])
        loiapplistet[[ii]][[hk]]=loiapp(Xt[[ii]][[hk]],
                                        Lt[Lt[, hk] > 0, hk],
                                        loiini[[ii]],
                                        Trposs[[ii]],
                                        bruit,
                                        trt$edge.length[hk],
                                        Tpst[[ii]][[hk]])
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
        NLt[sss2,hk],size=NL[sss2,hk]-1,prob=haut/(L[sss2,hk]*tr$edge.length[hk]),
        log = T
      ))
      corrb = corrb - sum(dbeta(Lt[sss3,hk], NLt[sss3,hk], 1, log = T))
      for (ii in 1:nch) {
        for (jj in c(i)){
          Xt[[ii]][[jj]] = sample(1:length(P[[ii]]),
                                  rpois(1, Lat[ii] * trt$edge.length[jj]),
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
        NLt[, jj] = rpois(nrow(Dat[[1]]), trt$edge.length[jj] * Rhot)
        Lt[, jj] = rbeta(nrow(Dat[[1]]), NLt[, jj], 1)
        corr=corr-sum(dpois(NLt[,jj], trt$edge.length[jj] * Rhot,log=T)) + sum(dpois(NL[,jj], tr$edge.length[jj] * Rho,log=T))
        corr=corr-sum(dbeta(Lt[Lt[,jj]>0,jj],  NLt[NLt[,jj]>0, jj], 1,log=T)) + sum(dbeta(L[L[,jj]>0,jj],  NL[NL[,jj]>0, jj], 1,log=T))
      }
      for (ii in 1:nch){
        for (jj in c(i)){
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
    corr=corr-sum(sapply(1:nch,function(z){dgamma(Lat[z],Prior$Prila[[z]][1],Prior$Prila[[z]][2],log=T)-
        dgamma(La[z],Prior$Prila[[z]][1],Prior$Prila[[z]][2],log=T)}))+
      dgamma(Rhot,Prior$Prirho[1],Prior$Prirho[2])-dgamma(Rho,Prior$Prirho[1],Prior$Prirho[2])+
      Prior$pribruit[[1]][[2]](bruitt)-Prior$pribruit[[1]][[2]](bruit)
    if (testarbre3(trt, trt$root[[2]], nnodes, feuilles) & !is.na(corrb+corr)) {
      return(list(trt, Xt, NLt, Lt, Tpst, corrb+corr, Mt,loiapplistet,0,0,rep(.25,nch),Lat,Rhot,bruitt))
    } else {
      return(list(tr, X, NL, L, Tps, 0,M,loiappliste,-1,-1,rep(-1,nch),La,Rho,bruit))
    }
  }