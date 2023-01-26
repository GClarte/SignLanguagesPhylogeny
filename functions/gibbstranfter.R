#on ajoute ou on retire une transformation, on passe une transformation au voisin

#X liste des transformations sur chaque branche,
#S liste des transf possibles
#renvoie 1 si ajout/retrait réussi 2 si ajout retrait raté, 3 si transfert réussi, 4 si transfert raté

gibbstransf <-
  function(tr,
           Dat,
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
           loiappliste,
           Tps) {
    Xc = X
    Mc = M
    Tpsc = Tps
    u = log(runif(1))
    v = sample(3:4,1,prob=c(5,1))#choix entre ajout, retrait, et prior
    #v=3
    prob=c(rep(1,length(tr$edge.length)-2),rep(10,2))
    if (v == 2)
      #retrait aléatoire
    {
      w = sample(1:length(tr$edge.length), 1,prob = prob)#choix de la branche
      I = tr$edge[w, 1]#noeud de début de branche
      J = tr$edge[w, 2]#noeud de fin
      N = length(X[[w]])
      if (N > 0) {
        qui = sample(1:N, 1)
        Xc[[w]] = Xc[[w]][-qui]
        Tpsc[[w]] = Tpsc[[w]][-qui]
        Mc[[w]] = transitionmatrix(Xc[[w]], Trposs, nph, bruit, tr$edge.length[[w]])
        path = cheminrac(tr, I, root)
        loiapplistec = loiappliste
        loiapplistec[[w]] = loiapp(Xc[[w]], L[L[, w] > 0, w], loiini, Trposs, bruit, tr$edge.length[w], Tpsc[[w]])
        listc = pruningeco(tr,
                           Xc,
                           Mc,
                           Dat,
                           nph,
                           path,
                           Lkldinter,
                           Trposs,
                           L,
                           bruit,
                           loiini,
                           loiapplistec)
        #listc=pruninglkldlistini(tr,Xc,Mc,L,Dat,nph,length(M),root,Trposs,loiini,bruit)
        if (ncogn>0){
          lkld1 = sum(sapply(root,function(xx){log(loiini %*% listc[[xx]][, 1:ncogn])-log(loiini %*%
                                                                                            Lkldinter[[xx]][, 1:ncogn])}))
        } else {
          lkld1=0
        }
        lkld2 = dpois(N - 1, tr$edge.length[w] * La, log = T) - dpois(N, tr$edge.length[w] *
                                                                        La, log = T)
        datestemp = c(0, Tps[[w]], 1)
        #calcul du jacobien
        #quand on enlève, on a juste à choisir laquelle, quand on ajoute, on choisit où, qui et quand
        q = log(pr[X[[w]][qui]] * (datestemp[qui + 2] - datestemp[qui])) +log(N)
      } else {
        return(list(X, M, 0, Lkldinter, loiappliste, Tps))
      }
    } else if (v == 1) {
      w = sample(1:length(tr$edge.length), 1,prob = prob)#choix de la branche
      k = sample(1:length(Trposs), 1, prob = pr)#choix de la transf
      I = tr$edge[w, 1]#noeud de début de branche
      J = tr$edge[w, 2]#noeud de fin
      if (length(Xc[[w]]) == 0) {
        Xc[[w]] = k
        Tpsc[[w]] = runif(1)
        ubi = .5
        datestemp = c(0, 1)
      } else {
        ubi = sample(0:length(Xc[[w]]), 1) + .5
        datestemp = c(0, Tpsc[[w]], 1)
        pp = 1:length(Xc[[w]])
        Xc[[w]] = c(Xc[[w]][pp < ubi], k, Xc[[w]][pp > ubi])
        Tpsc[[w]] = c(Tpsc[[w]][pp < ubi], runif(1, datestemp[ubi + .5], datestemp[ubi +
                                                                                     1.5]), Tpsc[[w]][pp > ubi])
      }
      Mc[[w]] = transitionmatrix(Xc[[w]], Trposs, nph, bruit, tr$edge.length[w])
      path = cheminrac(tr, I, root)
      loiapplistec = loiappliste
      loiapplistec[[w]] = loiapp(Xc[[w]], L[L[, w] > 0, w], loiini, Trposs, bruit, tr$edge.length[w], Tpsc[[w]])
      listc = pruningeco(tr,
                         Xc,
                         Mc,
                         Dat,
                         nph,
                         path,
                         Lkldinter,
                         Trposs,
                         L,
                         bruit,
                         loiini,
                         loiapplistec)
      if (ncogn>0){
        lkld1 = sum(sapply(root,function(xx){log(loiini %*% listc[[xx]][, 1:ncogn])-log(loiini %*%
                                                                                          Lkldinter[[xx]][, 1:ncogn])}))
        
      } else {
        lkld1=0
      }
      q = -log(pr[k] / ((-datestemp[ubi + .5] + datestemp[ubi + 1.5])))-log(length(Xc[[w]]))
      lkld2 = dpois(length(Xc[[w]]), tr$edge.length[w] * La, log = T) - dpois(length(X[[w]]), tr$edge.length[w] *
                                                                                La, log = T)
    } else if (v == 3) {
      w = sample(1:length(tr$edge.length), 1,prob = prob)#choix de la branche
      I = tr$edge[w, 1]#noeud de début de branche
      J = tr$edge[w, 2]#noeud de fin
      Xc[[w]] = sample(
        1:length(Trposs),
        rpois(1, tr$edge.length[w] * La),
        replace = T,
        prob = pr
      )
      Q = length(Xc[[w]])
      if (Q > 0) {
        Y = cumsum(rexp(Q + 1, 1))
        Tpsc[[w]] = Y[1:Q] / Y[Q + 1]
      }
      else{
        Tpsc[[w]] = numeric(0)
      }
      Mc[[w]] = transitionmatrix(Xc[[w]], Trposs, nph, bruit, tr$edge.length[w])
      path = cheminrac(tr, I, root)
      loiapplistec = loiappliste
      loiapplistec[[w]] = loiapp(Xc[[w]], L[L[, w] > 0, w], loiini, Trposs, bruit, tr$edge.length[w], Tpsc[[w]])
      listc = pruningeco(tr,
                         Xc,
                         Mc,
                         Dat,
                         nph,
                         path,
                         Lkldinter,
                         Trposs,
                         L,
                         bruit,
                         loiini,
                         loiapplistec)
      #listc=pruninglkldlistini(tr,Xc,Mc,L,Dat,nph,nnodes,root,Trposs,loiini,bruit)
      if (ncogn>0){
        lkld1 = sum(sapply(root,function(xx){log(loiini %*% listc[[xx]][, 1:ncogn])-log(loiini %*%
                                                                                          Lkldinter[[xx]][, 1:ncogn])}))
      } else {
        lkld1=0
      }
      q = 0
      lkld2 = 0
    } else if (v==4){
      loiapplistec = loiappliste
      for (ki in 1:length(Xc)){
        Xc[[ki]]=sample(
          1:length(Trposs),
          rpois(1, tr$edge.length[ki] * La),
          replace = T,
          prob = pr
        )
        Q = length(Xc[[ki]])
        if (Q > 0) {
          Y = cumsum(rexp(Q + 1, 1))
          Tpsc[[ki]] = Y[1:Q] / Y[Q + 1]
        }
        else{
          Tpsc[[ki]] = numeric(0)
        }
        Mc[[ki]] = transitionmatrix(Xc[[ki]], Trposs, nph, bruit, tr$edge.length[ki])
        loiapplistec[[ki]] = loiapp(Xc[[ki]], L[L[, ki] > 0, ki], loiini, Trposs, bruit, tr$edge.length[ki], Tpsc[[ki]])
      }
      listc = pruninglkldlistini(tr,
                         Xc,
                         Mc,
                         L,
                         Dat,
                         nph,
                         length(loiapplistec),
                         root,
                         Trposs,
                         loiini,
                         bruit,
                         loiapplistec)
      #listc=pruninglkldlistini(tr,Xc,Mc,L,Dat,nph,nnodes,root,Trposs,loiini,bruit)
      if (ncogn >0){
        lkld1 = sum(sapply(root,function(xx){log(loiini %*% listc[[xx]][, 1:ncogn])-log(loiini %*%
                                                                                          Lkldinter[[xx]][, 1:ncogn])}))
      } else {
        lkld1=0
      }
      q = 0
      lkld2 = 0
    }
    if (is.nan(q)) {
      q = -Inf
    }
    r = lkld1 + lkld2 + q
    #print(exp(r))
    if (is.nan(r)) {
      r = -Inf
    }
    if (u < r)
    {
      #print(v)
      return(list(Xc, Mc, v, listc, loiapplistec, Tpsc))
    } else {
      #print(-v)
      return(list(X, M, -v, Lkldinter, loiappliste, Tps))
    }
  }
