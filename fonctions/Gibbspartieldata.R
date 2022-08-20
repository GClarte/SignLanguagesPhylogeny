
gibbsbr2 <- function(Da, ncogn, Iter, Prior, Param, statini,fin)
{
  #updates statini using data Da, the meanings 1:ncogn, takes Iter[1] points on the chain with Iter[2] steps between each
  #Prior and Param are the same in SMCbruit
  #statini is the starting state
  #fin : if the noise is fixed
  #returns a list, the first element is the usefull information for all points, the second is the number of success of the proposals
  #the last is the full state at the end of the chain, this is the one that must be used for more steps, the last one counts how many transformations there were on the edges that has been changed by the topology move used for debug
  tra = statini$Tr #we start by computing and storing all the values needed from statini, Param and Prior
  nph = Param$nph
  nch = Param$nch
  Trpossliste = Param$Trpossliste
  if (fin){
    Prob = Param$Probfin
  } else {
    Prob = Param$Prob
  }
  Pri = Prior$Prila
  hyperbruit = Prior$pribruit
  hyperbeta = Prior$pribeta
  Prirho = Prior$Prirho
  tipprior=Prior$tipprior
  nnodes = max(tra$edge)
  agemax=Param$agemax
  commenttopo=matrix(0,nrow=2,ncol=6)
  hyperpbini = Prior$hyperpbini
  loiini = Param$loiini
  m = Iter[1]
  mm = Iter[2]
  bruit = statini$bruit
  taux = statini$taux
  tr = tra
  pr = statini$P#probabilities of those transformations
  res = list()
  ZZZut=matrix(nrow=5,ncol=0)
  Xt = statini$X#transformations
  quoisur=matrix(nrow=nch,ncol=0)
  Trposs = lapply(1:nch, function(y) {
    lapply(Trpossliste[[y]], function(x) {
      A = diag(rep(1, nph[y]))
      A[x[1], x[2]] = taux[y]
      A[x[1], x[1]] = 1 - taux[y]
      return(A)
    })
  })
  Mt = lapply(1:nch, function(z) { #list containing the transition matricies of each edge
    lapply(1:length(tr$edge.length), function(y) {
      transitionmatrix(Xt[[z]][[y]], Trposs[[z]], nph[z], bruit, tr$edge.length[y])
    })
  })
  la = statini$La
  rho = statini$Rho#paramètre d'apparition des mots
  L = statini$L
  NL = statini$NL
  priortree = Prior$priortree
  Tps = statini$Tps
  qquel=numeric()
  comment=matrix(0,ncol=7,nrow=2)
  loiapplistetot = lapply(1:nch, function(x) {#law of the apparitions
    lapply(1:length(tra$edge.length), function(y) {
      loiapp(Xt[[x]][[y]], L[L[, y] > 0, y], loiini[[x]], Trposs[[x]], bruit, tra$edge.length[y], Tps[[x]][[y]])
    })
  })
  lkldinter = lapply(1:nch, function(x) {#the matrix recursively defined on the tree that allows to compute the likelihood
    pruninglkldlistini(tr,
                       Xt[[x]],
                       Mt[[x]],
                       L,
                       Da[[x]],
                       nph[x],
                       nnodes,
                       tr$root[[2]],
                       Trposs[[x]],
                       loiini[[x]],
                       bruit,
                       loiapplistetot[[x]])
  })
  for (i in 1:m)
  {
    #idem pour la partie sampling
    for (l in 1:mm) {
      quel = sample(1:7, 1, prob = Prob)
      qquel=c(qquel,quel)
      if (quel == 1) {
        for (j in 1:nch) {
          #change the transformations on the edges for each character
          VV = gibbstransf(
            tr,
            Da[[j]],
            la[j],
            Xt[[j]],
            Mt[[j]],
            Trposs[[j]],
            nph[j],
            tr$root[[2]],
            loiini[[j]],
            bruit,
            lkldinter[[j]],
            L,
            pr[[j]],
            ncogn,
            loiapplistetot[[j]],
            Tps[[j]]
          )
          Xt[[j]] = VV[[1]]
          Mt[[j]] = VV[[2]]
          comment[1+(VV[[3]]<0),quel]=comment[1+(VV[[3]]<0),quel]+1
          lkldinter[[j]] = VV[[4]]
          Tps[[j]] = VV[[6]]
          loiapplistetot[[j]] = VV[[5]]
          #qquel=c(qquel,quel)
          #cat(l,quel,"\n",file="zut.txt",sep="\t",append=T)
        }
      } else if (quel == 2) {
        #fparameters
        for (j in 1:nch) {
          VV = gibbslambda(Xt[[j]], tr, Pri[[j]][1], Pri[[j]][2],Pri[[j]][c(3,4)])
          VVV = gibbsprobas(Xt[[j]], hyperpbini[[j]], length(Trposs[[j]]))
          la[j] = VV[[1]]
          pr[[j]] = VVV
          comment[1,quel]=comment[1,quel]+1
        }
        VV = tryCatch(#some try catch exists because of possible bugs when choosing the prior
          { gibbsrho(rho, NL, tr, Prirho)           
            },
          error=rho)
        rho = VV
        #qquel=c(qquel,quel)
      } else if (quel == 3) {
        #change the length of the edges
        VV = gibbstemps(
          tr,
          Da,
          la,
          Xt,
          Mt,
          Trposs,
          nph,
          tr$root[[2]],
          loiini,
          bruit,
          lkldinter,
          L,
          pr,
          ncogn,
          NL,
          rho,
          priortree,
          loiapplistetot,
          Tps,
          agemax
        )
        tr = VV[[1]]
        lkldinter = VV[[2]]
        loiapplistetot = VV[[4]]
        comment[1+(VV[[3]]<0),quel]=comment[1+(VV[[3]]<0),quel]+1
        #qquel=c(qquel,quel)
      } else if (quel == 4) {
        #change the apparitions
        VV = tryCatch(
          {gibbsL(
          tr,
          Mt,
          Xt,
          Da,
          Trposs,
          nph,
          nnodes,
          nch,
          tr$root[[2]],
          loiini,
          bruit,
          lkldinter,
          L,
          NL,
          rho,
          ncogn,
          loiapplistetot,
          Tps
        )},
        error=list(L,lkldinter,NL,loiapplistetot,-1))
        L = VV[[1]]
        lkldinter = VV[[2]]
        NL = VV[[3]]
        loiapplistetot = VV[[4]]
        comment[1+(VV[[5]]<0),quel]=comment[1+(VV[[5]]<0),quel]+1
        #qquel=c(qquel,quel)
      } else if (quel == 5) {
        #change in the topology
        # VV = tryCatch(
        #   {
        #     gibbstopopart2(
        #       tr,
        #       Da,
        #       la,
        #       Xt,
        #       Mt,
        #       Trposs,
        #       nph,
        #       tr$root,
        #       loiini,
        #       bruit,
        #       lkldinter,
        #       L,
        #       NL,
        #       rho,
        #       ncogn,
        #       priortree,
        #       loiapplistetot,
        #       Tps,
        #       agemax,
        #       pr,
        #       tipprior
        #     )
        #   },
        #   error=list(tr,lkldinter,0,-6,loiapplistetot,Xt,NL,L,Tps,Mt,rep(-3,nch)))
        VV = gibbstopopart2(
              tr,
              Da,
              la,
              Xt,
              Mt,
              Trposs,
              nph,
              tr$root,
              loiini,
              bruit,
              lkldinter,
              L,
              NL,
              rho,
              ncogn,
              priortree,
              loiapplistetot,
              Tps,
              agemax,
              pr,
              tipprior,
              Prior,
              fin
            )
        lkldinter = VV[[2]]
        tr = VV[[1]]
        loiapplistetot = VV[[5]]
        comment[1+(VV[[4]]<0),quel]=comment[1+(VV[[4]]<0),quel]+1
          commenttopo[1+(VV[[4]]<0),abs(VV[[4]])]=commenttopo[1+(VV[[4]]<0),abs(VV[[4]])]+1
        Xt=VV[[6]]
        NL=VV[[7]]
        L=VV[[8]]
        Tps=VV[[9]]
        Mt=VV[[10]]
        la=VV[[12]]
        rho=VV[[13]]
        bruit=VV[[14]]
        quoisur=cbind(quoisur,VV[[11]])
        ZZZut=cbind(ZZZut,VV[[3]])
        #qquel=c(qquel,quel)
      } else if (quel == 6 & !is.matrix(Prior$pribruit)) {
        #change the noise
        VV = tryCatch(
          { gibbsbruit(
          tr,
          Da,
          la,
          Xt,
          Mt,
          Trposs,
          nph,
          tr$root,
          loiini,
          bruit,
          lkldinter,
          L,
          NL,
          rho,
          ncogn,
          taux,
          nch,
          hyperbruit,
          loiapplistetot,
          Tps
        )
          },
        error=list(bruit,lkldinter,loiapplistetot,Mt,-7))
        bruit = VV[[1]]
        lkldinter = VV[[2]]
        loiapplistetot = VV[[3]]
        Mt=VV[[4]]
        comment[1+(VV[[5]]<0),quel]=comment[1+(VV[[5]]<0),quel]+1
      } else if (quel==6 & is.matrix(Prior$pribruit)){
        #not used anymore
        VV = gibbsbruitfauxpost(
          tr,
          Da,
          la,
          Xt,
          Mt,
          Trposs,
          nph,
          tr$root,
          loiini,
          bruit,
          lkldinter,
          L,
          NL,
          rho,
          ncogn,
          taux,
          nch,
          hyperbruit,
          loiapplistetot,
          Tps
        )
        bruit = VV[[1]]
        lkldinter = VV[[2]]
        loiapplistetot = VV[[3]]
        Mt=VV[[4]]
        #qquel=c(qquel,quel)
      } else if (quel == 7) {
        #change beta
        VV = gibbstaux(
          tr,
          Da,
          la,
          Xt,
          Mt,
          Trposs,
          Trpossliste,
          nph,
          tr$root,
          loiini,
          bruit,
          lkldinter,
          L,
          NL,
          rho,
          ncogn,
          taux,
          nch,
          hyperbeta,
          loiapplistetot,
          Tps
        )
        taux = VV[[1]]
        Trposs = VV[[2]]
        lkldinter = VV[[3]]
        loiapplistetot = VV[[4]]
        Mt=VV[[5]]
        #qquel=c(qquel,quel)
      }
    }
    #we recompute everything from time to time
    lkldinter = lapply(1:nch, function(x) {
      pruninglkldlistini(tr,
                         Xt[[x]],
                         Mt[[x]],
                         L,
                         Da[[x]],
                         nph[x],
                         nnodes,
                         tr$root[[2]],
                         Trposs[[x]],
                         loiini[[x]],
                         bruit,
                         loiapplistetot[[x]])
    })
    #we do not return everything for every steps
    res[[i]] = list(
      X = Xt,
      L = L,
      La = la,
      Tr = tr,
      NL = NL,
      P = pr,
      Rho = rho,
      taux = taux,
      bruit = bruit,
      Tps = Tps,
      weight=statini$weight
    )
  }
  staterenvoi = list(
    X = Xt,
    L = L,
    La = la,
    Lin=lkldinter,
    Tr = tr,
    Ltp = lkldtopotout(tr, Xt, la, NL, rho, L,pr),
    P = pr,
    Rho = rho,
    NL = NL,
    taux = taux,
    bruit = bruit,
    Tps = Tps,
    weight = statini$weight
  )
  return(list(res,comment,staterenvoi,commenttopo,ZZZut))
}
