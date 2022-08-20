# calcule la vraisemblance par pruning

#initialisation
pruninglkldlistini <-
  function(tr,
           xt,
           Mt,
           L,
           Dat,
           nph,
           nnodes,
           root,
           Trposs,
           loiini,
           br,
           loiappliste) {
    #initialisation de la liste qui contient pour chaque noeud la vraisemblance des données SACHANT la valeur en ce noeud (ie une matrice de taille nbétats*nbcognats)
    Lkldinter = list()
    for (i in 1:nnodes) {
      Lkldinter[[i]] = list()
    }
    for (j in root){
      Lkldinter=pruninglkldrec(
        tr,
        xt,
        Mt,
        L,
        Dat,
        j,
        nph,
        Lkldinter,
        Trposs,
        loiini,
        br,
        loiappliste
      )
    }
    return(
      Lkldinter
    )
  }
# RR250420 parlons de ... dans les fonctions


pruninglkldrec <-
  function(tr,
           xt,
           Mt,
           L,
           Dat,
           i,
           nph,
           listini,
           Trposs,
           loiini,
           br,
           loiappliste) {
    #remplit la liste récursivement
    listt = listini
    M1b = matrix(0, ncol = nrow(Dat), nrow = nph)
    M2b = matrix(0, ncol = nrow(Dat), nrow = nph)
    II=which(tr$edge[, 1] == i)
    JJ=tr$edge[II,2]
    if (length(which(tr$edge[, 1] == JJ[1])) < 1) {
      #si feuille on remplit juste le bout avec des 0 et 1
      M = matrix(0, ncol = nrow(Dat), nrow = nph)
      for (jj in 1:nrow(Dat)){
        if (is.na(Dat[jj,JJ[1]])){
          M[,jj]=1
        } else {
          M[Dat[jj,JJ[1]],jj]=1 
        }
      }
      M1b = nvelleproba(L[, II[1]],
                       xt[[II[1]]],
                       loiini,
                       M,
                       Mt[[II[1]]],
                       Trposs,
                       br,
                       tr$edge.length[II[1]],
                       loiappliste[[II[1]]])
    } else {
      listt = pruninglkldrec(tr,
                             xt,
                             Mt,
                             L,
                             Dat,
                             JJ[1],
                             nph,
                             listt,
                             Trposs,
                             loiini,
                             br,
                             loiappliste)
      M1b = nvelleproba(L[, II[1]],
                       xt[[II[1]]],
                       loiini,
                       listt[[JJ[1]]],
                       Mt[[II[1]]],
                       Trposs,
                       br,
                       tr$edge.length[II[1]],
                       loiappliste[[II[1]]])
      
      
    }
    if (length(which(tr$edge[, 1] == JJ[2])) < 1) {
      M = matrix(0, ncol = nrow(Dat), nrow = nph)
      #si feuille on remplit juste le bout avec des 0 et 1
      for (jj in 1:nrow(Dat)){
        if (is.na(Dat[jj,JJ[2]])){
          M[,jj]=1
        } else {
          M[Dat[jj,JJ[2]],jj]=1 
        }
      }
      M2b = nvelleproba(L[, II[2]],
                        xt[[II[2]]],
                        loiini,
                        M,
                        Mt[[II[2]]],
                        Trposs,
                        br,
                        tr$edge.length[II[2]],
                        loiappliste[[II[2]]])
    } else {
      listt = pruninglkldrec(tr,
                             xt,
                             Mt,
                             L,
                             Dat,
                             JJ[2],
                             nph,
                             listt,
                             Trposs,
                             loiini,
                             br,
                             loiappliste)
      M2b = nvelleproba(L[, II[2]],
                        xt[[II[2]]],
                        loiini,
                        listt[[JJ[2]]],
                        Mt[[II[2]]],
                        Trposs,
                        br,
                        tr$edge.length[II[2]],
                        loiappliste[[II[2]]])
      
      
    }
      listt[[i]] = (M1b * M2b)
      return(listt)
    
  }


pruningeco <-
  function(tr,
           xt,
           Mt,
           Dat,
           nph,
           noeuds,
           lkldint,
           Trposs,
           L,
           br,
           loiini,
           loiappliste) {
    #recalcule la lkld des noeuds dans noeuds
    #la liste est dans l'ordre : du plus proche ancêtre commun à changer vers la racine,
    #jamais aucune feuille dans noeuds
    listt = lkldint
    for (i in 1:length(noeuds)) {
      j = noeuds[i]
      II = which(tr$edge[, 1] == j)
      JJ = tr$edge[II, 2]
      if (length(which(tr$edge[, 1] == JJ[1])) < 1) {
        M = matrix(0, ncol = nrow(Dat), nrow = nph)
        
        #si feuille on remplit juste le bout avec des 0 et 1
        for (jj in 1:nrow(Dat)){
          if (is.na(Dat[jj,JJ[1]])){
            M[,jj]=1
          } else {
            M[Dat[jj,JJ[1]],jj]=1 
          }
        }
        M1b = nvelleproba(L[, II[1]],
                          xt[[II[1]]],
                          loiini,
                          M,
                          Mt[[II[1]]],
                          Trposs,
                          br,
                          tr$edge.length[II[1]],
                          loiappliste[[II[1]]])
      } else {
        M1b = nvelleproba(L[, II[1]],
                          xt[[II[1]]],
                          loiini,
                          listt[[JJ[1]]],
                          Mt[[II[1]]],
                          Trposs,
                          br,
                          tr$edge.length[II[1]],
                          loiappliste[[II[1]]])
        
        
      }
      if (length(which(tr$edge[, 1] == JJ[2])) < 1) {
        M = matrix(0, ncol = nrow(Dat), nrow = nph)
        
        #si feuille on remplit juste le bout avec des 0 et 1
        for (jj in 1:nrow(Dat)){
          if (is.na(Dat[jj,JJ[2]])){
            M[,jj]=1
          } else {
            M[Dat[jj,JJ[2]],jj]=1 
          }
        }
        M2b = nvelleproba(L[, II[2]],
                          xt[[II[2]]],
                          loiini,
                          M,
                          Mt[[II[2]]],
                          Trposs,
                          br,
                          tr$edge.length[II[2]],
                          loiappliste[[II[2]]])
      } else {
        M2b = nvelleproba(L[, II[2]],
                          xt[[II[2]]],
                          loiini,
                          listt[[JJ[2]]],
                          Mt[[II[2]]],
                          Trposs,
                          br,
                          tr$edge.length[II[2]],
                          loiappliste[[II[2]]])
        
        
      }
      listt[[j]] = (M1b * M2b)
    }
    return(listt)
  }


cheminrac <- function(tr, i, root) {
  #calcule le chemin jusqu'à la racine en partant d'un noeud, à utiliser avec pruningeco
  if (i %in% root) {
    return(i)
  } else {
    u = tr$edge[which(tr$edge[, 2] == i), 1]
    return(c(i, cheminrac(tr, u, root)))
  }
}



loiapp <- function(x,l,loiini,Trposs,bruit,Length,Tps){
  #calcule la loi d'une apparition sur une branche
  # loi de l'état d'un cognat en bas d'une branche, conditionnellement
  # à ce qu'il est apparu en l sachant les transformations x
  if (length(l)==0){#si pas d'apparition renvoie matrice vide
    return(matrix(ncol=0,nrow=length(loiini)))
  } else {#sinon
    Q=length(x)#nombre de transf
    k=length(loiini)#nombre d'états
    if(Q>0){
      #U=sort(runif(Q))
      r=sapply(l,function(tt){sum(Tps>tt)})#quelles transf se sont appliquées à quelles app
      M=matrix(NA,ncol=Q+1,nrow=length(loiini))
      M[,1]=loiini
      for (i in 2:(max(r+1,2))){#calcule une fois pour toutes les lois d'apparitions selon les transf appliquées
        M[,i]=M[,i-1]%*%Trposs[[x[i-1]]]#loi de ces apparitions
      }
      #rajoute le bruit
      return(M[,r+1]*matrix(exp(-bruit*(1-l)*Length),ncol=length(r),nrow=length(loiini),byrow=T) +matrix(1-exp(-bruit*(1-l)*Length),ncol=length(r),nrow=length(loiini),byrow=T)/k)
    } else {
      M=matrix(loiini,ncol=length(l),nrow=length(loiini))
      return(M*matrix(exp(-bruit*(1-l)*Length),ncol=length(l),nrow=length(loiini),byrow=T) +matrix(1-exp(-bruit*(1-l)*Length),ncol=length(l),nrow=length(loiini),byrow=T)/k)
    }
  }
}

lkldcogn <- function(lkldinter,nch,root,quelcogn,loiini){
  if (quelcogn==0){
    return(sum(sapply(root[[1]],function(xx){sapply(1:nch, function(z) {
      sum(log(loiini[[z]] %*% lkldinter[[z]][[xx]][, 1:quelcogn]))
    })})))
  } else if (quelcogn==ncol(lkldinter[[1]][[root[[2]]]])) {
    return(sum(sapply(root[[2]],function(xx){sapply(1:nch, function(z) {
    sum(log(loiini[[z]] %*% lkldinter[[z]][[xx]][, 1:quelcogn]))
  })})))
  } else {
    return(sum(sapply(root[[1]],function(xx){sapply(1:nch, function(z) {
      sum(log(loiini[[z]] %*% lkldinter[[z]][[xx]][, (quelcogn+1):ncol(lkldinter[[1]][[1]])]))
    })})) + 
      sum(sapply(root[[2]],function(xx){sapply(1:nch, function(z) {
        sum(log(loiini[[z]] %*% lkldinter[[z]][[xx]][, 1:quelcogn]))
      })})))
  }
}

lkldcognpart <- function(lkldinter,quich,root,quelcogn,loiini){
  #pareil que lkldcogn mais pour un ch particulier
  if (quelcogn==0){
    return(sum(sapply(root[[1]],function(xx){log(loiini[[quich]] %*% lkldinter[[quich]][[xx]][, 1:quelcogn])})))
  } else if (quelcogn==ncol(lkldinter[[1]][[1]])) {
    return(sum(sapply(root[[2]],function(xx){log(loiini[[quich]] %*% lkldinter[[quich]][[xx]][, 1:quelcogn])})))
  } else {
    return(sum(sapply(root[[1]],function(xx){log(loiini[[quich]] %*% lkldinter[[quich]][[xx]][, (quelcogn+1):ncol(lkldinter[[1]][[1]])])})) +
             sum(sapply(root[[2]],function(xx){log(loiini[[quich]] %*% lkldinter[[quich]][[xx]][, 1:quelcogn])})))
  } 
}



loiappold <- function(x,l,loiini,Trposs,bruit,Length){return(loiini)}

#loiapplistetot=lapply(1:nch,function(x){lapply(1:length(tr$edge.length),function(y){loiapp(Xt[[x]][[y]],Lt[,y],loiini[[x]],Trposs[[x]],bruit,tr$edge.length[y])})})

nvelleproba <- function(Lbranche,x,loiini,Matprob,Mtb,Trposs,bruit,Length,loiappmat){
  #Lbarnche qui contient les apparitions pour la branche
  #x transf de la branche
  #Matprob, lkld qu'on remonte
  #Mtb transf de la branche sans apparition
  nph=nrow(Matprob)
  qui1=Lbranche==0
  Matfin=Matprob
  Matfin[,qui1]=Mtb%*%Matprob[,qui1]
  if (sum(!qui1)>0){
  Matfin[,!qui1]=matrix(rep(colSums(loiappmat*Matprob[,!qui1]),nph),byrow = T,nrow=nph)}
  return(Matfin)
}

