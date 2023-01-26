#simule des transformations
transf <- function(tr,Trposs,la,pr){
  k=length(tr$edge.length)
  n=length(Trposs)
  LL=list()
  for(i in 1:k){
    p=rpois(1,tr$edge.length[i]*la)
    y=sample(1:n,p,replace=T,prob = pr)
    LL[[i]]=y
  }
  return(LL)
}

#fait évoluer u par x
#length longueur de la branche sur laquelle on évolue pour scale le bruit
evolbranche <- function(x,u,Trposs,nph,bruit,L,loiini,length)
{
  nph=length(loiini)
  if (length(x)>0){
    Q=length(x)
    v=u
    Trans=transitionmatrix(x,Trposs,nph,bruit,length)
    U=sort(runif(length(x)))
    R=matrix(NA,ncol=Q+1,nrow=length(loiini))
    R[,1]=loiini
    for (i in 2:(Q+1)){
      R[,i]=R[,i-1]%*%Trposs[[x[i-1]]]
    }
    for (i in 1:length(v)){
      if(L[i]==0){
        v[i]=sample(1:nph,1,prob=Trans[v[i],])
      } else {
        taux=exp(-bruit*(1-L[i])*length)
        v[i]=sample(1:nph,1,prob=(1-taux)/nph*rep(1,length(loiini)) + taux*R[,(Q+1-sum(U>L[i]))])
      }
    }
    return(v)
  } else {
    Q=length(x)
    v=u
    for (i in 1:length(v)){
      if(L[i]!=0){
        taux=exp(-bruit*(1-L[i])*length)
        v[i]=sample(1:nph,1,prob=(1-taux)/nph*rep(1,length(loiini)) + taux*loiini)
      }
    }
    return(v)
  }
}

#A matrice à faire descendre,
evolbrancheloi <- function(AA,x,Tps,Trposs,nph,bruit,L,loiini,l)
{
  A=AA
  nph=length(loiini)
  qui=which(L!=0)
  if (length(x)>0){
    Q=length(x)
    Trans=transitionmatrix(x,Trposs,nph,bruit,l)
    A[,-qui]=t(Trans)%*%A[,-qui]
  }
  loiappliste=loiapp(x, L[L!=0], loiini, Trposs, bruit, l, Tps)
  A[,qui]=loiappliste
  A=A*(exp(-l*bruit))+(1-exp(-l*bruit))*matrix(rep(loiini,ncol(A)),ncol=ncol(A))
  return(A)
}

evolloich <- function(tr,X,Tps,L,nph,bruit,loiini,tip,ncogn,Trposs){
  AA=nearestcommonancestor(tr,tip)
  temp=AA
  ou=c()
  while(temp!=tr$root){
    mm=which(tr$edge[,2]==temp)
    ou=c(mm,ou)
    temp=tr$edge[mm,1]
  }
  A=matrix(loiini,ncol=ncogn,nrow=length(loiini))
  for (u in ou){
    A=evolbrancheloi(A,X[[u]],Tps[[u]],Trposs,nph,bruit,L[,u],loiini,tr$edge.length[u])
  }
  return(A)
}

evolloitout <- function (V,Param){
  loiini=Param$loiini
  X=V$X
  Tps=V$Tps
  bruit=V$noise
  ncogn=nrow(V$L)
  L=V$L
  taux=V$beta
  nph=sapply(loiini,length)
  Trposs = lapply(1:Param$nch, function(y) {
    lapply(Param$Trpossliste[[y]], function(x) {
      A = diag(rep(1, nph[y]))
      A[x[1], x[2]] = taux[y]
      A[x[1], x[1]] = 1 - taux[y]
      return(A)
    })
  })
  return(lapply(1:length(Param$ReconstructClade),function(tip){
    lapply(1:length(Param$loiini),function(ch){
      evolloich(V$Tr,X[[ch]],Tps[[ch]],L,nph[ch],bruit,loiini[[ch]],
                Param$ReconstructClade[[tip]],ncogn,Trposs[[ch]])*V$OtherReconstruct[[tip]][[ch]]
    })
  }))
}

#evolution à patir d'un noeud pour un son
evolarbreaux <- function(tr,X,j,u,N,Trposs,nph,bruit,L,loiini)
{
  M=N
  M[,j]=u
  if (length(which(tr$edge[,1]==j))<1)
  {
    return(M)
  } else {
    a=which(tr$edge[,1]==j)
      u1=evolbranche(X[[a[1]]],u,Trposs,nph,bruit,L[,a[1]],loiini,tr$edge.length[a[1]])
      u2=evolbranche(X[[a[2]]],u,Trposs,nph,bruit,L[,a[2]],loiini,tr$edge.length[a[2]])
      M=evolarbreaux(tr,X,tr$edge[a[1],2],u1,M,Trposs,nph,bruit,L,loiini)
      M=evolarbreaux(tr,X,tr$edge[a[2],2],u2,M,Trposs,nph,bruit,L,loiini)
      return(M)
  }
}

evolarbre <- function(tr,X,ini,k,root,nnodes,Trposs,nph,bruit,L)
{
  a=matrix(0,ncol=nnodes,nrow=k)
  u=sample(1:length(ini),k,prob=ini,replace=T)
  a=evolarbreaux(tr,X,root,u,a,Trposs,nph,bruit,L,ini)
  return(a)
}
