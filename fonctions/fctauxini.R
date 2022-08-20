
transitionmatrix <- function(m,Trposs,k,bruit,length)
{
  #calcule la matrice des transformations sur une branche, avec un bruit
  n=length(m)
  M=diag(rep(1,k))
  if (n>0){
    for (i in 1:n)
    {
      M=M%*%Trposs[[m[i]]]
    }
  }
  taux=exp(-bruit*length)
  Z=matrix((1-taux)/k,ncol=k,nrow=k)
  M=M*taux + Z
  return(M)
}

vectadj <-function(tr)
{
  #calcule la "matrice d'adjacence", ie matrice de 4 colonnes, res[i,]=[fils droit, branche droite, fils gauche, branche gauche]
  #non utilisée
  nnode=max(tr$edge)
  res=matrix(0,ncol=4,nrow=nnode)
  for (i in 1:nnode)
  {
    u=which(tr$edge[,1]==i)
    if (length(u)>0){
      res[i,]=c(tr$edge[u[1],2],u[1],tr$edge[u[2],2],u[2])
    }
  }
  return(res)
}
