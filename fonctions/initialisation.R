
#Plus aucune de ces fonctions n'est utilisée


#reconstruit le chemin à partir de la matrice des voisins construite par matdist.
chemin <- function(V,i,j)
{
  if (V[i,j]==0)
  {
    return(j)
  } else {
    return(c(i,chemin(V,V[i,j],j)))
  }
}
#calcule l'état du noeud parentsachant l'état des enfants.
#i noeud de droite, j noeud de gauche
milieu <- function(i,j,S,n,AA)
{
  if (i==0 & j==0)
  {
    return(list(0,matrix(0,nrow=2,ncol=0),matrix(0,nrow=2,ncol=0)))
  } else if (i==0) {
    return(list(j,matrix(0,nrow=2,ncol=0),matrix(0,nrow=2,ncol=0)))
  } else if (j==0) {
    return(list(i,matrix(0,nrow=2,ncol=0),matrix(0,nrow=2,ncol=0)))
  } else {
    v=chemin(AA,i,j)
    k=ceiling(length(v)/2)
    A=list()
    A[[1]]=v[k]
    if (k>1)
    {
      A[[2]]=matrix(v[c(k:2,(k-1):1)],nrow=2,byrow = TRUE)
    } else {
      A[[2]]=matrix(0,nrow=2,ncol=0)
    }
    if (length(v)-k > 0) {
      A[[3]]=matrix(v[c(k:(length(v)-1),(k+1):length(v))],nrow=2,byrow=TRUE)
    } else {
      A[[3]]=matrix(0,nrow=2,ncol=0)
    }
    return(A)}
}

#fonction qui produit une liste de listes, chacune contenant les résultats de milieu.
#calcul ne dépednant que de S, faire avant tous les autres calculs.
paths <- function(S,k)
{
  AA=matdist(S,k)
  l=list()
  for (i in 1:k)
  {
    l[[i]]=list()
    for (j in 1:k)
    {
      l[[i]][[j]]=milieu(i,j,S,k,AA[[2]])
    }
  }
  return(l)
}

#fonction récursive
baseaux <- function(tr,d,adj,x,trf,path,M)
{
  u=adj[x,]
  X=trf
  N=M
  if (is.na(sum(u))){
    A1=baseaux(tr,d,adj,u[1],X,path,N)
    X=A1[[2]]
    N=A1[[3]]
    V=rep(0,length(A1[[1]]))
    V=A1[[1]]
    N[,x]=V
  }else  if (sum(u)==0)
  {
    N[,x]=d[,x]
    return(list(d[,x],trf,N))
  } else {
    A1=baseaux(tr,d,adj,u[1],X,path,N)
    X=A1[[2]]
    N=A1[[3]]
    A2=baseaux(tr,d,adj,u[3],X,path,N)
    X=A2[[2]]
    N=A2[[3]]
    V=rep(0,length(A1[[1]]))
    for (i in 1:nrow(d))
    {
      if (A1[[1]][i]==A2[[1]][i])
      {
        V[i]=A1[[1]][i]
      } else if (A1[[1]][i]==0){
        V[i]=A2[[1]][i]
      } else if (A2[[1]][i]==0) {
        V[i]=A1[[1]][i]
      } else {
        U=path[[A1[[1]][i]]][[A2[[1]][i]]]
        V[i]=U[[1]]
        X[[u[2]]]=fusion(X[[u[2]]],U[[2]])
        X[[u[4]]]=fusion(X[[u[4]]],U[[3]])
      }
    }
    N[,x]=V
  }
  return(list(V,X,N))
}

#renvoie un état d'initialisation correct
initXZ <- function(tr,d,adj,root,path,nnodes)
{
  l=list()
  k=nrow(tr$edge)
  M=matrix(NA,nrow=nrow(d),ncol=nnodes)
  for (i in 1:k)
  {
    l[[i]]=matrix(0,nrow=2,ncol=0)
  }
  AAAA=baseaux(tr,d,adj,root,l,path,M)
  return(list(AAAA[[2]],AAAA[[3]]))
}
