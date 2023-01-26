is.monophyletic.perso <- function(tr,tip){
  k=lapply(tip,function(x){ancetres(tr,x)})
  rootmax=sapply(k,function(x){x[1]})
  W=unique(rootmax)
  if (length(W)==1){
    ll=sapply(k,length)
    i=0
    while ((i)<min(ll) & all(sapply(k,function(x){x[i+1]})==k[[1]][i+1])){
      i=i+1
    }
    rr=k[[1]][i]
    dscd=desc(tr,rr)
  } else {
    dscd=c()
    for (j in 1:length(W)){
      dscd=c(dscd,desc(tr,W[j]))
    }
  }
  return(length(dscd)==length(tip))
}

ancetres <- function(tr,x){
  i=x
  u=c(x)
  if (is.list(tr$root)){
    rroott=tr$root[[2]]
  } else {
    rroott=tr$root
  }
  while (!i%in%rroott){
    i=tr$edge[which(tr$edge[,2]==i),1]
    u=c(i,u)
  }
  return(u)
}


ageconstraints <- function(tr,clade){
  AA=nearestcommonancestor(tr,clade[[1]])
  temp=AA
  l=0
  mm=which(tr$edge[,1]==temp)
  while(length(mm)>0){
    temp=tr$edge[mm[1],2]
    l=l+tr$edge.length[mm[1]]
    mm=which(tr$edge[,1]==temp)
  }
  return(l>clade[[2]][1] & l<clade[[2]][2])
}

desc <- function(tr,i){
  w=which(tr$edge[,1]==i)
  if (length(w)>0){
    return(c(desc(tr,tr$edge[w[1],2]),desc(tr,tr$edge[w[2],2])))
  } else {
    return(i)
  }
}
