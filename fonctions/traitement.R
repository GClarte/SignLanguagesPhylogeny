tousarbres = lapply(VV2[[1]], function(x){x$Tr})


qui=list(tousarbres[[1]]$edge)
cb=list(c(VV2[[2]][3,1]))
for (i in 2:5000){
  u=sapply(qui,function(x){identical(x,tousarbres[[i]]$edge)})
  if (any(u)){
    cb[[which(u==T)]]=c(cb[[which(u==T)]],VV2[[2]][3,i])
  } else {
    qui=c(qui,list(tousarbres[[i]]$edge))
    cb=c(cb,VV2[[2]][3,i])
  }
}

#pds=sapply(VV3[[1]], function(x){x$weight})

#WW=sapply(VV3[[1]],function(statet){lkldcogn(statet$Lin,Param$nch,statet$Tr$root,nrow(Dat[[1]]),Param$loiini)})

ggsave("lambdas_simu.pdf",height=5,width=7,unit="cm")

h=ggplot(data=data.frame(value=sapply(VV2[[1]],function(x){x[[4]]*Param$agemax})),aes(x=value))+geom_density()+
  theme(axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                  axis.ticks.y=element_blank()) + labs(x=expression(rho))#+geom_vline(xintercept = 0.0001^*Param$agemax)

h
ggsave("rho_simu.pdf",height=5,width=7,unit="cm")

Prt=sapply(VV2[[1]],function(x){x[[3]][[1]]})
lgt=sapply(VV2[[1]],function(x){TT=x$Tr;a=which(TT$edge[,1]==TT$root);return(TT$edge.length[a])})

#calcul de la lkld marginale

sum(apply(VV2[[3]][,-1],2,function(z){log(sum(exp(z-max(z))))+(max(z))}))


#pour le plot du bruit


ggsave("bruit_simu.pdf",height=5,width=7,unit="cm")

  npart=length(VV2[[1]])
  betat=t(sapply(VV2[[1]],function(x){x[[7]]}))
  betat=data.frame(value=c(betat[,1],betat[,2],betat[,3],betat[,4]),character=c(rep("1",npart),rep("2",npart),rep("3",npart),rep("4",npart)))
  
  h = ggplot(data=betat) + geom_density(aes(x=value,colour=character)) + theme(axis.title.y=element_blank(),
                                                                             axis.text.y=element_blank(),
                                                                             axis.ticks.y=element_blank()) + #geom_vline(xintercept = c(0.001,0.0001))  +
    labs(x=expression(beta))
  h
  ggsave("beta_simu.pdf",height=5,width=7,unit="cm")


#fait un dernier resampling pour voir
npart=ncol(VV2[[2]])
pds=sapply(VV2[[1]],function(x){x[[10]]})
JJJ=1:npart
if (length(table(pds))!=1){
  JJJ=sample(1:npart,replace = T,prob = sapply(pds,length))
}

multiphylo=lapply(JJJ,function(x){VV2[[1]][[x]][[1]]})
class(multiphylo) <- "multiPhylo"

writeNexus(multiphylo,"europe6.nex")
multiphylo=read.nexus("europe6.nex")
writeNexus(consensus.edges(multiphylo),"europe6_consensus.nex")
plot(consensus.edges(multiphylo))

writeNexus(tr,"vraisimu.nex")

#calcule la vraisemblance en fin de chaine

lkldtraitement <- function(v){
  nch=Param$nch
  Trposs = lapply(1:nch, function(y) {
    lapply(Param$Trpossliste[[y]], function(x) {
      A = diag(rep(1, Param$nph[y]))
      A[x[1], x[2]] = v[[7]][y]
      A[x[1], x[1]] = 1 - v[[7]][y]
      return(A)
    })
  })
  L = v[[9]]
  L[v[[9]] != 0] = rbeta(length(which(v[[9]] != 0)), L[v[[9]] != 0], 1)
  loiapplistetot = lapply(1:nch, function(x) {
    lapply(1:length(v[[1]]$edge.length), function(y) {
      loiapp(v[[5]][[x]][[y]], L[L[, y] > 0, y], loiini[[x]], Trposs[[x]], v[[8]], v[[1]]$edge.length[y], v[[6]][[x]][[y]])
    })
  })
  Mt = lapply(1:nch, function(z) {
    lapply(1:length(v[[1]]$edge.length), function(y) {
      transitionmatrix(v[[5]][[z]][[y]], Trposs[[z]], Param$nph[z], v[[8]], v[[1]]$edge.length[y])
    })
  })
  loiapplistetot = lapply(1:nch, function(x) {
    lapply(1:length(v[[1]]$edge.length), function(y) {
      loiapp(v[[5]][[x]][[y]], L[L[, y] > 0, y], Param$loiini[[x]], Trposs[[x]], v[[8]], v[[1]]$edge.length[y], v[[6]][[x]][[y]])
    })
  })
  lkldinter = lapply(1:nch, function(x) {
    pruninglkldlistini(v[[1]],
                       v[[5]][[x]],
                       Mt[[x]],
                       L,
                       Dat[[x]],
                       Param$nph[x],
                       1+length(Mt[[1]]),
                       v[[1]]$root,
                       Trposs[[x]],
                       Param$loiini[[x]],
                       bruit,
                       loiapplistetot[[x]])
  })
  tr=v[[1]]
  tr$root=list(c(1,2),v[[1]]$root)
  return(lkldtopotout(tr,v[[5]],v[[2]],v[[9]],v[[4]],L))
}

Lfin=sapply(VV2[[1]],lkldtraitement)

#pour plotter les probas des transformations

library(ggplot2)

qq=numeric()
qq2=numeric()
ord=order(rowMeans(Prt))
for (i in ord){
  qq=c(qq,Prt[i,])
  qq2=c(qq2,rep(paste(as.character(passages[[1]][[i]][1]),as.character(passages[[1]][[i]][2])),ncol(Prt)))
}

dataPrt <-data.frame(Freq <- qq,Transition <- qq2)

p<- ggplot(dataPrt,aes(x=Transition,y=Freq))+geom_violin()+ coord_flip()  + theme(axis.title.x=element_blank())
p
ggsave("transformations_ch1.pdf", height=10,width=7,unit="cm")


#pour regarder précisément les transfos ? utile ?
TRUC=list()
for (i in 1:5000){
  aaa=which(VV2[[1]][[i]]$Tr$edge[,1]==VV2[[1]][[i]]$Tr$root & VV2[[1]][[i]]$Tr$edge[,2]==5)
  if (length(aaa)==1){
    TRUC=c(TRUC,VV2[[1]][[i]]$X[[1]][[aaa]])
  }
}

#extrait les topologies juste avant la dernière fusion
topoav <- lapply(5002:10001,function(x){VV2[[4]][[x]]})


#cherche les feuilles à partir d'une racine
quileaf <- function(tr){
  lapply(tr$root,function(y){quidescend(tr,y)})
}

quidescend <- function(tr,i){
  u=which(tr$edge[,1]==i)
  if (length(u)==0){
    return(i)
  } else {
    return(c(quidescend(tr,tr$edge[u[1],2]),quidescend(tr,tr$edge[u[2],2])))
  }
}

UU=lapply(topoav,quileaf)
cbb=numeric()
for (i in 1:length(topoav)){
  AAA=quileaf(topoav[[i]])
  if (any(sapply(AAA,function(x){1 %in% x && 4 %in% x}))){
    cbb=c(cbb,VV2[[3]][3,i])
  }
}

topoav <- lapply(2:5001,function(x){VV2[[4]][[x]]})
cbb=numeric()
for (i in 1:length(topoav)){
  AAA=quileaf(topoav[[i]])
  if (any(sapply(AAA,function(x){(1 %in% x && 2 %in%x)}))){
    cbb=c(cbb,VV2[[3]][3,i])
  }
}


#pour les cognacy
#quelmeanings=1:100
trouvercheminaux <- function(tr,a,b,chemin){
  #trouve le chemin dans tr de a à b
  if (a==b){
    return(chemin)
  } else if (a == tr$root){
    w=which(tr$edge[,2]==b)
    return(trouvercheminaux(tr,a,tr$edge[w,1],c(chemin,w)))
  } else if(b==tr$root){
    w=which(tr$edge[,2]==a)
    return(trouvercheminaux(tr,tr$edge[w,1],b,c(chemin,w)))
  } else {
    w1=which(tr$edge[,2]==a)
    w2=which(tr$edge[,2]==b)
    return(trouvercheminaux(tr,tr$edge[w1,1],tr$edge[w2,1],c(chemin,w1,w2)))
  }
}

quelmeanings=1:100
tr$root=16

UU=sapply(1:2000, function(x){rowSums(VV2[[1]][[x]][[9]][,trouvercheminaux(VV2[[1]][[x]][[1]],1,5,numeric())])==0})
#vrai=rowSums(Lreelnb[,trouvercheminaux(tr,5,1,numeric())])==0
datacogn <- data.frame(Posterior=rowMeans(UU),Nb=1:nrow(Dat[[1]]),Meaning=quelmeanings)
datacogn$Meaning=factor(quelmeanings,levels=as.character(quelmeanings))
h=ggplot(data=datacogn,aes(x=Nb,y=Posterior))+geom_point() +geom_text(aes(label=Meaning),hjust=0, vjust=0)#
h=ggplot(data=datacogn,aes(x=Nb,y=Posterior,color=vrai))+geom_point() +geom_text(aes(label=Meaning),hjust=0, vjust=0)
ggsave("cognats_simu_5_1.pdf",width=17,height=10,unit="cm")

UU=sapply(1:2000, function(x){rowSums(VV2[[1]][[x]][[9]][,trouvercheminaux(VV2[[1]][[x]][[1]],10,7,numeric())])==0})
#vrai=rowSums(Lreelnb[,trouvercheminaux(tr,3,6,numeric())])==0
datacogn <- data.frame(Posterior=rowMeans(UU),Nb=1:nrow(Dat[[1]]),Meaning=quelmeanings)
datacogn$Meaning=factor(quelmeanings,levels=as.character(quelmeanings))
h=ggplot(data=datacogn,aes(x=Nb,y=Posterior))+geom_point() +geom_text(aes(label=Meaning),hjust=0, vjust=0)
#h=ggplot(data=datacogn,aes(x=Nb,y=Posterior,color=vrai))+geom_point() +geom_text(aes(label=Meaning),hjust=0, vjust=0)
ggsave("cognats_Pol_RSL.pdf",width=17,height=10,unit="cm")

UU=sapply(1:2000, function(x){rowSums(VV2[[1]][[x]][[9]][,trouvercheminaux(VV2[[1]][[x]][[1]],6,4,numeric())])==0})
#vrai=rowSums(Lreelnb[,trouvercheminaux(tr,3,6,numeric())])==0
datacogn <- data.frame(Posterior=rowMeans(UU),Nb=1:nrow(Dat[[1]]),Meaning=quelmeanings)
datacogn$Meaning=factor(quelmeanings,levels=as.character(quelmeanings))
h=ggplot(data=datacogn,aes(x=Nb,y=Posterior))+geom_point() +geom_text(aes(label=Meaning),hjust=0, vjust=0)
#h=ggplot(data=datacogn,aes(x=Nb,y=Posterior,color=vrai))+geom_point() +geom_text(aes(label=Meaning),hjust=0, vjust=0)
ggsave("cognats_LSF_PortSL.pdf",width=17,height=10,unit="cm")

UU=sapply(1:2000, function(x){rowSums(VV2[[1]][[x]][[9]][,trouvercheminaux(VV2[[1]][[x]][[1]],4,10,numeric())])==0})
#vrai=rowSums(Lreelnb[,trouvercheminaux(tr,3,6,numeric())])==0
datacogn <- data.frame(Posterior=rowMeans(UU),Nb=1:nrow(Dat[[1]]),Meaning=quelmeanings)
datacogn$Meaning=factor(quelmeanings,levels=as.character(quelmeanings))
h=ggplot(data=datacogn,aes(x=Nb,y=Posterior))+geom_point() +geom_text(aes(label=Meaning),hjust=0, vjust=0)
#h=ggplot(data=datacogn,aes(x=Nb,y=Posterior,color=vrai))+geom_point() +geom_text(aes(label=Meaning),hjust=0, vjust=0)
ggsave("cognats_NZSL_BSL.pdf",width=17,height=10,unit="cm")


Resla=sapply(1:2000,function(x){VV2[[1]][[x]][[2]]})
df <- data.frame(Lambda=c(Resla[1,],Resla[2,],Resla[3,],Resla[4,]),Trait=c(rep("1",2000),rep("2",2000),rep("3",2000),rep("4",2000)))
h1 <- ggplot()+geom_density(data=df,aes(x=Lambda,color=Trait))+geom_vline(xintercept = c(0.001,0.0001))#+facet_grid(.~Trait)
h1

ggsave("lambdavrai.pdf",height=5,width=7,unit="cm")

h2=ggplot()+geom_density()


tousPr1=sapply(1:2000,function(x){VV2[[1]][[x]][[3]][[1]]})
tousPr2=sapply(1:2000,function(x){VV2[[1]][[x]][[3]][[2]]})
tousPr3=sapply(1:2000,function(x){VV2[[1]][[x]][[3]][[3]]})


#généalogie
npart=ncol(VV2[[2]])
x1=numeric()
y1=numeric()
xend=numeric()
yend=numeric()
color=numeric()
actif=1:npart
for (i in ncol(VV2[[3]]):2){
  x1=c(x1,rep(i,npart))
  y1=c(y1,1:npart)
  xend=c(xend,rep(i-1,npart))
  yend=c(yend,VV2[[3]][,i])
  colortemp=rep(20/npart,npart)
  colortemp[actif]=1
  color=c(color,colortemp)
  actif=unique(VV2[[3]][actif,i])
}
geneal=data.frame(x1=x1,xend=xend,y1=y1,yend=yend,color=color)

b <- ggplot(geneal) + geom_segment(aes(x = x1, y = y1, xend = xend, yend = yend),alpha=color, data = geneal)+xlab("generation")+ylab("particle")
#b = b  + geom_vline(data=data.frame(pp=cumsum(sapply(Param$batches,length))),aes(xintercept=pp,color="red"))
b
    ggsave("geneal_simu.pdf",height=5,width=10, unit="cm")

#version juste généalogie
x1=numeric()
y1=numeric()
xend=numeric()
yend=numeric()
color=numeric()
actif=1:10000
for (i in ncol(VV2[[5]]):2){
  x1=c(x1,rep(i,length(actif)))
  y1=c(y1,actif)
  xend=c(xend,rep(i-1,length(actif)))
  yend=c(yend,VV2[[5]][actif,i])
  actif=unique(VV2[[5]][actif,i])
}

geneal=data.frame(x1=x1,xend=xend,y1=y1,yend=yend)

b <- ggplot(geneal) + geom_segment(aes(x = x1, y = y1, xend = xend, yend = yend), data = geneal)

#pour traitement des trucs en parallèle
ZZZ=lapply(VV2[[1]][[1]],function(x){x$Tr})
for (i in 2:length(VV2)){
 ZZZ=c(ZZZ,lapply(VV2[[i]][[1]],function(x){x$Tr})) 
}


ZZZ=lapply(ZZZ,function(x){y=x;y$root=y$root[[2]];return(y)})
class(ZZZ) <- "multiPhylo"

writeNexus(ZZZ,"MH_full_39_50_20k_twohands_handpart_poa2_priyule_30oct.nex")


#pour l'apparition des clades
sapply(VV2[[1]][[1]],function(x){is.monophyletic.perso(x$Tr,c(1,2))})


#pour le beta

U=sapply(VV2[[1]],function(x){x[[7]]})
  

#pour le traitement des trucs fin parallèle

multiphylo=lapply(UU[[1]],function(x){x$Tr})
class(multiphylo) <- "multiPhylo"

writeNexus(multiphylo,"bruit1seul_full_priBNZ_3kpts_1kpas_handshapesimple_handpart_PoA2_Mov2_priyule_Sat_Dec_26_21:59:35_2020.nex")
