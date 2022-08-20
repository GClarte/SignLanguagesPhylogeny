#gibbs uniquement sur les cognaits 1:ncogn
gibbsprior <- function(Iter, priortree, agemax, trini)
{
  #gc()
  #fait évoluer statini, données Da, mais n'utilise que Da[1:ncogn,]
  #renvoie une liste contenant une série d'états, Iter[1]=burn in, Iter[2] nombre de points, Iter[3] subsampling.
  tra = trini
  nnodes = max(tra$edge)
  feuilles=which(!(1:nnodes)%in%tr$edge[,1])
  m = Iter[1]
  mm = Iter[2]
  tr = tra
  res = list()
  for (i in 1:m)
  {
    #idem pour la partie sampling
    for (l in 1:mm) {
      quel = sample(1:2, 1)#on choisit quoi faire parmi les 7 possibilités
      if (quel == 1) {
        #change les temps dans l'arbre mais pas la topo
        VV = gibbstempsprior(
          tr,
          tr$root[[2]],
          priortree,
          agemax
        )
        tr = VV[[1]]
      } else if (quel == 2) {
        #change la topologie de l'arbre
        VV = gibbstopopart2prior(
          tr,
          tr$root,
          priortree,
          agemax,nnodes,feuilles
        )
        tr = VV[[1]]
      }
    }
    res[[i]] = tr
  }
  return(list(res))
}

