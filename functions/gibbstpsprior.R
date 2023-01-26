

gibbstempsprior <-  function(tr,
           root,
           priortree,
           agemax)
  {
    #modifie l'âge des branches
    #recommence entièrement le calcul de lkld.
    trc = modif(tr,agemax)
    u=log(runif(1))
    if (u < priortree(trc,agemax) - priortree(tr,agemax))
    {
      return(list(trc, 5))
    } else {
      return(list(tr, -5))
    }
  }



