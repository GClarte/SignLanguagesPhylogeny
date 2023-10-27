# Execute this R script to reproduce the analysis on the Asian Sign Languages.


J = read.csv("datasets/dataset.csv",
             header = T,
             stringsAsFactors = T)

source("functions/fctauxini.R")
source("functions/gibbsparam.R")
source("functions/gibbstps.R")
source("functions/gibbstranfter.R")
source("functions/initialisation.R")
source("functions/lkldpruningbis.R")
source("functions/simudata.R")
source("functions/Gibbspartieldata.R")
source("functions/gibbsrho.R")
source("functions/modiftopo.R")
source("functions/SMCforet.R")
source("functions/SMCforet_aux.R")
source("functions/testclades.R")

library("parallel")

# select characters to be included
qui = c(28, 29, 6, 10)
#qui=c(7,28,20,10)
#qui=c(7,28,20,6)

# select languages to be included
quelleslangues = c("Hongkong", "Chinese", "Japanese", "Taiwan", "AllMissing")

nlangues = length(quelleslangues)

rdirichlet <- function (n, alpha){
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x / as.vector(sm))
}

k = c()
for (i in 1:100) {
  if (sum(as.integer(J[which(J[, 4] %in% quelleslangues), 3]) == i) == nlangues &
      !any(is.na(sapply(1:3, function(z) {
        sum(as.integer(J[as.integer(J[, 1]) == i &
                         J[, 2] %in% quelleslangues, qui[z]]))
      })))) {
    k = c(k, i)
  }
}

# list of meanings to include
quelmeanings = 1:100

nch = length(qui)

Dtot = list()
for (i in 1:nch) {
  Dtot[[i]] = matrix(ncol = nlangues, nrow = length(quelmeanings))
  for (j in 1:length(quelleslangues)) {
    for (l in 1:length(quelmeanings)) {
      if (length(J[which(J[, 4] == quelleslangues[j] &
                         as.integer(J[, 3]) == quelmeanings[l]), qui[i]]) > 0) {
        Dtot[[i]][l, j] = as.integer(J[which(J[, 4] == quelleslangues[j] &
                                               as.integer(J[, 3]) == quelmeanings[l]), qui[i]])
      } else {
        Dtot[[i]][l, j] = NA
      }
    }
  }
}


toutestransf = function(n) {
  A = matrix(nrow = 2, ncol = 0)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        A = cbind(A, c(i, j))
      }
    }
  }
  return(A)
}

partopo = .001
priorcoal <- function(tr, agemax) {
  tip = Prior$tipprior
  if (!all(sapply(tipprior, 
                  function(tip) {is.monophyletic.perso(statet$Tr, tip)}
  ))) {
    return(-Inf)
  } else {
    L = sum(tr$edge.length)
    n = length(tr$tip.label)
    return(log((1 - exp(-partopo * L / 2)) ^ (n - 1) * (exp(partopo * t / 2))))
  }
}

etats = sapply(1:nch, function(x) {length(levels(J[, qui[x]]))} )

passages = lapply(etats, function(x) {
  S = toutestransf(x)
  return(split(S, col(S)))
})

priorbeta <- function(x, y) {
  if (length(x$root[[2]]) > 0) {
    0
  } else{
    dbeta((sum(x$edge.length) - 2 * y) / ((length(x$tip.label) - 2) * y), 1, 100, log = T)
  }
}

prioryule <- function(x, y) {
  tips = Prior$tipprior
  if (!all(sapply(tips, function(tip) {
    is.monophyletic.perso(x, tip)
  }))) {
    return(-Inf)
  } else {
    if (length(x$root[[2]] > 0)) {
      K = sapply(x$root[[2]], function(z) {
        y - hauteur(x, z)
      })
      return(-1 * length(x$tip.label) * log(sum(x$edge.length) + sum(K)))
    } else {
      return(-1 * length(x$tip.label) * log(sum(x$edge.length)))
    }
  }
}

priorunif <- function(tr, agemax) {
  tips = Prior$tipprior
  cladeage = Param$Cladeage
  if (!all(sapply(tips, function(tip) {
    is.monophyletic.perso(tr, tip)
  })) ||
  !all(sapply(cladeage, function(clade) {
    ageconstraints(tr, clade)
  }))) {
    return(-Inf)
  } else {
    hhh = hauteur(tr, tr$root[[2]])
    if (!(hhh > 150 & hhh < 1000)) {
      return(-Inf)
    } else {
      return(-(length(tr$tip.label) - 1) * log(hauteur(tr, tr$root[[2]])))
    }
  }
}

#paramètres

ncogntot = nrow(Dtot[[1]])

nch = length(etats)
Dat = Dtot
hyperpbini = lapply(passages, function(x) {
  rep(1, length(x))
})
#hyperpbini=list(2,10,2,20)

Pri = list(c(1, 500, 0, 1), c(1, 500, 0, 1), c(1, 500, 0, 1), 
           c(1, 500, 0, 1))

Prirho = c(1, 500)
loiini = lapply(etats, function(x) {rep(1 / x, x)})

pribruit = list(
  list(function() { runif(1, 10 ^ -5, 10 ^ -1)}, 
       function(x) {
         if (x > 10 ^ -5 & x < 10 ^ -1) {
           return(-x)
         } else {
           return(-Inf)
         }
       }),
  list(function(x) { rnorm(1, x, .001)},
       function(x, y) {dnorm(x, y, .001, log = T)}
  )
)


convcontraintes <- function(ages) {
  res = list()
  cpt = 1
  for (i in 1:(nrow(ages) - 1)) {
    for (j in (i + 1):nrow(ages)) {
      res[[cpt]] = list(c(ages[i, 1], ages[j, 1]), c(min(ages[i, 2], ages[j, 2]) -
                                                       10, Inf))
      cpt = cpt + 1
    }
  }
  return(res)
}

# List constraints on ages
Contraintesages = list(list(c(3, 4), c(60, 85)))

Prior = list(
  Prirho = Prirho,
  hyperpbini = hyperpbini,
  Prila = Pri,
  #priortree=prioryule,
  priortree = priorunif,
  #pribeta=c(.9,.9,.9,.9),
  pribeta = matrix(rep(c(10, 1), nch), ncol = nch, nrow = 2),
  #bruittemp=c(.1,.0005),
  bruittemp = c(.1, .05, .03, .02, .01, .009, .005, .001, .0005),
  pribruit = pribruit,
  #prior on noise
  tipprior = list()
)

Param = list(
  Trpossliste = passages,
  nch = nch,
  nph = etats,
  npart = 1500,
  npas = 2000,
  npasfin = 20000,
  Nmin = 750,
  loiini = loiini,
  Prob = c(1, .01, .1, .5, 1, 0, .01),
  Probfin = c(1, .01, .1, .5, 1, .1, .01),
  agemax = c(100, 1 / 2.5),
  tiplabel = quelleslangues,
  batches = c(lapply(1:(nlangues), function(x) {
    ncogntot
  })),
  ReconstructClade = list(),
  Cladeage = Contraintesages
)


# This will start the SMC sampler. This line will take a long time to execute.
# The last parameter is the number of cores used on the cluster, which should be adapted to your setup.
VV2 = SMCbruit(Dat, Param, Prior, 29)


save.image(paste("LSasia", date(), ".RData"), version = 2)
