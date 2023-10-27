# For any questions, please e-mail clarte $at$ ceremade $dot$ dauphiine $dot$ fr or clarteg $at$ gmail $dot$ com


# This file gives a template to adapt the Phylogenies from Matrices methodology to other data sets.
# To run it, several objects need to specified. They are listed with the TOFILL keyword.

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

#========================================================================================
#                             Preparing the Data
#=========================================================================================

#The dataset comes in the form of a CSV file, as presented in example_data.csv


#vector of name of the languages, in the order of the dataset,
#it also gives the name of the leaves of the tree, the number in the constraints
quelleslangues = c() # TOFILL

nlangues = length(quelleslangues)


# vector of the characters to include in the study, listed as column numbers in the data file
qui = c() # TOFILL


#vector of the indexes of the meanings
quelmeanings = c() # TOFILL

#number of characters in the dataset
nch = length(qui)


#The dataset must be transformed in order to have a list of length nch
#each element of the list is a matrix with nlanguages columns and a row for each meaning
#NA are possible in any of the matrices
Dat = list() # TOFILL

#this function creates all the pairwise transformations possible. This code does not need to be optimized.
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



#etats is a vector of size nch containing the number of values possible for each character
etats = c() # TOFILL

#passages is a list of length nch, for each character it gives the possible transformations
#for each character the possible transformation is a list of vectors of size 2 giving the starting and endind value of the transformation
#the following code transforms "etats" in order to have all pairwise transformations possible
passages = lapply(etats, function(x) {
  S = toutestransf(x)
  return(split(S, col(S)))
})

#uniform prior, always useful, we use this one
priorunif <- function(tr, agemax) {
  tips = Prior$tipprior
  cladeage = Param$Cladeage
  if (!all(sapply(tips, function(tip) {
    is.monophyletic.perso(tr, tip)
  })) ||
  !all(sapply(cladeage, function(clade) {
    ageconstraints(tr, clade)
  }))) {
    #tests for the constrains on clades and on internal ages
    return(-Inf)
  } else {
    hhh = hauteur(tr, tr$root[[2]])
    if (!(hhh > 200 &
          hhh < 1000)) {
      #add here the constrains on the root age.
      return(-Inf)
    } else {
      if (length(cladeage) == 0) {
        return(-length(tr$tip.label) * log(hauteur(tr, tr$root[[2]])))
      } else {
        return(-(length(tr$tip.label) - sum(sapply(cladeage, function(x) {
          length(x[[1]]) - 1
        }))) * log(hauteur(tr, tr$root[[2]])))
      }
    }
  }
}

#with gamma prior on the length of the tree, same as before, on simulated data it gives quite similar results
priorgamma <- function(tr, agemax) {
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
    if (!(hhh > 200 & hhh < 1000)) {
      return(-Inf)
    } else {
      if (length(cladeage) == 0) {
        return(-length(tr$tip.label) * log(hauteur(tr, tr$root[[2]])) +
                 dgamma(sum(tr$edge.length), 10, 1 / 90, log = T))
      } else {
        return(-(length(tr$tip.label) - sum(sapply(cladeage, function(x) {
          length(x[[1]]) - 1
        }))) * log(hauteur(tr, tr$root[[2]])) +
          dgamma(sum(tr$edge.length), 10, 1 / 90, log = T))
      }
    }
  }
}


#parameters of the algorithm

#number of meanings
ncogntot = nrow(Dtot[[1]])
#of characters
nch = length(etats)

#hyperparameter on the dirichlet distribution on the possible transformation, here uniform
hyperpbini = lapply(passages, function(x) {
  rep(1, length(x))
})
#prior on the transformation frequency, for each character, these are the parameter of a Gamma distribution
#beware that we use this prior ONLY because we have chosen a certain value for the root age
#in practice if the root age is around x, there should be in mean barplot(x*rgamma(1000,Pri[[i]][1],Pri[[i]][[2]]))
#transformations for character i
#here we propose this prior suitable for trees with root of age around a century, other parameters could be proposed,
#the general format is : a list of length nch, each element is a 4 sized vector, the two first are the parameters of the gamma prior
#and the last ones are the bound of the truncature.
Pri = lapply(1:nch, function(x) {
  c(1, 500, 0, 1)
})
#prior on Rho, also Gamma, idem to check the number of apparitions PER MEANING.
Prirho = c(1, 5000)
#initial distribution \pi_0 for each character, we choose a uniform distribution here
loiini = lapply(etats, function(x) {
  rep(1 / x, x)
})
#prior on \nu the noise parameter, it is a list givint the computation of the prior and the proposal
pribruit = list(list(function() { runif(1, 10 ^ -6, 10 ^ -2)}, 
                     function(x) {
                       #compute the value of the prior
                       if (x > 10 ^ -6 & x < 10 ^ -2) {
                         return(-log(x))
                       } else {
                         return(-Inf)
                       }
                     }),
                list(function(x) {
                  rnorm(1, x, .1)
                }, #an idoneous proposal
                function(x, y) {
                  dnorm(x, y, .1, log = T)
                }))


#we put everything in two nice lists
#in Prior all the prior informations
Prior = list(
  Prirho = Prirho,
  #prior on Rho
  hyperpbini = hyperpbini,
  #prior on the transf choice parameter
  Prila = Pri,
  #prior on the transf rate
  priortree = priorunif,
  #prior on the tree
  pribeta = matrix(rep(c(10, 1), nch), ncol = nch, nrow = 2),
  #prior on beta the acceptance rate, here Beta(10,1)
  #bruittemp=0.05-log(1:30)*(.05-.0005)/log(30),
  bruittemp = c(.1, .05, .03, .02, .01, .009, .008, .007, .006, .005, .001),
  #decreasing scheme for the noise, the final value must be in a non-zero area of the function pribruit
  pribruit = pribruit,
  #prior on noise
  tipprior = list()
) #here you can include several vectors indicating the indices of the leaves that must form a clade
#eg : tipprior=list(c(1,2,3),c(6,7))
#tipprior=list())
#
Param = list(
  Trpossliste = passages,
  #possible transformations
  nch = nch,
  nph = etats,
  npart = 3000,
  #number of particles. BEWARE, large value (over 4000 in my case) can cause memory overload, and R does not warn, see later
  npas = 5000,
  #number of step between each of the change in the noise value
  npasfin = 10000,
  #we finish with a very long MCMC
  Nmin = 2500,
  #if the ISS is smaller than Nmin, we resample, in practice we always resample
  loiini = loiini,
  Prob = c(1, .01, .1, .5, 1, 0, .01),
  #probability of the proposals in the MCMC when the noise is still fixes
  #in the order : transformations, parameters, length of the edges, apparitions, topology, noise, beta
  #the noise cannot change so Prob[7]=0
  Probfin = c(1, .01, .1, .5, 1, .1, .01),
  #when the noise can change Prob[7]!=0
  agemax = c(100, 1 / 3),
  #there are two possibilities, first you can fix it, giving a single number,
  #or you can lear it, in which case you must put two aruyments that will be used as parameter of a gamma proposal
  #just run hist(rgamma(1000,agemax[1],agemax[2])) and check the distribution of the rootage is close to what you want
  tiplabel = quelleslangues,
  #label of the tips IN THE SAME ORDER THAN IN THE DATASET
  batches = c(lapply(1:(nlangues), function(x) {
    ncogntot
  })),
  #don't change that, even if it's not used anymore
  ReconstructClade = list(),
  #list of vector of indices of the leaves whose common ancestor state must be reconstructed
  Cladeage = list()#list of constraints on clade age,
  #eg : ReconstructClade=list(c(10,9),c(6,2)),
  #Cladeage=list(list( c(10,9),c(10,200)),list( c(8,150),c(20,300)))#each element is a list whose first element is a vector of indices of leaves whose
  #common ancestor is constrained, and the second is the vector of lower and upper bound on the age
)

#4 parameters for the study. The three first are straightforward, the last one is the number of cores attributed on your cluster
VV2 = SMCbruit(Dat, Param, Prior, 29)

#The result should be a list of three elements the first one being of length npart
#If you only get a list of npart elements, there has been an error somewhere and this error will appear in the list somewhere
#if some elements of the lists contains NULL (ie are of length 0) there has been a memory overload, and you will have to decrease the number of particles

#save the image as it is intended to run on a cluster.
save.image(paste("test", date(), ".RData"), version = 2)

#the first element of the list is the result proper. Each element corresponds to a particle final state.
# For each particle you will find the following values :
#Tr = the tree
# La = lambda transformation rate
# Rho = apparition rate
# P the probability for each transformations
#X = the nature and order of the transformations
#Tps = the position of the transformations relatively to the edge
#beta = value of the acceptance probability of the transformations
#noise = the noise parameter
#NL = the number of apparition gor each meaning for each edge
#weight = the final weight of the particle

#the second element containts the weights of each particle at each resampling step
#the last one containts for each particle at each resampling step which particle was selected as ancestor.
