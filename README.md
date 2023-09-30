## Phylogenies with Matricial Datasets

This is the implementation of the numerical methods described in the paper "Computational phylogenetics reveal histories of sign languages" to sample from the posterior associated with the model we propose to deal with Matricial datasets.

# "fonctions"

This folder contains all the functions needed to run the analysis. The code might be difficult to read, but it is not needed to open the folder to run the analysis.

# "datasets"

Contains the datasets used for the experiments presented in the papers. The data consists in a csv file, where each line corresponds to a word, the first columns describe the meanings and language they belong to.

# "Results"

Contains the resulting samples of trees (EBNZ_1.nex, EBNZ_2.nex and Asia.nex), and the resulting consensus trees (EBNZ_1_annote.nex, EBNZ_2_annote.nex and Asia_annote.nex). We did not include the full output of the SMC as the weight would be too important.

# Reproducing the results of the papers

The R files at the root corresponds to the different experiments carried in the papers:
- AsieLS.R corresponds to the code for the study of the Asian dataset
- EuropeBNZLS.R corresponds to the code for the study of the European and British-New Zealand languages.

To reproduce the results, these R scripts needs to be started. We recommend using a large cluster as the running time is about a day with 40 cores, for this it is required to check the last lines of the code to change the number of cores used. The results can then be interpreted directly in R (for example for the plots of the parameters) or using additional softwares (requiring to save the tree files, that can be found in Results).

# Using the code

Everything is described in the GiveATry.R file, we did not include the data importation part, as this can be quite dependent on the datasets.
GiveAPlot.R plots all the parameters, and produces the .nex files needed for following studies.
