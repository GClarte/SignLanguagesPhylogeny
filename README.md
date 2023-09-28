## Phylogenies with Matricial Datasets

This is the implementation of the numerical methods described in the paper "Computational phylogenetics reveal histories of sign languages" to sample from the posterior associated with the model we propose to deal with Matricial datasets.

# "fonctions"

This folder contains all the functions needed to run the analysis. The code might be difficult to read, but it is not needed to open the folder to run the analysis.

# "datasets"

Contains the datasets used for the experiments presented in the papers.

# Reproducing the results of the papers

The R files at the root corresponds to the different experiments carried in the papers:
- AsieLS.R corresponds to the code for the study of the Asian dataset
- EuropeBNZLS.R corresponds to the code for the study of the European and British-New Zealand languages.

# Using the code

Everything is described in the GiveATry.R file, we did not include the data importation part, as this can be quite dependent on the datasets.
GiveAPlot.R plots all the parameters, and produces the .nex files needed for following studies.
