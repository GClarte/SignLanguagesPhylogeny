## Phylogenies with Matricial Datasets

This is the implementation of the numerical methods described in the paper "Computational phylogenetics reveal histories of sign languages". To reproduce the results of the paper see below, to run your own analysis, the best idea is to read the GiveATry.R file, that contains enough comments to describe each parameter choice and data information needed.

# fonctions

This folder contains all the functions needed to run the analysis. The code might be difficult to read, but it is not needed to open the folder to run the experiments. For additional changes to the code, you might want to contact Grégoire Clarté wherever he is at the moment where you are reading this file.

# datasets

Contains the datasets used for the experiments presented in the papers. The data consists in a csv file, where each line corresponds to a word, the first columns describe the meanings and language they belong to. Additional details can be found in the GiveATry.R file, please open this file.

# Results

Contains the resulting samples of trees (EBNZ_1.nex, EBNZ_2.nex and Asia.nex), and the resulting consensus trees (EBNZ_1_annote.nex, EBNZ_2_annote.nex and Asia_annote.nex). We did not include the full output of the SMC as the weight would be too important.

# Reproducing the results of the papers

The R files at the root corresponds to the different experiments carried in the papers:
- AsieLS.R corresponds to the code for the study of the Asian dataset
- EuropeBNZLS.R corresponds to the code for the study of the European and British-New Zealand languages.

To reproduce the results, these R scripts needs to be started. We recommend using a large cluster as the running time is about a day with 40 cores, for this it is required to check the last lines of the code to change the number of cores used. The results can then be interpreted directly in R (for example for the plots of the parameters) or using additional softwares (requiring to save the tree files, that can be found in Results).

# Using the code

The whole process: specification of the dataset, setting of the parameters, launching of the SMC, description of the output, is described in the GiveATry.R file. 
To get the different plots, please refer to the GiveAPlot.R file. It plots all the parameters and produces the .nex files needed for following phylogenetical analysis, plots will require additional softwares to build the consensus tree, annotate it, and plot the densitrees.
