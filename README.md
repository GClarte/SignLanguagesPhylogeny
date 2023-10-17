# Phylogenies with Matricial Datasets

This is the implementation of the numerical methods described in the paper "Computational phylogenetics reveal histories of sign languages". To reproduce the results of the paper see below, to run your own analysis, the best idea is to read the GiveATry.R file, that contains enough comments to describe each parameter choice and data information needed.

The code was working on R 4.2.1, the only library used is "parallel" to run on parallel CPU.

## functions

This folder contains all the core functions needed to run the analysis. It is not needed to open this file to replicate the experiments on this datasets or make small changes (values of parameters, use of another dataset, etc.), but for more substantial changes of the model it can be needed. In case of need, you can contact Grégoire Clarté.

## datasets

Contains the datasets used for the experiments presented in the papers. The data consists in a csv file, where each line corresponds to a word, the first columns describe the meanings and language they belong to. We also include an excel file containing exactly the same data.

## videos

Contains examples of the original videos which have been treated to produce the dataset presented in datasets, the folder also contains a video explaining the encoding process.

## Results

Contains the resulting samples of trees (EBNZ_1.nex, EBNZ_2.nex and Asia.nex), and the resulting consensus trees (EBNZ_1_annote.nex, EBNZ_2_annote.nex and Asia_annote.nex). We did not include the full output of the SMC as the weight would be too important.

## Reproducing the results of the papers

The R files at the root corresponds to the different experiments carried in the papers:
- AsieLS.R corresponds to the code for the study of the Asian dataset
- EuropeBNZLS.R corresponds to the code for the study of the European and British-New Zealand languages.


To reproduce the results, these R scripts needs to be started. We recommend using a large cluster as the running time is about a day with 40 cores, for this it is required to check the last lines of the code to change the number of cores used. The results can then be interpreted directly in R (for example for the plots of the parameters) or using additional softwares (requiring to save the tree files corresponding to the phlogenies. We stored our results in the Results folder).

Pay attention in the R files to the parameters used, especially the prior information on ages (written "Contraintesages"), the set of characters used ("qui") in the study, and the language included ("quelleslangues").

## Using the code

The whole process: formatting the dataset for the inference, setting of the parameters, launching of the SMC, description of the output, is described in the tutorial file GiveATry.R file. 
To get the different plots, please refer to the GiveAPlot.R file. It plots all the parameters and produces the .nex files needed for subsequent phylogenetical analysis, plots of the phylogenies will require additional softwares to build the consensus tree, annotate it, and plot the densitrees (we used
