# Phylogenies with Matricial Datasets

This is the implementation of the numerical methods described in the paper "Computational phylogenetics reveal histories of sign languages". Instructions are given below to reproduce all the results from the paper, and to re-use the methodology on other data sets.

The code requires R version 4.2.1, with the "parallel" library.

## functions

This folder contains all the core functions needed to run the analysis. Users wishing to reproduce the results from the paper, or to make small changes (parameters values, different data) will not need to open these files. These files may be of use to researchers wishing to make more substantial changes to the model. In case of need, you can contact Grégoire Clarté.

## datasets

This folder contains the datasets used for the experiments presented in the papers. The data consists in a csv file, where each line corresponds to a word in a language. The first columns give the meaning and the language; the other columns give all the characters used to encode the word. We also include an excel file containing exactly the same data.

## videos

Contains examples of the original videos which have been treated to produce the dataset presented in datasets. The folder also contains a video explaining the encoding process.

## Results

Contains the resulting samples of trees (EBNZ_1.nex, EBNZ_2.nex and Asia.nex), and the resulting consensus trees (EBNZ_1_annote.nex, EBNZ_2_annote.nex and Asia_annote.nex). We did not include the full output of the SMC as the files are too large; these files can be reproduced using the code below.

## Reproducing results from the paper

R scripts are available at the root of the folder. They correspond to different analyses:
- AsieLS.R corresponds to the code for the study of the Asian dataset
- EuropeBNZLS.R corresponds to the code for the study of the European sign languages (with New Zealand sign language).
- ToutesLS.R includes Asian, European and New Zealand sign languages.


To reproduce the results of the paper, execute the relevant R script. We recommend using a large cluster as the running time is about a day with 40 cores; the number of cores can be changed in the last lines of each script.

In all cases, you should then execute the GiveAPlot.R script, which post-processes the output and saves the trees in the Nexus format for interpretation in standard phylogenetic software.

Users wishing to make slight modifications to the analyses will presumably be particularly interested in changing the following parameter values:
- prior information on ages are set in the object "Contraintesages"
- the set of characters used is given by object "qui"
- the set of languages included is set by object "quelleslangues"


## Using the code on other data sets

The script GiveATry.R gives a template which can be adapted to apply the method to other data. To use this script, re-users will need to fill out certain parameters, which are marked explicitly. This script goes through the whole process: formatting the dataset for the inference, setting of the parameters, launching of the SMC, description of the output. 

Here too, you should then execute the GiveAPlot.R script. It plots all the parameters and produces the .nex files needed for subsequent phylogenetic analyses. The resulting .nex files can be fed to standard phylogenetic software.

## Thanks

The implementation of the Dirichlet distribution comes from the gtools package.
