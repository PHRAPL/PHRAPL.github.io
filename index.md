---
layout: default
title: Home
nav_order: 1
description: "Introduction"
permalink: /
---


#**<font color='#006579'>PH</font>ylogeographic** 
#**infe<font color='#006579'>R</font>ence** 
#**using**
#**<font color='#006579'>AP</font>proximated**
#**<font color='#006579'>L</font>ikelihood**

=======

<img src="https://github.com/ariadnamorales/phrapl-manual/blob/master/phrapl_logo.png?raw=true" width="200" height="200" />


## Why to use **<font color='#006579'>PHRAPL</font>**?
Phylogeographic research aims to understand the recent history of species. Over the last decades, researchers have increasingly incorporated demographic models in order to estimate parameters (i.e., divergence times, population sizes, and rates of migration and expansion) that can contribute phylogeographic inference. 
Typically, this is conducted via the use of software packages that contain specified models (n-island models or fixed topologies). Alternatively, simulation-based approaches allow researchers to customize models for the particular details of their system, and may be useful in testing preexisting biogeographic hypotheses. Because the demographic model is central to the analysis in either case, researchers may wish to assess the appropriateness of their model to the data. `PHRAPL` is designed to give such a tool to researchers. 

## How **<font color='#006579'>PHRAPL</font>** works?
`PHRAPL` simulates genealogies under a wide range of demographic models and compares the empirical genealogies to the simulated gene tree distributions. Demographic models that are probable given the data will contain many genealogies that match the estimated gene trees. Because the proportion of matching gene trees for a given model is equivalent to the probability of the data given the model and parameter values, we can use this value in an information theoretic framework to evaluate the relative weight of all models. This provides the researcher with an independent assessment of both the best model, given the data as well as the ability to calculate the model likelihoods of classes of models (e.g., n-island vs. isolation models).


### CODE
`PHRAPL` is written in `R`, but it uses `perl` and `ms` to perform simulations. The pre-CRAN (code under development) can be found in [github.](https://github.com/bomeara/phrapl)


### CITATION
- Jackson N, Morales AE, Carstens BC, O'Meara BC (2017) [PHRAPL: Phylogeographic Inference using Approximate likelihoods](https://academic.oup.com/sysbio/article/66/6/1045/2999288). Systematic Biology. 66:1045-1053.


### OTHER REFERENCES
- Jackson N, Carstens BC, Morales AE, O’Meara BC (2017) [Species delimitation with gene flow](https://academic.oup.com/sysbio/article/66/5/799/2726792?searchresult=1). Systematic Biology. 66:799-812.
- Morales AE, Jackson N, Dewey T, O’Meara BC, Carstens BC (2017) [Speciation with gene flow in North American Myotis bats](https://academic.oup.com/sysbio/article/66/3/440/2682289). Systematic Biology. 66:440-452.
- Carstens BC, Morales AE, Jackson N, O’Meara BC (2017) [Objective choice of Phylogeographic Models](https://www.sciencedirect.com/science/article/pii/S1055790317303160?via%3Dihub). Molecular Phylogenetics and Evolution. 116:136-140.



<sub>Template acknowledgement: <a href="https://github.com/pmarsceill/just-the-docs">Just the Docs</a>, a documentation theme for Jekyll.<sub>
