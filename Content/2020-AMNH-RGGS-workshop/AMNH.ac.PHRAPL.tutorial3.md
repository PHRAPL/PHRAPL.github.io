---
layout: default
title: Tutorial 3 - Species Delimitation
parent: AMNH-RGGS-workshop
nav_order: 3
---

# Tutorial 3: Testing species delimitation hypotheses
{: .no_toc }

If the question of interest pertains to the delimitation of putative species that are specified in the dataset, then for a given hypothesized tree, you will want to compare the fit of a "full tree" model for which all coalescence times are estimated (i.e., are non-zero) to "collapsed" models for which one or more coalescence times are set to be zero. To do this, rather than adding models to a `migrationArray`, these collapsed models should be specified when running a `GridSearch`. Thus, for a given set of models in a `migrationArray` that specify fully resolved trees, the collapsing of specific nodes can be specified using `GridSearch`s  `setCollapseZero` argument. This argument inputs a vector that defines which coalescence time parameters in a model one would like to be set to zero. So, for our toy dataset, setting `setCollapseZero = 1` when fitting an ((A,B),C) isolation-only model will effectively collapse populations A and B into a single population; setting `setCollapseZero = 1:2` will simulate a single panmictic population.

Below is an example for how to test the existance of 1, 2, or 3 species using the first three models run above using the toy dataset (see Tutorial 1). We have increased the number of nTrees to 10000 such there are enough simulated trees to calculate the log-likelihood for all the models. These runs will take a few minutes.

---
1. TOC
{:toc}
---

## Prepare data
First, load input data and set relevant parameters:

```
setwd("~/Desktop")
library(phrapl)
library(ape)
  
trees<-read.tree(paste(path.package("phrapl"),"/extdata/trees.tre",sep=""))
assignFile<-read.table(paste(path.package("phrapl"),"/extdata/cladeAssignments.txt",sep=""),
     header=TRUE,stringsAsFactors=FALSE)
     
popAssignments<-list(c(3,3,3))
assignmentsGlobal<-assignFile  
observedTrees<-trees  
popAssignments<-list(c(3,3,3))  
subsamplesPerGene<-10  
outgroup=TRUE  
outgroupPrune=TRUE  
  
observedTrees<-PrepSubsampling(assignmentsGlobal=assignmentsGlobal,observedTrees=observedTrees,
      popAssignments=popAssignments,subsamplesPerGene=subsamplesPerGene,outgroup=outgroup,
      outgroupPrune=outgroupPrune)
subsampleWeights.df<-GetPermutationWeightsAcrossSubsamples(popAssignments=popAssignments,
      observedTrees=observedTrees)
      
popVector<-popAssignments[[1]]  
maxK<-3  
maxMigrationK=1  
maxN0K=1  
maxGrowthK=0  
forceTree=TRUE  
forceSymmetricalMigration=TRUE  
migrationArray<-GenerateMigrationIndividuals(popVector=popVector,maxK=maxK,  
      maxMigrationK=maxMigrationK,maxN0K=maxN0K,maxGrowthK=maxGrowthK,
      forceTree=forceTree,forceSymmetricalMigration=forceSymmetricalMigration)

modelRange<-1:3  
nTrees<-10000  
```

## GridSearch

### Three species models
Then run and save analyses for three species models:

```
result<-GridSearch(migrationArray=migrationArray,modelRange=modelRange,
    popAssignments=popAssignments,observedTrees=observedTrees,
    subsampleWeights.df=subsampleWeights.df,subsamplesPerGene=subsamplesPerGene,
    nTrees=nTrees) 
save(list="result",file="phraplOutput_models1-3_3species.rda")  
```
### Two species models
Run and save analyses for two species models:

```
result<-GridSearch(migrationArray=migrationArray,modelRange=modelRange,
    popAssignments=popAssignments,observedTrees=observedTrees,
    subsampleWeights.df=subsampleWeights.df,subsamplesPerGene=subsamplesPerGene,
    nTrees=nTrees,setCollapseZero=1) 
save(list="result",file="phraplOutput_models1-3_2species.rda") 
```

### Single species models
And then run and save analyses for single species models:

```
result<-GridSearch(migrationArray=migrationArray,modelRange=modelRange,
    popAssignments=popAssignments,observedTrees=observedTrees,
    subsampleWeights.df=subsampleWeights.df,subsamplesPerGene=subsamplesPerGene,
    nTrees=nTrees,setCollapseZero=1:2) 
save(list="result",file="phraplOutput_models1-3_1species.rda") 
```

## Post-processing

### Concatenate results:


```
totalResults<-ConcatenateResults(migrationArray=migrationArray)  
```

(_real order_) | models | AIC | params.K | rank | dAIC | wAIC | params.vector | t1_1.2 | t2_1-2.3 | m1_1.2 | m1_1.3 | m1_2.1 | m1_3.1  | 
---| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---  | 
1 | 1 |  102.409023727  |        2 |     1 |    0.000 |  7.80744156813e-01 |              collapse_1 collapse_2 |  0.885953900419 |  4.28733673377 |              NA  |             NA |              NA   |            NA  | 
4 | 2 |  105.627484186   |       3 |     2 |    3.218 |  1.56217225852e-01 |  collapse_1 collapse_2 migration_1 |  1.429951012557 |  5.54272642869 |  0.137939540376 |              NA |  0.137939540376  |             NA | 
7 | 3 |  107.441537070  |        3 |     3 |    5.033 |  6.30386173283e-02 |  collapse_1 collapse_2 migration_1 |  0.656775270185 |  7.64571997328 |              NA |  0.109671500725 |              NA |  0.109671500725 | 
8 | 3 |  155.202041510  |        2 |     4 |   52.793 |  2.68320868533e-12  |            collapse_2 migration_1 |  0.000000000000 |  5.01247343036 |              NA |  1.539211074576 |              NA |  1.539211074576 | 
2 | 1 |  155.317868653  |        1 |     5 |   52.909 |  2.53200973488e-12  |                        collapse_2 |  0.000000000000 |  6.94471465031 |              NA |              NA |              NA |              NA | 
5 | 2 |  155.861451635  |        2  |    6 |   53.452 |  1.92998715801e-12  |            collapse_2 migration_1 |  0.000000000000 |  5.74580965294 |  1.544042397162 |              NA |  1.544042397162 |              NA | 
9 | 3 |  265.459991526  |        1  |    7 |  163.051 |  3.06502455529e-36  |                       migration_1 |  0.000000000000 |  0.00000000000 |              NA |  2.149987786716 |              NA |  2.149987786716 | 
6 | 2 |  292.400803632  |        1 |     8 |  189.992 |  4.32782946843e-42  |                       migration_1 |  0.000000000000 |  0.00000000000 |  0.460000005209 |              NA |  0.460000005209             NA | 
3 | 1 |  340.400803632  |        0  |    9 |  237.992 |  1.63381385280e-52  |                                   0.000000000000 |  0.00000000000 |              NA |              NA |              NA |              NA | 


### Model averages

And model average parameters

```
modelAverages<-CalculateModelAverages(totalResults, parmStartCol = 8, keep.na = TRUE)  
```

t1_1.2 | t2_1-2.3 | m1_1.2 | m1_1.3 | m1_2.1 | m1_3.1 | 
--- | --- | --- | --- | --- | --- | 
0.956488516171 | 4.69515806516 | 0.0215485323358 | 0.00691353977014 | 0.0215485323358 | 0.00691353977014 | 


## Calculate genealogical divergence index (_gdi_)
Finally, you can calculate the genealogical divergence index (gdi) between populations A and B using the model averaged parameter values of migration rate and divergence time. This index is a composite metric that estimates overall divergence (between 0 and 1) from the combined effects of genetic drift and gene flow (0 = panmixia; 1 = strong divergence).

```
gdi<-CalculateGdi(modelAverages$t1_1.2,modelAverages$m1_1.2)
```

method | x | n | mean | lower | upper | 
--- | --- | --- | --- | --- | --- | 
exact | 8827 | 10000 | 0.824137931034 | 0.81444013188 | 0.833500169879 | 

See [this paper](https://academic.oup.com/sysbio/article/66/5/799/2726792) for a reference on how to interpret this metric, and pay special attention to [Fig 6.](https://academic.oup.com/view-large/figure/95743527/syw117f6.tif)
_"Thus, as a rule of thumb, gdi values less than 0.2 suggest that a single species exists; gdi values above 0.7 suggest there are two species. Values in between indicate ambiguous delimitation (but of course, with values near 0.2 and 0.7 providing stronger or weaker evidence for a single species, respectively), which reflects the reality that there exists a speciation gray zone, where a definitive answer cannot easily be found."_