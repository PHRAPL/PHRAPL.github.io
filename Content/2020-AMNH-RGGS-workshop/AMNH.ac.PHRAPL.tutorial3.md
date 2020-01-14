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
library(phrapl)
  
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
### Model averages

And model average parameters

```
modelAverages<-CalculateModelAverages(totalResults)  
```

## Calculate genealogical divergence index (_gdi_)
Finally, you can calculate the genealogical divergence index (gdi) between populations A and B using the model averaged parameter values of migration rate and divergence time. This index is a composite metric that estimates overall divergence (between 0 and 1) from the combined effects of genetic drift and gene flow (0 = panmixia; 1 = strong divergence).

```
gdi<-CalculateGdi(modelAverages$t1_1.2,modelAverages$m1_1.2)
```
