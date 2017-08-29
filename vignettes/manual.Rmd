---
title: <center> An introduction to exploBATCH <center>
author: <center> Gift Nyamundanda, Pawan Poudel, Yatish Patil and Anguraj Sadanandam <center>
#date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{exploBTACH manual}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

<style>
body {
text-align: justify}
</style>

## Introduction

Explore batch effect, exploBATCH, is an R package for evaluating and correcting for batch effect in genomic studies. It contains two main functions, `findBATCH`, for performing formal statistical tests to detect batch effect, and `correctBATCH`, for removing batch effect in the data. Although here we provide a quick tour of exploBATCH, we recommend reading Nyamundanda $et$ $al$ (2017) for more details on exploBATCH.  

### Overview of exploBATCH
exploBATCH is based on the probabilistic version of principal component analysis which allows for covariates (PPCCA). The `findBATCH` function (within exploBATCH) employs PPCCA to statistically test if samples are distributed according to batches in the principal subspace by applying the following test statistic on each probabilistic principal component (pPC): 

$$ Δ_{bk} = \frac{\beta_{bk}}{\mbox{SE}(\beta_{bk})}, $$
where $\beta_{bk}$ is the regression coefficient that quantifies batch effect $b$ on the $k^{\mbox{th}}$ pPC, and SE is the corresponding standard error. exploBATCH correct for batch effect by calling `correctBATCH` function which subtracts batch effect on each affected pPC as follows,

$$ \underline{u}_{ck} = \underline{u}_{ak}  – \underline{x}_b \beta_{bk} $$

where $\underline{u}_{ak}$ and $\underline{u}_{ck}$ is a vector of scores of the $k^{th}$ pPC affected and corrected for batch effect, respectively, whilst $\underline{x}_b$ is the variable defining batches. The batch effect corrected data is recovered by predicting the observed data using PPCCA conditioning on the scores, **u** = (**uc**, **uu**), where **uc** and **uu** are scores of corrected and  uncorrected pPCs, respectively. 


### Availability and requirements
exploBATCH is freely available on github ([here](https://github.com/syspremed/exploBATCH)) and it is implemented within the R environment. It requires:

- Operating system(s): platform independent.
- The latest version of R: it can be downloaded from the R project webpage `http://www.r-
project.org/`.
- It depends on few other R packages: mclust, mvtnorm, sva, Rcpp, RcppArmadillo, rARPACK, parallel, foreach, doParallel, doMC, ggplot2, RColorBrewer, compiler and devtools.

### Installing exploBATCH

exploBATCH can be installed straight from our github account, syspremed, using devtools R package. Here are two lines of R code to install exploBATCH from syspremed on github.

```{r, eval=F, echo=T}
require(devtools)
install_github("syspremed/exploBATCH")
```

In order to run the two examples in Nyamundanda $et$ $al$ (2017) install the two R packages in syspremed with the two gene expression datasets, `Breast` and `Colon`.

```{r, eval=F, echo=T}
install_github("syspremed/exploBATCHbreast")
install_github("syspremed/exploBATCHcolon") 
```

### Running exploBATCH

The exploBATCH package can be run using a single interface `expBATCH` which calls `findBATCH` and `correctBATCH` to detect and correct for batch effect, respectively.  

             expBATCH(D, batchCL, Conf, mindim, maxdim, method, scale, SDselect)

- `D` is a matrix of expression data with features on the column.

- `batchCL` is a vector of the same sample size as `D` identifying which samples belong to which batch. 

- `Conf` is a vector of the same sample size as `D` identifying samples in different biological classes of interest.

- `mindim` is the minimum number of pPCs to be considered. It should be at least `2`.

- `maxdim` is the maximum number of pPCs to be considered. it should be at least equal to `mindim`.

- `method` is the method for batch correction, either `ppcca` or `combat`. 

- `scale` is the method for scaling the data, `none` no scaling, `unit` unit variance scaling and `pareto` square root of standard deviation scaling.

- `SDselect` is the standard deviation (SD) for filtering the number of genes to improve computational time.

### Output structure of exploBATCH

The results of exploBATCH are saved in the current working directory. A folder `exploBATCHresults` with seven sub-folders of exploBATCH output. Here is the folder structure
of the output of exploBATCH,

    exploBATCHresults: 
      + pcaBeforeBatchCorrection: results of applying PCA on data before batch correction. 
      + ppccaBeforeBatchCorrection: results of applying PPCCA on data before batch correction.
      + findBATCH: containts forest plots showing which pPC has batch effect.
      + correctBATCH: results of batch correction using correctBATCH.
      + assessCorrectBATCH: assessment of results of batch correction using correctBATCH.
      + correComBat: results of batch correction using ComBat.
      + assessComBat: assessment of results of batch correction using ComBat.

## Breast cancer dataset

To demostrate the application of exploBATCH we are applying it on the gene expression data of three batches (30, 22 and 18) of primary human breast tumors. This is the first example in Nyamundanda $et$ $al$ (2017) but for illustration purpose here we used `SDselect=2` to filter genes to 346 with SD>2 to improve computational time. 

Firstly, to run this example you need to load the exploBATCHbreast R package which contains the breast cancer dataset `Breast` and a vector identifying the three batches `batchBreast`.

```{r, eval=F, echo=T}
require(exploBATCH)
require(exploBATCHbreast)  # the package with breast cancer dataset
data(Breast)               # load the breast cancer data        
data(batchBreast)          # load the variable defining batches 
```

After loading the two packages you ready to run exploBATCH. You only need the function `expBATCH` to run exploBATCH,

```{r, eval=F, echo=T}
expBATCH(
  D=Breast,              # the breast cancer expression data matrix
  batchCL=batchBreast,   # the variable identifying the three batches
  Conf=NA,               # no biological variable of interest in this example
  mindim=2,              # the minimum number of pPCs
  maxdim=9,              # the maximum number of pPCs, we don't want this argument to be large,   
                         # otherwise it will slow down exploBATCH.
  method="ppcca",        # use correctBATCH as well as ComBat to correct batch effect
  SDselect=2             # set SD at 2 to reduce computational time.
  )
```

### pcaBeforeBatchCorrection 

This folder contains results of PCA on the breast cancer dataset. PCA results allow for exploring batch effect by simply visually inspecting PCA plots. 

<center><img src="/Users/gnyamundanda/exploBATCHresultsBreastData/pcaBeforeCorrection/pcaBeforeCorrection.png" height="300px" width="310px" /></center>

The clustering of breast tumors in the principal subspace was mainly driven by the three batches. This folder also contains plot of proportion of variation explained by each PC and  pairwise plots of the first ten PCs.

### ppccaBeforeBatchCorrection 

Firstly, exploBATCH has to determine the number of pPCs that can be used to represent this breast cancer dataset. The PPCCA model within exploBATCH employs the Bayesian information criterion (BIC) to identify the optimal number of pPCs. The higher the BIC value the better the model. The BIC plot below shows that the optimal number of pPCs is seven. All batch evaluation and correction for this dataset is based on these first seven pPCs.

<center><img src="/Users/gnyamundanda/exploBATCHresultsBreastData/ppccaBeforeCorrection/NumberOfpPCs.png" height="250px" width="250px" /></center>

This folder `ppccaBeforeBatchCorrection` contains BIC plot and other pairwise plots of pPCs from fitting PPCCA on this breast cancer dataset. 

### findBATCH

The function `findBATCH` employs Jack-knifing approach to quantify uncertainity associated with batch effect. This allows to assess statistical significance of batch effect on each pPC by building 95% confidence intervals (95% CIs) around this parameter of interest. 

<center><img src="/Users/gnyamundanda/exploBATCHresultsBreastData/findBATCH/BatchEffect.png" height="250px" width="250px" /></center>

The forest plot shows the estimated batch effect and the corresponding 95% CIs for each of the seven pPCs. There is significant batch effect on pPC1 and pPC2 since the corresponding 95% CIs don't include zero. The folder `findBATCH` contains this forest plot showing which pPCs are associated with batch effect.

### correctBATCH

After identifying which of the seven pPCs is associated with batch effect `correctBATCH` is employed to subtract this batch effect in the data. The two figures below show PCA plots before and after batch correction.

<center><img src="/Users/gnyamundanda/exploBATCHresultsBreastData/pcaBeforeCorrection/pcaBeforeCorrection.png" height="300px" width="310px" /><img src="/Users/gnyamundanda/exploBATCHresultsBreastData/correctBATCH/2017-08-03_pca_After_PPCCA_correction.png" height="300px" width="310px" /></center>

The folder `correctBATCH` contains breast cancer dataset which has been corrected for batch effect using `correctBATCH`.

### assessCorrectBATCH

After correcting for batch effect the function `findBATCH` is again employed as a formal statistical test to assess if batch effect has been properly corrected. Below are `findBATCH` plots before and after batch correction.

<center><img src="/Users/gnyamundanda/exploBATCHresultsBreastData/findBATCH/BatchEffect.png" height="250px" width="250px" /><img src="/Users/gnyamundanda/exploBATCHresultsBreastData/assessCorrectBATCH/BatchEffect.png" height="250px" width="250px" /></center>

The folder `assessCorrectBATCH` contains results of assessing batch effect after applying `correctBATCH`. 

### correctComBat

The exploBATCH package also allows the user to correct for batch effect using ComBat. Below is a comparison of batch correction using correctBATCH (left panel) and ComBat (right panel).

<center><img src="/Users/gnyamundanda/exploBATCHresultsBreastData/correctBATCH/2017-08-03_pca_After_PPCCA_correction.png" height="300px" width="310px" /><img src="/Users/gnyamundanda/exploBATCHresultsBreastData/correctComBat/pca_After_Combat_Correction.png" height="300px" width="310px" /></center>

The folder `correctComBat` contains breast cancer dataset which has been corrected for batch effect using ComBat.

### assessComBat

After correcting for batch effect using ComBat `findBATCH` can be used to assess if batch effect has been properly corrected. Below are `findBATCH` plots after batch correction using correctBATCH (left panel) and ComBat (right panel).

<center><img src="/Users/gnyamundanda/exploBATCHresultsBreastData/assessCorrectBATCH/BatchEffect.png" height="250px" width="250px" /><img src="/Users/gnyamundanda/exploBATCHresultsBreastData/assessComBat/BatchEffect.png" height="250px" width="250px" /></center>

The folder `assessComBat` contains results of assessing batch effect after applying ComBat. 

## Colon cancer dataset

This example briefly demonstrate that exploBATCH corrects for batch effect without distorting important biological structures in the data. The colon cancer dataset is a gene expression data consisting of two batches of 52 and 58 samples. Most (78%) of these samples are colorectal cancer (CRC) samples whilst the rest are normal.
The question is does batch correction using exploBATCH retains the normal/tumor biological effect in this dataset. 


This is the second example in Nyamundanda $et$ $al$ (2017) but for illustration purpose here we used `SDselect=1.2` to filter genes to 1549 with SD>1.2. To run this example you need to load the exploBATCHcolon R package, the colon dataset `Colon`, a vector identifying the three batches `batchBreast` and a variable defining the biological effect of interest `bioCL`.

```{r, eval=F, echo=T}
require(exploBATCH)
require(exploBATCHcolon)   
data(Colon)                       
data(batchColon)            
data(bioCL)               
```

After loading the two packages you can run exploBATCH on the colon cancer dataset as follows

```{r, eval=F, echo=T}
expBATCH(
  D=Colon,                 
  batchCL=batchColon,     
  Conf=bioCL,              # this is the variable defining the biological effect of interest
  mindim=2,                
  maxdim=11,                  
  method="ppcca",          
  SDselect=1.2             # set SD at 2 to reduce computational time.
  )
```

By inspecting the PCA plot on the top left pannel of the figure below, the samples cluster by batch denoted by the two colours. In addition, the CRC samples (dots) seem to be clustering away from the normal samples (squares). However, we can perform formal statistical tests using `findBATCH` to confirm if there is batch and normal/tumor biological effect. The BIC plot on the top right pannel shows that the first nine pPCs are enough to represent the variability in the colon cancer dataset. The two forest plots on the bottom of the
pannel show that batch and normal/tumor biological effects are intertwined in the first two pPCs. The challenge is to correct this batch effect without blurring the normal/tumor biological effect.

<center><img src="/Users/gnyamundanda/exploBATCHresultsColonData/pcaBeforeCorrection/pcaBeforeCorrection.png" height="280px" width="300px" /><img src="/Users/gnyamundanda/exploBATCHresultsColonData/ppccaBeforeCorrection/NumberOfpPCs.png" height="280px" width="270px" /><img src="/Users/gnyamundanda/exploBATCHresultsColonData/findBATCH/BatchEffect.png" height="280px" width="270px" /><img src="/Users/gnyamundanda/exploBATCHresultsColonData/findBATCH/BioEFFECT-.png" height="280px" width="270px" /></center>

The figures below show results of applying correctBATCH on the colon cancer dataset. The PCA plot highlight that the samples are no longer clustering by batch but the seperation between normal and tumor samples can still be seen. The forest plot from findBATCH confirms that none of the pPCs has batch effect but normal/tumor biological effect is retained in pPC 1 & 2. 

<center><img src="/Users/gnyamundanda/exploBATCHresultsColonData/correctBATCH/2017-08-03_pca_After_PPCCA_correction.png" height="260px" width="300px" /><img src="/Users/gnyamundanda/exploBATCHresultsColonData/assessComBat/BatchEffect.png" height="260px" width="240px" /><img src="/Users/gnyamundanda/exploBATCHresultsColonData/assessComBat/BioEFFECT-.png" height="260px" width="240px" /></center>

## Assessing runtime and memory usage 

Figure A below shows the runtime profiles and memory usage of exploBATCH applied to gene expression datasets with 70-breast cancer samples but varying the number of genes from 50 to 20 000. Figure B highlights the run time profiles of exploBATCH on several simulated data varying the number of samples and genes  

<center><img src="/Users/gnyamundanda/exploBATCHresultsColonData/FigureS5.png" height="260px" width="500px" /><center>

## Reference

- Nyamundanda, G., Poudel, P., Patil, Y. and Sadanandam, A., 2017. A Novel Statistical Method to Diagnose, Quantify and Correct Batch Effects in Genomic Studies.









