
title: “Semi-Supervised Elastic Net (ssenet)” author: “Amrit Singh”
date: “10 March, 2020”

``` r
knitr::opts_chunk$set(echo = TRUE, cache = FALSE, warning = FALSE, message = FALSE)

library(magrittr); 
library(knitr);
library(ssenet); ## devtools::install_github("singha53/ssenet")
```

    ## 
    ## Attaching package: 'ssenet'

    ## The following object is masked from 'package:stats':
    ## 
    ##     predict

``` r
library(ggplot2);
library(dplyr);
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(UpSetR);
```

# Analysis for the Abstract Submission to \#BIRSBioIntegration [Mathematical Frameworks for Integrative Analysis of Emerging Biological Data Types](https://www.birs.ca/events/2020/5-day-workshops/20w5197)

## Zhu *et al* 2018: [seqFISH paper](https://www.nature.com/articles/nbt.4260)

  - 43 genes (Supplementary Table 2) to map cell types in the seqFISH
    data:

### Step 1: mapping scRNASeq celltypes on seqFISH data

  - Randomly selected a subset of genes from the list of differentially
    expressed and applied a multiclass support vector machine; perform
    evaluated using cross-validation; 43 genes were used to map
    cell-types in the seqFISH data
  - Applied the SVM classification model to the bias-correct, quantile
    normalized seqFISH data to assign cell types.

<!-- end list -->

``` r
include_graphics("inst/extdata/suppTable2.png")
```

![](inst/extdata/suppTable2.png)<!-- -->

``` r
selectedGenes <- c("fbll1", "itpr2", "vps13c", "tnfrsf1b", "sox2",
  "hdx", "wrn", "sumf2", "vmn1r65", "rhob",
  "mrgprb1", "calb1", "pld1", "laptm5", "tbr1",
  "slc5a7", "abca9", "ankle1", "olr1", 
  "cecr2", "cpne5", "blzf1", "mertk",
  "nell1", "npy2r", "cdc5l", "slco1c1",
  "pax6", "cldn5", "cyp2j5", "mfge8",
  "col5a1", "bmpr1b", "rrm2", "gja1",
  "dcx", "spag6", "csf2rb2", "gda",
  "arhgef26", "slc4a8", "gm805", "omg")

plot(coord, col = mixOmics::color.mixo(as.numeric(seqfishLabels$V3)), pch = 21, 
  xlab = "x-coordinates", ylab = "y-coordinates")
points(coord, col = mixOmics::color.mixo(as.numeric(seqfishLabels$V3)), pch = 19)
```

![](README_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

### Step 2: a systemic approach to identify multicellular niche

  - HMRF (Hidden Markov Random Field): Systematically dissect the
    contributions of microenvironments on gene expression variation
  - Divide the visual cortex into domains with coherent gene expression
    patterns
  - HMRF enables the detection of spatial domains by systemically
    comparing the gene signature of each cell with its surrounding to
    search for coherent patterns
  - Domain state of each cell was influence by 2 sources: 1) gene
    expression pattern and 2) domain states of neighbouring cells
  - To enhance spatial domain detection 11 highly cell-specific genes
    were removed
  - HMRF revealed 9 spatial domains; four domains were in the outer
    layers of the cortex ( O1-O4), four domains were located on the
    inside of the cortex (I1A, I1B, I2, I3), domain IS was sporadically
    distributed across the inner layers of the cortex.

### Step 3: interactions between cell-type and spatial environment

  - The same cell-types expressed different genes depending on the
    domain state of the cell.
  - By visual inspection there were notable morphological variations
    near the boundary between different domains at multiple regions.

# Questions for the BIRSBiointegration workshop:

## 1\) Can scRNA-seq data be overlaid onto seqFISH for resolution enhancement?

## 2\) What is the minimal number of genes needed for data integration?

## 3\) Are there signatures of cellular co-localization or spatial coordinates in non-spatial scRNA-seq data?

### remove cells with little data

``` r
# set constants
M = 5;
iter = 5;
ncores = 5;
alpha = 1;
lambda_nfolds = 3;
family = "multinomial";
filter = "none";
max.iter = 50;
perc.full = 1;
thr.conf = 0.5;

## minimum number of samples required per cell-type class (required for hyperparameter tuning and cross-validation)
round(table(scrnaseqLabels$V1)/M/lambda_nfolds, 0)  # remove Oligodendrocyte.2
```

    ## 
    ##            Astrocyte     Endothelial Cell    GABA-ergic Neuron 
    ##                    3                    2                   51 
    ## Glutamatergic Neuron            Microglia    Oligodendrocyte.1 
    ##                   54                    1                    1 
    ##    Oligodendrocyte.2    Oligodendrocyte.3 
    ##                    0                    2

``` r
keepIndices <- which(scrnaseqLabels$V1 != "Oligodendrocyte.2")
xscrnaseq <- scrnaseq[, keepIndices]
yscrnaseq <- droplevels(scrnaseqLabels$V1[keepIndices])
```

## Apply Enet to scRNAseq data and apply to seqFISH to determine cell-type labels

``` r
fitEnet <- enet(xtrain = t(xscrnaseq), ytrain = yscrnaseq, alpha = alpha, lambda = NULL, lambda_nfolds = lambda_nfolds, family = "multinomial", filter = filter)
cvEnet <- predict(object = fitEnet, M = M, iter = iter, ncores = ncores)
cvEnet$perf
```

    ## # A tibble: 9 x 3
    ##   ErrName                Mean      SD
    ##   <chr>                 <dbl>   <dbl>
    ## 1 Astrocyte            0.0419 0.0104 
    ## 2 BER                  0.158  0.00478
    ## 3 Endothelial Cell     0.179  0.0378 
    ## 4 ER                   0.0578 0.00255
    ## 5 GABA-ergic Neuron    0.0520 0.00390
    ## 6 Glutamatergic Neuron 0.0409 0.00319
    ## 7 Microglia            0      0      
    ## 8 Oligodendrocyte.1    0.421  0      
    ## 9 Oligodendrocyte.3    0.368  0.0368

## Apply Semi-supervised Enet to scRNAseq+seqFISH data to determine cell-type labels

``` r
fitSSEnet <- ssenet(xtrain = t(cbind(xscrnaseq, seqfish)), 
  ytrain=factor(c(as.character(yscrnaseq), rep(NA, ncol(seqfish)))), 
  alpha = alpha, lambda = fitEnet$lambda, lambda_nfolds = lambda_nfolds, family = "multinomial", 
  filter = filter,
  max.iter = max.iter, perc.full = perc.full, thr.conf = thr.conf)
cvSSEnet <- predict(object = fitSSEnet, M = M, iter = iter, ncores = ncores)
cvSSEnet$perf
```

    ## # A tibble: 9 x 3
    ##   ErrName                Mean      SD
    ##   <chr>                 <dbl>   <dbl>
    ## 1 Astrocyte            0.0419 0.0195 
    ## 2 BER                  0.243  0.0263 
    ## 3 Endothelial Cell     0.283  0.0289 
    ## 4 ER                   0.0733 0.00518
    ## 5 GABA-ergic Neuron    0.0520 0.00150
    ## 6 Glutamatergic Neuron 0.0554 0.00603
    ## 7 Microglia            0.0455 0.0455 
    ## 8 Oligodendrocyte.1    0.653  0.109  
    ## 9 Oligodendrocyte.3    0.574  0.0421

## Compare supervised and semi-supervised Enet (Enet and SSEnet) performance using cross-validation

``` r
cvErr <- rbind(cvEnet$perf, cvSSEnet$perf) %>% 
  mutate(method = rep(c("Enet", "SSEnet"), each = nrow(cvEnet$perf))) %>% 
  mutate(ErrName = factor(ErrName, c(levels(yscrnaseq), "ER", "BER")))
pd <- position_dodge(0.5)
cvErr %>% 
  ggplot(aes(x = ErrName, y = Mean, color = method)) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size = 3) +
  theme_bw() +
  ylab("Average error rate (2x2 cross-validation)") +
  xlab("Cell-type, ER (error rate), BER (balanced error rate)") +
  customTheme(sizeStripFont = 15, xAngle = 40, hjust = 1, vjust = 1, 
    xSize = 10, ySize = 10, xAxisSize = 15, yAxisSize = 15) +
  ylab("Error") +
  xlab("Cell-type, ER (error rate), BER (balanced error rate)")
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

| ErrName              |      Mean |        SD | method |
| :------------------- | --------: | --------: | :----- |
| Astrocyte            | 0.0418605 | 0.0104003 | Enet   |
| BER                  | 0.1575556 | 0.0047798 | Enet   |
| Endothelial Cell     | 0.1793103 | 0.0377740 | Enet   |
| ER                   | 0.0577752 | 0.0025520 | Enet   |
| GABA-ergic Neuron    | 0.0520368 | 0.0038981 | Enet   |
| Glutamatergic Neuron | 0.0408867 | 0.0031877 | Enet   |
| Microglia            | 0.0000000 | 0.0000000 | Enet   |
| Oligodendrocyte.1    | 0.4210526 | 0.0000000 | Enet   |
| Oligodendrocyte.3    | 0.3677419 | 0.0367799 | Enet   |
| Astrocyte            | 0.0418605 | 0.0194572 | SSEnet |
| BER                  | 0.2434792 | 0.0262661 | SSEnet |
| Endothelial Cell     | 0.2827586 | 0.0288503 | SSEnet |
| ER                   | 0.0732673 | 0.0051831 | SSEnet |
| GABA-ergic Neuron    | 0.0520368 | 0.0014983 | SSEnet |
| Glutamatergic Neuron | 0.0554187 | 0.0060332 | SSEnet |
| Microglia            | 0.0454545 | 0.0454545 | SSEnet |
| Oligodendrocyte.1    | 0.6526316 | 0.1091392 | SSEnet |
| Oligodendrocyte.3    | 0.5741935 | 0.0420594 | SSEnet |

## Overlap between selected features with those used in the Nature paper (SVM)

``` r
panels = list(SVM = selectedGenes, Enet = fitEnet$enet.panel, SSEnet = fitSSEnet$enet.panel)

Input <- fromList(panels)
metadata <- data.frame(Methods=colnames(Input))
upset(Input, sets = colnames(Input))
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

## Abstract

### Your Name

> Amrit Singh

### Slack name on \#BIRSBioIntegration

> Amrit Singh

### Your Position

> trainee (post-doc)

### Name of supervisor

> Kim-Anh Le Cao/Bruce McManus

### Affiliation

> PROOF Centre of Excellence and The University of British Columbia

### Email

> <asingh@hli.ubc.ca>

### Co-authors

> none

### Which dataset(s) did you select for analysis?

> Spatial transcriptomics: seqFISH + scRNA-seq

### Why did you select this dataset(s) for analysis

> recommended by supervisor

### What integrative data analysis question have you addressed with the selected data and why?

> Can scRNA-seq data be overlaid onto seqFISH for resolution
> enhancement?

What are the advantages and performance of your approach? \> The
published approach trained a multiclass SVM on the scRNAseq data and
applied it to the seqFISH data to estimate the cell-types labels. My
approach uses a penalized regression method (glmnet) with a
semi-supervised appraoch in order to build a model using both the
scRNAseq+seqFISH data. This strategy uses a recursive approach that
invovles multiple rounds of training glmnet models using labeled data
(label and imputed) and predicting the cell-type labels of unlabeled
data. At each iteration, cell-type labels with high confidence
(probability \> 0.5) are retained for the next iteration, where a new
glmnet model is trained with the scRNAseq data and seqFISH data with
imputed cell-type labels with high confidence. This process is repeated
until all cell-types in the seqFISH data have been labeled or until 50
iterations have been reached (in order to reduce compute times). The
advantage of this approach is that more data in used for model training
such that the resulting model may generalize better to new data. The
performance of this appraoch was estimated using cross-validation, using
only the scRNAseq data as the test set.

### What were the specific challenges you have encountered so far?

> Compute times are significantly longer for the semi-supervised
> approach for model training. Thus, cross-validation takes even longer.
> The datasets are restricted to 113 genes and therefore the discovery
> space is very limited for the semi-supervised approach to learn
> classification rules that are superior to the supervised approach.

### How are you going to address those challenges?

> Cross-validation was parallelized such that each iteration of
> cross-validation was run on an independent cpu thread. If additional
> data is available for this study it may be better than the current
> results given that the genes are limited to those identified using the
> scRNAseq data only.

### Link to your preliminary code and results on a Github account (optional)

> <https://github.com/singha53/ssenet>

### Additional information you would like the organizers to know

> This is my first time looking at single cell data and this opportunity
> would expose me to knew methods, technologies and research in this
> field.

## References

1)  <https://github.com/mabelc/SSC>

# [Data files](https://github.com/BIRSBiointegration/Hackathon/tree/master/seqFISH)

``` r
# scrnaseq <- read.delim("inst/extdata/tasic_training_b2.txt", row.names = 1, header = FALSE)
# scrnaseqLabels <- read.delim("inst/extdata/tasic_labels.tsv", header = FALSE)
# seqfish <- read.delim("inst/extdata/seqfish_cortex_b2_testing.txt", row.names = 1, header = FALSE)
# seqfishLabels <- read.delim("inst/extdata/seqfish_labels.tsv", row.names = 1, header = FALSE)
# dim(scrnaseq); dim(scrnaseqLabels);
# dim(seqfish); dim(seqfishLabels);
# 
# coord <- read.delim("inst/extdata/fcortex.coordinates.txt", header = FALSE)
# coord <- lapply(1:nrow(coord), function(i){
#   as.numeric(strsplit(as.character(coord[i,]), " ")[[1]][-c(1, 2)])
# }) %>% 
#   do.call(rbind, .)
# dim(coord)
# 
# usethis::use_data(scrnaseq, overwrite = TRUE)
# usethis::use_data(scrnaseqLabels, overwrite = TRUE)
# usethis::use_data(seqfish, overwrite= TRUE)
# usethis::use_data(seqfishLabels, overwrite = TRUE)
# usethis::use_data(coord, overwrite = TRUE)
```
