Case study 2
================

-   [1 Index](#index)
-   [2 Running the script](#running-the-script)
-   [3 Input data](#input-data)
    -   [3.1 Biological context](#biological-context)
    -   [3.2 Summary](#summary)
    -   [3.3 Metadata](#metadata)
    -   [3.4 Data](#data)
-   [4 Input data](#input-data-1)
-   [5 Data quality control](#data-quality-control)
    -   [5.1 Accounting for missing
        values](#accounting-for-missing-values)
    -   [5.2 Accounting for unwanted
        variation](#accounting-for-unwanted-variation)
-   [6 Single omics analyses](#single-omics-analyses)
    -   [6.1 Parameter tuning](#parameter-tuning)
    -   [6.2 Running the analysis](#running-the-analysis)
    -   [6.3 Explanation of output](#explanation-of-output)
-   [7 Multi omics analyses](#multi-omics-analyses)
    -   [7.1 Parameter tuning](#parameter-tuning-1)
    -   [7.2 Running the analysis](#running-the-analysis-1)
    -   [7.3 Explanation of output](#explanation-of-output-1)
-   [8 Output data](#output-data)
-   [9 Acknowledgements](#acknowledgements)
-   [10 References](#references)

Copyright (c) 2020
<a href="https://orcid.org/0000-0002-9207-0385">Tyrone Chen
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>,
<a href="https://orcid.org/0000-0002-0827-866X">Melcy Philip
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>,
<a href="https://orcid.org/0000-0003-3923-1116">Kim-Anh Lê Cao
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>,
<a href="https://orcid.org/0000-0003-0181-6258">Sonika Tyagi
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>

Code in this package and git repository
<https://gitlab.com/tyagilab/sars-cov-2/> is provided under a [MIT
license](https://opensource.org/licenses/MIT). This documentation is
provided under a [CC-BY-3.0 AU
license](https://creativecommons.org/licenses/by/3.0/au/).

[Visit our lab website here.](https://bioinformaticslab.erc.monash.edu/)
Contact Sonika Tyagi at <sonika.tyagi@monash.edu>.

# 1 Index

-   [Introduction](introduction.html)
-   [Case study 1](case_study_1.html)
-   [Case study 2](case_study_2.html)

# 2 Running the script

Load the library.

``` r
library(multiomics)
```

A script to reproduce our analysis for case study 2 is shown. [You can
also download this
here](https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/src/case_study_2/example.sh).
This may take up to a few days to run.

    Rscript ../run_pipeline.R \
       ../../data/case_study_2/classes_diablo.tsv \
       --classes_secondary NA \
       --dropna_classes FALSE \
       --dropna_prop 0 \
       --data \
        ../../../data/case_study_2/data_lipidomics.tsv \
        ../../data/case_study_2/data_metabolomics.tsv \
        ../../data/case_study_2/data_proteomics.tsv \
        ../../data/case_study_2/data_transcriptomics_imputed_all.tsv \
       --data_names lipidome metabolome proteome transcriptome \
       --force_unique FALSE \
       --mappings NA \
       --ncpus 16 \
       --diablocomp 0 \
       --linkage 0.1 \
       --diablo_keepx 5 6 7 8 9 10 30 \
       --icomp 0 \
       --pcomp 10 \
       --plsdacomp 2 \
       --splsdacomp 2 \
       --splsda_keepx 5 6 7 8 9 10 30 \
       --dist_plsda centroids.dist \
       --dist_splsda centroids.dist \
       --dist_diablo centroids.dist \
       --contrib max \
       --outfile_dir ../../results/case_study_2/EXAMPLE/ \
       --rdata RData.RData \
       --plot Rplots.pdf \
       --args Rscript.sh

# 3 Input data

## 3.1 Biological context

## 3.2 Summary

The data used as input to this pipeline available in gitlab:

-   [biological class information]()
-   [lipidome]()
-   [metabolome]()
-   [proteome]()
-   [transcriptome]()

## 3.3 Metadata

The `classes_diablo.tsv` sample information file is loaded as a vector:

    > classes
      [1] "More severe" "More severe" "Less severe" "Less severe" "More severe"
      [6] "More severe" "Less severe" "More severe" "Less severe" "More severe"
     [11] "More severe" "Less severe" "Less severe" "Less severe" "Less severe"
     [16] "Less severe" "Less severe" "Less severe" "More severe" "More severe"
     [21] "More severe" "More severe" "Less severe" "More severe" "More severe"
     [26] "More severe" "More severe" "More severe" "More severe" "More severe"
     [31] "Less severe" "Less severe" "Less severe" "More severe" "Less severe"
     [36] "More severe" "More severe" "Less severe" "Less severe" "More severe"
     [41] "Less severe" "More severe" "More severe" "More severe" "Less severe"
     [46] "More severe" "More severe" "More severe" "Less severe" "More severe"
     [51] "Less severe" "More severe" "Less severe" "More severe" "Less severe"
     [56] "Less severe" "Less severe" "Less severe" "Less severe" "Less severe"
     [61] "Less severe" "More severe" "More severe" "More severe" "Less severe"
     [66] "More severe" "Less severe" "More severe" "Less severe" "Less severe"
     [71] "Less severe" "More severe" "Less severe" "Less severe" "Less severe"
     [76] "Less severe" "More severe" "More severe" "Less severe" "Less severe"
     [81] "More severe" "More severe" "Less severe" "More severe" "More severe"
     [86] "More severe" "Less severe" "Less severe" "Less severe" "More severe"
     [91] "More severe" "Less severe" "More severe" "More severe" "Less severe"
     [96] "Less severe" "Less severe" "More severe" "More severe" "More severe"

No repeated measurements were known to be carried out.

## 3.4 Data

The lipidome, metabolome, proteome and transcriptome data have 100
matched samples and an arbitrary number of features.

# 4 Input data

For reference, you can also access all these loaded data below and
results in the `RData` object. Data objects in both this walkthrough and
the `RData` object have identical names. The code blocks in this
walkthrough will reproduce the same data structures. Not all values may
be identical since some functions may be non-deterministic. This is true
even if a seed is specified, since `R >=3.6` these are not reproducible
across machines.

<details>
<summary>
Click to expand code block
</summary>

``` r
  url_rdata <- "https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/results/case_study_2/RData.RData"
  download.file(url_rdata, "RData.RData")
  load("RData.RData")
  ls()
  #>  [1] "argv"                "classes"             "data"               
  #>  [4] "data_imp"            "data_pca_multilevel" "data_plsda"         
  #>  [7] "data_splsda"         "diablo"              "dist_diablo"        
  #> [10] "dist_plsda"          "dist_splsda"         "linkage"            
  #> [13] "mappings"            "pca_impute"          "pca_withna"         
  #> [16] "pch"                 "perf_diablo"         "tuned_diablo"       
  #> [19] "tuned_splsda"        "url_rdata"
```

> **NOTE**: There are some differences in the R data object for case
> studies 1 and 2. Case study 1 was originally performed using an early
> version of the code. It has since been streamlined for case study 2,
> which contains more information. However, both analyses are
> reproducible and the user can if needed further investigate the
> internal data structures after loading them.

</details>

For case study 2, it is not recommended to run the pipeline directly
because of resource usage. This may take more than a day on 16 cpus.
Instead, we use the `RData` object created as a result and load data
from there directly to skip long-running steps. In those cases, we
provide the code for reference only.

Obtain the input data from the git repository.

<details>
<summary>
Click to expand code block
</summary>

``` r
  url_class <- "https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/data/case_study_2/classes_diablo.tsv"
  url_lipi <- "https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/data/case_study_2/data_lipidome.tsv"
  url_meta <- "https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/data/case_study_2/data_metabolome.tsv"
  url_prot <- "https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/data/case_study_2/data_proteome.tsv"
  url_tran <- "https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/data/case_study_2/data_transcriptome.tsv"
  
  urls <- c(url_class, url_lipi, url_meta, url_prot, url_tran)
  file_names <- sapply(strsplit(urls, "/"), tail, 1)
  mapply(function(x, y) download.file(x, y), urls, file_names, SIMPLIFY=FALSE)
  if (any(file.exists(file_names)) != TRUE) {stop("Files incorrectly downloaded!")}
```

</details>

Load the data and metadata with the following functions. All files must
be in the same order!

<details>
<summary>
Click to expand code block
</summary>

Using the pipeline:

``` r
  paths <- c(
    "data_lipidomics.tsv", "data_metabolomics.tsv",
    "data_proteomics.tsv", "data_transcriptomics.tsv"
  )
  classes <- parse_classes("classes_diablo.tsv")
  data <- lapply(paths, parse_data, missing_as=NA, rmna=TRUE)
  # these must be unique
  data_names <- c("lipidome", "metabolome", "proteome", "translatome")
  # design matrix specifying linkage across blocks (for multiomics)
  design <- create_design(data, 0.1)
```

Using the R data object:

``` r
  load("RData.RData")
  classes
  #>   [1] "More severe" "More severe" "Less severe" "Less severe" "More severe"
  #>   [6] "More severe" "Less severe" "More severe" "Less severe" "More severe"
  #>  [11] "More severe" "Less severe" "Less severe" "Less severe" "Less severe"
  #>  [16] "Less severe" "Less severe" "Less severe" "More severe" "More severe"
  #>  [21] "More severe" "More severe" "Less severe" "More severe" "More severe"
  #>  [26] "More severe" "More severe" "More severe" "More severe" "More severe"
  #>  [31] "Less severe" "Less severe" "Less severe" "More severe" "Less severe"
  #>  [36] "More severe" "More severe" "Less severe" "Less severe" "More severe"
  #>  [41] "Less severe" "More severe" "More severe" "More severe" "Less severe"
  #>  [46] "More severe" "More severe" "More severe" "Less severe" "More severe"
  #>  [51] "Less severe" "More severe" "Less severe" "More severe" "Less severe"
  #>  [56] "Less severe" "Less severe" "Less severe" "Less severe" "Less severe"
  #>  [61] "Less severe" "More severe" "More severe" "More severe" "Less severe"
  #>  [66] "More severe" "Less severe" "More severe" "Less severe" "Less severe"
  #>  [71] "Less severe" "More severe" "Less severe" "Less severe" "Less severe"
  #>  [76] "Less severe" "More severe" "More severe" "Less severe" "Less severe"
  #>  [81] "More severe" "More severe" "Less severe" "More severe" "More severe"
  #>  [86] "More severe" "Less severe" "Less severe" "Less severe" "More severe"
  #>  [91] "More severe" "Less severe" "More severe" "More severe" "Less severe"
  #>  [96] "Less severe" "Less severe" "More severe" "More severe" "More severe"
  lapply(data, dim)
  #> $lipidome
  #> [1]  100 3357
  #> 
  #> $metabolome
  #> [1] 100 150
  #> 
  #> $proteome
  #> [1] 100 517
  #> 
  #> $transcriptome
  #> [1]   100 13263
  data_names <- argv$data_names
  design <- create_design(data, 0.1)
```

</details>

# 5 Data quality control

Before analysing data further we perform standard checks for common
quality issues.

## 5.1 Accounting for missing values

We discovered a small proportion of missing values within the
transcriptome data. Missing values in data can bias analyses. There are
multiple ways to address this. In our case we use imputation.

<details>
<summary>
Click to expand code block
</summary>

``` r
  # the data object has already been cleaned, we load the raw data for reference
  url_tran <- "https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/data/case_study_2/data_transcriptomics.tsv"
  download.file(url_tran, "data_transcriptomics.tsv")
  
  unimputed_tran_path <- "data_transcriptomics.tsv"
  unimputed_tran <- read.table(unimputed_tran_path, sep="\t", header=TRUE, row.names=1)
  unimputed_tran[unimputed_tran == 0] <- NA
  
  na_prop_tran <- show_na_prop(unimputed_tran, "Transcriptome")
```

![](figs_2/check_na-1.png)<!-- -->
</details>
<details>
<summary>
Click to expand code block
</summary>

Using the pipeline:

``` r
  data <- lapply(data, remove_na_class, classes)
```

</details>

We corrected for the missing values in the transcruptome data (\~0.5% of
original data) by imputation with the NIPALS algorithm, effective on
data with &lt; 20% missing values. We considered that the proportion of
missing values are low and imputation would be effective. Note that
imputation can take some time, increasing with data size and component
count.

<details>
<summary>
Click to expand code block
</summary>

Using the pipeline:

``` r
  # this step is important, some functions use the names internally
  names(data) <- data_names
  data_imp <- impute_missing(data, rep(10, length(data)), outdir="./")
```

Using the R data object:

``` r
  # this data was pre-imputed and then loaded back in internally
  dim(data$transcriptome)
  #> [1]   100 13263
```

</details>

> **NOTE**: This step also writes out the imputed, un-replaced data.
> Since imputation can take some time, this is mainly for convenience.
> You can load the imputed data directly as input if it meets your
> requirements.

To test that imputation has not introduced significant technical
variation into the data, we observe the correlation between variates of
the principal components.

> **NOTE**: All PCAs are centered and scaled.

<details>
<summary>
Click to expand code block
</summary>

Using the pipeline (if data is imputed internally):

``` r
  pca_imputed <- mixOmics::pca(data$transcriptome, ncomp=10)
  pca_unimputed <- mixOmics::pca(unimputed_tran, ncomp=10)
  cor_imputed_unimputed_(pca_imputed, pca_unimputed, "Transcriptome")
  #> [1] "Plotting correlation between unimputed and imputed components"
```

![](figs_2/compare_na-1.png)<!-- -->

Using the R Data object:

To improve efficiency during code run, a imputed data file was generated
as output in the first run of the pipeline. This was then used in
subsequent runs of the pipeline as direct input. You can reproduce this
manually with the code in the `make_manuscript_figures.R` [script in our
gitlab repository](https://gitlab.com/tyagilab/sars-cov-2/).
</details>

In this case, there is a strong correlation between the variates on at
least the first 5 principal components corresponding to at least 50% of
the variation in the data.

## 5.2 Accounting for unwanted variation

The experimental design of this study contains no repeated measurements
on the same sample. There is also no known variation from batch effects.

# 6 Single omics analyses

We next apply the PLSDA (Partial Least Squares Discriminant Analysis)
and sPLSDA (sparse variant of PLSDA) method for each block of
single-omics data, and as before internally perform a multilevel
decomposition to account for the repeated measurements within each cell
culture.

## 6.1 Parameter tuning

To investigate the parameters best suited for the methods, leave-one-out
cross validation was performed. The number of components and features
selected were tuned internally with a function in the mixOmics package.

<details>
<summary>
Click to expand code block
</summary>

Using the pipeline:

``` r
  # this step can take some time
  tuned_splsda <- tune_splsda(
    data, classes, data_names, data.frame(pch), ncomp=4, nrepeat=10,
    logratio="none", test_keepX=c(10, 25, 50, 100), validation="loo", folds=10,
    dist="centroids.dist", cpus=6, progressBar=FALSE
  )
  splsda_keepx <- lapply(tuned_splsda, `[`, "choice.keepX")
  splsda_ncomp <- lapply(tuned_splsda, `[`, "choice.ncomp")
  
  print("Tuned splsda to use number of components:")
  splsda_ncomp <- lapply(splsda_ncomp, `[`, "ncomp")
  splsda_ncomp <- unlist(splsda_ncomp, recursive=FALSE)
  names(splsda_ncomp) <- data_names
  print(splsda_ncomp)
  
  print("Tuned the number of variables selected on each component to:")
  print(splsda_keepx)
  splsda_keepx <- unlist(splsda_keepx, recursive=FALSE)
  names(splsda_keepx) <- data_names
  print(splsda_keepx)
```

Using the R data object:

``` r
  names(tuned_splsda)
  #> [1] "lipidome"      "metabolome"    "proteome"      "transcriptome"
  
  # keep optimal number of features to keep
  splsda_keepx <- lapply(tuned_splsda, `[`, "choice.keepX")
  splsda_keepx <- unlist(splsda_keepx, recursive=FALSE)
  names(splsda_keepx) <- data_names
```

</details>

## 6.2 Running the analysis

With the tuned parameters, we run sPLSDA (subset of features).

<details>
<summary>
Click to expand code block
</summary>

Using the pipeline:

``` r
  # this step can take some time
  data_splsda <- classify_splsda(
    data=data, classes=classes, pch=pch, title=data_names,
    ncomp=2, keepX=splsda_keepx, contrib="max", outdir="./",
    mappings=NULL, dist="centroids.dist", bg=TRUE
  )
```

Using the R data object:

``` r
  lapply(data_splsda, names)
  #> $lipidome
  #> [1] "data_splsda" "perf_splsda"
  #> 
  #> $metabolome
  #> [1] "data_splsda" "perf_splsda"
  #> 
  #> $proteome
  #> [1] "data_splsda" "perf_splsda"
  #> 
  #> $transcriptome
  #> [1] "data_splsda" "perf_splsda"
```

</details>

We also run PLSDA (all features) for comparison.

<details>
<summary>
Click to expand code block
</summary>

Using the pipeline:

``` r
  # this step can take some time
  data_plsda <- classify_plsda(
    data=data, classes=classes, pch=pch, title=data_names,
    ncomp=2, contrib="max", outdir="./",
    mappings=NULL, dist="centroids.dist", bg=TRUE
  )
```

Using the R data object:

``` r
  lapply(data_plsda, names)
  #> $lipidome
  #> [1] "data_plsda" "perf_plsda"
  #> 
  #> $metabolome
  #> [1] "data_plsda" "perf_plsda"
  #> 
  #> $proteome
  #> [1] "data_plsda" "perf_plsda"
  #> 
  #> $transcriptome
  #> [1] "data_plsda" "perf_plsda"
```

</details>

These automatically generate a large series of plots. Figures are
ordered sequentially by each block of omics data in `data_name`, in this
case `lipidome`, `metabolome`, `proteome`, followed by `translatome`.
Some `txt` files containing feature loadings are also written to the
output directory.

## 6.3 Explanation of output

Detailed technical information on the individual figures types and
underlying methods are available at the [mixOmics
website](http://mixomics.org/). This walkthrough will explain the plots
in context of the biological system under study only. The plots are
explained in the order that they appear in this document.

### 6.3.1 Scatter plots

Plotting the first few components of the sPLSDA reveals a spectrum of
phenotypes from less severe to more severe, with a degree of overlap.

### 6.3.2 Classification error

> **NOTE**: Balanced error rate can be used in cases where classes are
> imbalanced.

To show accuracy of the method, the classification error rate across
multiple components for maximum, centroids and mahalanobis distance are
plotted. The centroids distance metric appeared to perform the best in
this case, with the lowest error rate after 3 components for proteome
and 4 for translatome.

### 6.3.3 Feature stability

To assess how stable the feature is across cross-validation, each of the
horizontal bars represent a feature. The height of the bar corresponds
to stability. We observe three subsets of highly stable features,
moderately stable features and less stable features.

### 6.3.4 ROC curves

As a supplementary layer of validation, ROC curves showing
classification accuracy are available, but we note that these have
limited applicability given the specific context of the method. The
underlying method already internally specifies the prediction cutoff to
achieve maximal sensitivity and specificity.

### 6.3.5 Arrow plots

To review agreement between matching data sets, we use an arrow plot.
Short arrows suggest a strong similarity, while long arrows indicate
dissimilarity. Two sample groups in the proteome data (Virus 6h, Virus
24h) appear to be dissimilar.

### 6.3.6 Clustered image maps

To investigate the relationship between samples, these plots can be
interpreted like heatmaps. They are coloured by biological class as well
as batch information. Individual components can be extracted for custom
visualisation if needed.

### 6.3.7 Variable loadings

To understand how much a feature contributes to the biological
classification, loading weights are shown. Direction of the bar
indicates the abundance of that feature in the data (left for less,
right for more). Bars are colour coded to biological class. Individual
components can be extracted for custom visualisation if needed.

Note that only a subset of these features are visualised. The full list
is exported into a tab-separated `txt` file, similar to that from a
`limma` differential expression analysis.

### 6.3.8 PLSDA

To supplement the sPLSDA, we also compared the performance of PLSDA (a
non sparse variant of sPLSDA keeping all features). This showed
similarities in patterns across the datasets. The plots are conceptually
identical, but several plots are excluded as they are not applicable
(feature selection, variable stability).

# 7 Multi omics analyses

Having assessed the major sources of variation and features of interest
contributing to biological conditions within the individual blocks of
omics data, we can use this information to guide our multi-omics
integration.

We applied a latent variable approach to identify a highly correlated
multi-omics signature. This analysis is carried out in a conceptually
similar way to the previous sPLSDA with similar parameter requirements,
except with multiple omics data blocks corrected for longitudinal study
effects specified as input. We illustrate the correlation between
features across these omics blocks with a circos plot.

## 7.1 Parameter tuning

To investigate the parameters best suited for the methods, leave-one-out
cross validation was performed. Similar to sPLSDA, the number of
components and features selected were tuned internally with a function
in the mixOmics package.

<details>
<summary>
Click to expand code block
</summary>

Using the pipeline:

``` r
  # tune number of components
  tuned_diablo <- tune_diablo_ncomp(data, classes, design, ncomp=2, cpus=6)
  print("Parameters with lowest error rate:")
  tuned_diablo <- tuned_diablo$choice.ncomp$WeightedVote["Overall.BER",]
  diablo_ncomp <- tuned_diablo[which.max(tuned_diablo)]
  print("Number of components:")
  print(diablo_ncomp)
  
  # tune keepx
  diablo_keepx <- c(5,10,12,14,16,18,20,30)
  diablo_keepx <- tune_diablo_keepx(
    data, classes, diablo_ncomp, design, diablo_keepx, cpus=6,
    dist="mahalanobis.dist", progressBar=FALSE
  )
  print("Diablo keepx:")
  print(diablo_keepx)
```

Using the R object:

``` r
  tuned_diablo
  #>         max.dist   centroids.dist mahalanobis.dist 
  #>                2                2                2
  diablo$keepX
  #> $lipidome
  #> [1] 7 7
  #> 
  #> $metabolome
  #> [1]  6 10
  #> 
  #> $proteome
  #> [1] 8 5
  #> 
  #> $transcriptome
  #> [1]  7 30
```

</details>

## 7.2 Running the analysis

With the tuned parameters, we run multi-block sPLSDA (DIABLO).

<details>
<summary>
Click to expand code block
</summary>

Using the pipeline:

``` r
  diablo_ncomp <- tuned_diablo[which.max(tuned_diablo)]
  # these values were precomputed, newer versions of code save these values for use
  diablo_keepx <- list(lipidome=c(7, 7)
                       metabolome=c(6, 10)
                       proteome=c(8, 5),
                       translatome=c(7, 30))
  
  # run the algorithm
  diablo <- run_diablo(data, classes, diablo_ncomp, design, diablo_keepx)
```

Using the R data object:

``` r
  diablo
  #> 
  #> Call:
  #>  block.splsda(X = data, Y = classes, ncomp = ncomp, keepX = keepx, design = design) 
  #> 
  #>  sGCCA with 2 components on block 1 named lipidome 
  #>  sGCCA with 2 components on block 2 named metabolome 
  #>  sGCCA with 2 components on block 3 named proteome 
  #>  sGCCA with 2 components on block 4 named transcriptome 
  #>  sGCCA with 2 components on the outcome Y
  #> 
  #>  Dimension of block 1 is  100 3357 
  #>  Dimension of block 2 is  100 150 
  #>  Dimension of block 3 is  100 517 
  #>  Dimension of block 4 is  100 13263 
  #>  Outcome Y has 2 levels 
  #> 
  #>  Selection of 7 7 variables on each of the sGCCA components on the block 1 
  #>  Selection of 6 10 variables on each of the sGCCA components on the block 2 
  #>  Selection of 8 5 variables on each of the sGCCA components on the block 3 
  #>  Selection of 7 30 variables on each of the sGCCA components on the block 4 
  #> 
  #>  Main numerical outputs: 
  #>  -------------------- 
  #>  loading vectors: see object$loadings 
  #>  variates: see object$variates 
  #>  variable names: see object$names 
  #> 
  #>  Functions to visualise samples: 
  #>  -------------------- 
  #>  plotIndiv, plotArrow, cimDiablo, plotDiablo 
  #> 
  #>  Functions to visualise variables: 
  #>  -------------------- 
  #>  plotVar, plotLoadings, network, circosPlot 
  #> 
  #>  Other functions: 
  #>  -------------------- 
  #>  selectVar, perf, auc
```

</details>

## 7.3 Explanation of output

Many of these plots can be interpreted in conceptually similar ways to
that of the single-omic sPLSDA above. However, some extra plots are
created to better illustrate correlation between datasets.

> **NOTE**: As usual, the full list of plots and results are available
> in the git repository.

### 7.3.1 Classification error

> **NOTE**: Balanced error rate can be used in cases where classes are
> imbalanced.

To show accuracy of the method, the classification error rate across
multiple components for maximum, centroids and mahalanobis distance are
plotted. The centroids distance metric appeared to perform the best in
this case, with the lowest error rate after 3 components for proteome
and 4 for translatome.

<details>
<summary>
Click to expand code block
</summary>

``` r
  # assess performance
  perf_diablo <- mixOmics::perf(
    diablo, validation="loo", nrepeats=10, auc=TRUE, cpus=6, progressBar=TRUE
  )
  #> 
  #> Performing repeated cross-validation...
  #>   |                                                                              |                                                                      |   0%
  plot(perf_diablo)
```

![](figs_2/perf_diablo-1.png)<!-- -->
</details>

### 7.3.2 Feature stability

To assess how stable the feature is across cross-validation, each of the
horizontal bars represent a feature. The height of the bar corresponds
to stability. We observe three subsets of highly stable features,
moderately stable features and less stable features.

<details>
<summary>
Click to expand code block
</summary>

``` r
  
  mapply(
    function(x, y) plot(
      x$comp1, type="h", main="Comp 1", las=2,
      ylab="Stability", xlab="Features", xaxt="n"),
    perf_diablo$features$stable$nrep1,
    c("Lipidome", "Metabolome", "Proteome", "Transcriptome")
  )
```

![](figs_2/perf_diablo_stability-1.png)<!-- -->![](figs_2/perf_diablo_stability-2.png)<!-- -->![](figs_2/perf_diablo_stability-3.png)<!-- -->![](figs_2/perf_diablo_stability-4.png)<!-- -->

      #> $lipidome
      #> NULL
      #> 
      #> $metabolome
      #> NULL
      #> 
      #> $proteome
      #> NULL
      #> 
      #> $transcriptome
      #> NULL

</details>

### 7.3.3 Correlation plot

A correlation score is provided per block of omics data. In this case
there is just one score as there are two blocks of omics data.

<details>
<summary>
Click to expand code block
</summary>

``` r
  # to keep the case study concise we use a custom function to output main plots
  mixOmics::plotDiablo(diablo, ncomp=1)
```

![](figs_2/diablo-1.png)<!-- -->
</details>

### 7.3.4 Scatter plots

Plotting the first few components of the multi-block sPLSDA (DIABLO)
reveals a spectrum of phenotypes from less severe to more severe, with a
degree of overlap.

<details>
<summary>
Click to expand code block
</summary>

``` r
  # to keep the case study concise we use a custom function to output main plots
  mixOmics::plotIndiv(
    diablo, ind.names=FALSE, legend=TRUE, title='DIABLO', ellipse=TRUE
  )
```

![](figs_2/perf_scatter-1.png)<!-- -->
</details>

### 7.3.5 ROC curves

As a supplementary layer of validation, ROC curves showing
classification accuracy are available, but we note that these have
limited applicability given the specific context of the method. The
underlying method already internally specifies the prediction cutoff to
achieve maximal sensitivity and specificity.

> **NOTE**: No ROC curves were generated for this particular case study.

### 7.3.6 Arrow plots

To review agreement between matching data sets, we use an arrow plot.
Short arrows suggest a strong similarity, while long arrows indicate
dissimilarity. There appears to be a spectrum of phenotypes from less
severe to more severe, with a degree of overlap.

<details>
<summary>
Click to expand code block
</summary>

``` r
  # to keep the case study concise we use a custom function to output main plots
  mixOmics::plotArrow(diablo, ind.names=FALSE, legend=TRUE, title='DIABLO')
```

![](figs_2/perf_arrow-1.png)<!-- -->
</details>

### 7.3.7 Correlation circle plots

The correlation circle plots highlight the contribution of each variable
to each component. A strong correlation between variables is indicated
by clusters of points.

<details>
<summary>
Click to expand code block
</summary>

``` r
  # to keep the case study concise we use a custom function to output main plots
  mixOmics::plotVar(diablo, style='graphics', legend=TRUE, comp=c(1,2),
    title="DIABLO 1/2", var.names=FALSE
  )
```

![](figs_2/perf_correlation-1.png)<!-- -->
</details>

### 7.3.8 Clustered image maps

To investigate the relationship between samples, these plots can be
interpreted like heatmaps. They are coloured by biological class as well
as batch information. Individual components can be extracted for custom
visualisation if needed.

<details>
<summary>
Click to expand code block
</summary>

``` r
  # to keep the case study concise we use a custom function to output main plots
  mixOmics::cimDiablo(diablo, size.legend=0.5, col.names=FALSE)
```

![](figs_2/cim-1.png)<!-- -->
</details>

### 7.3.9 Variable loadings

To understand how much a feature contributes to the biological
classification, loading weights are shown. Direction of the bar
indicates the abundance of that feature in the data (left for less,
right for more). Bars are colour coded to biological class. Individual
components can be extracted for custom visualisation if needed.

Note that only a subset of these features are visualised. The full list
is exported into a tab-separated `txt` file, similar to that from a
`limma` differential expression analysis.

<details>
<summary>
Click to expand code block
</summary>

``` r
  # to keep the case study concise we use a custom function to output main plots
  mixOmics::plotLoadings(diablo, contrib="max", comp=1, max.name.length=6,
    method='median', ndisplay=20, size.name=0.6,
    size.legend=0.6, title=paste("DIABLO max loadings")
  )
```

![](figs_2/loadings-1.png)<!-- -->
</details>

### 7.3.10 Circos plot

To visualise correlations between different blocks of omics data, a
circos plot is generated. The blue block represents lipidome data, the
green block represents metabolome data, the red block represents
proteome data and the orange block represents transcriptome data. Each
point on the circle is a single feature. Lines linking features show
correlations between features that pass the user-specified correlation
threshold, in this case 0.8. Red lines indicate positive correlation and
blue lines negative correlation.

<details>
<summary>
Click to expand code block
</summary>

``` r
  # to keep the case study concise we use a custom function to output main plots
  mixOmics::circosPlot(
    diablo, cutoff=0.8, line=FALSE, size.legend=0.5, size.variables=0.001
  )
```

![](figs_2/circos-1.png)<!-- -->
</details>

### 7.3.11 Network plot

Similar to the circos plot, a network plot can be generated. It is
compatible with cytoscape and can be exported as a `gml` file.

<details>
<summary>
Click to expand code block
</summary>

``` r
  # to keep the case study concise we use a custom function to output main plots
  mixOmics::network(
    diablo, blocks=c(1,2), color.node=c('darkorchid','lightgreen'), cutoff=0.4,
    col.names=FALSE, row.names=FALSE
  )
```

![](figs_2/network-1.png)<!-- -->
</details>

# 8 Output data

Files output by the pipeline include:

-   a `pdf` file of all plots generated by the pipeline
-   tab-separated `txt` files containing feature contribution weights to
    each biological class
-   tab-separated `txt` file containing correlations between each omics
    data block

[A `RData` object with all input and output is available in the git
repository.](https://gitlab.com/tyagilab/sars-cov-2/-/blob/master/results/case_study_2/RData.RData)
This is not included directly in the `multiomics` package because of
size constraints, and includes data from four omics datasets.

# 9 Acknowledgements

[Please refer to introduction.](introduction.html)

# 10 References

[Please refer to introduction.](introduction.html)
