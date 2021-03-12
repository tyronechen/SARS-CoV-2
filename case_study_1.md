Case study 1
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

-   [Introduction](introduction.md)
-   [Case study 1](case_study_1.md)
-   [Case study 2](case_study_2.md)

# 2 Running the script

Load the library.

``` r
library(multiomics)
```

A script to reproduce our analysis for case study 1 is shown. [You can
also download this
here](https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/src/case_study_1/example.sh).
This may take up to a few hours to run.

    Rscript ../run_pipeline.R \
       --classes ../../data/case_study_1/classes_diablo.txt  \
       --classes_secondary ../../data/case_study_1/pch.txt \
       --dropna_classes TRUE \
       --dropna_prop 0 \
       --data ../../data/case_study_1/diablo_proteome.txt ../../data/case_study_1/diablo_translatome.txt \
       --data_names proteome translatome \
       --mappings ../../data/case_study_1/proteome_mapfile.txt ../../data/case_study_1/translatome_mapfile.txt \
       --ncpus 6 \
       --diablocomp 0 \
       --diablo_keepx 5 10 12 14 16 18 20 30 \
       --icomp 24 \
       --pcomp 10 \
       --plsdacomp 4 \
       --splsdacomp 4 \
       --splsda_keepx 10 25 50 100 \
       --dist_plsda centroids.dist \
       --dist_splsda centroids.dist \
       --dist_diablo mahalanobis.dist \
       --contrib max \
       --outfile_dir ../../results/case_study_1/EXAMPLE \
       --rdata RData.RData \
       --plot Rplots.pdf \
       --args Rscript.sh

Alternatively, pass the command line arguments as a json file. **This will override any arguments specified on the command line.** [An example file is available here](https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/src/case_study_2/args.json).

    Rscript ../run_pipeline.R --json args.json


# 3 Input data

## 3.1 Biological context

This dataset contains two omics data: proteome and translatome. There
are three biological replicates for 8 sample types. Cell cultures
(uninfected/infected with SARS-CoV-2) were resampled over four
timepoints. The original publication with the source data is here:

-   Bojkova, D., Klann, K., Koch, B. et al. Proteomics of
    SARS-CoV-2-infected host cells reveals therapy targets. *Nature*
    **583,** 469–472 (2020). <https://doi.org/10.1038/s41586-020-2332-7>

## 3.2 Summary

The data used as input to this pipeline available in gitlab:

-   [biological class
    information](https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/data/case_study_1/classes_diablo.txt)
-   [longitudinal measurement
    information](https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/data/case_study_1/pch.txt)
-   [proteome
    mappings](https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/data/case_study_1/proteome_mapfile.txt)
-   [translatome
    mappings](https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/data/case_study_1/translatome_mapfile.txt)
-   [proteome](https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/data/case_study_1/diablo_proteome.txt)
-   [translatome](https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/data/case_study_1/diablo_translatome.txt)

## 3.3 Metadata

The `classes_diablo.txt` sample information file is loaded as a vector:

    > classes
     [1] "Control_2h"  "Control_2h"  "Control_2h"  "Control_6h"  "Control_6h"
     [6] "Control_6h"  "Control_10h" "Control_10h" "Control_10h" "Control_24h"
    [11] "Control_24h" "Control_24h" "Virus_2h"    "Virus_2h"    "Virus_2h"
    [16] "Virus_6h"    "Virus_6h"    "Virus_6h"    "Virus_10h"   "Virus_10h"
    [21] "Virus_10h"   "Virus_24h"   "Virus_24h"   "Virus_24h"

The `pch.txt` repeated measurements information file is loaded as a
vector:

    > pch
     [1] 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3

## 3.4 Data

The proteome and translatome data have 24 matched samples and an
arbitrary number of features.

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
  url_rdata <- "https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/results/case_study_1/RData.RData"
  download.file(url_rdata, "RData.RData")
  load("RData.RData")
  ls()
  #>  [1] "argv"                "classes"             "data"               
  #>  [4] "data_imp"            "data_pca_multilevel" "data_plsda"         
  #>  [7] "data_splsda"         "diablo"              "diablo_all"         
  #> [10] "dist_diablo"         "dist_splsda"         "mappings"           
  #> [13] "pca_impute"          "pca_withna"          "pch"                
  #> [16] "tuned_diablo"        "tuned_splsda"        "url_rdata"
```

> **NOTE**: There are some differences in the R data object for case
> studies 1 and 2. Case study 1 was originally performed using an early
> version of the code. It has since been streamlined for case study 2,
> which contains more information. However, both analyses are
> reproducible and the user can if needed further investigate the
> internal data structures after loading them.

</details>

For case study 1, you can run the pipeline directly but it will take
some time (up to a few hours to run from end to end). In this
walkthrough, we use the `RData` object created as a result and load data
from there directly to skip long-running steps.

Obtain the input data from the git repository.

<details>
<summary>
Click to expand code block
</summary>

``` r
  url_class <- "https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/data/case_study_1/classes_diablo.txt"
  url_pch <- "https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/data/case_study_1/pch.txt"
  url_prot_map <- "https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/data/case_study_1/proteome_mapfile.txt"
  url_tran_map <- "https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/data/case_study_1/translatome_mapfile.txt"
  url_prot <- "https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/data/case_study_1/diablo_proteome.txt"
  url_tran <- "https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/data/case_study_1/diablo_translatome.txt"

  urls <- c(url_class, url_pch, url_prot_map, url_tran_map, url_prot, url_tran)
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
  paths <- c("diablo_proteome.txt", "diablo_translatome.txt")
  mappings <- c("proteome_mapfile.txt", "translatome_mapfile.txt")
  classes <- parse_classes("classes_diablo.txt")
  pch <- parse_classes("pch.txt")
  data <- lapply(paths, parse_data, missing_as=NA, rmna=TRUE)
  mapped <- lapply(mappings, parse_mappings)
  # these must be unique
  data_names <- c("proteome", "translatome")
  # design matrix specifying linkage across blocks (for multiomics)
  design <- create_design(data, 0.1)
```

Using the R data object:

``` r
  load("RData.RData")
  classes
  #>  [1] "Control_2h"  "Control_2h"  "Control_2h"  "Control_6h"  "Control_6h"
  #>  [6] "Control_6h"  "Control_10h" "Control_10h" "Control_10h" "Control_24h"
  #> [11] "Control_24h" "Control_24h" "Virus_2h"    "Virus_2h"    "Virus_2h"   
  #> [16] "Virus_6h"    "Virus_6h"    "Virus_6h"    "Virus_10h"   "Virus_10h"  
  #> [21] "Virus_10h"   "Virus_24h"   "Virus_24h"   "Virus_24h"
  lapply(data, dim)
  #> $proteome
  #> [1]   24 6380
  #>
  #> $translatome
  #> [1]   24 1595
  data_names <- c("proteome", "translatome")
  design <- create_design(data, 0.1)
  pch
  #>  [1] 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3
```

</details>

# 5 Data quality control

Before analysing data further we perform standard checks for common
quality issues.

## 5.1 Accounting for missing values

We discovered a high proportion of missing values within the translatome
data. Missing values in data can bias analyses. There are multiple ways
to address this. In our case we use a mixture of filtering and
imputation.

<details>
<summary>
Click to expand code block
</summary>

``` r
  # the data object has already been cleaned, we load the raw data for reference
  unimputed_prot_path <- "diablo_proteome.txt"
  unimputed_prot <- read.table(unimputed_prot_path, sep="\t", header=TRUE, row.names=1)
  unimputed_tran_path <- "diablo_translatome.txt"
  unimputed_tran <- read.table(unimputed_tran_path, sep="\t", header=TRUE, row.names=1)
  unimputed_prot[unimputed_prot == 0] <- NA
  unimputed_tran[unimputed_tran == 0] <- NA
  na_prop_prot <- show_na_prop(unimputed_prot, "Proteome")
```

![](figs_1/check_na-1.png)<!-- -->

``` r
  na_prop_tran <- show_na_prop(unimputed_tran, "Translatome")
```

![](figs_1/check_na-2.png)<!-- -->

``` r

  # need to drop col where all missing
  unimputed_prot <- remove_na_class(unimputed_prot, classes, missing_as=NA)
  #> [1]   24 6381
  #> [1] "Dropping features where at least one class is NA"
  unimputed_tran <- remove_na_class(unimputed_tran, classes, missing_as=NA)
  #> [1]   24 2712
  #> [1] "Dropping features where at least one class is NA"
```

</details>

We corrected for the missing values in the translatome data (\~47% of
original data) by a mixture of filtering and imputation. We considered
that filtering alone would be too aggressive and imputation alone would
be ineffective. Filtering was performed by dropping all features which
were not represented across each biological sample group.

<details>
<summary>
Click to expand code block
</summary>

``` r
  data <- lapply(data, remove_na_class, classes)
  #> [1]   24 6380
  #> [1] "Dropping features where at least one class is NA"
  #> [1]   24 1595
  #> [1] "Dropping features where at least one class is NA"
```

</details>

This reduced the quantity of missing values to \~17% of the original
data. An imputation was performed with the NIPALS algorithm, which is
effective on data with &lt; 20% missing values. Note that imputation can
take some time, increasing with data size and component count.

<details>
<summary>
Click to expand code block
</summary>

Using the pipeline:

``` r
  # this step is important, some functions use the names internally
  names(data) <- data_names
  data_imp <- impute_missing(data, rep(24, length(data)), outdir="./")

  # only replace the missing values and preserve the original values
  data <- replace_missing(data, data_imp)
```

Using the R data object:

``` r
  lapply(data_imp, dim)
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

``` r
  pca_unimputed_prot <- mixOmics::pca(unimputed_prot, ncomp=10)
  pca_unimputed_tran <- mixOmics::pca(unimputed_tran, ncomp=10)
  pca_imputed_prot <- mixOmics::pca(data$proteome, ncomp=10)
  pca_imputed_tran <- mixOmics::pca(data$translatome, ncomp=10)
  cor_imputed_unimputed_(pca_imputed_prot, pca_unimputed_prot, "Proteome")
  #> [1] "Plotting correlation between unimputed and imputed components"
```

![](figs_1/compare_na-1.png)<!-- -->

``` r
  cor_imputed_unimputed_(pca_imputed_tran, pca_unimputed_tran, "Translatome")
  #> [1] "Plotting correlation between unimputed and imputed components"
```

![](figs_1/compare_na-2.png)<!-- --> &gt; **NOTE**: Imputation may take
some time, go eat lunch

</details>

In both cases, there is a strong correlation between the variates on at
least the first 5 principal components corresponding to at least 50% of
the variation in the data.

## 5.2 Accounting for unwanted variation

We observed a “sample effect” in the data in the above PCA, which is
likely caused by the longitudinal study design, where sets of cell
cultures were resampled over a time series.

<details>
<summary>
Click to expand code block
</summary>

``` r
  data_pca_multilevel <- plot_pca_multilevel(
    data, classes, pch=pch, ncomp=10,
    title=paste("With NA.")
  )
  #> [1] "Removing 0 variance columns from data..."
  #> [1] "Plotting PCA multilevel component contribution..."
```

![](figs_1/sample_effect-1.png)<!-- -->![](figs_1/sample_effect-2.png)<!-- -->

      #> [1] "Plotting PCA multilevel..."

![](figs_1/sample_effect-3.png)<!-- -->![](figs_1/sample_effect-4.png)<!-- -->![](figs_1/sample_effect-5.png)<!-- -->![](figs_1/sample_effect-6.png)<!-- -->![](figs_1/sample_effect-7.png)<!-- -->![](figs_1/sample_effect-8.png)<!-- -->

``` r
  data_imp <- replace_missing(data, data_imp)
  names(data_imp) <- data_names
  data_pca_multilevel <- plot_pca_multilevel(
    data_imp, classes, pch=pch, ncomp=10,
    title=paste("Imputed.")
  )
  #> [1] "Removing 0 variance columns from data..."
  #> [1] "Plotting PCA multilevel component contribution..."
```

![](figs_1/sample_effect-9.png)<!-- -->![](figs_1/sample_effect-10.png)<!-- -->

      #> [1] "Plotting PCA multilevel..."

![](figs_1/sample_effect-11.png)<!-- -->![](figs_1/sample_effect-12.png)<!-- -->![](figs_1/sample_effect-13.png)<!-- -->![](figs_1/sample_effect-14.png)<!-- -->![](figs_1/sample_effect-15.png)<!-- -->![](figs_1/sample_effect-16.png)<!-- -->
</details>

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
    data_imp, classes, data_names, data.frame(pch), ncomp=4, nrepeat=10,
    logratio="none", test_keepX=c(5,6,7,8,9,10,30), validation="loo", folds=10,
    dist="centroids.dist", cpus=6, progressBar=FALSE
  )
```

![](figs_1/tune_splsda-1.png)<!-- -->![](figs_1/tune_splsda-2.png)<!-- -->

``` r
  splsda_keepx <- lapply(tuned_splsda, `[`, "choice.keepX")
  splsda_ncomp <- lapply(tuned_splsda, `[`, "choice.ncomp")

  print("Tuned splsda to use number of components:")
  #> [1] "Tuned splsda to use number of components:"
  splsda_ncomp <- lapply(splsda_ncomp, `[`, "ncomp")
  splsda_ncomp <- unlist(splsda_ncomp, recursive=FALSE)
  names(splsda_ncomp) <- data_names
  print(splsda_ncomp)
  #> $proteome
  #> NULL
  #>
  #> $translatome
  #> NULL

  print("Tuned the number of variables selected on each component to:")
  #> [1] "Tuned the number of variables selected on each component to:"
  print(splsda_keepx)
  #> $proteome
  #> $proteome$choice.keepX
  #> comp1 comp2 comp3 comp4
  #>    30    30    30    30
  #>
  #>
  #> $translatome
  #> $translatome$choice.keepX
  #> comp1 comp2 comp3 comp4
  #>     5     5    30     6
  splsda_keepx <- unlist(splsda_keepx, recursive=FALSE)
  names(splsda_keepx) <- data_names
  print(splsda_keepx)
  #> $proteome
  #> comp1 comp2 comp3 comp4
  #>    30    30    30    30
  #>
  #> $translatome
  #> comp1 comp2 comp3 comp4
  #>     5     5    30     6
```

Using the R data object:

``` r
  names(tuned_splsda)

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

``` r
  # this step can take some time
  data_splsda <- classify_splsda(
    data=data_imp, classes=classes, pch=pch, title=data_names,
    ncomp=4, keepX=splsda_keepx, contrib="max", outdir="./",
    mappings=NULL, dist="centroids.dist", bg=TRUE
  )
  #> [1] "splsda components:"
  #> [1] 4
  #> [1] "number of variables on each component:"
  #> comp1 comp2 comp3 comp4
  #>    30    30    30    30
```

![](figs_1/splsda_result-1.png)<!-- -->![](figs_1/splsda_result-2.png)<!-- -->![](figs_1/splsda_result-3.png)<!-- -->![](figs_1/splsda_result-4.png)<!-- -->

      #> [1] "Getting performance metrics"
      #> [1] "Plotting error rates..."
      #>
      #> comp 1
      #>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
      #> comp 2
      #>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
      #> comp 3
      #>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
      #> comp 4
      #>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
      #> $overall
      #>        max.dist centroids.dist mahalanobis.dist
      #> comp1 0.7500000     0.58333333        0.5833333
      #> comp2 0.6250000     0.25000000        0.5000000
      #> comp3 0.4166667     0.12500000        0.2083333
      #> comp4 0.3333333     0.08333333        0.2500000
      #>
      #> $BER
      #>        max.dist centroids.dist mahalanobis.dist
      #> comp1 0.7500000     0.58333333        0.5833333
      #> comp2 0.6250000     0.25000000        0.5000000
      #> comp3 0.4166667     0.12500000        0.2083333
      #> comp4 0.3333333     0.08333333        0.2500000

![](figs_1/splsda_result-5.png)<!-- -->

      #> [1] "Plotting stability of sPLSDA..."

![](figs_1/splsda_result-6.png)<!-- -->![](figs_1/splsda_result-7.png)<!-- -->![](figs_1/splsda_result-8.png)<!-- -->![](figs_1/splsda_result-9.png)<!-- -->![](figs_1/splsda_result-10.png)<!-- -->![](figs_1/splsda_result-11.png)<!-- -->![](figs_1/splsda_result-12.png)<!-- -->

      #> [1] "Plotting arrow plot..."

![](figs_1/splsda_result-13.png)<!-- -->

      #> [1] "Getting loadings and plotting clustered image maps"

![](figs_1/splsda_result-14.png)<!-- -->![](figs_1/splsda_result-15.png)<!-- -->![](figs_1/splsda_result-16.png)<!-- -->![](figs_1/splsda_result-17.png)<!-- -->

      #> [1] "Writing sPLSDA loadings to:"
      #> [1] ".//proteome_1_sPLSDA_max.txt"
      #> [1] ".//proteome_1_sPLSDA_min.txt"

![](figs_1/splsda_result-18.png)<!-- -->![](figs_1/splsda_result-19.png)<!-- -->![](figs_1/splsda_result-20.png)<!-- -->

      #> [1] "Writing sPLSDA loadings to:"
      #> [1] ".//proteome_2_sPLSDA_max.txt"
      #> [1] ".//proteome_2_sPLSDA_min.txt"

![](figs_1/splsda_result-21.png)<!-- -->![](figs_1/splsda_result-22.png)<!-- -->![](figs_1/splsda_result-23.png)<!-- -->

      #> [1] "Writing sPLSDA loadings to:"
      #> [1] ".//proteome_3_sPLSDA_max.txt"
      #> [1] ".//proteome_3_sPLSDA_min.txt"

![](figs_1/splsda_result-24.png)<!-- -->![](figs_1/splsda_result-25.png)<!-- -->![](figs_1/splsda_result-26.png)<!-- -->

      #> [1] "Writing sPLSDA loadings to:"
      #> [1] ".//proteome_4_sPLSDA_max.txt"
      #> [1] ".//proteome_4_sPLSDA_min.txt"
      #> [1] "splsda components:"
      #> [1] 4
      #> [1] "number of variables on each component:"
      #> comp1 comp2 comp3 comp4
      #>     5     5    30     6

![](figs_1/splsda_result-27.png)<!-- -->![](figs_1/splsda_result-28.png)<!-- -->![](figs_1/splsda_result-29.png)<!-- -->![](figs_1/splsda_result-30.png)<!-- -->

      #> [1] "Getting performance metrics"
      #> [1] "Plotting error rates..."
      #>
      #> comp 1
      #>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
      #> comp 2
      #>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
      #> comp 3
      #>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
      #> comp 4
      #>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
      #> $overall
      #>        max.dist centroids.dist mahalanobis.dist
      #> comp1 0.7500000     0.50000000       0.50000000
      #> comp2 0.5833333     0.41666667       0.45833333
      #> comp3 0.4583333     0.20833333       0.33333333
      #> comp4 0.2916667     0.08333333       0.08333333
      #>
      #> $BER
      #>        max.dist centroids.dist mahalanobis.dist
      #> comp1 0.7500000     0.50000000       0.50000000
      #> comp2 0.5833333     0.41666667       0.45833333
      #> comp3 0.4583333     0.20833333       0.33333333
      #> comp4 0.2916667     0.08333333       0.08333333

![](figs_1/splsda_result-31.png)<!-- -->

      #> [1] "Plotting stability of sPLSDA..."

![](figs_1/splsda_result-32.png)<!-- -->![](figs_1/splsda_result-33.png)<!-- -->![](figs_1/splsda_result-34.png)<!-- -->![](figs_1/splsda_result-35.png)<!-- -->![](figs_1/splsda_result-36.png)<!-- -->![](figs_1/splsda_result-37.png)<!-- -->![](figs_1/splsda_result-38.png)<!-- -->

      #> [1] "Plotting arrow plot..."

![](figs_1/splsda_result-39.png)<!-- -->

      #> [1] "Getting loadings and plotting clustered image maps"

![](figs_1/splsda_result-40.png)<!-- -->![](figs_1/splsda_result-41.png)<!-- -->![](figs_1/splsda_result-42.png)<!-- -->![](figs_1/splsda_result-43.png)<!-- -->

      #> [1] "Writing sPLSDA loadings to:"
      #> [1] ".//translatome_1_sPLSDA_max.txt"
      #> [1] ".//translatome_1_sPLSDA_min.txt"

![](figs_1/splsda_result-44.png)<!-- -->![](figs_1/splsda_result-45.png)<!-- -->![](figs_1/splsda_result-46.png)<!-- -->

      #> [1] "Writing sPLSDA loadings to:"
      #> [1] ".//translatome_2_sPLSDA_max.txt"
      #> [1] ".//translatome_2_sPLSDA_min.txt"

![](figs_1/splsda_result-47.png)<!-- -->![](figs_1/splsda_result-48.png)<!-- -->![](figs_1/splsda_result-49.png)<!-- -->

      #> [1] "Writing sPLSDA loadings to:"
      #> [1] ".//translatome_3_sPLSDA_max.txt"
      #> [1] ".//translatome_3_sPLSDA_min.txt"

![](figs_1/splsda_result-50.png)<!-- -->![](figs_1/splsda_result-51.png)<!-- -->![](figs_1/splsda_result-52.png)<!-- -->

      #> [1] "Writing sPLSDA loadings to:"
      #> [1] ".//translatome_4_sPLSDA_max.txt"
      #> [1] ".//translatome_4_sPLSDA_min.txt"

</details>

We also run PLSDA (all features) for comparison.

<details>
<summary>
Click to expand code block
</summary>

``` r
  # this step can take some time
  data_plsda <- classify_plsda(
    data=data_imp, classes=classes, pch=pch, title=data_names,
    ncomp=4, contrib="max", outdir="./",
    mappings=NULL, dist="centroids.dist", bg=TRUE
  )
  #> [1] "Plotting multi level partial least squares discriminant analysis"
```

![](figs_1/plsda_result-1.png)<!-- -->![](figs_1/plsda_result-2.png)<!-- -->![](figs_1/plsda_result-3.png)<!-- -->![](figs_1/plsda_result-4.png)<!-- -->

      #> [1] "Getting performance metrics"
      #> [1] "Plotting error rates..."
      #>
      #> comp 1
      #>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
      #> comp 2
      #>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
      #> comp 3
      #>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
      #> comp 4
      #>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
      #> $overall
      #>        max.dist centroids.dist mahalanobis.dist
      #> comp1 0.7500000      0.6250000        0.6250000
      #> comp2 0.6250000      0.4583333        0.3750000
      #> comp3 0.4166667      0.2500000        0.2500000
      #> comp4 0.2083333      0.1666667        0.1666667
      #>
      #> $BER
      #>        max.dist centroids.dist mahalanobis.dist
      #> comp1 0.7500000      0.6250000        0.6250000
      #> comp2 0.6250000      0.4583333        0.3750000
      #> comp3 0.4166667      0.2500000        0.2500000
      #> comp4 0.2083333      0.1666667        0.1666667

![](figs_1/plsda_result-5.png)<!-- -->![](figs_1/plsda_result-6.png)<!-- -->![](figs_1/plsda_result-7.png)<!-- -->![](figs_1/plsda_result-8.png)<!-- -->![](figs_1/plsda_result-9.png)<!-- -->

      #> [1] "Plotting arrow plot..."

![](figs_1/plsda_result-10.png)<!-- -->

      #> [1] "Getting loadings and plotting clustered image maps"

![](figs_1/plsda_result-11.png)<!-- -->![](figs_1/plsda_result-12.png)<!-- -->![](figs_1/plsda_result-13.png)<!-- -->![](figs_1/plsda_result-14.png)<!-- -->

      #> [1] "Writing PLSDA loadings to:"
      #> [1] ".//proteome_1_PLSDA_max.txt"
      #> [1] ".//proteome_1_PLSDA_min.txt"

![](figs_1/plsda_result-15.png)<!-- -->![](figs_1/plsda_result-16.png)<!-- -->![](figs_1/plsda_result-17.png)<!-- -->

      #> [1] "Writing PLSDA loadings to:"
      #> [1] ".//proteome_2_PLSDA_max.txt"
      #> [1] ".//proteome_2_PLSDA_min.txt"

![](figs_1/plsda_result-18.png)<!-- -->![](figs_1/plsda_result-19.png)<!-- -->![](figs_1/plsda_result-20.png)<!-- -->

      #> [1] "Writing PLSDA loadings to:"
      #> [1] ".//proteome_3_PLSDA_max.txt"
      #> [1] ".//proteome_3_PLSDA_min.txt"

![](figs_1/plsda_result-21.png)<!-- -->![](figs_1/plsda_result-22.png)<!-- -->![](figs_1/plsda_result-23.png)<!-- -->

      #> [1] "Writing PLSDA loadings to:"
      #> [1] ".//proteome_4_PLSDA_max.txt"
      #> [1] ".//proteome_4_PLSDA_min.txt"
      #> [1] "Plotting multi level partial least squares discriminant analysis"

![](figs_1/plsda_result-24.png)<!-- -->![](figs_1/plsda_result-25.png)<!-- -->![](figs_1/plsda_result-26.png)<!-- -->![](figs_1/plsda_result-27.png)<!-- -->

      #> [1] "Getting performance metrics"
      #> [1] "Plotting error rates..."
      #>
      #> comp 1
      #>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
      #> comp 2
      #>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
      #> comp 3
      #>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
      #> comp 4
      #>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
      #> $overall
      #>        max.dist centroids.dist mahalanobis.dist
      #> comp1 0.7500000      0.6250000        0.6250000
      #> comp2 0.6250000      0.3333333        0.5833333
      #> comp3 0.3333333      0.2083333        0.4583333
      #> comp4 0.1250000      0.1250000        0.4166667
      #>
      #> $BER
      #>        max.dist centroids.dist mahalanobis.dist
      #> comp1 0.7500000      0.6250000        0.6250000
      #> comp2 0.6250000      0.3333333        0.5833333
      #> comp3 0.3333333      0.2083333        0.4583333
      #> comp4 0.1250000      0.1250000        0.4166667

![](figs_1/plsda_result-28.png)<!-- -->![](figs_1/plsda_result-29.png)<!-- -->![](figs_1/plsda_result-30.png)<!-- -->![](figs_1/plsda_result-31.png)<!-- -->![](figs_1/plsda_result-32.png)<!-- -->

      #> [1] "Plotting arrow plot..."

![](figs_1/plsda_result-33.png)<!-- -->

      #> [1] "Getting loadings and plotting clustered image maps"

![](figs_1/plsda_result-34.png)<!-- -->![](figs_1/plsda_result-35.png)<!-- -->![](figs_1/plsda_result-36.png)<!-- -->![](figs_1/plsda_result-37.png)<!-- -->

      #> [1] "Writing PLSDA loadings to:"
      #> [1] ".//translatome_1_PLSDA_max.txt"
      #> [1] ".//translatome_1_PLSDA_min.txt"

![](figs_1/plsda_result-38.png)<!-- -->![](figs_1/plsda_result-39.png)<!-- -->![](figs_1/plsda_result-40.png)<!-- -->

      #> [1] "Writing PLSDA loadings to:"
      #> [1] ".//translatome_2_PLSDA_max.txt"
      #> [1] ".//translatome_2_PLSDA_min.txt"

![](figs_1/plsda_result-41.png)<!-- -->![](figs_1/plsda_result-42.png)<!-- -->![](figs_1/plsda_result-43.png)<!-- -->

      #> [1] "Writing PLSDA loadings to:"
      #> [1] ".//translatome_3_PLSDA_max.txt"
      #> [1] ".//translatome_3_PLSDA_min.txt"

![](figs_1/plsda_result-44.png)<!-- -->![](figs_1/plsda_result-45.png)<!-- -->![](figs_1/plsda_result-46.png)<!-- -->

      #> [1] "Writing PLSDA loadings to:"
      #> [1] ".//translatome_4_PLSDA_max.txt"
      #> [1] ".//translatome_4_PLSDA_min.txt"

</details>

These automatically generate a large series of plots. Figures are
ordered sequentially by each block of omics data in `data_name`, in this
case `proteome` followed by `translatome`. Some `txt` files containing
feature loadings are also written to the output directory.

## 6.3 Explanation of output

Detailed technical information on the individual figures types and
underlying methods are available at the [mixOmics
website](http://mixomics.org/). This walkthrough will explain the plots
in context of the biological system under study only. The plots are
explained in the order that they appear in this document.

### 6.3.1 Scatter plots

Plotting the first few components of the sPLSDA reveals several distinct
sample groups, with the main distinction in both omics data blocks as
the difference between late stage infected samples with all other
samples. This matches the independent observations of the laboratory
which originally generated the data. There are multiple secondary
distinctions between data groups in both omics data blocks, mostly
between groups of timepoints. It is also interesting to note that some
sub-groups of data include both infected and uninfected samples.

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
  tuned_diablo <- tune_diablo_ncomp(data_imp, classes, design, ncomp=8, cpus=6)
  print("Parameters with lowest error rate:")
  tuned_diablo <- tuned_diablo$choice.ncomp$WeightedVote["Overall.BER",]
  diablo_ncomp <- tuned_diablo[which.max(tuned_diablo)]
  print("Number of components:")
  print(diablo_ncomp)

  # tune keepx
  diablo_keepx <- c(5,10,12,14,16,18,20,30)
  diablo_keepx <- tune_diablo_keepx(
    data_imp, classes, diablo_ncomp, design, diablo_keepx, cpus=6,
    dist="mahalanobis.dist", progressBar=FALSE
  )
  print("Diablo keepx:")
  print(diablo_keepx)
```

Using the R object:

``` r
  tuned_diablo
  #>         max.dist   centroids.dist mahalanobis.dist
  #>                8                1                8
  diablo$keepX
  #> $proteome
  #> [1] 14 14 30  5 10 10 10  5
  #>
  #> $translatome
  #> [1] 10  5 20  5  5  5 10  5
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
  diablo_keepx <- list(proteome=c(14, 14, 30, 5, 10, 10, 10, 5),
                       translatome=c(10, 5, 20, 5, 5, 5, 10, 5))

  data_imp <- force_unique_blocks(data_imp)

  # run the algorithm
  diablo <- run_diablo(data_imp, classes, diablo_ncomp, design, diablo_keepx)
```

Using the R data object:

``` r
  diablo
  #>
  #> Call:
  #>  block.splsda(X = data, Y = classes, ncomp = ncomp, keepX = keepx, design = design)
  #>
  #>  sGCCA with 8 components on block 1 named proteome
  #>  sGCCA with 8 components on block 2 named translatome
  #>  sGCCA with 8 components on the outcome Y
  #>
  #>  Dimension of block 1 is  24 6380
  #>  Dimension of block 2 is  24 1595
  #>  Outcome Y has 8 levels
  #>
  #>  Selection of 14 14 30 5 10 10 10 5 variables on each of the sGCCA components on the block 1
  #>  Selection of 10 5 20 5 5 5 10 5 variables on each of the sGCCA components on the block 2
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

![](figs_1/perf_diablo-1.png)<!-- -->
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
  print("Proteome")
  #> [1] "Proteome"
  plot(perf_diablo$features$stable$nrep1$proteome$comp1,
    type="h", main="Comp 1", las=2, ylab="Stability", xlab="Features", xaxt='n'
  )
```

![](figs_1/perf_diablo_stability-1.png)<!-- -->

``` r
  print("Translatome")
  #> [1] "Translatome"
  plot(perf_diablo$features$stable$nrep1$translatome$comp1,
    type="h", main="Comp 1", las=2, ylab="Stability", xlab="Features", xaxt='n'
  )
```

![](figs_1/perf_diablo_stability-2.png)<!-- -->
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

![](figs_1/diablo-1.png)<!-- -->
</details>

### 7.3.4 Scatter plots

Plotting the first few components of the multi-block sPLSDA (DIABLO)
reveals several distinct sample groups, with the main distinction in
both omics data blocks as the difference between late stage infected
samples with all other samples. This matches the independent
observations of the laboratory which originally generated the data.
There are multiple secondary distinctions between data groups in both
omics data blocks, mostly between groups of timepoints. It is also
interesting to note that some sub-groups of data include both infected
and uninfected samples.

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

![](figs_1/perf_scatter-1.png)<!-- -->
</details>

### 7.3.5 ROC curves

As a supplementary layer of validation, ROC curves showing
classification accuracy are available, but we note that these have
limited applicability given the specific context of the method. The
underlying method already internally specifies the prediction cutoff to
achieve maximal sensitivity and specificity.

<details>
<summary>
Click to expand code block
</summary>

``` r
  # to keep the case study concise we use a custom function to output main plots
  sink("/dev/null")
  mixOmics::auroc(diablo, roc.comp=1)
```

![](figs_1/perf_auroc-1.png)<!-- -->

      #> $proteome
      #> $proteome$comp1
      #>                            AUC  p-value
      #> Control_10h vs Other(s) 0.3968 0.570500
      #> Control_24h vs Other(s) 0.9683 0.010020
      #> Control_2h vs Other(s)  0.5238 0.895800
      #> Control_6h vs Other(s)  0.6667 0.359400
      #> Virus_10h vs Other(s)   0.2222 0.126600
      #> Virus_24h vs Other(s)   1.0000 0.005968
      #> Virus_2h vs Other(s)    0.4762 0.895800
      #> Virus_6h vs Other(s)    0.7460 0.176100
      #>
      #> $proteome$comp2
      #>                            AUC  p-value
      #> Control_10h vs Other(s) 0.8571 0.049530
      #> Control_24h vs Other(s) 1.0000 0.005968
      #> Control_2h vs Other(s)  0.7619 0.149800
      #> Control_6h vs Other(s)  0.7619 0.149800
      #> Virus_10h vs Other(s)   0.7143 0.238600
      #> Virus_24h vs Other(s)   1.0000 0.005968
      #> Virus_2h vs Other(s)    0.9524 0.012860
      #> Virus_6h vs Other(s)    0.6667 0.359400
      #>
      #> $proteome$comp3
      #>                            AUC  p-value
      #> Control_10h vs Other(s) 1.0000 0.005968
      #> Control_24h vs Other(s) 1.0000 0.005968
      #> Control_2h vs Other(s)  0.9524 0.012860
      #> Control_6h vs Other(s)  0.9048 0.026030
      #> Virus_10h vs Other(s)   0.8571 0.049530
      #> Virus_24h vs Other(s)   1.0000 0.005968
      #> Virus_2h vs Other(s)    0.9206 0.020720
      #> Virus_6h vs Other(s)    0.7460 0.176100
      #>
      #> $proteome$comp4
      #>                            AUC  p-value
      #> Control_10h vs Other(s) 0.9524 0.012860
      #> Control_24h vs Other(s) 1.0000 0.005968
      #> Control_2h vs Other(s)  0.9683 0.010020
      #> Control_6h vs Other(s)  1.0000 0.005968
      #> Virus_10h vs Other(s)   0.9206 0.020720
      #> Virus_24h vs Other(s)   1.0000 0.005968
      #> Virus_2h vs Other(s)    0.9206 0.020720
      #> Virus_6h vs Other(s)    0.8413 0.060560
      #>
      #> $proteome$comp5
      #>                           AUC  p-value
      #> Control_10h vs Other(s) 1.000 0.005968
      #> Control_24h vs Other(s) 1.000 0.005968
      #> Control_2h vs Other(s)  1.000 0.005968
      #> Control_6h vs Other(s)  1.000 0.005968
      #> Virus_10h vs Other(s)   1.000 0.005968
      #> Virus_24h vs Other(s)   1.000 0.005968
      #> Virus_2h vs Other(s)    1.000 0.005968
      #> Virus_6h vs Other(s)    0.873 0.040240
      #>
      #> $proteome$comp6
      #>                            AUC  p-value
      #> Control_10h vs Other(s) 1.0000 0.005968
      #> Control_24h vs Other(s) 1.0000 0.005968
      #> Control_2h vs Other(s)  1.0000 0.005968
      #> Control_6h vs Other(s)  1.0000 0.005968
      #> Virus_10h vs Other(s)   1.0000 0.005968
      #> Virus_24h vs Other(s)   1.0000 0.005968
      #> Virus_2h vs Other(s)    1.0000 0.005968
      #> Virus_6h vs Other(s)    0.8889 0.032470
      #>
      #> $proteome$comp7
      #>                         AUC  p-value
      #> Control_10h vs Other(s)   1 0.005968
      #> Control_24h vs Other(s)   1 0.005968
      #> Control_2h vs Other(s)    1 0.005968
      #> Control_6h vs Other(s)    1 0.005968
      #> Virus_10h vs Other(s)     1 0.005968
      #> Virus_24h vs Other(s)     1 0.005968
      #> Virus_2h vs Other(s)      1 0.005968
      #> Virus_6h vs Other(s)      1 0.005968
      #>
      #> $proteome$comp8
      #>                         AUC  p-value
      #> Control_10h vs Other(s)   1 0.005968
      #> Control_24h vs Other(s)   1 0.005968
      #> Control_2h vs Other(s)    1 0.005968
      #> Control_6h vs Other(s)    1 0.005968
      #> Virus_10h vs Other(s)     1 0.005968
      #> Virus_24h vs Other(s)     1 0.005968
      #> Virus_2h vs Other(s)      1 0.005968
      #> Virus_6h vs Other(s)      1 0.005968
      #>
      #>
      #> $translatome
      #> $translatome$comp1
      #>                            AUC  p-value
      #> Control_10h vs Other(s) 0.1746 0.073550
      #> Control_24h vs Other(s) 0.7302 0.205600
      #> Control_2h vs Other(s)  0.5556 0.760000
      #> Control_6h vs Other(s)  0.6984 0.275200
      #> Virus_10h vs Other(s)   0.4127 0.631200
      #> Virus_24h vs Other(s)   1.0000 0.005968
      #> Virus_2h vs Other(s)    0.5079 0.965200
      #> Virus_6h vs Other(s)    0.9206 0.020720
      #>
      #> $translatome$comp2
      #>                            AUC  p-value
      #> Control_10h vs Other(s) 0.8571 0.049530
      #> Control_24h vs Other(s) 1.0000 0.005968
      #> Control_2h vs Other(s)  0.6508 0.407000
      #> Control_6h vs Other(s)  0.7619 0.149800
      #> Virus_10h vs Other(s)   0.6032 0.570500
      #> Virus_24h vs Other(s)   1.0000 0.005968
      #> Virus_2h vs Other(s)    0.9048 0.026030
      #> Virus_6h vs Other(s)    0.7302 0.205600
      #>
      #> $translatome$comp3
      #>                            AUC  p-value
      #> Control_10h vs Other(s) 1.0000 0.005968
      #> Control_24h vs Other(s) 1.0000 0.005968
      #> Control_2h vs Other(s)  0.9048 0.026030
      #> Control_6h vs Other(s)  0.8095 0.088740
      #> Virus_10h vs Other(s)   0.8413 0.060560
      #> Virus_24h vs Other(s)   1.0000 0.005968
      #> Virus_2h vs Other(s)    0.9524 0.012860
      #> Virus_6h vs Other(s)    0.7302 0.205600
      #>
      #> $translatome$comp4
      #>                            AUC  p-value
      #> Control_10h vs Other(s) 0.7937 0.106400
      #> Control_24h vs Other(s) 1.0000 0.005968
      #> Control_2h vs Other(s)  0.9524 0.012860
      #> Control_6h vs Other(s)  1.0000 0.005968
      #> Virus_10h vs Other(s)   0.9206 0.020720
      #> Virus_24h vs Other(s)   1.0000 0.005968
      #> Virus_2h vs Other(s)    0.9206 0.020720
      #> Virus_6h vs Other(s)    0.5556 0.760000
      #>
      #> $translatome$comp5
      #>                            AUC  p-value
      #> Control_10h vs Other(s) 0.8730 0.040240
      #> Control_24h vs Other(s) 0.9841 0.007762
      #> Control_2h vs Other(s)  1.0000 0.005968
      #> Control_6h vs Other(s)  1.0000 0.005968
      #> Virus_10h vs Other(s)   1.0000 0.005968
      #> Virus_24h vs Other(s)   1.0000 0.005968
      #> Virus_2h vs Other(s)    0.8571 0.049530
      #> Virus_6h vs Other(s)    0.7143 0.238600
      #>
      #> $translatome$comp6
      #>                            AUC  p-value
      #> Control_10h vs Other(s) 0.9683 0.010020
      #> Control_24h vs Other(s) 0.9841 0.007762
      #> Control_2h vs Other(s)  1.0000 0.005968
      #> Control_6h vs Other(s)  1.0000 0.005968
      #> Virus_10h vs Other(s)   1.0000 0.005968
      #> Virus_24h vs Other(s)   1.0000 0.005968
      #> Virus_2h vs Other(s)    1.0000 0.005968
      #> Virus_6h vs Other(s)    0.7143 0.238600
      #>
      #> $translatome$comp7
      #>                            AUC  p-value
      #> Control_10h vs Other(s) 0.9683 0.010020
      #> Control_24h vs Other(s) 0.9841 0.007762
      #> Control_2h vs Other(s)  1.0000 0.005968
      #> Control_6h vs Other(s)  1.0000 0.005968
      #> Virus_10h vs Other(s)   1.0000 0.005968
      #> Virus_24h vs Other(s)   1.0000 0.005968
      #> Virus_2h vs Other(s)    1.0000 0.005968
      #> Virus_6h vs Other(s)    1.0000 0.005968
      #>
      #> $translatome$comp8
      #>                            AUC  p-value
      #> Control_10h vs Other(s) 0.9841 0.007762
      #> Control_24h vs Other(s) 0.9841 0.007762
      #> Control_2h vs Other(s)  1.0000 0.005968
      #> Control_6h vs Other(s)  1.0000 0.005968
      #> Virus_10h vs Other(s)   1.0000 0.005968
      #> Virus_24h vs Other(s)   1.0000 0.005968
      #> Virus_2h vs Other(s)    1.0000 0.005968
      #> Virus_6h vs Other(s)    1.0000 0.005968
      sink()

</details>

### 7.3.6 Arrow plots

To review agreement between matching data sets, we use an arrow plot.
Short arrows suggest a strong similarity, while long arrows indicate
dissimilarity. Two sample groups in the proteome data (Virus 6h, Virus
24h) appear to be dissimilar.

<details>
<summary>
Click to expand code block
</summary>

``` r
  # to keep the case study concise we use a custom function to output main plots
  mixOmics::plotArrow(diablo, ind.names=FALSE, legend=TRUE, title='DIABLO')
```

![](figs_1/perf_arrow-1.png)<!-- -->
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

![](figs_1/perf_correlation-1.png)<!-- -->
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

![](figs_1/cim-1.png)<!-- -->
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

![](figs_1/loadings-1.png)<!-- -->
</details>

### 7.3.10 Circos plot

To visualise correlations between different blocks of omics data, a
circos plot is generated. The blue block represents proteome data and
the green block represents translatome data. Each point on the circle is
a single feature. Lines linking features show correlations between
features that pass the user-specified correlation threshold, in this
case 0.95. Red lines indicate positive correlation and blue lines
negative correlation.

<details>
<summary>
Click to expand code block
</summary>

``` r
  # to keep the case study concise we use a custom function to output main plots
  mixOmics::circosPlot(
    diablo, cutoff=0.95, line=FALSE, size.legend=0.5, size.variables=0.001
  )
```

![](figs_1/circos-1.png)<!-- -->
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

![](figs_1/network-1.png)<!-- -->
</details>

# 8 Output data

Files output by the pipeline include:

-   a `pdf` file of all plots generated by the pipeline
-   tab-separated `txt` files containing feature contribution weights to
    each biological class
-   tab-separated `txt` file containing correlations between each omics
    data block

[A `RData` object with all input and output is available in the git
repository.](https://gitlab.com/tyagilab/sars-cov-2/-/blob/master/results/case_study_1/RData.RData)
This is not included directly in the `multiomics` package because of
size constraints, and includes data from two omics datasets.

# 9 Acknowledgements

[Please refer to introduction.](introduction.html)

# 10 References

[Please refer to introduction.](introduction.html)
