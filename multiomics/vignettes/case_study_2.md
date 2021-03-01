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
    -   [6.3 Diagnostic plots](#diagnostic-plots)
    -   [6.4 Results](#results)
-   [7 Multi omics analyses](#multi-omics-analyses)
    -   [7.1 Parameter tuning](#parameter-tuning-1)
    -   [7.2 Diagnostic plots](#diagnostic-plots-1)
    -   [7.3 Results](#results-1)
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
```

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
  url_lipi <- "https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/data/case_study_2/proteome_mapfile.txt"
  url_meta <- "https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/data/case_study_2/translatome_mapfile.txt"
  url_prot <- "https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/data/case_study_2/diablo_proteome.txt"
  url_tran <- "https://gitlab.com/tyagilab/sars-cov-2/-/raw/master/data/case_study_2/diablo_translatome.txt"
  
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

> **NOTE**: There are some differences in the R data object for case
> studies 1 and 2. Case study 1 was originally performed using an early
> version of the code. It has since been streamlined for case study 2,
> which contains more information. However, both analyses are
> reproducible and the user can if needed further investigate the
> internal data structures after loading them.

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

Using the pipeline:

``` r
  # you may need this line if NA is represented as 0 in your data
  # in this specific case study the data has already been cleaned
  data <- mapply(function(x) zero_to_na(x), data, SIMPLIFY=FALSE)
  missing <- lapply(data, count_missing)
  #> [1] "Percentage of missing values in data:"
  #> [1] 0
  #> [1] "Percentage of missing values in data:"
  #> [1] 0
  #> [1] "Percentage of missing values in data:"
  #> [1] 0
  #> [1] "Percentage of missing values in data:"
  #> [1] 0
  mapply(function(x, y) show_na_prop(x, y), data, data_names, SIMPLIFY=FALSE)
```

![](figs_2/check_na-1.png)<!-- -->![](figs_2/check_na-2.png)<!-- -->![](figs_2/check_na-3.png)<!-- -->![](figs_2/check_na-4.png)<!-- -->

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
  missing <- lapply(data, count_missing)
  pca_withna <- plot_pca_single(
    data, classes, ncomp=10,
    title=paste("With NA")
  )
  pca_impute <- plot_pca_single(
    data_imp, classes, ncomp=10,
    title=paste("Imputed")
  )
  heatmaps <- cor_imputed_unimputed(pca_withna, pca_impute, data_names)
```

Using the R Data object:

To improve efficiency during code run, a imputed data file was generated
as output in the first run of the pipeline. This was then used in
subsequent runs of the pipeline as direct input. You can reproduce this
manually with the code in the `make_manuscript_figures.R` [script in our
gitlab repository](https://gitlab.com/tyagilab/sars-cov-2/).
</details>

In both cases, there is a strong correlation between the variates on at
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
  
  # keep number of components or use user specified
  splsda_ncomp <- lapply(tuned_splsda, `[`, "choice.ncomp")
  splsda_ncomp <- unlist(splsda_ncomp, recursive=FALSE)
  
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
    data=data_imp, classes=classes, pch=pch, title=data_names,
    ncomp=splsda_ncomp, keepX=splsda_keepx, contrib="max", outdir="./",
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
    data=data_imp, classes=classes, pch=pch, title=data_names,
    ncomp=4, contrib="max", outdir="./",
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

## 6.3 Diagnostic plots

## 6.4 Results

# 7 Multi omics analyses

## 7.1 Parameter tuning

## 7.2 Diagnostic plots

## 7.3 Results

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
