Introduction
================

-   [1 Index](#index)
-   [2 Installation](#installation)
    -   [2.1 Quick](#quick)
    -   [2.2 Manual](#manual)
-   [3 Usage](#usage)
-   [4 Pipeline minimum input data](#pipeline-minimum-input-data)
    -   [4.1 Input data](#input-data)
    -   [4.2 Input files](#input-files)
-   [5 Pipeline output data](#pipeline-output-data)
-   [6 Acknowledgements](#acknowledgements)
-   [7 References](#references)
    -   [7.1 Data source in package and case
        studies](#data-source-in-package-and-case-studies)
    -   [7.2 Methods used in package](#methods-used-in-package)
    -   [7.3 R package compilation](#r-package-compilation)

Copyright (c) 2020
<a href="https://orcid.org/0000-0002-9207-0385">Tyrone Chen
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>,
<a href="https://orcid.org/0000-0002-4146-2848">Al J Abadi
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>,
<a href="https://orcid.org/0000-0003-3923-1116">Kim-Anh Lê Cao
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>,
<a href="https://orcid.org/0000-0003-0181-6258">Sonika Tyagi
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>

Code in this package and git repository
<https://github.com/tyronechen/SARS-CoV-2/> is provided under a [MIT
license](https://opensource.org/licenses/MIT). This documentation is
provided under a [CC-BY-3.0 AU
license](https://creativecommons.org/licenses/by/3.0/au/).

[Visit our lab website here.](https://bioinformaticslab.erc.monash.edu/)
Contact Sonika Tyagi at <sonika.tyagi@monash.edu>.

# 1 Index

> **NOTE**: The pipeline API has changed since the original publication. To reproduce the results in the original COVID-19 paper for case study 1 and 2, please use the specific version of the pipeline available on zenodo only. Please refer to case study 3 for latest usage.

-   [Introduction](introduction.md)
-   [Case study 1](case_study_1.md)
-   [Case study 2](case_study_2.md)
-   [Case study 3](case_study_3.md)

# 2 Installation

## 2.1 Quick

You can install this directly as an R package from gitlab. Note that you may get errors if you don't have `libgit2` and `freetype` libraries installed (these are not R packages).

> **NOTE**: This version of the pipeline is compatible with `R` version `4.2.3`.

```
    install.packages("devtools")
    library("devtools")
    install_git("https://github.com/tyronechen/SARS-CoV-2.git", subdir="multiomics", build_vignettes=FALSE, INSTALL_opts="--no-multiarch")
```

The actual script used to run the pipeline is not directly callable but
provided as a separate script.

```
    # this will show you the path to the script
    system.file("scripts", "run_pipeline.R", package="multiomics")
```

## 2.2 Manual (for developers or if above doesnt work)

Alternatively, clone the git repository with:

```
    git clone "https://github.com/tyronechen/SARS-CoV-2.git"
```

### 2.2.1 Install dependencies

[With `conda`](https://bioconda.github.io/user/install.html):

```
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

    install_me="r-argparser r-brio r-colorspace r-diffobj r-dplyr r-ellipsis \
     r-farver r-ggplot2 r-ggrepel r-igraph r-isoband r-matrixStats r-mixOmics \
     r-parallel r-plyr r-rARPACK r-Rcpp r-RcppEigen r-reshape2 r-RSpectra \
     r-stringi r-testthat r-tibble r-tidyr r-utf8 r-vctrs r-zeallot \
     bioconductor-biocparallel"

    conda create -n my_environment install ${install_me}
```

You can also install dependencies in `R` directly:

```
    install_me <- c(
      "argparser", "brio", "colorspace", "diffobj", "dplyr", "ellipsis", "farver",
      "ggplot2", "ggrepel", "igraph", "isoband", "matrixStats", "mixOmics",
      "parallel", "plyr", "rARPACK", "Rcpp", "RcppEigen", "reshape2", "RSpectra",
      "stringi", "testthat", "tibble", "tidyr", "utf8", "vctrs", "zeallot",
      "BiocParallel")
    sapply(install_me, install.packages)
```

If you run into any issues with the manual install, please double check the library versions against `multiomics/DESCRIPTION`.

# 3 Usage

Load the library.

``` r
library(multiomics)
```

If you installed this pipeline as an `R` package, an example pipeline
script is included. You can find it by running this command in your `R`
environment:

```
system.file("scripts", "run_pipeline.R", package="multiomics")
```

Otherwise, you can find a copy of this script in the public git
repository:
<https://github.com/tyronechen/SARS-CoV-2/blob/master/src/run_pipeline.R>

To inspect the arguments to the script, run this command:

```
Rscript run_pipeline.R --help
```

A minimal script to run the pipeline is below. [You can also download
this
here](https://github.com/tyronechen/SARS-CoV-2/blob/master/src/case_study_3/example.sh).
This example may take a few hours to run fully.

Data is provided as part of the `multiomics` package and not directly as files. Extract it first with this:

```
Rscript -e 'library(multiomics); data(BPH2819); names(BPH2819); export <- function(name, data) {write.table(data.frame(data), paste(name, ".tsv", sep=""), quote=FALSE, sep="\t", row.names=TRUE, col.names=NA)}; mapply(export, names(BPH2819), BPH2819, SIMPLIFY=FALSE)'
```

Four files will be generated in the current working directory, where classes contains sample information and remaining files contain corresponding omics data:

```
classes.tsv
metabolome.tsv
proteome.tsv
transcriptome.tsv
```

Then run the multiomics pipeline on the data:

```
Rscript run_pipeline.R \
  --classes classes.tsv \
  --data metabolome.tsv \
         proteome.tsv \
         transcriptome.tsv \
  --data_names metabolome proteome transcriptome \
  --ncpus 2 \
  --icomp 12 \
  --pcomp 10 \
  --plsdacomp 2 \
  --splsdacomp 2 \
  --diablocomp 2 \
  --dist_plsda "centroids.dist" \
  --dist_splsda "centroids.dist" \
  --dist_diablo "centroids.dist" \
  --cross_val "Mfold" \
  --cross_val_folds 5 \
  --cross_val_nrepeat 50 \
  --corr_cutoff 0.1 \
  --outfile_dir BPH2819 \
  --contrib "max" \
  --progress_bar
```       

# 4 Pipeline minimum input data

## 4.1 Input data

The minimum input data needed is a file of classes and at least two
files of quantitative omics data. Tab separated data is expected by
default.

A small example subset of test data is included in the package for
reference. In this test case, data has already been log2 transformed and missing values filled in with imputation.

``` r
data(BPH2819)
names(BPH2819)
#> [1] "classes"       "metabolome"    "proteome"      "transcriptome"
```

The class information is available as a vector:

``` r
BPH2819$classes
#> [1] "RPMI" "RPMI" "RPMI" "RPMI" "RPMI" "RPMI"
#> [6] "Sera" "Sera" "Sera" "Sera" "Sera" "Sera"
```

Each of the three omics data blocks have 12 matched samples and an
arbitrary number of features.

``` r
sapply(BPH2819, dim)
#> $classes
#> NULL
#> $metabolome
#> [1]  12 153
#> $proteome
#> [1]   12 1451
#> $transcriptome
#> [1]   12 2771
```

``` r
BPH2819$metabolome[,1:3]
#>         X3.Aminoglutaric.acid HMDB0000005 HMDB0000008
#> RPMI_0             -1.7814083   -9.103010   -3.471373
#> RPMI_1             -1.9108074   -5.401229   -3.488496
#> RPMI_2             -1.5458964  -10.898804   -2.845025
#> RPMI_3             -2.1842312   -9.563557   -1.232155
#> RPMI_4             -1.3106881   -4.755440   -1.723564
#> RPMI_5             -0.9600247   -4.771127   -1.403044
#> Sera_6             -0.8764074   -6.507606   -2.884537
#> Sera_7             -1.4139388  -11.175670   -2.861640
#> Sera_8             -3.7537269  -10.883382   -1.238028
#> Sera_9             -2.6902848  -10.744718   -2.635249
#> Sera_10            -3.3605788  -10.439710   -1.774845
#> Sera_11            -2.9362071   -5.850829   -1.139523
```

Important notes on input data:

1.  Class information and the **sample order in each omics dataset must
    be identical**.
2.  Ideally data should be already preprocessed and missing values should be   below 20%.
3.  Feature names in each omics dataset may be truncated. Too long names
    cause issues with visualisation.
4.  `R` silently replaces all non-alphanumeric characters in feature
    names with `.`.

To work around (3) and (4), you can rename your feature names to a short
alphanumeric ID in your files, and remap them back later.

## 4.2 Input files

If you did not install the `R` package, you can obtain these example files from
github:

-   [classes](https://github.com/tyronechen/SARS-CoV-2/blob/master/multiomics/data/classes.tsv)
-   [metabolome](https://github.com/tyronechen/SARS-CoV-2/blob/master/multiomics/data/metabolome.tsv)
-   [proteome](https://github.com/tyronechen/SARS-CoV-2/blob/master/multiomics/data/proteome.tsv)
-   [transcriptome](https://github.com/tyronechen/SARS-CoV-2/blob/master/multiomics/data/transcriptome.tsv)

# 5 Pipeline output data

Files output by the pipeline include:

-   a `pdf` file of all plots generated by the pipeline
-   tab-separated `txt` files containing feature contribution weights to
    each biological class
-   tab-separated `txt` file containing correlations between each omics
    data block

[A `RData` object with all input and output is available in the git
repository.](https://github.com/tyronechen/SARS-CoV-2/tree/master/results/case_study_3/data.RData)
This is not included directly in the `multiomics` package because of
size constraints, and includes data from three omics datasets.

# 6 Acknowledgements

[We thank David A. Matthews](https://orcid.org/0000-0003-4611-8795) for
helpful discussions and feedback. [We thank Yashpal
Ramakrishnaiah](https://orcid.org/0000-0002-2213-8348) for performing an
extended analysis of the primary data. [We thank Melcy
Philip](https://orcid.org/0000-0002-0827-866X) for performing downstream
analysis of the data. This work was supported by the [MASSIVE HPC
facility](www.massive.org.au) and the authors thank the HPC team at
Monash eResearch Centre for their continuous personnel support. This R
package was compiled referring to information from blog posts or books
by [Hilary
Parker](https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/),
[Fong Chun
Chan](https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html),
[Karl Broman](https://kbroman.org/pkg_primer/pages/data.html), [Yihui
Xie, J. J. Allaire, Garrett
Grolemund](https://bookdown.org/yihui/rmarkdown/) as well as [Jenny
Bryan and Hadley Wickham](https://r-pkgs.org/). [We acknowledge and pay
respects to the Elders and Traditional Owners of the land on which our 4
Australian campuses
stand](https://www.monash.edu/indigenous-australians/about-us/recognising-traditional-owners).

# 7 References

## 7.1 Data source in package and case studies

-   [Bojkova, D., Klann, K., Koch, B. et al. Proteomics of
    SARS-CoV-2-infected host cells reveals therapy targets. Nature 583,
    469–472 (2020).
    https://doi.org/10.1038/s41586-020-2332-7](https://doi.org/10.1038/s41586-020-2332-7)
-   [Overmyer, K.A., Shishkova, E., Miller, I.J., Balnis, J., Bernstein,
    M.N., Peters-Clarke, T.M., Meyer, J.G., Quan, Q., Muehlbauer, L.K.,
    Trujillo, E.A. and He, Y., 2021. Large-scale multi-omic analysis of
    COVID-19 severity. Cell
    systems.](https://dx.doi.org/10.1016%2Fj.cels.2020.10.003)
-   [Mu, A., Klare, W.P., Baines, S.L. et al. Integrative omics identifies conserved and pathogen-specific responses of sepsis-causing bacteria. Nat Commun 14, 1530 (2023)](https://doi.org/10.1038/s41467-023-37200-w)

## 7.2 Methods used in package

-   [Amrit Singh, Casey P Shannon, Benoît Gautier, Florian Rohart,
    Michaël Vacher, Scott J Tebbutt, Kim-Anh Lê Cao, DIABLO: an
    integrative approach for identifying key molecular drivers from
    multi-omics assays, Bioinformatics, Volume 35, Issue 17, 1 September
    2019, Pages 3055–3062,
    https://doi.org/10.1093/bioinformatics/bty1054](https://doi.org/10.1093/bioinformatics/bty1054)
-   Butte, A. J., Tamayo, P., Slonim, D., Golub, T. R. and Kohane, I. S.
    (2000). Discovering functional relationships between RNA expression
    and chemotherapeutic susceptibility using relevance networks.
    *Proceedings of the National Academy of Sciences of the USA* *97*,
    12182-12186.
-   Chavent, Marie and Patouille, Brigitte (2003). Calcul des
    coefficients de regression et du PRESS en regression PLS1. *Modulad
    n*, *30* 1-11.
-   Eisen, M. B., Spellman, P. T., Brown, P. O. and Botstein, D. (1998).
    Cluster analysis and display of genome-wide expression patterns.
    *Proceeding of the National Academy of Sciences of the USA* *95*,
    14863-14868.
-   Eslami, A., Qannari, E. M., Kohler, A., and Bougeard, S. (2013).
    Multi-group PLS Regression: Application to Epidemiology. In New
    Perspectives in Partial Least Squares and Related Methods, pages
    243-255. Springer.
-   González I., Lê Cao K.A., Davis M.J., Déjean S. (2012). Visualising
    associations between paired ‘omics’ data sets. *BioData Mining*;
    *5*(1)
-   [H.M. Blalock, A. Aganbegian, F.M. Borodkin, Raymond Boudon,
    Vittorio Capecchi. Path Models with Latent Variables: The NIPALS
    Approach. In International Perspectives on Mathematical and
    Statistical Modeling (1975).
    https://doi.org/10.1016/B978-0-12-103950-9.50017-4](https://doi.org/10.1016/B978-0-12-103950-9.50017-4)
-   Lê Cao, K.-A., Martin, P.G.P., Robert-Granie, C. and Besse, P.
    (2009). Sparse canonical methods for biological data integration:
    application to a cross-platform study. *BMC Bioinformatics* *10*:34
-   [Liquet, B., Cao, K.L., Hocini, H. et al. A novel approach for
    biomarker selection and the integration of repeated measures
    experiments from two assays. BMC Bioinformatics 13, 325 (2012).
    https://doi.org/10.1186/1471-2105-13-325](https://doi.org/10.1186/1471-2105-13-325)
-   Mevik, B.-H., Cederkvist, H. R. (2004). Mean Squared Error of
    Prediction (MSEP) Estimates for Principal Component Regression (PCR)
    and Partial Least Squares Regression (PLSR). *Journal of
    Chemometrics* *18*(9), 422-429.
-   Moriyama, M., Hoshida, Y., Otsuka, M., Nishimura, S., Kato, N.,
    Goto, T., Taniguchi, H., Shiratori, Y., Seki, N. and Omata, M.
    (2003). Relevance Network between Chemosensitivity and Transcriptome
    in Human Hepatoma Cells. *Molecular Cancer Therapeutics* *2*,
    199-205.
-   [Rohart F, Gautier B, Singh A, Lê Cao KA (2017) mixOmics: An R
    package for ‘omics feature selection and multiple data integration.
    PLOS Computational Biology 13(11): e1005752.
    https://doi.org/10.1371/journal.pcbi.1005752](https://doi.org/10.1371/journal.pcbi.1005752)
-   Weinstein, J. N., Myers, T. G., O’Connor, P. M., Friend, S. H.,
    Fornace Jr., A. J., Kohn, K. W., Fojo, T., Bates, S. E.,
    Rubinstein, L. V., Anderson, N. L., Buolamwini, J. K., van Osdol, W.
    W., Monks, A. P., Scudiero, D. A., Sausville, E. A., Zaharevitz, D.
    W., Bunow, B., Viswanadhan, V. N., Johnson, G. S., Wittes, R. E. and
    Paull, K. D. (1997). An information-intensive approach to the
    molecular pharmacology of cancer. *Science* *275*, 343-349.

## 7.3 R package compilation

-   <https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html>

-   <https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/>

-   <https://kbroman.org/pkg_primer/pages/data.html>

-   Wickham, H., 2015. R packages: organize, test, document, and share
    your code. " O’Reilly Media, Inc.", <https://r-pkgs.org/>.

-   Xie Y (2021). knitr: A General-Purpose Package for Dynamic Report
    Generation in R. R package version 1.31, <https://yihui.org/knitr/>.

-   Xie Y (2015). Dynamic Documents with R and knitr, 2nd edition.
    Chapman and Hall/CRC, Boca Raton, Florida. ISBN 978-1498716963,
    <https://yihui.org/knitr/>.

-   Xie Y (2014). “knitr: A Comprehensive Tool for Reproducible Research
    in R.” In Stodden V, Leisch F, Peng RD (eds.), Implementing
    Reproducible Computational Research. Chapman and Hall/CRC. ISBN
    978-1466561595,
    <http://www.crcpress.com/product/isbn/9781466561595>.
