#' BPH2819 sepsis study
#'
#' Data from a sepsis study. 12 samples of Staphylococcus aureus strain BPH2819.
#' Divided into two categories, bacteria grown on RPMI media and from human sera.
#' Three omics datasets are featured - metabolome, proteome, transcriptome.
#' Original data had a small amount of missing values.
#' Missing values were filled in with data imputation.
#' For convenience and reproducibility this data object contains the imputed data.
#' Original unimputed data is available for reference.
#'
#' @docType data
#'
#' @usage data(BPH2819)
#'
#' @format A list of named dataframes containing class information and omics data.
#'
#' @keywords datasets
#'
#' @references Mu, A., Klare, W.P., Baines, S.L. et al. Integrative omics identifies conserved and pathogen-specific responses of sepsis-causing bacteria. Nat Commun 14, 1530 (2023). https://doi.org/10.1038/s41467-023-37200-w
#' (\href{https://doi.org/10.1038/s41467-023-37200-w}{DOI})
#'
#' @examples
#' data(BPH2819)
"BPH2819"
