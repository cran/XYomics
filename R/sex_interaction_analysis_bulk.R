#' Perform Sex-Phenotype Interaction Analysis for Bulk Data (Interaction Term)
#'
#' @description
#' This function performs a formal interaction analysis on bulk expression data
#' to identify genes whose expression is significantly modulated by the interaction
#' between sex and a given phenotype/condition. It uses a linear model with a
#' multiplicative interaction term (`phenotype * sex`).
#'
#' @param x A numeric matrix of expression data (features x samples).
#' @param phenotype A character or factor vector indicating the condition for each sample.
#' @param gender A character or factor vector indicating the sex for each sample.
#' @param phenotype_labels Character vector. Labels for phenotype groups (default: c("WT", "TG")).
#' @param sex_labels Character vector. Labels for sexes (default: c("F", "M")).
#'
#' @return A data frame with differential expression statistics for the
#'         interaction term, including logFC, t-statistic, P-value, and
#'         adjusted P-value.
#'
#' @details This function constructs a design matrix that includes a formal interaction term
#'          between the phenotype and sex (e.g., `~ phenotype * sex`). It then uses `limma`
#'          to test for genes where the effect of the phenotype differs significantly
#'          between sexes. This is a statistically rigorous approach to identify
#'          sex-modulated genes.
#
#' @import limma
#' @importFrom stats model.matrix
#' @importFrom utils tail
#' @export
sex_interaction_analysis_bulk <- function(x, phenotype, gender,
                                          phenotype_labels = c("WT", "TG"),
                                          sex_labels = c("F", "M")) {
  # Input validation
  if (ncol(x) != length(phenotype) || length(phenotype) != length(gender)) {
    stop("The number of samples in x, phenotype, and gender must be equal.")
  }

  # Ensure factors have the correct levels
  phenotype <- factor(phenotype, levels = phenotype_labels)
  gender <- factor(gender, levels = sex_labels)

  # Create design matrix with interaction term
  design <- model.matrix(~ phenotype * gender)

  # Fit linear model
  fit <- lmFit(x, design)
  fit <- eBayes(fit)

  # The interaction term is the last coefficient in this model design
  interaction_coef <- tail(colnames(design), 1)
  
  results <- topTable(fit, coef = interaction_coef, number = Inf, sort.by = "P")

  return(results)
}
