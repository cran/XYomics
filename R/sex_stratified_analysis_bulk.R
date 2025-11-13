#' Perform differential expression analysis within each sex
#'
#' This function identifies differentially expressed genes between conditions separately for each sex
#' using a linear modeling approach.
#'
#' @param x A numeric matrix of expression data (features × samples).
#' @param phenotype A vector indicating condition labels for each sample.
#' @param gender A vector indicating gender for each sample. Labels must start with "f" (female) and "m" (male).
#' @param analysis_type Character. Type of analysis to perform: "dimorphic" (difference in differences),
#'        "female" (female condition effect), or "male" (male condition effect). Default is "dimorphic".
#'
#' @return A data frame with differential expression statistics including logFC, AveExpr, t-statistic,
#'         P-value, and adjusted P-value.
#'
#' @details This function performs differential expression analysis within each sex separately.
#'          For male analysis, it compares conditions within males.
#'          For female analysis, it compares conditions within females.
#'          For dimorphic analysis, it tests for difference in condition effects between sexes.
#'          Note: To identify truly sex-specific genes, use the output of this function as input
#'          for identify_sex_specific_genes().
#'
#' @importFrom stats model.matrix
#'
#' @export
sex_stratified_analysis_bulk <- function(x, phenotype, gender, analysis_type = c("male", "female")) {
  # Input validation
  if (ncol(x) != length(phenotype)) {
    stop("The number of columns in x must match the length of phenotype vector")
  }
  if (length(phenotype) != length(gender)) {
    stop("The length of phenotype and gender vectors must be equal")
  }
  
  # Validate and standardize gender labels
  genderlabs <- sort(unique(substr(gender, 1, 1)))
  if (!(all(c("f", "m") %in% genderlabs))) {
    stop("Gender labels must start with 'f' for female and 'm' for male")
  }
  
  # Match argument and select appropriate samples
  analysis_type <- match.arg(analysis_type)
  if (analysis_type == "male") {
    selected_samples <- which(substr(gender, 1, 1) == "m")
  } else {
    selected_samples <- which(substr(gender, 1, 1) == "f")
  }
  
  # Subset data for selected sex
  x_subset <- x[, selected_samples]
  phenotype_subset <- phenotype[selected_samples]
  
  # Create design matrix for the selected sex
  # Convert phenotype to factor ensuring first level is reference
  phenotype_factor <- factor(phenotype_subset, levels = unique(phenotype_subset))
  design <- model.matrix(~ phenotype_factor)
  
  # Fit linear model
  fit <- lmFit(x_subset, design)
  fit <- eBayes(fit)
  
  # Get results
  results <- topTable(fit, coef = 2, number = Inf)
  return(results)
}


#' Identify sex-specific and sex-dimorphic genes
#'
#' This function identifies truly sex-specific and sex-dimorphic genes by analyzing
#' differential expression results from both sexes.
#'
#' @param male_results Data frame of differential expression results for males (from differential_expression).
#' @param female_results Data frame of differential expression results for females (from differential_expression).
#' @param target_fdr Numeric. FDR threshold for significant differential expression (default: 0.05).
#' @param exclude_fdr Numeric. FDR threshold for excluding effects in the opposite sex (default: 0.5).
#'
#' @return A data frame with identified genes categorized as:
#'         - male-specific: significant in males, not significant in females
#'         - female-specific: significant in females, not significant in males
#'         - sex-dimorphic: significant in both sexes with opposite effects
#'         - sex-shared: significant in both sexes with same direction
#'         Including columns for gene IDs, logFC values, and FDR values for both sexes.
#'
#' @details This function implements a two-step approach to identify sex-specific effects:
#'          1. Identifies genes significantly affected in one sex (target_fdr)
#'          2. Confirms lack of effect in the other sex (exclude_fdr)
#'          Additionally identifies genes with opposite (dimorphic) or same (shared) effects in both sexes.
#'
#' @export
identify_sex_specific_genes <- function(male_results, female_results, 
                                      target_fdr = 0.05, exclude_fdr = 0.5) {
  # Identify male-specific genes
  male_deg <- rownames(male_results)[which(male_results$adj.P.Val < target_fdr)]
  male_specific <- male_deg[which(female_results[match(male_deg, 
                                rownames(female_results)),]$adj.P.Val > exclude_fdr)]
  
  # Identify female-specific genes
  female_deg <- rownames(female_results)[which(female_results$adj.P.Val < target_fdr)]
  female_specific <- female_deg[which(male_results[match(female_deg, 
                                  rownames(male_results)),]$adj.P.Val > exclude_fdr)]
  
  # Identify shared and dimorphic genes
  shared_genes <- intersect(male_deg, female_deg)
  dimorphic_genes <- NULL
  shared_same_dir <- NULL
  
  if(length(shared_genes) > 0) {
    # Compare direction of effects
    male_fc <- male_results[match(shared_genes, rownames(male_results)),]$logFC
    female_fc <- female_results[match(shared_genes, rownames(female_results)),]$logFC
    
    # Different directions = dimorphic
    dimorphic_genes <- shared_genes[sign(male_fc) != sign(female_fc)]
    # Same direction = shared
    shared_same_dir <- shared_genes[sign(male_fc) == sign(female_fc)]
  }
  
  # Create results data frame
  result <- data.frame(
    "gene_type" = c(rep("male-specific", length(male_specific)),
                    rep("female-specific", length(female_specific)),
                    rep("sex-dimorphic", length(dimorphic_genes)),
                    rep("sex-shared", length(shared_same_dir))),
    "gene_id" = c(male_specific, female_specific, dimorphic_genes, shared_same_dir),
    stringsAsFactors = FALSE
  )
  
  # Add logFC and FDR values for both sexes
  result$male_logFC <- male_results[match(result$gene_id, rownames(male_results)),]$logFC
  result$female_logFC <- female_results[match(result$gene_id, rownames(female_results)),]$logFC
  result$male_FDR <- male_results[match(result$gene_id, rownames(male_results)),]$adj.P.Val
  result$female_FDR <- female_results[match(result$gene_id, rownames(female_results)),]$adj.P.Val
  
  return(result)
}
