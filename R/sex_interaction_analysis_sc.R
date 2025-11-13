#' Perform Sex-Phenotype Interaction Analysis for Single-Cell Data
#'
#' @description Performs differential difference analysis for a given cell type
#' to identify genes modulated by sex-phenotype interactions using limma.
#'
#' @param seurat_obj A Seurat object.
#' @param target_cell_type Character. Cell type to analyze.
#' @param sex_col Character. Column name for sex (default "sex").
#' @param phenotype_col Character. Column name for phenotype (default "status").
#' @param celltype_col Character. Column name for cell type (default "cell_type").
#' @param min_logfc Numeric. Minimum absolute log fold change (default 0.25).
#' @param fdr_threshold Numeric. FDR threshold for significance (default 0.05).
#' @param sex_labels Character vector of sex labels (default c("F","M")).
#' @param phenotype_labels Character vector of phenotype groups (default c("WT","TG")).
#'
#' @return A list with complete DE results, significant results, and summary statistics.
#'
#'
#' @importFrom Seurat DefaultAssay GetAssayData
#' @importFrom edgeR DGEList filterByExpr calcNormFactors
#' @importFrom limma voom lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom stats model.matrix
#' @importFrom Seurat DefaultAssay
#' @import SeuratObject
#' 
#' @export
sex_interaction_analysis_sc <- function(seurat_obj,
                                        target_cell_type,
                                        sex_col = "sex",
                                        phenotype_col = "status",
                                        celltype_col = "cell_type",
                                        min_logfc = 0.25,
                                        fdr_threshold = 0.05,
                                        sex_labels = c("F", "M"),
                                        phenotype_labels = c("WT", "TG")) {

  if (!all(c(sex_col, phenotype_col, celltype_col) %in% names(seurat_obj@meta.data))) {
    stop("Required metadata columns not found in Seurat object")
  }
  
  if (!(target_cell_type %in% unique(seurat_obj@meta.data[[celltype_col]]))) {
    stop(paste0("Target cell type '", target_cell_type, "' not found in '", celltype_col, "' column."))
  }
  

  # Subset metadata for the target cell type
  cell_subset_metadata <- seurat_obj@meta.data[seurat_obj@meta.data[[celltype_col]] == target_cell_type, ]
  
  # Ensure factors have correct levels
  cell_subset_metadata[[sex_col]] <- factor(cell_subset_metadata[[sex_col]], levels = sex_labels)
  cell_subset_metadata[[phenotype_col]] <- factor(cell_subset_metadata[[phenotype_col]], levels = phenotype_labels)

  # Create a combined group factor for differential difference analysis
  cell_subset_metadata$Group <- factor(paste(cell_subset_metadata[[phenotype_col]], cell_subset_metadata[[sex_col]], sep = "_"))

  # Ensure all four groups are present for interaction model
  required_groups <- c(
    paste0(phenotype_labels[1], "_", sex_labels[1]),
    paste0(phenotype_labels[2], "_", sex_labels[1]),
    paste0(phenotype_labels[1], "_", sex_labels[2]),
    paste0(phenotype_labels[2], "_", sex_labels[2])
  )

  if (!all(required_groups %in% levels(cell_subset_metadata$Group))) {
    missing_groups <- required_groups[!required_groups %in% levels(cell_subset_metadata$Group)]
    stop(paste("Insufficient groups for interaction analysis in cell type", target_cell_type, ". Missing:", paste(missing_groups, collapse = ", ")))
  }

  # Subset Seurat object for the target cell type
  cell_subset_seurat <- subset(seurat_obj, cells = rownames(cell_subset_metadata))

  # Define design matrix
  design <- model.matrix(~ 0 + Group, data = cell_subset_metadata)
  colnames(design) <- gsub("Group", "", colnames(design)) # Clean up column names

  DefaultAssay(cell_subset_seurat) <- "RNA"
  dge <- edgeR::DGEList(counts = as.matrix(Seurat::GetAssayData(cell_subset_seurat, slot = "counts")))
  keep <- edgeR::filterByExpr(dge, design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge, method = "TMMwsp")

  vm <- limma::voom(dge, design, plot = FALSE)
  fit <- limma::lmFit(vm, design)

  # Define contrasts for differential difference
  # Effect of phenotype in female: (TG_F - WT_F)
  # Effect of phenotype in male: (TG_M - WT_M)
  # Interaction: (TG_M - WT_M) - (TG_F - WT_F)
  
  # Dynamically create contrast strings
  group_wt_f <- paste0(phenotype_labels[1], "_", sex_labels[1])
  group_tg_f <- paste0(phenotype_labels[2], "_", sex_labels[1])
  group_wt_m <- paste0(phenotype_labels[1], "_", sex_labels[2])
  group_tg_m <- paste0(phenotype_labels[2], "_", sex_labels[2])

  contrast_female_effect <- paste0(group_tg_f, " - ", group_wt_f)
  contrast_male_effect <- paste0(group_tg_m, " - ", group_wt_m)
  
  # Construct the full contrast string for makeContrasts
  interaction_contrast_definition <- paste0("Interaction = (", contrast_male_effect, ") - (", contrast_female_effect, ")")

  # Construct the complete makeContrasts call as a string
  make_contrasts_call_str <- paste0("limma::makeContrasts(", interaction_contrast_definition, ", levels = design)")

  # Evaluate the string as an R expression
  contrasts_matrix <- eval(parse(text = make_contrasts_call_str))

  fit2 <- limma::contrasts.fit(fit, contrasts_matrix)
  fit2 <- limma::eBayes(fit2)

  res <- limma::topTable(fit2, coef = "Interaction", number = Inf, adjust.method = "BH")
  res$cell_type <- target_cell_type # Use target_cell_type here
  res$gene <- rownames(res)
  # Return results in a list structure similar to the original
  return(list(
    all_results = list(res), # Wrap in list to maintain consistency with original output structure
    sig_results = list(subset(res, res$adj.P.Val < fdr_threshold & abs(res$logFC) > min_logfc)),
    summary_stats = data.frame(
      cell_type = target_cell_type,
      n_total_genes = nrow(res),
      n_sig_genes = nrow(subset(res, res$adj.P.Val < fdr_threshold & abs(res$logFC) > min_logfc)),
      stringsAsFactors = FALSE
    )
  ))
}
