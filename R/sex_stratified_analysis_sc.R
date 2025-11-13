#' Compute sex-specific differentially expressed genes (DEGs)
#'
#' @description Identifies differentially expressed genes (DEGs) separately for male and female samples
#' within different cell types using the Seurat package. Compares gene expression between control
#' and perturbed groups in each sex.
#'
#' @param seurat_obj Seurat object containing the single-cell data.
#' @param sex_column Character. Column name in metadata for sex (default "sex").
#' @param phenotype_column Character. Column name in metadata for phenotype (default "status").
#' @param celltype_column Character. Column name in metadata for cell type (default "cell_type").
#' @param sex_labels_vector Character vector of sex labels (default c("F","M")).
#' @param phenotype_labels_vector Character vector of phenotype groups (default c("WT","TG")).
#' @param method Character. Statistical test to use for differential expression (default "wilcox").
#' @param min_logfc Numeric. Minimum absolute log fold change threshold (default 0.25).
#'
#' @return A list with male and female DEGs results.
#'
#'
#' @importFrom Seurat DefaultAssay GetAssayData FindMarkers Idents
#' @importFrom Seurat Idents "Idents<-"
#' @importFrom Seurat DefaultAssay "DefaultAssay<-"
#' @importFrom edgeR DGEList filterByExpr calcNormFactors
#' @importFrom limma voom lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom stats p.adjust
#' @import SeuratObject
#'
#' @export
sex_stratified_analysis_sc <- function(seurat_obj,
                                      sex_column = "sex",
                                      phenotype_column = "status",
                                      celltype_column = "cell_type",
                                      sex_labels_vector = c("F", "M"),
                                      min_logfc = 0.25,
                                      phenotype_labels_vector = c("WT", "TG"),
                                      method="wilcox") {

  # Input Validation
  required_columns <- c(sex_column, phenotype_column, celltype_column)
  if (!all(required_columns %in% names(seurat_obj@meta.data))) {
    missing_cols <- required_columns[!required_columns %in% names(seurat_obj@meta.data)]
    stop(paste("Required metadata columns not found in Seurat object:", paste(missing_cols, collapse = ", ")))
  }

  # Prepare metadata
  metadata <- seurat_obj@meta.data
  metadata[[sex_column]] <- factor(metadata[[sex_column]], levels = sex_labels_vector)
  metadata[[phenotype_column]] <- factor(metadata[[phenotype_column]], levels = phenotype_labels_vector)
  metadata[[celltype_column]] <- as.factor(metadata[[celltype_column]])

  # Store results
  male_results <- list()
  female_results <- list()
  # Loop over each cell type
  cell_types <- unique(metadata[[celltype_column]])
  for (cell in cell_types) {
    message(paste("Processing cell type:", cell))
    Idents(seurat_obj) <- celltype_column
    cell_type_subset <- subset(seurat_obj, idents = cell)

    for (sex in sex_labels_vector) {
      Idents(cell_type_subset) <- sex_column
      sex_subset <- subset(cell_type_subset, idents=sex)

      if (length(which(sex_subset[[phenotype_column]] == phenotype_labels_vector[1])) < 5 ||
          length(which(sex_subset[[phenotype_column]] == phenotype_labels_vector[2])) < 5) {
        next
      }

      tryCatch({
        markers <- FindMarkers(
          object = sex_subset,
          group.by = phenotype_column,
          ident.1 = phenotype_labels_vector[2],
          ident.2 = phenotype_labels_vector[1],
          only.pos=FALSE,
          logfc.threshold = min_logfc,
          test.use = method
        )
        markers$p_val_adj <- p.adjust(markers$p_val, method = "BH")
        markers$cell_type <- cell
        markers$sex <- sex
        markers$gene <- rownames(markers)
        markers$direction <- ifelse(markers$avg_log2FC < 0, "down", "up")

        if (sex == "M") {
          male_results[[cell]] <- markers
        } else {
          female_results[[cell]] <- markers
        }
      }, error = function(e) {
        message(paste("Error processing cell type", cell, "and sex", sex, ":", e$message))
      })
    }
  }

  return(list(
    male_DEGs = male_results,
    female_DEGs = female_results
  ))
}


#' Compute sex-specific differentially expressed genes (DEGs) per category
#' @description
#' Identifies male-specific, female-specific, sex-dimorphic, and sex-neutral DEGs from differential expression results.
#'
#' @param male_degs Data frame containing male differential expression results from one specific cell-type or bulk dataset.
#' @param female_degs Data frame containing female differential expression results from one specific cell-type or bulk dataset.
#' @param target_fdr Numeric. FDR threshold for significance.
#' @param exclude_pval Numeric. P-value threshold for excluding genes in opposite sex.
#' @param min_abs_logfc Numeric. Minimum absolute log2 fold change threshold.
#'
#' @return Data frame containing categorized DEGs with associated statistics.
#'
#' @importFrom utils download.file
#' @importFrom stats na.omit
#' @importFrom clusterProfiler bitr
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#'
#' @export
#'

categorize_sex_sc <- function(male_degs,
                               female_degs,
                               target_fdr = 0.05,
                               exclude_pval = 0.5,
                               min_abs_logfc = 0.25) {

  # Helper function to extract logFC and FDR values
  get_deg_stats <- function(genes, deg_data) {
    matched_idx <- match(genes, rownames(deg_data))
    list(
      logFC = deg_data[matched_idx, "avg_log2FC"],
      fdr = deg_data[matched_idx, "p_val_adj"]
    )
  }

  # Identify male-specific genes
  male_candidates <- rownames(male_degs)[which(male_degs$p_val_adj < target_fdr)]
  male_specific <- male_candidates[which(female_degs[match(male_candidates, rownames(female_degs)),]$p_val > exclude_pval)]
  male_specific <- male_specific[which(abs(male_degs[match(male_specific, rownames(male_degs)),]$avg_log2FC) > min_abs_logfc)]


  # Identify female-specific genes
  female_candidates <- rownames(female_degs)[which(female_degs$p_val_adj < target_fdr)]
  female_specific <- female_candidates[which(male_degs[match(female_candidates, rownames(male_degs)),]$p_val > exclude_pval)]
  female_specific <- female_specific[which(abs(female_degs[match(female_specific, rownames(female_degs)),]$avg_log2FC) > min_abs_logfc)]


  # Identify shared and dimorphic genes
  dimorphic_genes <- NULL
  shared_genes <- NULL
  intersecting_genes <- intersect(male_candidates, female_candidates)

  if (length(intersecting_genes) > 0) {
    logfc_male <- male_degs[match(intersecting_genes, rownames(male_degs)),]$avg_log2FC
    logfc_female <- female_degs[match(intersecting_genes, rownames(female_degs)),]$avg_log2FC

    dimorphic_indices <- intersect(which(sign(logfc_male) != sign(logfc_female)),
                                   intersect(which(abs(logfc_male) > min_abs_logfc),
                                             which(abs(logfc_female) > min_abs_logfc)))
    dimorphic_genes <- intersecting_genes[dimorphic_indices]

    shared_indices <- intersect(which(sign(logfc_male) == sign(logfc_female)),
                                intersect(which(abs(logfc_male) > min_abs_logfc),
                                          which(abs(logfc_female) > min_abs_logfc)))
    shared_genes <- intersecting_genes[shared_indices]
  }

  # Compile results into a data frame
  all_categories <- list(
    "male-specific" = male_specific,
    "female-specific" = female_specific,
    "sex-dimorphic" = dimorphic_genes,
    "sex-neutral" = shared_genes
  )

  results <- data.frame(
    DEG_Type = rep(names(all_categories), sapply(all_categories, length)),
    Gene_Symbols = unlist(all_categories),
    stringsAsFactors = FALSE
  )

  # Add male statistics
  male_stats <- get_deg_stats(results$Gene_Symbols, male_degs)
  results$Male_avg_logFC <- male_stats$logFC
  results$Male_FDR <- male_stats$fdr

  # Add female statistics
  female_stats <- get_deg_stats(results$Gene_Symbols, female_degs)
  results$Female_avg_logFC <- female_stats$logFC
  results$Female_FDR <- female_stats$fdr

  return(results)
}
