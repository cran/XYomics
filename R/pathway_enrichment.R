#' Improved Pathway Enrichment Analysis
#'
#' Performs pathway enrichment analysis on a set of sex-biased genes using clusterProfiler.
#'
#' @param gene_list A character vector of gene identifiers.
#' @param enrichment_db Character string specifying the database for enrichment.
#'   Options include "KEGG", "GO", and "Reactome". Default is "KEGG".
#' @param organism Character string specifying the organism code (e.g., "hsa" for human).
#' @param org_db database of the organism (e.g: "")
#' @param pvalueCutoff Numeric. P-value cutoff for enrichment (default: 0.05).
#' @param qvalueCutoff Numeric. Q-value cutoff for enrichment (default: 0.2).
#'
#' @return An enrichment result object.
#'
#' @export
improved_pathway_enrichment <- function(gene_list, 
                                        enrichment_db = "KEGG", 
                                        organism = "hsa", 
                                        org_db = org.Hs.eg.db::org.Hs.eg.db,
                                        pvalueCutoff = 0.05, 
                                        qvalueCutoff = 0.2) {

  
  enrichment_result <- NULL
  db <- toupper(enrichment_db)
  
  if (db == "KEGG") {
    enrichment_result <- clusterProfiler::enrichKEGG(gene = gene_list, 
                                                     organism = organism, 
                                                     pvalueCutoff = pvalueCutoff, 
                                                     qvalueCutoff = qvalueCutoff)
  } else if (db == "GO") {
    # Note: For GO enrichment, an OrgDb object must be available (e.g.,  for human)
    enrichment_result <- enrichGO(gene = gene_list, 
                                                   OrgDb = org_db, 
                                                   keyType = "ENTREZID", 
                                                   ont = "ALL", 
                                                   pvalueCutoff = pvalueCutoff, 
                                                   qvalueCutoff = qvalueCutoff)
  } else if (db == "REACTOME") {
    # Adjust organism code for Reactome
    reactome_organism <- ifelse(organism == "hsa", "hsapiens", organism)
    
    # Perform Reactome pathway enrichment analysis
    enrichment_result <- tryCatch(
      gost(query=gene_list, sources="REAC", organism=reactome_organism,
           user_threshold=pvalueCutoff, evcodes=TRUE, correction_method="fdr")$result,
      error = function(e) NULL
    )
    
    if (!is.null(enrichment_result) && nrow(enrichment_result) > 0) {
      enrichment_result$p.adjust <- enrichment_result$p_value
      enrichment_result$qvlaue <- enrichment_result$p_value 
      enrichment_result$GeneRatio <- enrichment_result$intersection_size / enrichment_result$term_size
      enrichment_result$category <- enrichment_result$term_name
      enrichment_result$geneID <- gsub(",", "/", enrichment_result$intersection)
      cols_to_keep <- c("category", "p_value", "p.adjust", "qvlaue", "term_size", "term_id", "term_name", "intersection_size" , "intersection", "GeneRatio", "geneID")
      enrichment_result <- enrichment_result[, cols_to_keep, drop = FALSE]
    } else {
      enrichment_result <- data.frame(matrix(ncol = 11, nrow = 0))
      colnames(enrichment_result) <- c("category", "pvalue", "p.adjust", "qvlaue", "term_size", "term_id", "term_name", "intersection_size" , "intersection", "GeneRatio", "geneID")
    }
    colnames(enrichment_result)[colnames(enrichment_result) == "term_name"] <- "Description"
    colnames(enrichment_result)[colnames(enrichment_result) == "intersection_size"] <- "Count"
    colnames(enrichment_result)[colnames(enrichment_result) == "p_value"] <- "pvalue"
    colnames(enrichment_result)[colnames(enrichment_result) == "term_id"] <- "ID"
  } else {
    stop("Unsupported enrichment database. Please choose 'KEGG', 'GO', or 'Reactome'.")
  }
  
  return(enrichment_result)
}
