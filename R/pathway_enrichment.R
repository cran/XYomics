#' Improved Pathway Enrichment Analysis
#'
#' Performs pathway enrichment analysis on a set of sex-biased genes using clusterProfiler.
#'
#' @param gene_list A character vector of gene identifiers.
#' @param enrichment_db Character string specifying the database for enrichment.
#'   Options include "KEGG", "GO", and "Reactome". Default is "KEGG".
#' @param organism Character string specifying the organism code (e.g., "hsa" for human).
#' @param org_db database of the organism (e.g: "org.Hs.eg.db")
#' @param pvalueCutoff Numeric. P-value cutoff for enrichment (default: 0.05).
#' @param qvalueCutoff Numeric. Q-value cutoff for enrichment (default: 0.2).
#'
#' @return An enrichment result object.
#'
#' @export
improved_pathway_enrichment <- function(gene_list, 
                                        enrichment_db = "KEGG", 
                                        organism = "hsa", 
                                        org_db = org.Hs.eg.db,
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
    # Note: For GO enrichment, an OrgDb object must be available (e.g., org.Hs.eg.db for human)
    enrichment_result <- enrichGO(gene = gene_list, 
                                                   OrgDb = org_db, 
                                                   keyType = "ENTREZID", 
                                                   ont = "ALL", 
                                                   pvalueCutoff = pvalueCutoff, 
                                                   qvalueCutoff = qvalueCutoff)
  } else if (db == "REACTOME") {
    enrichment_result <- enrichPathway(gene = gene_list, 
                                                        organism = organism, 
                                                        pvalueCutoff = pvalueCutoff, 
                                                        qvalueCutoff = qvalueCutoff)
  } else {
    stop("Unsupported enrichment database. Please choose 'KEGG', 'GO', or 'Reactome'.")
  }
  
  return(enrichment_result)
}
