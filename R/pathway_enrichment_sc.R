#' Perform Pathway Enrichment Analysis for Pre-Categorized Differentially Expressed Genes (DEGs)
#' 
#' @description
#' This function performs pathway enrichment analysis for differentially expressed genes (DEGs),
#' which are already categorized into different types (e.g., Dimorphic, Neutral, Sex-specific) via
#' the `categorize_sex_sc` function. The function analyzes their enrichment in KEGG, GO, or Reactome pathways.
#' 
#' @param DEGs_category Data frame containing gene symbols and their corresponding DEG types.
#'                      Must include columns 'DEG_Type' (DEGs categories) and 'Gene_Symbols'.
#' @param enrichment_db Character string specifying the enrichment database to use:
#'                      "KEGG", "GO", or "REACTOME" (default: "KEGG").
#' @param organism Character string representing the organism code. For KEGG enrichment,
#'                 use "hsa" (default). For Reactome enrichment, use "human".
#' @param org_db databse of the organism (e.g: Org.Hs.eg.db)
#' @param pvalueCutoff Numeric value specifying the p-value cutoff for statistical significance
#'                     (default: 0.05).
#' @param qvalueCutoff Numeric value specifying the q-value cutoff for multiple testing correction
#'                     (default: 0.2).
#' 
#' @return A named list of enriched pathways for each DEG category, structured as a data frame.
#' 
#' @details
#' - The input DEGs are already categorized by the `categorize_sex_sc` function.
#' - For GO enrichment, an appropriate OrgDb object (e.g., org.Hs.eg.db for humans) must be available.
#' - For KEGG and Reactome enrichment, gene symbols are first converted to ENTREZ IDs.
#' - Requires the 'clusterProfiler' package for enrichment analysis.
#' - Ensures appropriate error handling for missing genes or database issues.
#' @importFrom clusterProfiler enrichKEGG enrichGO
#' @importFrom ReactomePA enrichPathway
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#'
#' @export

categorized_enrich_sc <- function(DEGs_category, enrichment_db = "KEGG", 
                                organism = "hsa", 
                                org_db = org.Hs.eg.db,
                                pvalueCutoff = 0.05, 
                                qvalueCutoff = 0.2){
  
  
  db <- toupper(enrichment_db)
  
  # Validate input: Check required columns in DEGs_category
  if (!all(c("DEG_Type", "Gene_Symbols") %in% colnames(DEGs_category))) {
    stop("DEGs_category must contain columns 'DEG_Type' and 'Gene_Symbols'.")
  }
  
  # Standardize DEG type names (replace '-' with '_')
  DEGs_category$DEG_Type <- gsub("-", "_", DEGs_category$DEG_Type)
  
  # Identify unique DEG categories
  all_categories <- unique(DEGs_category$DEG_Type)
  results <- list()
  
  # Iterate over each DEG category and perform enrichment analysis
  for (type in all_categories){
    
    # Extract gene symbols for the current DEG category
    degs <- DEGs_category[DEGs_category[["DEG_Type"]] == type, "Gene_Symbols"] %>%
      na.omit()
    
    if (length(degs) > 0) {
      
      if (db == "GO") {
        # GO enrichment analysis
        enrichment_result <- enrichGO(gene = degs, 
                                                  OrgDb = org_db, 
                                                  keyType = "SYMBOL", 
                                                  ont = "ALL", 
                                                  pvalueCutoff = pvalueCutoff, 
                                                  qvalueCutoff = qvalueCutoff)
      } else {
        # Convert gene symbols to ENTREZ IDs
        entrez_ids <- tryCatch({
          bitr(degs, fromType="SYMBOL", toType="ENTREZID", OrgDb=org_db)
        }, error = function(e) {
          warning(paste("ID conversion failed for category:", type))
          return(NULL)
        })
        
        if (is.null(entrez_ids) || nrow(entrez_ids) == 0) {
          warning(paste("No valid ENTREZ IDs found for category:", type))
          next
        }
        
        if (db == "KEGG") {
          # Perform KEGG enrichment analysis
          enrichment_result <- clusterProfiler::enrichKEGG(gene = entrez_ids$ENTREZID, 
                                                           organism = organism, 
                                                           pvalueCutoff = pvalueCutoff, 
                                                           qvalueCutoff = qvalueCutoff)
        } else if (db == "REACTOME") {
          # Adjust organism code for Reactome
          reactome_organism <- ifelse(organism == "hsa", "human", organism)
          
          # Perform Reactome pathway enrichment analysis
          enrichment_result <- enrichPathway(gene = entrez_ids$ENTREZID, 
                                                         organism = reactome_organism, 
                                                         pvalueCutoff = pvalueCutoff, 
                                                         qvalueCutoff = qvalueCutoff)
        } else {
          stop("Unsupported enrichment database. Please choose 'KEGG', 'GO', or 'REACTOME'.")
        }
      }
      
      # Convert enrichment results to a data frame
      
      results[[type]] <- as.data.frame(enrichment_result)
      
    }
  }
  
  # Return the list of enriched pathways for each DEG category
  return(results)
}


#' Convert Data Frame to enrichResult
#'
#' Converts a data frame containing enrichment results into a clusterProfiler enrichResult object.
#' Assumes the data frame has columns: ID, geneID, pvalue, and optionally p.adjust.
#'
#' @param df Data frame containing enrichment results.
#' @param pvalueCutoff Numeric. P-value cutoff for the enrichment object (default: 0.1).
#' @param pAdjustMethod Character string specifying the p-value adjustment method (default: "BH").
#'
#' @return An enrichResult object compatible with clusterProfiler plotting functions.
#' @importFrom methods new
#' @export


convertdf2enr <- function(df, 
                          pvalueCutoff = 0.1, 
                          pAdjustMethod = "BH") {
  
  if (!"p.adjust" %in% colnames(df)) {
    df$p.adjust <- p.adjust(df$pvalue, method = pAdjustMethod)
  }
  
  geneSets <- strsplit(as.character(df[["geneID"]]), "[/]")
  names(geneSets) <- df[["ID"]]
  
  gene <- unique(unlist(geneSets))
  universe <- gene
  
  x <- new("enrichResult",
           result = df,
           pvalueCutoff = pvalueCutoff,
           pAdjustMethod = pAdjustMethod,
           gene = as.character(gene),
           universe = universe,
           geneSets = geneSets,
           organism = "UNKNOWN",
           keytype = "UNKNOWN",
           ontology = "UNKNOWN",
           readable = FALSE)
  x
}