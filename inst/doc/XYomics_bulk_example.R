params <-
list(report_date = structure(20461, class = "Date"))

## ----simulate-data, message=FALSE, warning=FALSE------------------------------
library(XYomics)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

simulate_omics_data <- function(n_samples = 100, n_genes = 500) {
  set.seed(123)
  all_entrez <- keys(org.Hs.eg.db, keytype = "ENTREZID")
  selected_entrez <- sample(all_entrez, n_genes)
  
  sex <- rep(c("m", "f"), each = n_samples / 2)
  phenotype <- rep(rep(c("Control", "Disease"), each = n_samples / 4), 2)
  
  expression_data <- matrix(rnorm(n_samples * n_genes, mean = 10, sd = 0.1), 
                            nrow = n_genes, ncol = n_samples)
  
  male_disease <- which(sex == "m" & phenotype == "Disease")
  female_disease <- which(sex == "f" & phenotype == "Disease")
  
  expression_data[1:25, male_disease] <- expression_data[1:25, male_disease] + 8
  expression_data[26:50, female_disease] <- expression_data[26:50, female_disease] + 8
  expression_data[51:75, male_disease] <- expression_data[51:75, male_disease] + 6
  expression_data[51:75, female_disease] <- expression_data[51:75, female_disease] - 6
  
  colnames(expression_data) <- paste0("Sample", 1:n_samples)
  rownames(expression_data) <- selected_entrez
  
  return(list(expression_data = expression_data, phenotype = phenotype, sex = sex))
}

simulated_data <- simulate_omics_data()
omics_data <- simulated_data$expression_data
phenotype <- simulated_data$phenotype
sex <- simulated_data$sex

cat("Dimensions of expression data:", dim(omics_data), "\n")

## ----stratified-analysis------------------------------------------------------
male_deg <- sex_stratified_analysis_bulk(
  x = omics_data,
  phenotype = phenotype,
  gender = sex,
  analysis_type = "male"
)

female_deg <- sex_stratified_analysis_bulk(
  x = omics_data,
  phenotype = phenotype,
  gender = sex,
  analysis_type = "female"
)

cat("Top differentially expressed genes in males (stratified analysis):\n")
print(head(male_deg))

## ----interaction-term-analysis-bulk-------------------------------------------
interaction_term_results_bulk <- sex_interaction_analysis_bulk(
  x = omics_data,
  phenotype = phenotype,
  gender = sex,
  phenotype_labels = c("Control", "Disease"),
  sex_labels = c("f", "m")
)

cat("Top genes from Interaction Term analysis (bulk):\n")
print(head(interaction_term_results_bulk))

## ----sex-specific-analysis----------------------------------------------------
sex_specific_results <- identify_sex_specific_genes(
  male_results = male_deg,
  female_results = female_deg,
  target_fdr = 0.05,
  exclude_fdr = 0.5
)

cat("Number of genes in each category (from stratified analysis):\n")
print(table(sex_specific_results$gene_type))

## ----visualization------------------------------------------------------------
# Male-specific gene
top_male_gene <- sex_specific_results %>%
  filter(gene_type == "male-specific") %>%
  arrange(male_FDR) %>%
  pull(gene_id) %>%
  head(1)

generate_boxplot(
  x = omics_data,
  index = top_male_gene,
  phenotype = phenotype,
  gender = sex,
  title = paste("Expression of Top Male-Specific Gene:", top_male_gene)
)

# Female-specific gene
top_female_gene <- sex_specific_results %>%
  filter(gene_type == "female-specific") %>%
  arrange(female_FDR) %>%
  pull(gene_id) %>%
  head(1)

generate_boxplot(
  x = omics_data,
  index = top_female_gene,
  phenotype = phenotype,
  gender = sex,
  title = paste("Expression of Top Female-Specific Gene:", top_female_gene)
)

# Sex-dimorphic gene
top_dimorphic_gene <- sex_specific_results %>%
  filter(gene_type == "sex-dimorphic") %>%
  arrange(male_FDR) %>%
  pull(gene_id) %>%
  head(1)

generate_boxplot(
  x = omics_data,
  index = top_dimorphic_gene,
  phenotype = phenotype,
  gender = sex,
  title = paste("Expression of Top Sex-Dimorphic Gene:", top_dimorphic_gene)
)

## ----pathway-analysis---------------------------------------------------------
# GO analysis for female-specific genes
female_specific_genes <- sex_specific_results %>%
  filter(gene_type == "female-specific") %>%
  pull(gene_id)

cat("\nPerforming GO enrichment analysis for female-specific genes...\n")
go_results <- improved_pathway_enrichment(
  gene_list = female_specific_genes,
  enrichment_db = "GO"
)
if (!is.null(go_results) && nrow(go_results) > 0) {
  cat("Top GO Terms:\n")
  print(head(as.data.frame(go_results)))
} else {
  cat("No significant GO terms found.\n")
}

# GO analysis for male-specific genes
male_specific_genes <- sex_specific_results %>%
  filter(gene_type == "male-specific") %>%
  pull(gene_id)

cat("\nPerforming GO pathway analysis for male-specific genes...\n")
go_results_male <- improved_pathway_enrichment(
  gene_list = male_specific_genes,
  enrichment_db = "GO"
)
if (!is.null(go_results_male) && nrow(go_results_male) > 0) {
  cat("Top GO Pathways:\n")
  print(head(as.data.frame(go_results_male)))
} else {
  cat("No significant GO pathways found.\n")
}

## ----get-string-network-------------------------------------------------------
# In a real analysis, you would run:
# string_network <- get_string_network(organism = "9606", score_threshold = 700)

# For this vignette, we create a small dummy network
set.seed(123)
dummy_nodes <- unique(c(male_specific_genes, female_specific_genes))
if(length(dummy_nodes) > 100) dummy_nodes <- sample(dummy_nodes, 100)
dummy_edges <- data.frame(
  from = sample(dummy_nodes, 50, replace = TRUE),
  to = sample(dummy_nodes, 50, replace = TRUE)
)
string_network <- igraph::graph_from_data_frame(dummy_edges, directed = FALSE)
string_network <- igraph::simplify(string_network)
cat("Using a dummy network for demonstration purposes.\n")

## ----construct-pcsf-network---------------------------------------------------
# Use male DE results to define prizes
prizes <- -log10(male_deg$P.Value)
names(prizes) <- rownames(male_deg)
prizes <- prizes[is.finite(prizes)]

# Construct the network
ppi <- construct_ppi_pcsf(
  g = string_network,
  prizes = prizes,
  w = 2,
  b = 1,
  mu = 5e-04
)

cat("Constructed PPI with", igraph::vcount(ppi), "nodes and", igraph::ecount(ppi), "edges.\n")

## ----visualize-network--------------------------------------------------------
visualize_network(
  g = ppi,
  female_res = female_deg,
  male_res = male_deg,
  vertex.size = 8,
  vertex.label.cex = 0.7,
  main = "Sex-Specific Protein-Protein Interaction Network"
)

