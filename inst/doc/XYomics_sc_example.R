params <-
list(report_date = structure(20461, class = "Date"))

## ----simulate-data, message=FALSE, warning=FALSE------------------------------
library(Seurat)
library(org.Hs.eg.db)
library(stringr)
library(clusterProfiler)
library(igraph)
library(XYomics)
library(ggrepel)
library(ggraph)
library(dplyr)
library(tidyr)

set.seed(123)
gene_sample <- keys(org.Hs.eg.db, keytype = "SYMBOL") %>% 
  as.data.frame() %>% 
  setNames("genes") %>%
  filter(!str_starts(genes, "LOC")) 

gene_sample <- sample(gene_sample$genes, 300, replace = FALSE)

matrices <- lapply(1:10, function(i) {
  batch_effect <- rgamma(300, shape = 2, scale = 1.5)
  m <- matrix(
    ifelse(runif(300*100) < 0.1, 0, rnbinom(300*100, size = 1/0.5, mu = 10 * batch_effect)),
    nrow = 300,
    dimnames = list(gene_sample, paste("Sample", i, ": Cell", 1:100))
  )
  CreateSeuratObject(counts = m, meta.data = data.frame(sample = paste("Sample", i)))
})

sim.seurat <- Reduce(merge, matrices)
sim.seurat@meta.data$sample <- sapply(strsplit(rownames(sim.seurat@meta.data), ":"), "[", 1)

## ----preprocessing, message=FALSE, warning=FALSE, results='hide'--------------
sim.seurat <- NormalizeData(sim.seurat)
sim.seurat <- FindVariableFeatures(sim.seurat, selection.method = "vst", nfeatures = 2000)
sim.seurat <- ScaleData(sim.seurat)
sim.seurat <- RunPCA(sim.seurat)
sim.seurat <- FindNeighbors(sim.seurat, dims = 1:10)
sim.seurat <- FindClusters(sim.seurat)
sim.seurat <- RunUMAP(sim.seurat, dims = 1:10)
DimPlot(sim.seurat)

## ----annotate-data, warning=FALSE---------------------------------------------
# Assign mock cell types
cellTypes <- c("cell type 1", "cell type 2", "cell type 3", "cell type 4", "cell type 5")
sim.seurat@meta.data$cell_type <- sample(cellTypes, nrow(sim.seurat@meta.data), replace = TRUE)
Idents(sim.seurat) <- "cell_type"

# Assign mock status (WT/TG) and sex (M/F)
samples <- sim.seurat@meta.data$sample
sim.seurat@meta.data$status <- ifelse(grepl("1|3|5|8|9", samples), "TG", "WT")
sim.seurat@meta.data$sex <- ifelse(grepl("1|2|4|8|10", samples), "M", "F")


## ----sex-stratified-analysis, message=FALSE, warning=FALSE--------------------
# Run for all cell types
sim.seurat <- JoinLayers(sim.seurat) #this is for the simulated object but not always required for your own object
results <- sex_stratified_analysis_sc(sim.seurat, min_logfc = -Inf)

sex_degs <- lapply(levels(as.factor(sim.seurat@meta.data$cell_type)), function(cell) {
  list(
    male = results$male_DEGs[[cell]],
    female = results$female_DEGs[[cell]]
  )
})

names(sex_degs) <- levels(as.factor(sim.seurat@meta.data$cell_type))
result_categories <- lapply(sex_degs, function(degs) categorize_sex_sc(degs$male, degs$female))

# Example for one cell type
result_one <- result_categories$`cell type 1`
cat("\nTop categorized DEGs for 'cell type 1' (from stratified analysis):\n")
head(result_one)

## ----sex-interaction-analysis, message=FALSE, warning=FALSE-------------------
# Example for one cell type (e.g., "cell type 1")
target_cell <- "cell type 1"
interaction_results_one_cell <- sex_interaction_analysis_sc(sim.seurat, target_cell_type = target_cell)
cat(paste0("\nSummary of Sex-Phenotype Interaction analysis results for '", target_cell, "':\n"))
print(interaction_results_one_cell$summary_stats)
cat(paste0("\nTop genes from Sex-Phenotype Interaction analysis for '", target_cell, "':\n"))
print(head(interaction_results_one_cell$all_results[[1]]))

## ----dot-plots, message=FALSE, warning=FALSE, fig.width=10, fig.height=8------
# Get top gene from each category for 'cell type 1'
top_male_gene <- result_one %>% filter(DEG_Type == "male-specific") %>% top_n(1, -Male_FDR) %>% pull(Gene_Symbols)
top_female_gene <- result_one %>% filter(DEG_Type == "female-specific") %>% top_n(1, -Female_FDR) %>% pull(Gene_Symbols)
top_dimorphic_gene <- result_one %>% filter(DEG_Type == "sex-dimorphic") %>% top_n(1, -Male_FDR) %>% pull(Gene_Symbols)

top_genes <- c(top_male_gene, top_female_gene, top_dimorphic_gene)

# Create dot plot
if (length(top_genes) > 0) {
  sim.seurat$group_plot <- paste(sim.seurat$cell_type, sim.seurat$sex, sim.seurat$status, sep = "_")
  Idents(sim.seurat) <- "group_plot"
  # Seurat v5 requires specifying assay for DotPlot
  DotPlot(sim.seurat, features = top_genes, cols = c("blue", "red"), assay = "RNA") +
    coord_flip() +
    labs(title = "Expression of Top Sex-Specific and Dimorphic Genes") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}

## ----pathway-analysis, message=FALSE, warning=FALSE---------------------------
# Run for all cell types
pathway_category <- lapply(result_categories, function(cell) categorized_enrich_sc(cell, enrichment_db = "REACTOME"))
cat("\nTop pathway results for 'cell type 1':\n")
print(head(pathway_category$`cell type 1`))

# Run for one specific cell type
pathway_category_one <- categorized_enrich_sc(result_one)
cat("\nTop pathway results for 'cell type 1' (re-run separately):\n")
print(head(pathway_category_one))

## ----build-network, warning=FALSE, message=FALSE------------------------------
# Fetch STRING network (can be replaced with a custom network)
# g <- get_string_network(organism = "9606", score_threshold = 900)

# Load a pre-existing network from a file
g <- readRDS(system.file("extdata", "string_example_network.rds", package = "XYomics")) # This should load the 'g' variable


## ----network-analysis, warning=FALSE, message=FALSE---------------------------

### Run for all cell type


network_results <- list()

for (cell_type in names(result_categories)) {
  
  # Extract DEG results for the current cell type
  cell_results <- result_categories[[cell_type]]
  
  # Filter DEGs by type
  male_specific       <- cell_results[cell_results$DEG_Type == "male-specific", ]
  female_specific     <- cell_results[cell_results$DEG_Type == "female-specific", ]
  sex_dimorphic       <- cell_results[cell_results$DEG_Type == "sex-dimorphic", ]
  sex_neutral         <- cell_results[cell_results$DEG_Type == "sex-neutral", ]
  
  # Convert to log-transformed prizes
  male_prizes <- -log10(male_specific$Male_FDR)
  names(male_prizes) <- male_specific$Gene_Symbols
  
  female_prizes <- -log10(female_specific$Female_FDR)
  names(female_prizes) <- female_specific$Gene_Symbols
  
  dimorphic_prizes <- -log10((sex_dimorphic$Male_FDR + sex_dimorphic$Female_FDR) / 2)
  names(dimorphic_prizes) <- sex_dimorphic$Gene_Symbols
  
  neutral_prizes <- -log10((sex_neutral$Male_FDR + sex_neutral$Female_FDR) / 2)
  names(neutral_prizes) <- sex_neutral$Gene_Symbols
  
  # Construct ppis
  male_network     <- construct_ppi_pcsf(g = g, prizes = male_prizes)
  female_network   <- construct_ppi_pcsf(g = g, prizes = female_prizes)
  dimorphic_network <- construct_ppi_pcsf(g = g, prizes = dimorphic_prizes)
  neutral_network   <- construct_ppi_pcsf(g = g, prizes = neutral_prizes)
  
  # Store all results per cell type
  network_results[[cell_type]] <- list(
    male_network      = male_network,
    female_network    = female_network,
    dimorphic_network = dimorphic_network,
    neutral_network   = neutral_network
  )
}

network_results


# Example for one specific cell type

# Create prizes from dimorphic DEGs in 'cell type 1'
dimorphic_specific <- result_one[result_one$DEG_Type == "sex-dimorphic", ]
dimorphic_prizes <- -log10((dimorphic_specific$Male_FDR + dimorphic_specific$Female_FDR) / 2)
names(dimorphic_prizes) <-dimorphic_specific$Gene_Symbols

# Construct the PCSF subnetwork
dimorphic_network <- construct_ppi_pcsf(g = g, prizes = dimorphic_prizes)

## ----visualize-network, warning=FALSE, message=FALSE--------------------------
if (!is.null(dimorphic_network) && igraph::vcount(dimorphic_network) > 0) {
  #Generate network visualization
  plot_network(dimorphic_network, "Cell type 1", "sex-dimorphic", result_one)
} else {
  cat("No network could be constructed for this category.")
}

## ----report, eval=FALSE-------------------------------------------------------
# # This command generates a comprehensive HTML report
# network_results  <- lapply(network_results , function(cell) {
#   Filter(Negate(is.null), cell)
# })
# 
# 
# generate_cat_report(result_categories, pathway_category, network_results)
# 
# 

