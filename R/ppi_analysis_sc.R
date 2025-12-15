#' Construct Protein-protein interaction Network using Prize-Collecting Steiner Forest
#'
#' Constructs a condition-specific gene regulatory network based on differential
#' expression results using the PCSF algorithm.
#'
#' @param g An igraph object representing the base network.
#' @param prizes A named numeric vector of gene scores (prizes). Names must match vertex names in g.
#' @param w Numeric. Edge cost scaling weight. Default is 2.
#' @param b Numeric. Balance between prizes and edge costs. Default is 1.
#' @param mu Numeric. Trade-off parameter for sparsity. Default is 5e-04.
#' @param seed Integer. Random seed. Default is 1.
#' @param min_nodes Integer. Minimum number of nodes in subnetwork. Default is 1.
#' @return An igraph object representing the extracted subnetwork. Returns NULL invisibly 
#'   if no prize genes are present, the subnetwork is too small, or the PCSF algorithm fails.
#' @importFrom igraph V E E<- V<- graph_from_data_frame simplify vcount ecount neighbors induced_subgraph edge_attr_names is_igraph
#'
#' @return An igraph object representing the extracted subnetwork. Returns NULL invisibly if no prize genes are present, the subnetwork is too small, or the PCSF algorithm fails
#'
#' @export
construct_ppi_pcsf <- function(g, prizes, w = 2, b = 1, mu = 5e-04, seed = 1, min_nodes = 1) {
  set.seed(seed)
  
  common_genes <- intersect(names(prizes), V(g)$name)
  if (length(common_genes) == 0) {
    warning("None of the prize genes are present in the network.")
    return(invisible())
  }
  
  prize_nodes <- which(V(g)$name %in% common_genes)
  neighbors <- unique(unlist(lapply(prize_nodes, function(x) neighbors(g, x, mode = "all"))))
  nodes_to_keep <- unique(c(prize_nodes, neighbors))
  
  sub_g <- induced_subgraph(g, nodes_to_keep)
  
  if (vcount(sub_g) < min_nodes) {
    warning("Subnetwork too small: ", vcount(sub_g), " node(s). Skipping PCSF.")
    return(invisible())
  }
  
  prizes_sub <- rep(0, vcount(sub_g))
  names(prizes_sub) <- V(sub_g)$name
  prizes_sub[common_genes] <- prizes[common_genes]
  
  if (!"weight" %in% edge_attr_names(sub_g)) {
    E(sub_g)$weight <- 1
  }
  
  result_network <- tryCatch({
      result <- PCSF(
        ppi = sub_g,
        terminals = prizes_sub,
        w = w,
        b = b,
        mu = mu
      )

   
    
    if (!is.null(result) && is_igraph(result)) {
      V(result)$prize <- prizes_sub[V(result)$name]
      return(result)
    } else {
      warning("PCSF did not return a valid network.")
      return(invisible())
    }
  }, error = function(e) {
    warning("PCSF failed: ", e$message)
    return(invisible())
  })
  
  return(result_network)
}


#' Download and Process STRING Protein-Protein Interaction Network
#'
#' Downloads and processes the STRING protein-protein interaction network,
#' converting it to a simplified igraph object. The function downloads the
#' network from STRING database, filters interactions by confidence score,
#' converts STRING IDs to ENTREZ IDs, and returns the largest connected
#' component as an undirected graph.
#'
#' @param organism Character string specifying the NCBI taxonomy identifier.
#'        Default is "9606" (Homo sapiens).
#' @param score_threshold Numeric value between 0 and 1000 specifying the minimum
#'        combined score threshold for including interactions. Default is 700.
#' @param use_default it will return the default network (9606 and score of 700)
#'
#' @return An igraph object representing the largest connected component of
#'         the filtered STRING network, with the following properties:
#'         \itemize{
#'           \item Undirected edges
#'           \item No self-loops
#'           \item No multiple edges
#'           \item Edge weights (1000 - combined_score)
#'           \item Vertex names as ENTREZ IDs
#'         }
#'
#' @details The function performs the following steps:
#' \enumerate{
#'   \item Downloads protein interactions from STRING database
#'   \item Filters interactions based on combined score
#'   \item Downloads and processes STRING ID to ENTREZ ID mappings
#'   \item Creates an igraph object with filtered interactions
#'   \item Removes self-loops and multiple edges
#'   \item Extracts the largest connected component
#' }
#'
#'
#' @importFrom data.table fread
#' @importFrom igraph V E E<- V<- graph_from_data_frame simplify vcount ecount neighbors induced_subgraph edge_attr_names is_igraph graph.data.frame components
#' @importFrom utils download.file
#' @importFrom stats na.omit
#'
#' @export
get_string_network <- function(organism = "9606", score_threshold = 700, use_default = TRUE) {
  
  if (use_default) {
    return(readRDS(system.file("extdata", "STRING_9606_700.rds", package = "XYomics")))
  }
  
  temp_dir <- tempdir()
  temp_links <- ""
  temp_alias <- ""
  tryCatch({
    
    string_url <- paste0("https://stringdb-static.org/download/protein.links.v11.5/", 
                         organism, ".protein.links.v11.5.txt.gz")
    temp_links <- file.path(temp_dir, "string_links.txt.gz")
    download.file(string_url, temp_links, mode = "wb", method = "curl")
    if (endsWith(tolower(temp_links), ".gz")) {
      con <- gzfile(temp_links, "rt")
      close(con)
      interactions <- data.table::fread(temp_links, header = TRUE, sep = " ", data.table = FALSE)
    }
    interactions <- interactions[interactions$combined_score >= score_threshold, ]
    
    
    alias_url <- paste0("https://stringdb-static.org/download/protein.aliases.v11.5/", 
                        organism, ".protein.aliases.v11.5.txt.gz")
    temp_alias <- file.path(temp_dir, "string_alias.txt.gz")
    download.file(alias_url, temp_alias, mode = "wb", method = "curl")
    
    
    if (endsWith(tolower(temp_alias), ".gz")) {
      con <- gzfile(temp_alias, "rt")
      close(con)
      string_to_entrez <- data.table::fread(temp_alias, header = T, quote = "",
                                            sep = "\t", data.table = FALSE)
      colnames(string_to_entrez) <- c("string_protein_id", "alias_name", "source")
    }
    
    
    
    grep("Symbol",string_to_entrez$source)
    entrez_sources <- c("Ensembl_EntrezGene")
    
    string_to_entrez <- string_to_entrez[string_to_entrez$source %in% entrez_sources, ]
    
    id_map <- string_to_entrez$alias_name
    names(id_map) <- gsub(paste0(organism, "\\."), "", string_to_entrez$string_protein_id)
    
    protein1 <- gsub(paste0(organism, "\\."), "", interactions$protein1)
    
    protein2 <- gsub(paste0(organism, "\\."), "", interactions$protein2)
    
    edges <- data.frame(
      from = id_map[protein1],
      to = id_map[protein2],
      weight = interactions$combined_score/1000,
      stringsAsFactors = FALSE
    )
    
    
    edges <- na.omit(edges)
    
    if (nrow(edges) >1){
      g <- graph.data.frame(edges, directed = FALSE)
      g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb="mean")
      
      components <- components(g)
      largest_component <- which.max(components$csize)
      g <- induced_subgraph(g, 
                            V(g)[components$membership == largest_component])
      
      
      
      return(g)
    }
    
    
  }, finally = {
    unlink(temp_links)
    unlink(temp_alias)
  })
}



#' Plot a Condition-Specific protein-protein interaction network with DEG Annotations
#'
#' Visualizes a gene regulatory or protein–protein interaction network for a given 
#' cell type and differential expression group. Nodes are sized and colored by degree, 
#' and key hub genes are optionally annotated with their barplots of log fold-changes 
#' across sexes.
#'
#' @param g An `igraph` object representing the gene or protein interaction network.
#' @param cell_type Character string. The cell type label used in the plot title.
#' @param DEG_type Character string. The differential expression category to visualize (e.g., `"sex-dimorphic"`).
#' @param result_categories A `data.frame` or tibble containing at least the columns:
#'   `"DEG_Type"`, `"Gene_Symbols"`, `"Male_avg_logFC"`, and `"Female_avg_logFC"`.
#'
#' @return A `ggplot` object representing the visualized network.
#' 
#' @importFrom igraph V E E<- V<- graph_from_data_frame simplify vcount ecount neighbors induced_subgraph edge_attr_names is_igraph degree
#' @import ggraph
#' @import ggplot2
#' @import ggrepel
#' @import dplyr
#' @import tidyr
#' @import scales
#' @import grid
#'
#'
#' @export
#'
#'
plot_network <- function(g, cell_type, DEG_type, result_categories) {
  
  # Adjust node size scaling based on graph size
  min_size <- if (vcount(g) > 30) 1 else 3
  
  # Extract the largest connected component
  comps <- igraph::clusters(g, mode = "weak")
  main_id <- which.max(comps$csize)
  g <- igraph::induced_subgraph(g, vids = V(g)[comps$membership == main_id])
  
  # Compute node degrees
  degree_vals <- degree(g, mode = "all")
  V(g)$degree <- degree_vals
  
  # Identify top hub nodes
  top_n <- min(20, vcount(g))
  top_nodes <- names(sort(degree_vals, decreasing = TRUE))[1:top_n]
  
  # Layout for visualization
  layout <- create_layout(g, layout = "fr")
  
  # Label only top nodes (or all if graph is small)
  layout$label <- if (vcount(g) < 20) {
    layout$name
  } else {
    ifelse(layout$name %in% top_nodes, layout$name, "")
  }
  
  DEG_Type <- NULL
  Gene_Symbols <- NULL
  
  logFC <- NULL
  Female_avg_logFC <- NULL
  Male_avg_logFC <- NULL
  Sex <- ""
  x <- y <- label <- NULL
  # Base network plot
  gplot_net <- ggraph(layout) +
    geom_edge_link(aes(edge_alpha = 0.5), color = "grey", show.legend = FALSE) +
    geom_node_point(aes(size = degree, fill = degree), shape = 21, color = "black", stroke = 0.5,  show.legend = c(size = FALSE, fill = TRUE)) +
    scale_size(range = c(min_size, 10)) +
    scale_fill_gradient(low = "azure", high = "darkcyan",  labels = number_format(accuracy = 1)) +
    geom_label_repel(
      aes(x = x, y = y, label = label),
      color = "black", size = 3, fill = "white", box.padding = 0.5,
      point.padding = 0.5, max.overlaps = 1000
    ) +
    theme_void() +
    labs(
      title = paste("Protein-protein interaction network ", cell_type, DEG_type),
      fill = "Node Degree"
    ) +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 14)) 
  
  # Prepare DEG data for barplots
  df <- result_categories %>% filter(DEG_Type==DEG_type) %>%
    filter(Gene_Symbols %in% V(g)$name) %>%
    dplyr::select(Gene_Symbols, Male_avg_logFC, Female_avg_logFC)
  df <- df %>% pivot_longer( cols = ends_with("FC"),
                             names_to = "Sex",
                             values_to = "logFC",)
  
  # Create a named list of mini barplots indexed by gene
  barplots <- lapply(1:length(V(g)$name), function(i) {
    ggplot(df[df$Gene_Symbols==V(g)$name[i],], aes(Sex, logFC, fill = Sex)) +
      geom_col(show.legend = FALSE) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +  # zero line
      theme_void() +
      theme(
        plot.margin = margin(0, 0, 0, 0),
        panel.background = element_rect(fill = "transparent", color = NA)
      )
  })
  
  # Scale for mini-plots relative to graph size
  plot_scale <- min(vcount(g), 100) * 0.015
  
  # Embed mini barplots for hub nodes that have DEG data
  for (i in seq_len(nrow(layout))) {
    if (layout$name[i] %in% intersect(top_nodes, df$Gene_Symbols)){
      gplot_net <- gplot_net + annotation_custom(
        grob = ggplotGrob(barplots[[i]]),
        xmin = layout$x[i] - plot_scale ,
        xmax = layout$x[i] + plot_scale ,
        ymin = layout$y[i] -  plot_scale ,
        ymax = layout$y[i] +  plot_scale 
      )
    }
    
  }
  return(gplot_net)
}
