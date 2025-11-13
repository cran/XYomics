#' Generate Boxplots for Expression Data
#'
#' Creates boxplots to visualize expression differences across conditions and genders.
#'
#' @param x Expression data matrix.
#' @param index Numeric vector indicating which features (rows) to plot.
#' @param phenotype Vector of phenotype labels.
#' @param gender Vector of gender labels.
#' @param title Title for the plot.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#'
#' @return A boxplot is generated.
#' 
#' @importFrom graphics boxplot
#' 
#' @export
generate_boxplot <- function(x, index, phenotype, gender, 
                             title = "Expression Boxplot", 
                             xlab = "Conditions", 
                             ylab = "Expression Level") {
  boxplot(x[index, ] ~ factor(paste(phenotype, gender, sep = "_")),
          main = title, xlab = xlab, ylab = ylab,
          col = c("lightblue", "darkblue", "red", "darkred"))
}

#' Visualize Gene Regulatory Network with Pie Charts
#'
#' Plots a network with nodes represented by pie charts that display male and female effects.
#'
#' @param g An igraph network object.
#' @param female_res Differential expression results for females.
#' @param male_res Differential expression results for males.
#' @param vertex.size Size of the network nodes.
#' @param vertex.label.cex Text size for vertex labels.
#' @param ... Additional graphical parameters.
#'
#' @return The modified igraph object with visualization attributes.
#' 
#'
#' @export
visualize_network <- function(g, female_res, male_res, vertex.size = 5, vertex.label.cex = 0.8, ...) {
  # Map log fold changes from results to network vertices
  female_fc <- female_res$logFC
  male_fc <- male_res$logFC
  names(female_fc) <- rownames(female_res)
  names(male_fc) <- rownames(male_res)
  
  map_f <- match(V(g)$name, names(female_fc))
  map_m <- match(V(g)$name, names(male_fc))
  
  node_f <- ifelse(is.na(map_f), 0, female_fc[map_f])
  node_m <- ifelse(is.na(map_m), 0, male_fc[map_m])
  
  # Create color gradients (simple example)
  node_colors <- ifelse(node_f > node_m, "pink", "lightblue")
  
  plot(g, vertex.size = vertex.size, vertex.label.cex = vertex.label.cex,
       vertex.color = node_colors, ...)
  return(g)
}

#' Generate a Comprehensive Analysis Report
#'
#' Creates an integrated HTML report combining differential expression results,
#' enrichment analyses (GO, KEGG, GSEA), and gene regulatory network (GRN) data.
#' Uses a parameterized R Markdown template for rendering.
#'
#' @param de_results Data frame or list with differential expression results.
#' @param enrichment_results List of enrichment results (e.g., BP, MF, KEGG, GSEA).
#' @param grn_object An igraph object of the gene regulatory network.
#' @param output_file Output report name (default: "analysis_report.html").
#' @param template_path Path to the R Markdown template. If NULL, uses the built-in template.
#' @param params_list Named list of extra parameters passed to the R Markdown report.
#' @param quiet Logical; if TRUE (default), rendering is quiet.
#'
#' @return Character string with the path to the rendered report.
#'
#'
#' @export
generate_report <- function(de_results,
                            enrichment_results,
                            grn_object,
                            output_file = "analysis_report.html",
                            template_path = NULL,
                            params_list = list(),
                            quiet = TRUE) {

  
  # If no template is supplied, use the built-in one from the package.
  # Use built-in template if none is supplied
  if (is.null(template_path)) {
    template_path <- system.file("extdata", "Template_report_bulk.Rmd", package = "XYomics")
  }
  
  
  # Prepare the list of parameters to pass to the R Markdown report.
  # These parameters should be referenced in your template Rmd file.
  report_params <- c(list(
    de_results = de_results,
    enrichment_results = enrichment_results,
    grn_object = grn_object,
    report_date = Sys.Date()
  ), params_list)
  
  # Render the R Markdown report into the specified output file.
  rendered_report <- rmarkdown::render(
    input = template_path,
    output_file = output_file,
    params = report_params,
    envir = new.env(parent = globalenv()),
    quiet = quiet
  )
  
  return(rendered_report)
}