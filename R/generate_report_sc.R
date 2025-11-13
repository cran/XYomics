#' Generate a Comprehensive Analysis Report
#'
#' This function creates an integrated report that combines key analysis outputs,
#' @param results_cat A data frame or list containing differential expression results.
#' @param enrichment_cat A list with enrichment objects (e.g., BP, MF, KEGG, and optionally GSEA results).
#' @param grn_object An igraph object representing the gene regulatory network (e.g., from PCSF analysis).
#' @param output_file Character. The desired name (and optionally path) for the rendered report (default: "analysis_report.html").
#' @param output_dir Character. Output directory to save the report to. 
#' @param template_path Character. Path to the R Markdown template file. If \code{NULL}, the function uses the built-in template located in \code{inst/rmd/template_report.Rmd}.
#' @param quiet Logical. If \code{TRUE} (default), rendering will be quiet.
#'
#' @return A character string with the path to the rendered report.
#'
#'
#' @export
generate_cat_report <- function(results_cat = results_cat,
                                enrichment_cat =  results_cat,
                                grn_object =  grn_object,
                                output_file = "cat_analysis_report.html",
                                output_dir = tempdir(), 
                                template_path = NULL,
                                quiet = TRUE) {

  # Load the RDS files (lists of data frames)
  de_results <- results_cat
  enrichment_results <- enrichment_cat
  grn_object <-  grn_object 
  

  output_file <- file.path(output_dir, output_file) 
  
  message("Generating report: ", output_file)

  if (is.null(template_path)) {
    template_path <- system.file("extdata", "Template_report.Rmd", package = "XYomics")
    # Fallback to dev path
    if (template_path == "") {
      template_path <- "Template_report.Rmd"
    }
  }
  # Create parameter list
  report_params <- list(
    de_results = de_results,
    enrichment_results = enrichment_results,
    grn_object = grn_object,
    report_date = Sys.Date()
  )
  
  # Render report
  rmarkdown::render(
    input = template_path,
    output_file = output_file,
    params = report_params,
    envir = new.env(parent = globalenv()),
    quiet = quiet
  )
  
  return(invisible(output_file))
}