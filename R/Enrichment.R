#' Functional enrichment analysis for analyzing the DE results of different normalization methods and biologically interpreting the results
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param de_res data table resulting of run_DE
#' @param comparison String of comparison (must be a valid comparison saved in de_res)
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in de_res)
#' @param id_column String specifying the column of the rowData of the SummarizedExperiment object which includes the gene names
#' @param organism Organism name (gprofiler parameter)
#' @param source Data source to use (gprofiler parameter, example: KEGG)
#' @param signif_thr Significance threshold
#'
#' @return list of ggplot objects
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' data(tuberculosis_TMT_de_res)
#' plot_intersection_enrichment(tuberculosis_TMT_se, tuberculosis_TMT_de_res,
#'                 ain = c("IRS_on_Median", "IRS_on_RobNorm", "RobNorm"),
#'                 comparison = "Rx-TBL", id_column = "Gene.Names",
#'                 organism = "hsapiens", source = "GO:BP", signif_thr = 0.05)
#'
plot_intersection_enrichment <- function(se, de_res, comparison, ain = NULL, id_column = "Gene.Names", organism = "hsapiens", source = "KEGG", signif_thr = 0.05){
  # Check parameters
  tmp <- check_plot_DE_parameters(de_res, ain, c(comparison))
  de_res <- tmp[[1]]
  ain <- tmp[[2]]
  comparison <- tmp[[3]]
  
  # Check id_column
  if(!id_column %in% colnames(de_res)){
    # if id_column not in rowData
    stop(paste0(id_column, " not in DE results!"))
  }
  
  # Prepare DE results
  de_res <- de_res[de_res$Change != "No Change",]
  de_res <- de_res[de_res$Assay %in% ain,]

  de_res <- de_res[de_res$Comparison %in% c(comparison),]
  dt <- de_res
  
  assay_ordering <- levels(de_res$Assay)
  if(is.null(assay_ordering)){
    assay_ordering <- unique(de_res$Assay)
  }
  
  # Prepare queries for gProfiler
  queries <- lapply(ain, function(method) {
    dt <- dt[dt$Assay == method, ]
    query <- dt[[id_column]]
    query <- unique(query[query != ""])
    # split by ;
    query <- unlist(strsplit(query, ";"))
    query
  })
  names(queries) <- ain
  
  # Run gProfiler
  gres <- gprofiler2::gost(queries, organism = organism, 
                           sources = c(source), user_threshold = signif_thr)
  gres <- gres$result
  gres <- gres[, c("query", "term_name", "source")]
  gres$present <- 1
  dt <- data.table::dcast(data.table::as.data.table(gres), 
                          ... ~ query, value.var = "present")
  dt[is.na(dt)] <- 0
  
  
  # Prepare results for clustering and Jaccard
  cluster_dt <- dt
  terms <- cluster_dt$term_name
  cluster_dt$term_name <- NULL
  cluster_dt$source <- NULL
  cluster_dt <- as.data.frame(cluster_dt)
  rownames(cluster_dt) <- terms

  # Clustering for heatmap
  tryCatch({
    dist.mat <- vegan::vegdist(cluster_dt, method = "jaccard")
    clust.res <- stats::hclust(dist.mat)
    ordering_terms <- rownames(cluster_dt)[clust.res$order]
    
    dist.mat <- vegan::vegdist(t(cluster_dt), method = "jaccard")
    clust.res <- stats::hclust(dist.mat)
    ordering_assays <- rownames(t(cluster_dt))[clust.res$order]
  }, error = function(e){
    message("Error in clustering: ", e)
    message("Rows and columns in enrichment heatmap not clustered!")
    ordering_terms <- rownames(cluster_dt)
    ordering_assays <- colnames(cluster_dt)
  })
  
  # Jaccard distance
  jaccard <- function(x, y){
    M.11 <- sum(x == 1 & y == 1)
    M.10 <- sum(x == 1 & y == 0)
    M.01 <- sum(x == 0 & y == 1)
    return(M.11 / (M.11 + M.10 + M.01))
  }
  
  assays_without_terms <- assay_ordering[!assay_ordering %in% colnames(dt)]
  # add tables with 0 for assays without terms
  for(assay in assays_without_terms){
    cluster_dt[, assay] <- 0
  }
  
  cluster_dt <- cluster_dt[, assay_ordering]
  
  m <- matrix(data = NA, nrow = ncol(cluster_dt), ncol = ncol(cluster_dt))
  for(i in seq_len(ncol(cluster_dt))){
    for(j in seq_len(ncol(cluster_dt))){
      col1 <- colnames(cluster_dt)[i]
      col2 <- colnames(cluster_dt)[j]
      if(col1 == col2){
        m[i, j] <- 1
      } else if(i > j){
        m[i, j] <- jaccard(cluster_dt[, col1], cluster_dt[, col2])
      }
    }
  }
  colnames(m) <- colnames(cluster_dt)
  rownames(m) <- colnames(cluster_dt)
  
  melted_m <- reshape2::melt(m, measure.vars = colnames(m), na.rm = TRUE)
  
  # Jaccard heatmap
  melted_m$Var1 <- factor(melted_m$Var1, levels = assay_ordering)
  melted_m$Var2 <- factor(melted_m$Var2, levels = assay_ordering)
  
  jaccard_heatmap <- ggplot2::ggplot(melted_m, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
    ggplot2::geom_tile(color = "black") +
    ggplot2::geom_text(ggplot2::aes(label = round(value,digits = 2), color = value > 0.5)) +
    ggplot2::scale_fill_gradient(low = "white", high = "#0072B2", limits = c(0,1)) +
    ggplot2::scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
    ggplot2::labs(x = "Normalization Method", y = "Normalization Method", fill = "Jaccard Similarity") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust =0.5))
  
  # Prepare data for intersection heatmap
  dt[dt == 1] <- "Yes"
  dt[dt == 0] <- "No"
  melted_dt <- data.table::melt(data.table::as.data.table(dt), 
                                value.name = "Present", variable.name = "Assay", 
                                measure.vars = colnames(dt)[!colnames(dt) %in% 
                                                              c("term_name", "source")])
  
  melted_dt$term_name <- factor(melted_dt$term_name, levels = ordering_terms)
  melted_dt$Assay <- factor(melted_dt$Assay, levels = ordering_assays)
  enrich_heatmap <- ggplot2::ggplot(melted_dt, ggplot2::aes(x = get("Assay"), 
                                                            y = get("term_name"), fill = get("Present"))) + ggplot2::geom_tile(color = "white") + 
    ggplot2::scale_fill_manual(name = "Significant", 
                               values = c(No = "grey80", Yes = "#D55E00")) + 
    ggplot2::labs(x = "Normalization Method", y = "Terms") +  ggplot2::theme_bw() +
    
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                       vjust = 0.5,
                                                       hjust=1))
  return(list("enrichment_term_plot" = enrich_heatmap, "jaccard_intersection_plot" = jaccard_heatmap))
}