




#' load_se_from_files
#'
#' Create a SummarizedExperiment object (dataset_se) from files.
#'
#' This convenience function makes a SummarizedExperiment object in a form that
#' should work for celaref functions. Specifically, that means it will have an
#' 'ID' feild for genes (view with \code{rowData(dataset_se)}), and both
#' 'cell_sample' and 'group' feild for cells (view with
#' \code{colData(dataset_se)}). See parameters for detail.
#' Additionally, the counts will be an integer matrix (not a
#' sparse matrix), and the \emph{group} feild (but not \emph{cell_sample}
#' or \emph{ID}) will be a factor.
#'
#' Note that data will be subsetted to cells present in both the counts matrix
#' and cell info files, this is handy for loading subsets of cells.
#' However, if \bold{gene_info_file} is defined, all genes must match exactly.
#'
#' @param counts_file A tab-separated file of a matrix of read counts for each gene
#' (row) and each cell (column). First column should be gene names,
#' and top row cell ids.
#'
#' @param cell_info_file Tab-separated text file of cell information.
#' Columns must have names. If there is a column labelled
#' \emph{cell_sample}, that will be used as the unique cell identifiers. If not,
#' the first column is assumed to be cell identifiers, and will be copied to a
#' new feild labelled \emph{cell_sample}.
#' Similarly - the clusters of these cells should be listed in one column -
#' which can be called 'group' (case-sensitive) or specified with
#' \bold{group_col_name}. \emph{Minimal file format: <cell_sample> <group>.}
#'
#' @param gene_info_file Optional tab-separated text file of gene information.
#' Columns must have names. If there is a column labelled
#' \emph{ID}, that will be used as the gene identifiers (they must be unique!).
#' If not, the first column is assumed to be a gene identifier, and will be copied to a
#' new feild labelled \emph{ID}. Must match all rownames in counts_matrix file.
#' If omitted, ID wll be generated from counts file rownames. Default=NA
#'
#' @param group_col_name Name of the column in \bold{cell_info_file} containing
#' the cluster/group that each cell belongs to. Case-sensitive. Default='group'
#'
#'
#' @return A SummarisedExperiment object containing the count data, cell info
#' and gene info.
#'
#' @examples
#'
#' # Use demo data files
#' counts_filepath    <- system.file("extdata", "sim_query_counts.tab",    package = "celaref")
#' cell_info_filepath <- system.file("extdata", "sim_query_cell_info.tab", package = "celaref")
#' gene_info_filepath <- system.file("extdata", "sim_query_gene_info.tab", package = "celaref")
#'
#' demo_se <- load_se_from_files(counts_file=counts_filepath, cell_info_file=cell_info_filepath)
#'
#' # Same counts, but with some extra info on the genes stored
#' demo_se <- load_se_from_files(counts_file=counts_filepath, cell_info_file=cell_info_filepath, gene_info_filepath=gene_info_filepath )
#'
#' @seealso \code{\link[RangedSummarizedExperiment-class]{SummarizedExperiment}} For general doco on the SummarizedExperiment objects.
#'
#'
#'@export
load_se_from_files <- function(counts_file, cell_info_file, gene_info_file = NA, group_col_name="group") {

  counts_matrix   <- as.matrix(read.table(counts_file, row.names=1, header=TRUE, sep = "\t", stringsAsFactors = FALSE ))

  cell_info_table <- read.table(cell_info_file, header=TRUE, sep = "\t", stringsAsFactors = FALSE )
  # If there's no cell_sample, make the first 'cell_sample'
  if (! "cell_sample" %in% colnames(cell_info_table)) {
    cell_info_table <- cbind.data.frame(cell_sample=cell_info_table[,1], cell_info_table, stringsAsFactors = FALSE )
  }

  # Check for 'group' anything from group_col_name will be copied into group
  if (! group_col_name %in% colnames(cell_info_table)) {
    stop( paste0("Couldn't find group/cluster column ", group_col_name ," in cell_info_table ",cell_info_table) )
  }
  if (group_col_name != "group") {cell_info_table$group <- cell_info_table[,group_col_name]}


  # Only keep common cells, match order
  cells <- intersect(cell_info_table$cell_sample, colnames(counts_matrix))
  if (length(cells) <= 1) { stop("Couldn't find cells in common between counts matrix (col names) and cell_info_file (first col) ") }
  if (length(cells) != nrow(cell_info_table) || length(cells) != ncol(counts_matrix) ) {
    message(paste0("Not all cells were listed in both counts matrix and cell_info_file. Is this expected? Keeping the ", length(cells), " in common"))
  }
  cell_info_table <- cell_info_table[match(cells, cell_info_table$cell_sample),]
  counts_matrix <- counts_matrix[,cells]


  # NB factorising group after removal of unmatched cells
  cell_info_table$group <- factor(cell_info_table$group)


  # Make summarised experiment.
  # With or without gene info file.
  dataset_se <- NA
  if (is.na(gene_info_file)) {
    dataset_se  <- SummarizedExperiment(counts_matrix,
                                        colData=DataFrame(cell_info_table))
    rowData(dataset_se)$ID <- rownames(assay(dataset_se))
  }
  else {

    # If there's no ID col, make the first 'ID'
    gene_info_table <- read.table(gene_info_file, header=TRUE, sep = "\t", stringsAsFactors = FALSE )
    if (! "ID" %in% colnames(gene_info_table)) {
      gene_info_table <- cbind.data.frame("ID"=gene_info_table[,1], gene_info_table)
    }

    # Cells might not, but genes should be matching.
    genes <- intersect(gene_info_table$ID, rownames(counts_matrix))
    if (length(genes) != nrow(gene_info_table) || length(genes) != nrow(counts_matrix))
    { stop("Gene IDs did not match between ID feild of gene_info_file (or first column), and row names of counts matrix ") }

    # Create a summarised experiment object.
    dataset_se  <- SummarizedExperiment(counts_matrix,
                                        colData=DataFrame(cell_info_table),
                                        rowData=DataFrame(gene_info_table))
  }

  return(dataset_se)
}





















#' trim_small_groups_and_low_expression_genes
#'
#' Filter and return a SummarizedExperiment object (dataset_se) by several
#' metrics:
#' \itemize{
#'   \item Cells with at least \bold{min_lib_size} total reads.
#'   \item Genes expressed in at least \bold{min_detected_by_min_samples} cells, at a threshold of
#'   \bold{min_reads_in_sample} per cell.
#'   \item Remove entire groups (clusters) of cells where there are fewer than
#'   \bold{min_group_membership} cells in that group.
#' }
#'
#' If it hasn't been done already, it is highly reccomended to use this function
#' to filter out genes with no/low total counts (especially in single cell data,
#' there can be many) - without expression they are not useful and may reduce
#' statistical power.
#'
#' Likewise, very small groups (<5 cells) are unlikely to give useful
#' results with this method. And cells with abnormally small library sizes may
#' not be desireable.
#'
#'
#' Of course 'reasonable' thresholds for filtering cells/genes are subjective.
#' Defaults are moderately sensible starting points.
#'
#' @param dataset_se Summarised experiment object containing count data. Also
#' requires 'ID' and 'group' to be set within the cell information
#' (see \code{colData()})
#' @param min_lib_size Minimum library size. Cells with fewer than this many reads removed. Default = 1000
#' @param min_reads_in_sample Require this many reads to consider a gene detected in a sample. Default = 1
#' @param min_detected_by_min_samples Keep genes detected in this many samples. May change with experiment size. Default = 5
#' @param min_group_membership Throw out groups/clusters with fewer than this many cells. May change with experiment size. Default = 5
#'
#' @return A filtered dataset_se, ready for use.
#'
#' @examples
#' demo_query_se.filtered <- trim_small_groups_and_low_expression_genes(demo_query_se)
#'
#'
#'@export
trim_small_groups_and_low_expression_genes <- function(dataset_se,
                                                       min_lib_size=1000, min_group_membership=5,
                                                       min_reads_in_sample=1, min_detected_by_min_samples=5
                                                       ) {

    ## Filter by min lib size, num samples detected in
    dataset_se <- dataset_se[,Matrix::colSums(assay(dataset_se))>=min_lib_size ]
    dataset_se <- dataset_se[Matrix::rowSums(assay(dataset_se) >= min_reads_in_sample) >=  min_detected_by_min_samples, ]

    ## Less than a certain number of cells in a group, discard the group, and its cells.
    # NB: also removes 'NA' group entries.
    cell_group_sizes <- table(dataset_se$group)
    groups_to_keep   <- names(cell_group_sizes)[cell_group_sizes >= min_group_membership]
    dataset_se       <- dataset_se[,dataset_se$group %in% groups_to_keep]
    dataset_se$group <- droplevels(dataset_se$group)

    return(dataset_se)
}
