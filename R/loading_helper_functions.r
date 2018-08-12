#' load_se_from_tables
#'
#' Create a SummarizedExperiment object (dataset_se) from a count matrix, cell 
#' information and optionally gene information.
#'
#' This function makes a SummarizedExperiment object in a form that
#' should work for celaref functions. Specifically, that means it will have an
#' 'ID' feild for genes (view with \code{rowData(dataset_se)}), and both
#' 'cell_sample' and 'group' feild for cells (view with
#' \code{colData(dataset_se)}). See parameters for detail.
#' Additionally, the counts will be an integer matrix (not a
#' sparse matrix), and the \emph{group} feild (but not \emph{cell_sample}
#' or \emph{ID}) will be a factor.
#'
#' Note that data will be subsetted to cells present in both the counts matrix
#' and cell info, this is handy for loading subsets of cells.
#' However, if \bold{gene_info_file} is defined, all genes must match exactly.
#'
#' The \code{load_se_from_files} form of this function will run the same checks,
#' but will read everything from files in one go. The \code{load_se_from_tables}
#' form is perhaps more useful when the annotations need to be modified (e.g. 
#' programmatically adding a different gene identifier, renaming groups, 
#' removing unwanted samples). 
#' 
#' Note that the SummarizedExperiment object can also be created without using
#' these functions, it just needs the \emph{cell_sample}, \emph{ID} and
#' \emph{group} feilds as described above. Since sometimes it might be easier
#' to add these to an existing \emph{SummarizedExperiment} from upstream
#' analyses.
#'
#'
#' @param counts_matrix A tab-separated matrix of read counts for each gene
#' (row) and each cell (column). Columns and rows should be named.
#' 
#' @param cell_info_table Table of cell information. If there is a column labelled
#' \emph{cell_sample}, that will be used as the unique cell identifiers. If not,
#' the first column is assumed to be cell identifiers, and will be copied to a
#' new feild labelled \emph{cell_sample}.
#' Similarly - the clusters of these cells should be listed in one column -
#' which can be called 'group' (case-sensitive) or specified with
#' \bold{group_col_name}. \emph{Minimal data format: <cell_sample> <group>}
#' 
#' @param gene_info_table Optional table of gene information. If there is a
#' column labelled
#' \emph{ID}, that will be used as the gene identifiers (they must be unique!).
#' If not, the first column is assumed to be a gene identifier, and will be copied to a
#' new feild labelled \emph{ID}. Must match all rownames in \bold{counts_matrix}.
#' If omitted, ID wll be generated from the rownames of counts_matrix. Default=NA
#' 
#' @param group_col_name Name of the column in \bold{cell_info_table} containing
#' the cluster/group that each cell belongs to. Case-sensitive. Default='group'
#' 
#' @param cell_col_name Name of the column in \bold{cell_info_table} containing
#' a cell id. Ignored if \emph{cell_sample} column is already present. If omitted, 
#' (and no \emph{cell_sample} column) will use first column.
#' Case-sensitive. Default=NA
#'
#' @return A SummarisedExperiment object containing the count data, cell info
#' and gene info.
#'
#' @examples
#'
#' # From data frames (or a matrix for counts) :
#' demo_se <- load_se_from_tables(counts_matrix=demo_counts_matrix, 
#'                                cell_info_table=demo_cell_info_table)
#' demo_se <- load_se_from_tables(counts_matrix=demo_counts_matrix, 
#'                                cell_info_table=demo_cell_info_table, 
#'                                gene_info_table=demo_gene_info_table)
#'
#' # Or from data files : 
#' counts_filepath    <- system.file("extdata", "sim_query_counts.tab",    package = "celaref")
#' cell_info_filepath <- system.file("extdata", "sim_query_cell_info.tab", package = "celaref")
#' gene_info_filepath <- system.file("extdata", "sim_query_gene_info.tab", package = "celaref")
#'
#' demo_se <- load_se_from_files(counts_file=counts_filepath, cell_info_file=cell_info_filepath)
#' demo_se <- load_se_from_files(counts_file=counts_filepath, cell_info_file=cell_info_filepath, 
#'                               gene_info_file=gene_info_filepath )
#'
#' @seealso \href{https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html}{SummarizedExperiment} For general doco on the SummarizedExperiment objects.
#'
#' @family Data-loading functions
#'
#' @import SummarizedExperiment
#' 
#' @export   
load_se_from_tables <- function(
   counts_matrix, cell_info_table, gene_info_table = NA, group_col_name="group",
   cell_col_name=NA 
) {
   
   # If there's no cell_sample, and no cell_col_name, make the first 'cell_sample', else use cell_col_name.
   if (! "cell_sample" %in% colnames(cell_info_table)) {
      if (is.na(cell_col_name)) {
         cell_info_table <- cbind.data.frame(cell_sample=cell_info_table[,1], 
                                             cell_info_table, 
                                             stringsAsFactors = FALSE )
      }
      else {
         stopifnot(cell_col_name %in% colnames(cell_info_table))
         cell_info_table <- cbind.data.frame(cell_sample=cell_info_table[,cell_col_name], 
                                             cell_info_table, 
                                             stringsAsFactors = FALSE )
      }
   }
   
   # Check for 'group' anything from group_col_name will be copied into group
   if (! group_col_name %in% colnames(cell_info_table)) {
      stop( paste0("Couldn't find group/cluster column ", group_col_name ,
                   " in cell_info_table ",cell_info_table) )
   }
   if (group_col_name != "group") {
      cell_info_table$group <- cell_info_table[,group_col_name]
   }
   
   
   # Only keep common cells, match order
   cells <- intersect(cell_info_table$cell_sample, colnames(counts_matrix))
   if (length(cells) <= 1) { 
      stop(paste("Couldn't find cells in common between counts matrix",
                 "(col names) and cell_info_file (cell_sample column,",
                 "first col or specified as cell_col_name)")) 
   }
   if (   length(cells) != nrow(cell_info_table) 
       || length(cells) != ncol(counts_matrix)   ) {
      message(paste0(
         "Not all cells were listed in both counts matrix and cell_info_file. ",
         "Is this expected? Keeping the ", length(cells), " in common"))
   }
   cell_info_table <-cell_info_table[match(cells, cell_info_table$cell_sample),]
   counts_matrix   <-counts_matrix[,cells]
   
   
   # NB factorising group after removal of unmatched cells
   cell_info_table$group <- factor(cell_info_table$group)
   
   

   
   # Make summarised experiment.
   # With or without gene info file.
   dataset_se <- NA
   if (all(is.na(gene_info_table))) {
      dataset_se  <- SummarizedExperiment(
                           counts_matrix,
                           colData=base::as.data.frame(cell_info_table))
      rowData(dataset_se)$ID <- rownames(assay(dataset_se))
   }
   else {
      
      # If there's no ID col, make the first 'ID'
      if (! "ID" %in% colnames(gene_info_table)) {
         gene_info_table <- cbind.data.frame(
            "ID"=as.character(gene_info_table[,1]), 
            gene_info_table, 
            stringsAsFactors=FALSE)
      }
      
      # Cells might not, but genes should be matching.
      genes     <- intersect(gene_info_table$ID, rownames(counts_matrix))
      num_genes <- length(genes)
      if (    num_genes != nrow(gene_info_table) 
           || num_genes != nrow(counts_matrix)   ) { 
         stop( paste("Gene IDs did not match between ID feild of",
                     " gene_info_file (or first column), and row names of ",
                     "counts matrix ")) 
      }
      
      # Create a summarised experiment object.
      dataset_se  <- SummarizedExperiment(
                        counts_matrix,
                        colData=S4Vectors::DataFrame(cell_info_table),
                        rowData=S4Vectors::DataFrame(gene_info_table))
   }
   
   return(dataset_se)
}




#' load_se_from_files
#'
#' \code{load_se_from_files} is a wrapper for \code{load_se_from_tables} that
#' will read in tables from specified files. 
#' 
#' @param counts_file A tab-separated file of a matrix of read counts. As per 
#' \bold{counts_matrix}. First column should be gene ID, and top row cell ids.
#'
#' @param cell_info_file Tab-separated text file of cell information, as per
#' \bold{cell_info_table}. Columns must have names. 
#'
#' @param gene_info_file Optional tab-separated text file of gene information, 
#' as per \bold{gene_info_file}. Columns must have names. Default=NA
#' 
#' @family Data loading functions
#'
#' @describeIn load_se_from_tables To read from files
#' 
#' @import SummarizedExperiment
#' 
#' @export   
load_se_from_files <- function(
   counts_file, cell_info_file, gene_info_file = NA, group_col_name="group", 
   cell_col_name=NA 
) {
   
   counts_matrix   <- as.matrix(utils::read.table(
      counts_file, row.names=1, header=TRUE, sep = "\t", 
      stringsAsFactors = FALSE, check.names=FALSE ))
   
   cell_info_table <- utils::read.table(cell_info_file, header=TRUE, 
                                        sep = "\t", stringsAsFactors = FALSE )
   
   # Read gene Info table, if specified
   gene_info_table <- NA
   if (! is.na(gene_info_file)) {
      gene_info_table <- utils::read.table(gene_info_file, 
                                           header=TRUE, 
                                           sep = "\t", 
                                           stringsAsFactors = FALSE )
   }
   
   return(load_se_from_tables(counts_matrix, cell_info_table, gene_info_table, 
                              group_col_name, cell_col_name = cell_col_name) )
}







   

#' load_dataset_10Xdata
#'
#' Convenience function to create a SummarizedExperiment object (dataset_se) 
#' from a the output of 10X cell ranger pipeline run. 
#' 
#' 
#' This function makes a SummarizedExperiment object in a form that
#' should work for celaref functions. Specifically, that means it will have an
#' 'ID' feild for genes (view with \code{rowData(dataset_se)}), and both
#' 'cell_sample' and 'group' feild for cells (view with
#' \code{colData(dataset_se)}). See parameters for detail.
#' Additionally, the counts will be an integer matrix (not a
#' sparse matrix), and the \emph{group} feild (but not \emph{cell_sample}
#' or \emph{ID}) will be a factor.
#' 
#' The clustering information can be read from whichever cluster is specified,
#' usually there will be several choices. 
#' 
#' This funciton is designed to work with output of version 2.0.1 of the 
#' cellRanger pipeline, may not work with others (will not work for 1.x).
#'     
#' @param dataset_path Path to the directory of 10X data, as generated by the 
#' cellRanger pipeline (versions 2.1.0 and 2.0.1). The directory should have 
#' subdirecotires \emph{analysis}, \emph{filtered_gene_bc_matrices} and
#' \emph{raw_gene_bc_matrices} (only the first 2 are read).
#' @param dataset_genome The genome that the reads were aligned against, 
#' e.g. GRCh38.  Check for this as a directory name under the 
#' \emph{filtered_gene_bc_matrices} subdirectory if unsure.
#' @param clustering_set The 10X cellRanger pipeline produces several different 
#' cluster definitions per dataset. Specify which one to use e.g. 
#' kmeans_10_clusters Find them as directory names under 
#' \emph{analysis/clustering/}
#' @param gene_id_cols_10X Vector of the names of the columns in the gene 
#' description file (\emph{filtered_gene_bc_matrices/GRCh38/genes.csv}). The 
#' first element of this will become the ID. 
#' Default = c("ensembl_ID","GeneSymbol")
#' @param id_to_use Column from \bold{gene_id_cols_10X} that defines the gene 
#' identifier to use as 'ID' in the returned SummarisedExperiment object.
#' Many-to-one relationships betwen the assumed unique first element of 
#' \bold{gene_id_cols_10X} and \bold{id_to_use} will be handled gracefully by 
#' \code{\link[celaref]{convert_se_gene_ids}}. 
#' Defaults to first element of \bold{gene_id_cols_10X}
#' 
#' @return A SummarisedExperiment object containing the count data, cell info
#' and gene info.
#'
#' @examples
#' example_10X_dir <- system.file("extdata", "sim_cr_dataset", package = "celaref")
#' dataset_se <- load_dataset_10Xdata(example_10X_dir, dataset_genome="GRCh38", 
#'     clustering_set="kmeans_4_clusters", gene_id_cols_10X=c("gene")) 
#' 
#' \dontrun{
#' dataset_se <- load_dataset_10Xdata('~/path/to/data/10X_pbmc4k', dataset_genome="GRCh38", 
#'     clustering_set="kmeans_7_clusters") 
#' } 
#'
#' @seealso \href{https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html}{SummarizedExperiment} 
#' For general doco on the SummarizedExperiment objects.
#' @seealso \code{\link[celaref]{convert_se_gene_ids}} describes method for 
#' converting IDs.
#' 
#' @family Data loading functions
#'
#' @import SummarizedExperiment
#' 
#' @export
load_dataset_10Xdata <- function(
   dataset_path, dataset_genome, clustering_set, 
   gene_id_cols_10X =c("ensembl_ID","GeneSymbol"), 
   id_to_use = gene_id_cols_10X[1] 
) {

   matrix_file <- file.path(dataset_path,"filtered_gene_bc_matrices",
                            dataset_genome,"matrix.mtx")
   cells_file  <- file.path(dataset_path,"filtered_gene_bc_matrices",
                            dataset_genome,"barcodes.tsv")
   genes_file  <- file.path(dataset_path,"filtered_gene_bc_matrices",
                            dataset_genome,"genes.tsv")
   
   #.../10X_pbmc5pExpr/analysis/clustering/kmeans_5_clusters/clusters.csv
   clustering_file <- file.path(dataset_path,"analysis","clustering",
                                clustering_set,"clusters.csv")
   clustering_table <- readr::read_csv(clustering_file, col_types=readr::cols())
   colnames(clustering_table) <- c('cell_sample', 'group')
   clustering_table$group <- factor(clustering_table$group)
   
   # gene info
   # Start with the first id (assumed uniq!), but change to specified after.
   genes_table    <- readr::read_tsv(genes_file, 
                                     col_names = gene_id_cols_10X, 
                                     col_types = readr::cols())
   genes_table$ID <- dplyr::pull(genes_table, gene_id_cols_10X[1])
   
   
   # Not there's no lables here, but <genes> rows and <cells> columns
   filtered_matrix <- as.matrix(Matrix::readMM(matrix_file)) # from Matrix
   storage.mode(filtered_matrix ) <- "integer"
   order_of_cells <- scan(cells_file, what=character())
   colnames(filtered_matrix) <- order_of_cells 
   rownames(filtered_matrix) <- genes_table$ID    
   
   # Create a summarised experiment objct.
   dataset_se  <- SummarizedExperiment(filtered_matrix, 
                                       colData=clustering_table,
                                       rowData=genes_table)
   
   # Optionally change id (handles m:1)
   rowData(dataset_se)$total_count <- Matrix::rowSums(assay(dataset_se))
   if (id_to_use != gene_id_cols_10X[1] ) { 
      dataset_se <- convert_se_gene_ids(dataset_se, 
                                        new_id=id_to_use, 
                                        eval_col='total_count')
   }
   
   return(dataset_se)
} 












#' convert_se_gene_ids
#'
#' Change the gene IDs in in the supplied datatset_se object to some other id 
#' already present in the gene info (as seen with \code{rowData()})
#'
#' @param dataset_se Summarised experiment object containing count data. Also
#' requires 'ID' and 'group' to be set within the cell information
#' (see \code{colData()})
#' @param new_id  A column within the feature information (view 
#' \code{colData(dataset_se)})) of the \bold{dataset_se}, which will become
#' the new ID column. Non-uniqueness of this column is handled gracefully! 
#' Any \emph{NAs} will be dropped.
#' @param eval_col Which column to use to break ties of duplicate \bold{new_id}.
#' Must be a column within the feature information (view 
#' \code{colData(dataset_se)})) of the \bold{dataset_se}. Total reads per gene
#' feature is a good choice.
#' @param find_max If false, this will choose the minimal \bold{eval_col} 
#' instead of max. Default = TRUE
#'
#' @return A modified dataset_se - ID will now be \bold{new_id}, and unique. 
#' It will have fewer genes if old ID to new ID was not a 1:1 mapping. 
#' The selected genes will be according to the eval col max (or min). 
#' \emph{should} pick the alphabetical first on ties, but could change. 
#'
#' @examples
#' 
#' # The demo dataset doesn't have other names, so make some up 
#' # (don't do this)
#' dataset_se <- demo_ref_se
#' rowData(dataset_se)$dummyname <- toupper(rowData(dataset_se)$ID)
#'
#' # If not already present, define a column to evaluate, typically total reads/gene.
#' rowData(dataset_se)$total_count <- rowSums(assay(dataset_se))
#' 
#' dataset_se <- convert_se_gene_ids(dataset_se, new_id='dummyname', eval_col='total_count') 
#'
#' @seealso \href{https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html}{SummarizedExperiment} 
#' For general doco on the SummarizedExperiment objects.
#' @seealso \code{\link[celaref]{load_se_from_files}} For reading data from flat 
#' files (not 10X cellRanger output)
#'
#' @import SummarizedExperiment
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
convert_se_gene_ids <- function(dataset_se, new_id, eval_col, find_max=TRUE) {
   
   old_id = "ID" 
   if (! all(c(old_id, new_id, eval_col) %in% colnames(rowData(dataset_se))) ) {
      stop(paste("Can't find all of", c(old_id, new_id, eval_col), 
                 "in rowData(dataset_se) colnames"))
   } 
   
   row_data_df <- BiocGenerics::as.data.frame(rowData(dataset_se))[,c(old_id,new_id, eval_col)]
   colnames(row_data_df) <- c("old_lab", "new_lab", "eval_lab")
   row_data_df <- row_data_df[!is.na(row_data_df$new_lab),] 
   if (find_max) {
      row_data_df <- row_data_df %>% 
         dplyr::arrange(.data$new_lab, 
                        dplyr::desc(.data$eval_lab), 
                        .data$old_lab)
   }else { #min
      row_data_df <- row_data_df %>% 
         dplyr::arrange(.data$new_lab, 
                        .data$eval_lab, 
                        .data$old_lab)
   }
   row_data_unique <- row_data_df %>%
      dplyr::group_by(.data$new_lab) %>% 
      dplyr::slice(1)
   
   # Subset to just those representative old ids, and give the unique new id 
   dataset_se<- dataset_se[  row_data_unique$old_lab , ]
   rowData(dataset_se)$ID <- rowData(dataset_se)[[new_id]] # overwrite the ID
   rownames(dataset_se)   <- rowData(dataset_se)[["ID"]]
   
   return(dataset_se)
} 







#' trim_small_groups_and_low_expression_genes
#'
#' Filter and return a SummarizedExperiment object (dataset_se) by several
#' metrics:
#' \itemize{
#'   \item Cells with at least \bold{min_lib_size} total reads.
#'   \item Genes expressed in at least \bold{min_detected_by_min_samples} cells, 
#'   at a threshold of \bold{min_reads_in_sample} per cell.
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
#' @param min_lib_size Minimum library size. Cells with fewer than this many 
#' reads removed. Default = 1000
#' @param min_reads_in_sample Require this many reads to consider a gene 
#' detected in a sample. Default = 1
#' @param min_detected_by_min_samples Keep genes detected in this many samples. 
#' May change with experiment size. Default = 5
#' @param min_group_membership Throw out groups/clusters with fewer than this 
#' many cells. May change with experiment size. Default = 5
#'
#' @return A filtered dataset_se, ready for use.
#'
#' @examples
#' 
#' demo_query_se.trimmed  <- 
#'    trim_small_groups_and_low_expression_genes(demo_query_se)
#' demo_query_se.trimmed2 <- 
#'    trim_small_groups_and_low_expression_genes(demo_ref_se, 
#'                                               min_group_membership = 10)
#'
#' @import SummarizedExperiment
#' 
#' @export
trim_small_groups_and_low_expression_genes <- function(
   dataset_se, min_lib_size=1000, min_group_membership=5,
   min_reads_in_sample=1, min_detected_by_min_samples=5 
) {
   
   ## Filter by min lib size, num samples detected in
   samples_per_gene <- Matrix::rowSums(assay(dataset_se) >= min_reads_in_sample)
   dataset_se <- dataset_se[,Matrix::colSums(assay(dataset_se))>=min_lib_size ]
   dataset_se <- dataset_se[ samples_per_gene >=  min_detected_by_min_samples, ]
   
   ## Less than a certain number of cells in a group, 
   # discard the group, and its cells.
   # NB: also removes 'NA' group entries.
   cell_group_sizes <- table(dataset_se$group)
   groups_to_keep   <- names(cell_group_sizes)[cell_group_sizes >= min_group_membership]
   dataset_se       <- dataset_se[,dataset_se$group %in% groups_to_keep]
   dataset_se$group <- droplevels(dataset_se$group)
   
   return(dataset_se)
}
