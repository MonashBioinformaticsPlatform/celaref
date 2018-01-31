



#' trim_small_groups_and_low_expression_genes
#'
#' \code{trim_small_groups_and_low_expression_genes} will filter and return a 
#' SummarizedExperiment object (dataset_se) by several metrics:
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
#' @param min_lib_size Minimum library size. Cells with fewer than this many reads removed. Default= 1000
#' @param min_reads_in_sample Require this many reads to consider a gene detected in a sample. Default = 1
#' @param min_detected_by_min_samples Keep genes detected in this many samples. May change with experiment size. Default=5
#' @param min_group_membership Throw out groups/clusters with fewer than this many cells. May change with experiment size. Default=5
#' 
#' @return A filtered dataset_se, ready for use.
#' 
#' @examples
#' \dontrun{
#' de_table.mousebrain <- contrast_each_group_to_the_rest(dataset_se.mousebrain)
#' }
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
