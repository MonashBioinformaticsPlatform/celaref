
#' make_ranking_violin_plot
#'
#' Plot a panel of violin plots showing the distribution of the 'top' genes of 
#' each of query group, across the reference dataset.
#'
#' In the plot output, each panel correponsds to a different group/cluster in 
#' the query experiment. The x-axis has the groups in the reference dataset. 
#' The y-axis is the rescaled rank of each 'top' gene from the query group, 
#' within each reference group.
#' 
#' Only the 'top' genes for each query group are plotted, forming the violin
#' plots - each individual gene is shown as a tickmark. Some groups have few 
#' top genes, and so their uncertanty can be seen on this plot. 
#' 
#' The thick black lines reprenset the median gene rescaled ranking for each 
#' query group / reference group combination. Having this fall above the dotted 
#' median threshold marker is a quick indication of potential similarity. 
#' A complete lack of similarity would have a median rank around 0.5. Median rankings 
#' much less than 0.5 are common though (an 'anti-cell-groupA' signature), 
#' because genes overrepresented in one group in an experiment, are likely to 
#' be relatively 'underrepresented' in the other groups. Taken to an extreme, 
#' if there are only two reference groups, they'll be complete opposites.
#' 
#'
#' @param de_table.marked The output of \code{\link{get_the_up_genes_for_all_possible_groups}} for the contrast of interest.
#' @param median_rank_threshold Where to draw horizontal threshold line as a rough indicator of similarity. Number betweeen 0 (top-rank) and 1. Default 0.25.
#' @param log10trans  Plot on a log scale? Useful for distinishing multiple similar, yet distinct cell type that bunch at top of plot. Default=FALSE.
#' 
#' @return  A ggplot object.
#'
#' @examples
#'
#' de_table.demo_ref   <- contrast_each_group_to_the_rest(demo_ref_se,    
#'                                                        dataset_name="demo_ref",
#'                                                        num_cores=2)
#' de_table.demo_query <- contrast_each_group_to_the_rest(demo_query_se,  
#'                                                        dataset_name="demo_query",
#'                                                        num_cores=2)
#' de_table.marked.query_vs_ref <- 
#'      get_the_up_genes_for_all_possible_groups( de_table.test=de_table.demo_query, 
#'                                                de_table.ref=de_table.demo_ref, 
#'                                                'demo_query')
#' make_ranking_violin_plot(de_table.marked.query_vs_ref)
#'
#' @seealso  \code{\link{get_the_up_genes_for_all_possible_groups}} To make the input data.
#'
#' 
#'@export
make_ranking_violin_plot <- function(de_table.marked, median_rank_threshold=0.25, log10trans=FALSE) {
   
   if (log10trans) { 
      de_table.marked$rescaled_rank <- log10(de_table.marked$rescaled_rank) #happily, it'll never be 0
      median_rank_threshold <- log10(median_rank_threshold)
   }
   
   p <- ggplot2::ggplot(de_table.marked, ggplot2::aes_string(y='rescaled_rank', x='group', fill='group')) +
      ggplot2::geom_violin(ggplot2::aes_string(colour='group')) +
      ggplot2::geom_point(alpha=0.5, size=3, pch='-', show.legend = FALSE) +
      ggplot2::geom_hline(yintercept = median_rank_threshold, lty=2) + 
      ggplot2::scale_y_reverse() +
      ggplot2::ylab("Test geneset rank in reference cluster") + ggplot2::xlab("") +  
      ggplot2::stat_summary(fun.y = stats::median, fun.ymin = stats::median, fun.ymax = stats::median, geom = "crossbar", col="black", show.legend = FALSE) +
      ggplot2::theme_bw() +    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust=0.5)) +
      ggplot2::facet_wrap(~test_group)
   return (p)
}


