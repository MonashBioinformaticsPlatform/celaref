


#' make_ref_similarity_names_for_groups
#'
#' Given the comparison between the two dataset, construct some sensible labels 
#' for the groups in the query group.
#' 
#' This function aims to report a) the top most similar reference group, if 
#' there's a clear frontrunner, b) A list of multiple similar groups if they 
#' have similar similarity, or c) 'No similarity', if there is none.
#' 
#' Each group is named according to the following rules. Testing for significant 
#' (smaller) differences with a one-directional Mann-Whitney U test on their rescaled ranks:
#' \enumerate{
#'   \item The first (as ranked by median rescaled rank) reference group is 
#'   significantly more similar than the next: Report \emph{first only}.
#'   \item When comparing differences betwen groups stepwise ranked by 
#'   median rescaled rank - no group is significantly different to its 
#'   neighbour: Report \emph{no similarity}
#'   \item There's no significant differences in the stepwise comparisons 
#'   of the first N reference groups - but there is a significant 
#'   difference later on : Report \emph{multiple group similarity}
#' }
#' 
#' The similarity is formatted into a group label. Where there are 
#' multiple similar groups, they're listed from most to least similar by their 
#' median ranks.
#'  
#' For instance, a query dataset of clusters c1, c2 and c3 againsts a cell-type
#' labelled reference datatset might get names like:
#' E.g.
#' \itemize{
#'   \item c1:macrophage
#'   \item c2:endotheial|mesodermal
#'   \item c3:no_similarity
#' }
#' 
#' 
#'   
#'
#' @param de_table.ref.marked The output of \code{\link{get_the_up_genes_for_all_possible_groups}} for the contrast of interest.
#' @param the_test_dataset A short meaningful name for the experiment. (Should match \emph{test_dataset} column in \bold{de_table.marked})
#' @param the_ref_dataset A short meaningful name for the experiment. (Should match \emph{dataset} column in \bold{de_table.marked})
#' @param pval Differences between the rescaled ranking distribution of 'top'
#' genes on different reference groups are tested with a Mann-Whitney U test. If 
#' they are \emph{significantly different}, only the top group(s) are reported. 
#' ie. A more stringent \bold{pval} is more likely to report multiple 
#' reference groups. 
#' This parameter unlikely to need change. Default = 0.01. 
#'
#' @return A table of automagically-generated labels for each query group, given 
#' their similarity to reference groups. Short_lab and long_lab are just 
#' different versions of the name.
#' 
#'
#' @examples
#' 
#' de_table.demo_query <- contrast_each_group_to_the_rest(demo_query_se, "demo_query", num_cores=2)
#' de_table.demo_ref   <- contrast_each_group_to_the_rest(demo_ref_se,   "demo_ref", num_cores=2)
#' de_table.marked.query_vs_ref <- get_the_up_genes_for_all_possible_groups(
#'      de_table.demo_query, de_table.demo_ref, 'demo_query')
#'
# make_ref_similarity_names_for_groups(de_table.marked.query_vs_ref,
#                               the_test_dataset="demo_query",
#                               the_ref_dataset="demo_ref")
#'
#' @seealso \code{\link[celaref]{get_the_up_genes_for_all_possible_groups}} To prepare the \bold{de_table.ref.marked} input.
#' 
#' @importFrom magrittr %>%
#' @export
make_ref_similarity_names_for_groups <- function(de_table.ref.marked, the_test_dataset, 
                                                 the_ref_dataset,  pval=0.01){
   
   test_groups <- base::unique(de_table.ref.marked$test_group)
   
   
   #NB: mwtest table won't include anything that doesn't have any similarity
   mwtest_res_table <- dplyr::bind_rows(lapply(FUN=get_ranking_and_test_results, X=test_groups, 
                                               de_table.ref.marked=de_table.ref.marked, 
                                               the_test_dataset=the_test_dataset, 
                                               the_ref_dataset = the_ref_dataset,
                                               pval=pval))
   
   # Just for when there are labells to be made
   labels_by_similarity_table <- mwtest_res_table %>% 
      dplyr::filter(.data$in_name==TRUE) %>%
      dplyr::arrange(.data$test_group, .data$grouprank ) %>%
      dplyr::group_by(.data$test_group) %>%
      dplyr::mutate(shortlab             = paste0(.data$test_group, ":", paste0(.data$group, collapse = "|")), 
                    pval                 = S4Vectors::tail(.data$pval_to_next, n=1),
                    pval_sim_to_all_rest = .data$pval_sig_sim_to_rest[1],
                    pval_top_to_rest     = .data$pval_top_to_rest[1],
                    longlab              = paste0(.data$test_group, ":", paste0(.data$group, collapse = "|"),":",.data$ref_dataset,"_similarity")  ) %>%
      dplyr::select(.data$test_group, .data$shortlab,  
                    .data$pval, .data$pval_sim_to_all_rest, .data$pval_top_to_rest, .data$longlab) %>% base::unique()
   
   # Add no-similarity groups back in.
   if (any ( ! test_groups %in% labels_by_similarity_table$test_group)) {
      # Add no-similarity groups back in
      # Pval2next will only be top to second 
      labels_by_similarity_table.nomatch <- mwtest_res_table %>% 
         dplyr::filter( ! .data$test_group %in% labels_by_similarity_table$test_group)  %>%
         dplyr::arrange(.data$test_group, .data$grouprank ) %>%
         dplyr::group_by(.data$test_group) %>%
         dplyr::mutate(shortlab  = paste0(.data$test_group, ":no_similarity"), 
                       pval                  = .data$pval_to_next[1],
                       pval_sim_to_all_rest  = .data$pval_sig_sim_to_rest[1],
                       pval_top_to_rest      = .data$pval_top_to_rest[1],
                       longlab=paste0(.data$test_group, ":",the_ref_dataset,"_no_similarity")) %>%
         dplyr::select(.data$test_group, .data$shortlab, 
                       .data$pval, .data$pval_sim_to_all_rest,  .data$pval_top_to_rest, .data$longlab) %>% base::unique()
      
      labels_by_similarity_table <- dplyr::bind_rows(labels_by_similarity_table, labels_by_similarity_table.nomatch)
   }
   
   return(labels_by_similarity_table)
}






#' get_ranking_and_test_results
#'
#' Internal function to get reference group similarity contrasts for an individual query qroup. 
#' 
#' For use by \bold{make_ref_similarity_names_for_groups}, see that function for parameter details.
#' This function just runs this for a single query group \bold{the_test_group}
#'
#' @param de_table.ref.marked see \link[celaref]{make_ref_similarity_names_for_groups}
#' @param the_test_group The group to calculate the stats on.
#' @param the_test_dataset see \link[celaref]{make_ref_similarity_names_for_groups}
#' @param the_ref_dataset see \link[celaref]{make_ref_similarity_names_for_groups}
#' @param pval see \link[celaref]{make_ref_similarity_names_for_groups}
#' @param median_rank_threshold see \link[celaref]{make_ref_similarity_names_for_groups}
#'  
#' @seealso \code{\link[celaref]{make_ref_similarity_names_for_groups}} which calls this. 
#' @importFrom magrittr %>%
#' @importFrom rlang .data
get_ranking_and_test_results <- function (de_table.ref.marked, the_test_group, the_test_dataset, the_ref_dataset, pval=0.01) {
   
   # Should only be one test datast, and one reference dataset in de_table.ref_marked. Not needed, but forced here for paranoia reasons.
   de_table.ref.marked <- de_table.ref.marked %>% 
      dplyr::filter(.data$test_dataset==the_test_dataset, .data$dataset==the_ref_dataset )
   
   
   # Get the average rankings in order.
   # For this test_group
   rankstat_table <- get_rankstat_table(de_table.ref.marked, the_test_group)
   #group                median_rank mean_rank     n grouprank
   #<fct>                      <dbl>     <dbl> <int>     <int>
   #1 microglia                0.00314    0.0574    60         1
   #2 interneurons             0.354      0.519     60         2
   #3 oligodendrocytes         0.757      0.616     60         3
   last_rank <- base::nrow(rankstat_table)
   
   # Sanity checks    
   if (base::nrow(rankstat_table) ==1) { rlang::warn("Only 1 reference group. Something odd is happening.") }
   
   
   ## Run 'stepped' contrasts.
   # from this ranking, Pair first with second, secont to third, e.t.c.
   # But last entry is the second last to last ranks.
   pairset.step <- dplyr::bind_cols(groupA=as.character(rankstat_table$group[c(1:(last_rank-1))]),
                                    groupB=as.character(rankstat_table$group[c(2:(last_rank))])) 
   
   mwtest_res_table.step <- pairset.step %>% dplyr::rowwise() %>% 
      dplyr::do (run_pair_test_stats(de_table.ref.marked=de_table.ref.marked, the_test_group, .data$groupA, .data$groupB))
   
   # Run checks
   rankstat_table$group <- as.character(rankstat_table$group)
   ranking_and_mwtest_results <- base::merge(x=rankstat_table, y=mwtest_res_table.step, 
                                             by.x="group", by.y="groupB", all.x=TRUE)                        
   ranking_and_mwtest_results <- ranking_and_mwtest_results %>% 
      dplyr::arrange(.data$grouprank) %>%
      dplyr::rename(up_contrast_group_step=.data$groupA) %>% 
      tibble::add_column( test_group=the_test_group, 
                          test_dataset=the_test_dataset, 
                          ref_dataset=the_ref_dataset, .before=1)
   
   # SHuffly the pval to be a pval-to-next not p-val-from-previous 
   # (want for included names - last becomes NA instead of first)
   ranking_and_mwtest_results$pval_to_next <- c(ranking_and_mwtest_results$pval[2:last_rank], NA)
   
   
   # Where's the first jump that is sigificantly different? If any.
   # Note that first pval is NA because test is vs higher ranked group.
   # I.e min will be 2. 
   ranking_and_mwtest_results$in_name <- FALSE
   if (any(ranking_and_mwtest_results$pval[-1] < pval) ) {
      first_sig_different_grouprank <- min(which(ranking_and_mwtest_results$pval <= pval))
      ranking_and_mwtest_results[1:(first_sig_different_grouprank-1), 'in_name'] <- TRUE
   }
   
   
   
   
   # Compare top group (sig or not) to all else
   # Probably remove this test.
   top_group <- ranking_and_mwtest_results$group[1]
   ranks.top_group <- de_table.ref.marked %>% 
      dplyr::filter( test_group == the_test_group, group==top_group ) %>% 
      dplyr::pull(rescaled_rank)
   ranks.not_top_group <- de_table.ref.marked %>% 
      dplyr::filter( test_group == the_test_group, group!=top_group ) %>% 
      dplyr::pull(rescaled_rank)
   pval.top2rest <- suppressMessages(suppressWarnings(stats::wilcox.test( alternative = "less",  
                                        x = ranks.top_group,
                                        y = ranks.not_top_group )$p.value))
   
   
   
   # And, what's the difference of first (or up to first sig diff) to the 
   # rest of the groups compbind?
   # For 6 groups at 60 gehes thats: 60 vs 60*5, or 60*2 vs 60*4.
   # If there are no significant differences to draw the line at - report NA.
   pval.sigsimilar2test <- NA
   if (any(ranking_and_mwtest_results$in_name)) {
      
      sigsim_groups <- ranking_and_mwtest_results$group[ranking_and_mwtest_results$in_name]
      
      ranks.sigsim_groups <- de_table.ref.marked %>% 
         dplyr::filter( test_group == the_test_group, group %in% sigsim_groups ) %>% 
         dplyr::pull(rescaled_rank)
      
      ranks.not_sigsim_group <- de_table.ref.marked %>% 
         dplyr::filter( test_group == the_test_group, ! (group %in% top_group) ) %>% 
         dplyr::pull(rescaled_rank)
      
      pval.sigsimilar2test <- suppressMessages(suppressWarnings(stats::wilcox.test( alternative = "less",
                                                  x = ranks.sigsim_groups,
                                                  y = ranks.not_sigsim_group)$p.value))
      
   }
   
   # query group level results, so repeat over whole table.
   ranking_and_mwtest_results$pval_top_to_rest     <- pval.top2rest
   ranking_and_mwtest_results$pval_sig_sim_to_rest <- pval.sigsimilar2test
   
   
   return(ranking_and_mwtest_results)
}




#' get_rankstat_table 
#' 
#' Summarise the comparison of the specified query group against in the 
#' comparison in \bold{de_table.ref.marked} - number of 'top' genes and their 
#' median rank in each of the reference groups, with reference group rankings.
#'
#' @param de_table.ref.marked The output of \code{\link{get_the_up_genes_for_all_possible_groups}} for the contrast of interest.
#' @param the_test_group Name of query group to test
#'
#' @return A tibble of query group name (test_group), number of 'top' genes (n), 
#' reference dataset group (group) with its ranking (grouprank) and the median 
#' (rescaled 0..1) ranking of 'top' genes (median_rank).
#'
#' @examples
#'
#' de_table.demo_query <- contrast_each_group_to_the_rest(demo_query_se, 
#'                            "demo_query", num_cores=2)
#' de_table.demo_ref   <- contrast_each_group_to_the_rest(demo_ref_se,   
#'                            "demo_ref", num_cores=2)
#' de_table.marked.query_vs_ref <- get_the_up_genes_for_all_possible_groups(
#'     de_table.demo_query, de_table.demo_ref, 'demo_query')
#'
#' get_rankstat_table(de_table.marked.query_vs_ref, "Group3")
#'
#' @seealso \code{\link[celaref]{get_the_up_genes_for_all_possible_groups}} To prepare the \bold{de_table.ref.marked} input.
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr n
get_rankstat_table <- function(de_table.ref.marked, the_test_group){
   
   # Get the average rankings in order.
   # For this test_group
   rankstat_table <- de_table.ref.marked %>%
      dplyr::group_by(.data$group) %>%
      dplyr::filter(.data$test_group == the_test_group) %>% #only test group
      dplyr::summarise(median_rank=stats::median(.data$rescaled_rank),
                       n=n() ) %>%
      dplyr::arrange(.data$median_rank)
   rankstat_table$grouprank <- base::as.integer(base::rownames(rankstat_table))
   return(rankstat_table)
   
   # NB:
   #n() is part of dplyr, but yeilds error if specified as 
   #dplyr::n 'Evaluation error: This function 
   #should not be called directly.'
   
}




#' run_pair_test_stats
#'
#' Internal function to compare the distribution of a query datasets 'top' 
#' genes between two different reference datasete groups with a 
#' Mann–Whitney U test. One directional test if groupA median < group B.
#' 
#' For use by make_ref_similarity_names_for_groups
#'
#' @param de_table.ref.marked The output of \code{\link{get_the_up_genes_for_all_possible_groups}} for the contrast of interest.
#' @param the_test_group Name of the test group in query dataset.
#' @param groupA One of the reference group names
#' @param groupB Another of the reference group names
#'
#' @return A tibble of wilcox / man-whitneyU test results for this contrast.
#'
#' @seealso  \code{\link[celaref]{make_ref_similarity_names_for_groups}} 
#'
#' @importFrom magrittr %>%
run_pair_test_stats <- function(de_table.ref.marked, the_test_group, groupA, groupB) {
   
   groupA_ranks <- de_table.ref.marked %>% 
      dplyr::filter(.data$test_group == the_test_group, .data$group==groupA ) %>%
      dplyr::pull(.data$rescaled_rank)
   
   groupB_ranks <- de_table.ref.marked %>% 
      dplyr::filter(.data$test_group == the_test_group, .data$group==groupB ) %>%
      dplyr::pull(.data$rescaled_rank)
   
   
   if (stats::median(groupA_ranks) > stats::median(groupB_ranks)) {
      stop("Running test comparision in wrong direction (B < A) this shouldn't happen")
   }
   
   suppressMessages(suppressWarnings(
      mwtest_res <- stats::wilcox.test(groupA_ranks,groupB_ranks, alternative = "less",  
                                       exact=FALSE, conf.int=FALSE) ))
   
   return( dplyr::bind_cols("groupA"     = groupA,
                            "groupB"     = groupB,
                            "meandiff"   = base::mean(groupA_ranks)    - base::mean(groupB_ranks),
                            "mediandiff" = stats::median(groupA_ranks) - stats::median(groupB_ranks),
                            #"CI95"      = as.numeric(mwtest_res$conf.int[2]),  # Mabye useful later, but don't bother now.
                            "pval"       = base::as.numeric(mwtest_res$p.value)))
   
}




