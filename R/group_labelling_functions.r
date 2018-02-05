



#' make_ref_similarity_names_for_groups_ks
#'
#' Given the comparison between the two dataset, construct some sensible labels 
#' for the groups in the query group.
#' 
#' This function aims to report a) the top most similar reference group, if 
#' there's a clear frontrunner, b) A list of multiple similar groups if they 
#' have similar similarity, or c) 'No similarity', if there is none.
#' 
#' Each group is named according to the following rules:
#' \enumerate{
#'   \item There is \emph{one} reference group under the \bold{median_rank_threshold}: Use that one.
#'   \item There are \emph{multiple} reference group under the
#'   \bold{median_rank_threshold}: Run a Kolmogorov-Smirnov test on the 
#'   \bold{rescaled_rank} distribution of query 'top' genes in the 2 most similar reference groups.
#'   \itemize{
#'      \item If there is a signifiant (<= \bold{ks_pval}) \emph{difference between 
#'      the two distributions} - report the first only. Because they're similar,
#'      but different enough to tell apart.
#'      \item If the distributions are \emph{not different}, report both. If there are 
#'      more groups under the \bold{median_rank_threshold} start stepping down 
#'      group by group to see if there's a difference somewhere (e.g. 3 groups 
#'      under threshold, but its possible to return only the top 2.).
#'   }
#'   \item There is nothing under the \bold{median_rank_threshold} : No simiarity. 
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
#' @param de_table.marked The output of \code{\link{get_the_up_genes_for_all_possible_groups}} for the contrast of interest.
#' @param the_test_dataset A short meaningful name for the experiment. (Should match \emph{test_dataset} column in \bold{de_table.marked})
#' @param the_ref_dataset A short meaningful name for the experiment. (Should match \emph{dataset} column in \bold{de_table.marked})
#' @param ks_pval When multiple reference groups pass the median rank threshold, 
#' the difference between them is testes with a Kolmogorov-Smirnov test - if 
#' they are \emph{significantly different}, only the top group(s) are reported. 
#' ie. A more stringent \bold{ks_pval} is more likely to report multiple 
#' reference groups. 
#' Note that comparisons with small numbers of 'top' genes will never be 
#' 'significantly' different - these should be taken with a grain of salt 
#' anyway. This parameter unlikely to need change. Default = 0.01.
#' #' @param median_rank_threshold A threshold for similarity. Number betweeen 0 
#' (top-rank) and 1.  If the median ranking of the 'top' query genes is 
#' below this, consider possible similarity. 
#' Consider reducing for very similar cell clusters. Note that 0.5 is 
#' theoretically random (but do occur, see detail), so values larger than that
#' are pointless. Default 0.25.
#'  
#'
#' @return A table of automagically-generated labels for each query group, given 
#' their similarity to reference groups. Short_lab and long_lab are just 
#' different versions of the name.
#' 
#'
#' @examples
#' 
#' de_table.demo_query <- contrast_each_group_to_the_rest(demo_query_se, "demo_query")
#' de_table.demo_ref   <- contrast_each_group_to_the_rest(demo_ref_se,   "demo_ref")
#' de_table.marked.query_vs_ref <- get_the_up_genes_for_all_possible_groups(de_table.demo_query, de_table.demo_ref, 'demo_query')
#'
#' make_ref_similarity_names_for_groups_ksfunction(de_table.marked.query_vs_ref,
#'                               the_test_dataset="demo_query",
#'                               the_ref_dataset="demo_ref")
#'
#' @seealso \code{\link[celaref]{get_the_up_genes_for_all_possible_groups}} To prepare the \bold{de_table.ref.marked} input.
#' 
#' @importFrom magrittr %>%
#'@export
make_ref_similarity_names_for_groups_ks <- function(de_table.ref.marked, the_test_dataset, the_ref_dataset,
                                                         ks_pval=0.01, median_rank_threshold=0.25){
   
   test_groups <- base::unique(de_table.ref.marked$test_group)
   
   #rankstat_table <- get_rankstat_table(de_table.ref.marked, the_test_group)
   #NB: ks_res table won't include anything that doesn't have any similarity
   ks_res_table <- dplyr::bind_rows(lapply(FUN=get_ranking_and_ks_test_results, X=test_groups, 
                                    de_table.ref.marked=de_table.ref.marked, 
                                    the_test_dataset=the_test_dataset, the_ref_dataset = the_ref_dataset,
                                    ks_pval=ks_pval, median_rank_threshold=median_rank_threshold))
   
   
   
   # Just for when there are labells to be made
   labels_by_similarity_table <- ks_res_table %>% 
      dplyr::filter(in_name==TRUE) %>%
      dplyr::arrange(test_group, grouprank ) %>%
      dplyr::group_by(test_group) %>%
      dplyr::mutate(shortlab  = paste0(test_group, ":", paste0(group, collapse = "|")), 
                    longlab  = paste0(test_group, ":", paste0(group, collapse = "|"),":",ref_dataset,"_similarity")  ) %>%
      dplyr::select(test_group, shortlab, longlab) %>% base::unique()
   
   
   # Add no-similarity groups back in.
   if (any ( ! test_groups %in% labels_by_similarity_table$test_group)) {
      test_groups.nomatch <- test_groups[! test_groups %in% labels_by_similarity_table$test_group]
      
      labels_by_similarity_table.nomatch <- 
         dplyr::bind_cols("test_group"=test_groups.nomatch, 
                      shortlab=paste0(test_groups.nomatch, ":no_similarity"),
                      longlab=paste0(test_groups.nomatch, ":",the_ref_dataset,"_no_similarity"))
      labels_by_similarity_table <- dplyr::bind_rows(labels_by_similarity_table, labels_by_similarity_table.nomatch)
   }
   
   return(labels_by_similarity_table)
}






#' get_ranking_and_ks_test_results
#'
#' Internal function to get reference group similarity contrasts for an individual query qroup. 
#' 
#' For use by \bold{make_ref_similarity_names_for_groups_ks}, see that function for parameter details.
#' 
#' @param the_test_group Additional parameter for which query group is being tested.
#'
#' @seealso \code{\link[celaref]{make_ref_similarity_names_for_groups_ks}} which calls this. 
#' @importFrom magrittr %>%
get_ranking_and_ks_test_results <- function (de_table.ref.marked, the_test_group, the_test_dataset, the_ref_dataset, ks_pval=0.01,  median_rank_threshold=0.25) {
   
   # Anything that each stepwise ks test is NOT significant (with top contrast over threshold) :
   # Second might not be.
   # AAAA |  
   #  BBBB|B
   #    CC|C
   # <- top(0)
   # A+B+C
   
   # AAAAA|AA
   #     B|BBB
   # A+B    
   
   # Or anything that is ks test NOT sig with Top hit
   # (does this happen? AN easy test anyway.)
   # AAAAAAAA |
   #    BB    |
   #      CC  | 
   
   # Should only be one test datast, and one reference dataset in de_table.ref_marked. Not needed, but forced here for paranoia reasons.
   de_table.ref.marked <- de_table.ref.marked %>% dplyr::filter(test_dataset==the_test_dataset, dataset==the_ref_dataset )
   
   
   # Get the average rankings in order.
   # For this test_group
   rankstat_table <- get_rankstat_table(de_table.ref.marked, the_test_group)
   #                             group  median_rank     n grouprank
   #                        <fctr>        <dbl> <int>     <chr>
   #1                      c3-Tcells60 0.0002933584     9         1
   #2                     c7-NKcells10 0.0266956114     9         2
   #3 c6-Monocytes20_Tcells10_Bcells10 0.0543299695     9         3
   
   ## Which reference groups are above threshold?
   if (rankstat_table$median_rank[1] > median_rank_threshold) {
      return(NULL)     # If top under the theshold, can't say anything,
   }
   last_above_median_rank <- max(which(rankstat_table$median_rank <= median_rank_threshold))
   
   # Sanity checks    
   if (nrow(rankstat_table) ==1) { warn("Only 1 reference group. Something odd is happening.") }
   if (last_above_median_rank == nrow(rankstat_table)) { 
      stop("All references above median rank?? Something odd is happening. Can't handle.") 
   }
   
   
   ## Run 'stepped' contrasts.
   # from this ranking, Pair first with second, secont to third, e.t.c.
   # But only go down to groupA meeting meadian rank thresohold (and compare to first not meeting threshold)
   pairset.step <- dplyr::bind_cols(groupA=as.character(rankstat_table$group[c(1:last_above_median_rank)]),
                                    groupB=as.character(rankstat_table$group[c(2:(last_above_median_rank+1))])) 
   #pairset.step$contrast_rank <- 1:last_above_median_rank
   #pairset.step$contrast_run <- "step"
   
   ks_res_table.step <- pairset.step %>% dplyr::rowwise() %>% 
      dplyr::do (run_pair_ks_stats(de_table.ref.marked=de_table.ref.marked, the_test_group, .$groupA, .$groupB))
   
   
   
   
   ## Run 'first' contrasts
   # everything (even ns) vs top ranked group
   num_groups    <- base::nrow(rankstat_table)
   pairset.first <- dplyr::bind_cols(groupA=rep(as.character(rankstat_table$group[1]), num_groups-1),
                                     groupB=as.character(rankstat_table$group[c(2:(num_groups))])) 
   ks_res_table.first <- pairset.first %>% dplyr::rowwise() %>% 
      dplyr::do (run_pair_ks_stats(de_table.ref.marked=de_table.ref.marked, the_test_group, .$groupA, .$groupB))
   
   
   
   # Run checks
   rankstat_table$group <- as.character(rankstat_table$group)
   ranking_and_ks_results <- merge(x=rankstat_table,y=ks_res_table.first, 
                                   by.x="group", by.y="groupB", all.x=TRUE)
   ranking_and_ks_results <- merge(x=ranking_and_ks_results, y=ks_res_table.step, 
                                   by.x="group", by.y="groupB",suffixes=c("_first","_step"), all.x=TRUE)                        
   ranking_and_ks_results <- ranking_and_ks_results %>% dplyr::arrange(grouprank) %>%
      dplyr::rename(top_contrast_group_first=groupA_first) %>%  dplyr::rename(up_contrast_group_step=groupA_step) %>% 
      tibble::add_column( test_group=the_test_group, test_dataset=the_test_dataset, ref_dataset=the_ref_dataset,
                  .before=1)
   
   
   #                             group  median_rank n grouprank top_contrast_group_first   D_first   pval_first
   #1                      c3-Tcells60 0.0002933584 9         1                     <NA>        NA           NA
   #2                     c7-NKcells10 0.0266956114 9         2              c3-Tcells60 1.0000000 0.0001234098
   #            up_contrast_group_step    D_step    pval_step
   #1                             <NA>        NA           NA
   #2                      c3-Tcells60 1.0000000 0.0001234098
   
   
   
   # assume a non-sig differnece is similar (ie maybe drown from same distribution)
   ranking_and_ks_results$similar_to_first <- ifelse (!is.na(ranking_and_ks_results$pval_first), ! ranking_and_ks_results$pval_first <= ks_pval, FALSE)
   
   # Step by step where is a~b, b~c, c~d ... stop being *potentially* from same ditribution
   # NB: upper part of each conrast meets threshold already, there are NAs thereafter because we can ignore those.
   similar_to_one_step_up <- c(TRUE, !(ranking_and_ks_results$pval_step[-1] <= ks_pval)) #top doesn't get contrast :. T
   # Where is the last step that is similar to the one above it? Usually 1.
   last_step_not_na <- max(which(! is.na(similar_to_one_step_up)))  # Last that isn't NA.
   step_dropoff <- last_step_not_na 
   if (any(!similar_to_one_step_up[1:last_step_not_na]) )  { # Any false/not similar
      step_dropoff <- min(which(!similar_to_one_step_up[1:last_step_not_na]))-1
   } 
   ranking_and_ks_results$similar_via_steps <- FALSE
   ranking_and_ks_results$similar_via_steps[1:step_dropoff] <- TRUE
   
   # Either of those criteria, or just being the top hit - want to include in a label
   ranking_and_ks_results$in_name     <- ranking_and_ks_results$similar_to_first | ranking_and_ks_results$similar_via_steps 
   ranking_and_ks_results$in_name[1] <- TRUE # Sorted and, already know first was above threshold
   
   # median threshold now a hard limit. 
   ranking_and_ks_results$in_name[ranking_and_ks_results$median_rank > median_rank_threshold] <- FALSE
   
   return(ranking_and_ks_results)
}




#' get_rankstat_table 
#' 
#' Summarise the comparison of the specified query group against in the 
#' comparison in \bold{de_table.ref.marked} - number of 'top' genes and their 
#' median rank in each of the reference groups, with reference group rankings.
#'
#' @param de_table.marked The output of \code{\link{get_the_up_genes_for_all_possible_groups}} for the contrast of interest.
#' @param the_test_group Name of query group to test
#'
#' @return A tibble of query group name (test_group), number of 'top' genes (n), 
#' reference dataset group (group) with its ranking (grouprank) and the median 
#' (rescaled 0..1) ranking of 'top' genes (median_rank).
#'
#' @examples
#'
#' de_table.demo_query <- contrast_each_group_to_the_rest(demo_query_se, "demo_query")
#' de_table.demo_ref   <- contrast_each_group_to_the_rest(demo_ref_se,   "demo_ref")
#' de_table.marked.query_vs_ref <- get_the_up_genes_for_all_possible_groups(de_table.demo_query, de_table.demo_ref, 'demo_query')
#'
#' get_rankstat_table(de_table.marked.query_vs_ref, "Group3")
#'
#' @seealso \code{\link[celaref]{get_the_up_genes_for_all_possible_groups}} To prepare the \bold{de_table.ref.marked} input.
#'
#' @export
#' @importFrom magrittr %>%
get_rankstat_table <- function(de_table.ref.marked, the_test_group){
   
   # Get the average rankings in order.
   # For this test_group
   rankstat_table <- de_table.ref.marked %>%
      dplyr::group_by(group) %>%
      dplyr::filter(test_group == the_test_group) %>% #only test group
      dplyr::summarise(median_rank=median(rescaled_rank),
                n=n() ) %>% #n() is part of dplyr, but yeilds error if specified as dplyr::n 'Evaluation error: This function should not be called directly.'
      dplyr::arrange(median_rank)
   rankstat_table$grouprank <- base::rownames(rankstat_table)
   return(rankstat_table)
}




#' run_pair_ks_stats
#'
#' Internal function to compare the distribution of a query datasets 'top' 
#' genes between two different reference datasete groups with a KS test.
#'
#' For use by make_ref_similarity_names_for_groups_ks
#'
#' @param de_table.marked The output of \code{\link{get_the_up_genes_for_all_possible_groups}} for the contrast of interest.
#' @param the_test_group Name of the test group in query dataset.
#' @param groupA One of the reference group names
#' @param groupB Another of the reference group names
#'
#' @return A tibble of KS test results for this contrast (pval, D value).
#'
#' @seealso  \code{\link[celaref]{make_ref_similarity_names_for_groups_ks}} 
#'
#' @importFrom magrittr %>%
run_pair_ks_stats <- function(de_table.ref.marked, the_test_group, groupA, groupB) {
   
   groupA_ranks <- de_table.ref.marked %>% 
      dplyr::filter(test_group == the_test_group, group==groupA ) %>%
      dplyr::pull(rescaled_rank)
   
   groupB_ranks <- de_table.ref.marked %>% 
      dplyr::filter(test_group == the_test_group, group==groupB ) %>%
      dplyr::pull(rescaled_rank)
   
   # even specificed exact==FALSE, I still get warning about presence of ties.
   suppressMessages(suppressWarnings( ks.res <- ks.test(groupA_ranks,groupB_ranks, alternative = "greater", exact=FALSE)))
   
   return( dplyr::bind_cols("groupA"=groupA,"groupB"=groupB,
                     "D"=as.numeric(ks.res$statistic),  "pval"=as.numeric(ks.res$p.value)))
}






