


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
#' @param de_table.ref.marked The output of \code{\link{get_the_up_genes_for_all_possible_groups}} for the contrast of interest.
#' @param the_test_dataset A short meaningful name for the experiment. (Should match \emph{test_dataset} column in \bold{de_table.marked})
#' @param the_ref_dataset A short meaningful name for the experiment. (Should match \emph{dataset} column in \bold{de_table.marked})
#' @param pval Differences between the rescaled ranking distribution of 'top'
#' genes on different reference groups are tested with a Mann-Whitney U test. If 
#' they are \emph{significantly different}, only the top group(s) are reported. 
#' It isn't a simple cutoff threshold as it can change the number of similar groups reported. 
#' ie. A more stringent \bold{pval} is more likely to decide that groups are similar -
#' which would result in multiple group reporting, or no similarity at all.
#' Unlikely that this parameter will need to change. Default = 0.01. 
#' @param num_steps After ranking reference groups according to median 'top' gene
#' ranking, how many adjacent pairs to test for differences. 
#' Set to 1 to only compare each group to the next, or NA to perform an 
#' all-vs-all comparison. 
#' Setting too low may means it is possible to miss groups with some similarity to 
#' the reported matches (\emph{similar_non_match} column)).
#' Too high (or NA) with a large number of reference groups could be slow. 
#' Default = 5.
#' @return A table of automagically-generated labels for each query group, given 
#' their similarity to reference groups. ADD OUTPUT DESCRIPTION HERE XXXX.
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
                                                 the_ref_dataset,  pval=0.01, num_steps=5){
   
   
   # Run the contrasts, as many as requested.
   test_groups <- base::unique(de_table.ref.marked$test_group)
   mwtest_res_table <- dplyr::bind_rows(base::lapply(FUN=get_ranking_and_test_results, X=test_groups, 
                                                     de_table.ref.marked=de_table.ref.marked, 
                                                     the_test_dataset=the_test_dataset, 
                                                     the_ref_dataset = the_ref_dataset,
                                                     num_steps=num_steps,
                                                     pval=pval))
   # make labels
   labels_by_similarity_table <- dplyr::bind_rows(
      base::lapply(FUN=make_ref_similarity_names_for_group, test_groups,
                   mwtest_res_table=mwtest_res_table, de_table.ref.marked=de_table.ref.marked,
                   the_test_dataset=the_test_dataset, the_ref_dataset=the_ref_dataset,
                   the_pval=pval))
   
   return(labels_by_similarity_table)
}




#' make_ref_similarity_names_for_group
#'
#' Internal function, called by make_ref_similarity_names_for_groups 
#' for each group.
#'
#'
#' @param the_test_group Query group to make name for
#' @param mwtest_res_table Mann-whitney test results as constructed 
#' in \code{\link[celaref]{make_ref_similarity_names_for_groups}}
#' @param de_table.ref.marked The output of \code{\link{get_the_up_genes_for_all_possible_groups}} for the contrast of interest.
#' @param the_test_dataset A short meaningful name for the experiment. (Should match \emph{test_dataset} column in \bold{de_table.marked})
#' @param the_ref_dataset A short meaningful name for the experiment. (Should match \emph{dataset} column in \bold{de_table.marked})
#' @param the_pval pval as per \code{\link[celaref]{make_ref_similarity_names_for_groups}}
#' @return A tibble with just one group's labelling information, as per 
#' \code{\link[celaref]{make_ref_similarity_names_for_groups}}
#' 
#'
#' @seealso \code{\link[celaref]{make_ref_similarity_names_for_groups}} Only place that uses this function, details there.
#'
make_ref_similarity_names_for_group <- function(the_test_group, mwtest_res_table, de_table.ref.marked, the_test_dataset, 
                                                the_ref_dataset, the_pval){
   
   matches_col            = NA
   pval_col               = NA
   sim_outside_col        = NA
   sim_outside_detail_col = NA
   differences_within_col = NA
   stepped_pvals_col      = NA
   
   # Just this group
   mwtest_res_table.this <- mwtest_res_table %>% dplyr::filter(.data$test_group == the_test_group)
   
   #test_group test_dataset      ref_dataset             group median_rank n grouprank        next_group meandiff mediandiff       pval
   #1     Group1 a_demo_query a_demo_reference     Weird subtype        0.03 5         1          Exciting   -0.722      -0.93 0.01078587
   #2     Group1 a_demo_query a_demo_reference     Weird subtype        0.03 5         1             dunno   -0.340      -0.25 0.12457632
   #3     Group1 a_demo_query a_demo_reference     Weird subtype        0.03 5         1 Mystery cell type   -0.172      -0.16 0.08660856
   #4     Group1 a_demo_query a_demo_reference Mystery cell type        0.19 5         2     Weird subtype    0.172       0.16 0.94196303
   #5     Group1 a_demo_query a_demo_reference Mystery cell type        0.19 5         2          Exciting   -0.550      -0.77 0.03745645
   
   
   # First, get the stepped pvals 
   # Irrespective of any similarities/significance - this is information 
   # for when theres 'no similarity'.
   #dunno:0.000264,Weird subtype:0.226,Exciting:0.00705,Mystery cell type:NA
   mwtest_res_table.this.stepped <- mwtest_res_table.this %>% 
      dplyr::filter(.data$step == 1) %>% 
      dplyr::arrange(.data$grouprank) 
   stepped_pvals_col=paste0(mwtest_res_table.this.stepped$group, ":", base::signif(mwtest_res_table.this.stepped$pval, 3),
                            collapse=",")
   stepped_pvals_col=paste0(stepped_pvals_col, ",",
                            utils::tail(mwtest_res_table.this.stepped$next_group, 1),":NA")
   
   
   # Any sig in adjacent?
   # Still how a match is defined - will warn when assumptions are violated.
   mwtest_res_table.last_of_top <- mwtest_res_table.this.stepped  %>% 
      dplyr::filter(.data$pval <= the_pval) %>%  #siginifcant
      dplyr::slice(1) #top, is already sorted
   last_of_match <- ifelse(base::nrow(mwtest_res_table.last_of_top) == 0, NA, dplyr::pull(mwtest_res_table.last_of_top, .data$grouprank)[1])
   
   if (!is.na(last_of_match)) {
      
      mwtest_res_table.this$inmatch     <- mwtest_res_table.this$grouprank <= last_of_match
      mwtest_res_table.this$nextinmatch <- mwtest_res_table.this$grouprank + mwtest_res_table.this$step <= last_of_match
      matches <- mwtest_res_table.this %>% dplyr::filter(.data$inmatch==TRUE) %>% dplyr::pull(.data$group) %>% base::unique()
      
      # define the shortname
      matches_col = paste0(the_test_group,":",paste(matches, collapse="|"))
      pval_col    = base::signif(dplyr::pull(mwtest_res_table.last_of_top, .data$pval), 3)
      
      
      # Any significant differences within match (from higher to lower)
      if (length(matches) > 1) {
         
         de_table.ref.marked.thismatch <- de_table.ref.marked %>%
            dplyr::filter(.data$test_group==the_test_group, .data$group %in% matches) 
         
         # Do an all vs all test, just within the matches.
         # Care if sig in either direction
         mwtest_res_table.inmatches <- get_ranking_and_test_results(
            the_test_group = the_test_group,
            de_table.ref.marked = de_table.ref.marked.thismatch , 
            the_test_dataset = the_test_dataset, 
            the_ref_dataset = the_ref_dataset,
            num_steps=NA,
            pval=the_pval)
         
         differences_within_col = ""
         diffs_in_match <- mwtest_res_table.inmatches %>% dplyr::filter(.data$pval <= the_pval)
         if (base::nrow(diffs_in_match) > 0 ) {
            differences_within_col = paste0(diffs_in_match$group, " > ",diffs_in_match$next_group, 
                                            " (p=",base::signif(diffs_in_match$pval,3),")", 
                                            collapse="|")
         }

         
      }
      
      # Any non-signifiant difference between match and outside of it.
      # These are borderline maybe matches that aren't reported as matches 
      # because there's a significant difference above them.
      sim_outside_match <- mwtest_res_table.this %>% 
         dplyr::filter(.data$inmatch == TRUE, .data$nextinmatch == FALSE) %>%
         dplyr::filter(.data$pval > the_pval)
      
      if (base::nrow(sim_outside_match) > 0 ) {
         sim_outside_col        = paste0(base::unique(sim_outside_match$next_group), collapse="|")
         sim_outside_detail_col = paste0(sim_outside_match$group, " > ",sim_outside_match$next_group, 
                                         " (p=",base::signif(sim_outside_match$pval,3),",n.s)", 
                                         collapse="|")
      }
      
   } else {
      # No similarity
      matches_col = paste0(the_test_group,":No similarily")
   }
   
   labels_by_similarity_table <- tibble::tibble(
      test_group        = the_test_group,
      shortlab          = matches_col,
      pval              = pval_col,
      stepped_pvals     = stepped_pvals_col,
      similar_non_match = sim_outside_col,
      similar_non_match_detail = sim_outside_detail_col,
      differences_within       = differences_within_col ) 
   
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
#' @param num_steps see \link[celaref]{make_ref_similarity_names_for_groups}
#'  
#' @seealso \code{\link[celaref]{make_ref_similarity_names_for_groups}} which calls this. 
#' @importFrom magrittr %>%
#' @importFrom rlang .data
get_ranking_and_test_results <- function (de_table.ref.marked, the_test_group, the_test_dataset, the_ref_dataset, num_steps, pval=0.01) {
   
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
   
   
   ## Make a table of the contrasts to test
   # All vs all, or some subset.
   pairset <- NA
   if (is.na(num_steps) || num_steps == 0 ) { # All vs all.
      pairset <- base::expand.grid('groupArank'=1:last_rank, 'groupBrank'=1:last_rank) 
      pairset <- tibble::as.tibble(pairset) %>% dplyr::filter(.data$groupArank!=.data$groupBrank) # no self contrast
   } else { 
      # Shortcut table - only Group to Group+Nth rank check
      # for stepped contrats num_steps=1.
      # from n+1 to n+num_steps or end 
      get_nstep_ranks <- function(the_rank, num_steps, last_rank) {
         if (the_rank >= last_rank) {return (data.frame())} #empty, not a NA
         nstep_rank_vals <- (the_rank+1):min(the_rank + num_steps, last_rank)
         return ( 
            dplyr::bind_cols('groupArank'=as.integer(rep(the_rank, times=base::length(nstep_rank_vals))),
                             'groupBrank'=nstep_rank_vals)
         ) 
      }
      pairset <- dplyr::bind_rows(base::lapply(FUN=get_nstep_ranks, rankstat_table$grouprank, num_steps=num_steps, last_rank=last_rank))
   }
   
   
   # ranks back to names
   rank2group <- stats::setNames(rankstat_table$group, rankstat_table$grouprank)
   pairset$groupA <- rank2group[pairset$groupArank]
   pairset$groupB <- rank2group[pairset$groupBrank]  
   pairset$step <- pairset$groupBrank - pairset$groupArank
   
   # calculate diferences
   mwtest_res_table <- pairset %>% dplyr::rowwise() %>% 
      dplyr::do (run_pair_test_stats(de_table.ref.marked=de_table.ref.marked, 
                                     the_test_group, .data$groupA, .data$groupB,
                                     enforceAgtB=FALSE))
   
   # Collect information
   rankstat_table$group <- as.character(rankstat_table$group)
   ranking_and_mwtest_results <- base::merge(x=rankstat_table, y=mwtest_res_table, 
                                             by.x="group", by.y="groupA", all.x=TRUE) 
   
   # pretty names    
   ranking_and_mwtest_results <- ranking_and_mwtest_results %>% 
      dplyr::arrange(.data$grouprank) %>%
      dplyr::rename(next_group=.data$groupB) %>% 
      tibble::add_column( test_group   = the_test_group, 
                          test_dataset = the_test_dataset, 
                          ref_dataset  = the_ref_dataset, .before=1)
   
   # Add step size 
   ranking_and_mwtest_results <- base::merge(x=ranking_and_mwtest_results, 
                                             y=pairset %>% dplyr::select(.data$groupA, .data$groupB, .data$step),
                                             by.x=c("group", "next_group"),
                                             by.y=c("groupA", "groupB"), all=TRUE)
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
#' @param enforceAgtB Do a one tailed test of A 'less' B (more similar)? Or two-tailed. Default = TRUE.
#'
#' @return A tibble of wilcox / man-whitneyU test results for this contrast.
#'
#' @seealso  \code{\link[celaref]{make_ref_similarity_names_for_groups}} 
#'
#' @importFrom magrittr %>%
run_pair_test_stats <- function(de_table.ref.marked, the_test_group, groupA, groupB, enforceAgtB=TRUE) {
   
   groupA_ranks <- de_table.ref.marked %>% 
      dplyr::filter(.data$test_group == the_test_group, .data$group==groupA ) %>%
      dplyr::pull(.data$rescaled_rank)
   
   groupB_ranks <- de_table.ref.marked %>% 
      dplyr::filter(.data$test_group == the_test_group, .data$group==groupB ) %>%
      dplyr::pull(.data$rescaled_rank)
   
   
   if (enforceAgtB) {
      if (stats::median(groupA_ranks) > stats::median(groupB_ranks)) {
         stop("Running test comparision in wrong direction (B < A) this shouldn't happen")
      }
   }
   
   suppressMessages(suppressWarnings(
      mwtest_res <- stats::wilcox.test(groupA_ranks,groupB_ranks, alternative = "less",  
                                       exact=FALSE, conf.int=FALSE) ))
   
   return( dplyr::bind_cols("groupA"     = groupA,
                            "groupB"     = groupB,
                            "meandiff"   = base::mean(groupA_ranks)    - base::mean(groupB_ranks),
                            "mediandiff" = stats::median(groupA_ranks) - stats::median(groupB_ranks),
                            "pval"       = base::as.numeric(mwtest_res$p.value)))
   
}




