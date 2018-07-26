


#' make_ref_similarity_names
#'
#' Construct some sensible labels or the groups/clusters in the query dataset, 
#' based on similarity the reference dataset.
#' 
#' This function aims to report a) the top most similar reference group, if 
#' there's a clear frontrunner, b) A list of multiple similar groups if they 
#' have similar similarity, or c) 'No similarity', if there is none.
#' 
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
#' There are some further heuristic caveats:
#' \enumerate{
#'   \item The distribution of top genes in the last (or only) match group is tested 
#'   versus a theroetical random distribution around 0.5 (as reported in 
#'   \emph{pval_vs_random} column). If the distribution is not significantly above
#'   random  (It is possible in edge cases where there is a skewed dataset and no/few matches),
#'   \emph{no similarity} is reported. The significnat \emph{pval} column is left intact.
#'   \item The comparison is repeated reciprocally - reference groups vs the 
#'   query groups. This helps sensitivity of heterogenous query groups - and investigating 
#'   the reciprocal matches can be informative in these cases.
#'   If a query group doens't 'match' a reference group, but the reference 
#'   group does match that query group - it is reported in the group label in brackets.
#'   e.g. \emph{c1:th_lymphocytes(tc_lympocytes)}. 
#'   Its even possible if there was no match (and pval = NA) e.g. emph{c2:(tc_lymphocytes)}
#' }
#' 
#' 
#' 
#' The similarity is formatted into a group label. Where there are 
#' multiple similar groups, they're listed from most to least similar by their 
#' median ranks.
#'  
#' For instance, a query dataset of clusters c1, c2, c3 and c4 againsts a cell-type
#' labelled reference datatset might get names like:
#' E.g.
#' \itemize{
#'   \item c1:macrophage
#'   \item c2:endotheial|mesodermal
#'   \item c3:no_similarity
#'   \item c4:mesodermal(endothelial)
#' }
#'
#' Function \code{make_ref_similarity_names} is a convenience wrapper function for 
#' \code{make_ref_similarity_names_from_marked}. It accepts two SummarisedExpression
#' objects to compare and handles running
#' \code{\link[celaref]{get_the_up_genes_for_all_possible_groups}}. 
#' Sister function \code{make_ref_similarity_names_from_marked} may (rarely) be of use 
#' if the \bold{de_table.marked} object has already been created, 
#' or if reciprocal tests are not wanted. 
#' 
#' @param de_table.test A differential expression table of the query experiment,
#' as generated from \code{\link[celaref]{contrast_each_group_to_the_rest}}
#' @param de_table.ref A differential expression table of the reference dataset,
#' as generated from \code{\link[celaref]{contrast_each_group_to_the_rest}}
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
#' 
#' @return A table of automagically-generated labels for each query group, given 
#' their similarity to reference groups. 
#' 
#' The columns in this table:
#' \itemize{
#'   \item \bold{test_group} : Query group e.g. "c1"
#'   \item \bold{shortlab} : The cluster label described above e.g. "c1:macrophage"
#'   \item \bold{pval} : If there is a similarity flagged, this is the P-value from 
#'   a Mann-Whitney U test from the last 'matched' group to the adjacent  
#'   'non-matched' group. Ie. If only one label in shortlab, this will be the first 
#'   of the stepped_pvals, if there are 2, it will be the second. If there is 'no_similarity'
#'   this will be NA (Because there is no confidence in what is the most appropriate 
#'   of the all non-significant stepped pvalues.).
#'   \item \bold{stepped_pvals} : P-values from Mann-Whitney U tests across adjacent
#'   pairs of reference groups ordered from most to least similar (ascending median rank).
#'   ie. 1st-2nd most similar first, 2nd-3rd, 3rd-4th e.t.c. The last value will always 
#'   be NA (no more reference group).
#'   e.g. refA:8.44e-10,refB:2.37e-06,refC:0.000818,refD:0.435,refE:0.245,refF:NA
#'   \item \bold{pval_to_random} : P-value of test of median rank (of last 
#'   matched reference group) < random, from binomial test on top gene ranks (<0.5). 
#'   \item \bold{matches} : List of all reference groups that 'match', 
#'   as described, except it also includes (rare) examples where pval_to_random is not significant.
#'   "|" delimited.
#'   \item \bold{reciprocal_matches} : List of all reference groups that flagged 
#'   test group as a match when directon of comparison is reversed.
#'   (significant pval and pval_to_random). "|" delimited.
#'   \item \bold{similar_non_match}: This column lists any reference groups outside 
#'   of shortlab that are not signifcantly different to a reported match group. 
#'   Limited by \emph{num_steps}, and will never find anything if num_steps==1.
#'   "|" delimited. Usually NA.
#'   \item \bold{similar_non_match_detail} : P-values for any details about 
#'   similar_non_match results. These p-values will always be non-significant.
#'   E.g. "A > C (p=0.0214,n.s)". "|" delimited. Usually NA.
#'   \item \bold{differences_within} :  This feild lists any pairs of reference 
#'   groups in shortlab that are significantly different. "|" delimited. Usually NA.
#' }
#' 
#' 
#'
#' @examples
#' 
#' de_table.demo_query <- contrast_each_group_to_the_rest(demo_query_se, "demo_query")
#' de_table.demo_ref   <- contrast_each_group_to_the_rest(demo_ref_se,   "demo_ref")
#' 
#' make_ref_similarity_names(de_table.demo_query, de_table.demo_ref)
#' make_ref_similarity_names(de_table.demo_query, de_table.demo_ref, num_steps=3)
#' make_ref_similarity_names(de_table.demo_query, de_table.demo_ref, num_steps=NA)
#'
#' @seealso \code{\link[celaref]{contrast_each_group_to_the_rest}} For preparing de_table input
#' 
#' 
#' @importFrom magrittr %>%
#' @export
make_ref_similarity_names <- function(de_table.test, de_table.ref, pval=0.01, num_steps=5) {
   
   # Read one dataset name from test and query. 
   test_dataset <- unique(de_table.test$dataset)
   ref_dataset  <- unique(de_table.ref$dataset)   
   if (length(test_dataset) >1  || length(ref_dataset) > 1) {stop("Either test or reference includes more than one dataset in dataset column. Trim down or use make_ref_similarity_names_using_marked")}
   
   # Prepare the marked and reciprocal marked tables
   de_table.ref.marked   <- get_the_up_genes_for_all_possible_groups(de_table.test, de_table.ref)
   de_table.recip.marked <- get_the_up_genes_for_all_possible_groups(de_table.ref, de_table.test)
   
   return(make_ref_similarity_names_using_marked(de_table.ref.marked, de_table.recip.marked,
                                                 pval=pval, num_steps=num_steps))
   
}



#' make_ref_similarity_names_using_marked 
#'
#' This is a more low level/customisable version of \code{\link{make_ref_similarity_names}},
#' (usually be all that is needed).  
#' Suitable for rare cases to reuse an existing \bold{de_table.ref.marked} object.
#' Or use a \bold{de_table.ref.marked} table with more than one dataset present 
#' (discoraged). Or to skip the reciprocal comparison step. 
#' 
#' @param de_table.ref.marked The output of \code{\link{get_the_up_genes_for_all_possible_groups}} for the contrast of interest.
#' @param de_table.recip.marked Optional. The (reciprocal) output of \code{\link{get_the_up_genes_for_all_possible_groups}}
#' with the test and reference datasets swapped. 
#' If omitted a reciprocal test will not be done. Default = NA.
#' @param the_test_dataset Optional. A short meaningful name for the experiment. 
#' (Should match \emph{test_dataset} column in \bold{de_table.marked}). 
#' Only needed in a table of more than one dataset. Default = NA.
#' @param the_ref_dataset Optional. A short meaningful name for the experiment. 
#' (Should match \emph{dataset} column in \bold{de_table.marked}). 
#' Only needed in a table of more than one dataset. Default = NA.
#' 
#' @examples
#'
#' de_table.demo_query <- contrast_each_group_to_the_rest(demo_query_se, "demo_query")
#' de_table.demo_ref   <- contrast_each_group_to_the_rest(demo_ref_se,   "demo_ref")
#' 
#' de_table.marked.query_vs_ref <- get_the_up_genes_for_all_possible_groups(
#'      de_table.demo_query, de_table.demo_ref) 
#' de_table.marked.reiprocal <- get_the_up_genes_for_all_possible_groups(
#'      de_table.demo_ref, de_table.demo_query)
#'      
#'
#' make_ref_similarity_names_using_marked(de_table.marked.query_vs_ref, de_table.marked.reiprocal)
#' make_ref_similarity_names_using_marked(de_table.marked.query_vs_ref)
#' 
#'
#' @seealso \code{\link[celaref]{get_the_up_genes_for_all_possible_groups}} To prepare the \bold{de_table.ref.marked} input.
#' 
#' @describeIn make_ref_similarity_names Construct some sensible cluster labels, but using a premade marked table.  
#'
#' @export
make_ref_similarity_names_using_marked <- function(de_table.ref.marked, de_table.recip.marked=NA,
                                                   the_test_dataset=NA,the_ref_dataset=NA, 
                                                   pval=0.01, num_steps=5){
   
   # Sanity checks with nicer errors
   # datatset designations may be omoitted, IF only one dataset in the table (normal)
   ref_datasets_present  <- unique(de_table.ref.marked$dataset)
   test_datasets_present <- unique(de_table.ref.marked$test_dataset)   
   if (is.na(the_ref_dataset)) {
      if (base::length(ref_datasets_present) != 1 ) {stop("More than one reference datatset in de_table.ref.marked, specify with the_ref_dataset")}
      the_ref_dataset=ref_datasets_present[1]
   } 
   if (is.na(the_test_dataset)) {
      if (base::length(test_datasets_present) != 1 ) {stop("More than one reference datatset in de_table.ref.marked, specify with the_test_dataset")}
      the_test_dataset=test_datasets_present[1]
   } 
   
   # Also it has to be present if defined!
   if ( ! the_ref_dataset  %in% ref_datasets_present) {stop("Cannot find specified ref dataset in de_table.ref.marked")}
   if ( ! the_test_dataset %in% test_datasets_present) {stop("Cannot find specified ref dataset in de_table.ref.marked")}
   
   # Likewise, the reciprocal tests.
   have_recip <- typeof(de_table.recip.marked) == "list" #NA is logical
   if (have_recip) {
      if ( ! the_ref_dataset  %in% unique(de_table.recip.marked$test_dataset)) {stop("Cannot find specified ref dataset as query dataset in de_table.recip.marked")}
      if ( ! the_test_dataset %in% unique(de_table.recip.marked$dataset)) {stop("Cannot find specified test dataset as ref dataset in de_table.recip.marked")}
   }
   
   
   
   # Run the contrasts, as many as requested.
   test_groups <- base::unique(de_table.ref.marked$test_group)
   mwtest_res_table <- dplyr::bind_rows(base::lapply(FUN=get_ranking_and_test_results, X=test_groups, 
                                                     de_table.ref.marked=de_table.ref.marked, 
                                                     the_test_dataset=the_test_dataset, 
                                                     the_ref_dataset = the_ref_dataset,
                                                     num_steps=num_steps,
                                                     pval=pval))
   
   
   reciprocal_matches <- NA
   if (have_recip) {
      recip_test_groups <- base::unique(de_table.recip.marked$test_group) 
      mwtest_res_table.recip <- dplyr::bind_rows(base::lapply(FUN=get_ranking_and_test_results, X= recip_test_groups, 
                                                              de_table.ref.marked=de_table.recip.marked, 
                                                              the_test_dataset= the_ref_dataset, 
                                                              the_ref_dataset = the_test_dataset,
                                                              num_steps=1, # Simplistic check.
                                                              pval=pval))
      
      # Reciprocal matches
      reciprocal_matches <- get_reciprocal_matches(mwtest_res_table.recip, de_table.recip.marked, the_pval=pval) 
   }
   
   ## Make labels
   labels_by_similarity_table <- dplyr::bind_rows(
      base::lapply(FUN=make_ref_similarity_names_for_group, test_groups,
                   mwtest_res_table=mwtest_res_table, de_table.ref.marked=de_table.ref.marked,
                   reciprocal_matches,
                   the_test_dataset=the_test_dataset, the_ref_dataset=the_ref_dataset,
                   the_pval=pval))
   
   return(labels_by_similarity_table)
}









#' make_ref_similarity_names_for_group
#'
#' Internal function, called by make_ref_similarity_names_using_marked 
#' for each group.
#'
#'
#' @param the_test_group Query group to make name for
#' @param mwtest_res_table Mann-whitney test results as constructed 
#' in \code{\link[celaref]{make_ref_similarity_names_using_marked}}
#' @param de_table.ref.marked The output of \code{\link{get_the_up_genes_for_all_possible_groups}} for the contrast of interest.
#' @param reciprocal_matches Simplified table of reciprocal matches prepared within \code{\link{make_ref_similarity_names_using_marked}}. 
#' If omitted no reciprocal matching done. Default = NA.
#' @param the_test_dataset A short meaningful name for the experiment. (Should match \emph{test_dataset} column in \bold{de_table.marked})
#' @param the_ref_dataset A short meaningful name for the experiment. (Should match \emph{dataset} column in \bold{de_table.marked})
#' @param the_pval pval as per \code{\link[celaref]{make_ref_similarity_names_using_marked}}
#' @return A tibble with just one group's labelling information, as per 
#' \code{\link[celaref]{make_ref_similarity_names_using_marked}}
#' 
#'
#' @seealso \code{\link[celaref]{make_ref_similarity_names_using_marked}} 
#' Only place that uses this function, details there.
make_ref_similarity_names_for_group<- function(the_test_group, 
                                                    mwtest_res_table, de_table.ref.marked, 
                                                    reciprocal_matches = NA,
                                                    the_test_dataset, the_ref_dataset, the_pval){
   
   short_lab_col          = NA
   matches_col            = NA
   pval_col               = NA
   not_random_pval_col    = NA
   sim_outside_col        = NA
   sim_outside_detail_col = NA
   differences_within_col = NA
   stepped_pvals_col      = NA
   reciprocal_matches_col = NA
   
   matches               <- NULL
   reciprocal_match_list <- NULL
   
   # Just this group
   mwtest_res_table.this <- mwtest_res_table %>% dplyr::filter(.data$test_group == the_test_group)
   
   # Get a bunch of matches (stepped), or NA
   mwtest_res_table.matches <- get_matched_stepped_mwtest_res_table(mwtest_res_table.this, the_pval)
   
   
   # Reciprocal matches (simple)
   if (any(! is.na(reciprocal_matches))) {
      reciprocal_match_list <- reciprocal_matches$ref_group[reciprocal_matches$test_group == the_test_group]
   }
   
   
   # Much of the other checks only apply when a match is defined. 
   if (any(!is.na(mwtest_res_table.matches))) { 
      
      # This table is ordered by grouprank, starting at 1
      last_of_match <- nrow(mwtest_res_table.matches)
      matches       <- mwtest_res_table.matches$group 
      mwtest_res_table.this$inmatch     <- mwtest_res_table.this$grouprank <= last_of_match
      mwtest_res_table.this$nextinmatch <- mwtest_res_table.this$grouprank + mwtest_res_table.this$step <= last_of_match
      
      ## Pval from last of match group to first of non-match. 
      pval_col              <- base::signif(mwtest_res_table.matches$pval[last_of_match], 3) #the First sig diff 
      
      ## Pval to a random distribution of ranks.
      not_random_pval_col   <- get_vs_random_pval(de_table.ref.marked, the_group=mwtest_res_table.matches$group[last_of_match], the_test_group)
      
      
      ## Any non-signifiant difference between match and outside of it?
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
      
      ## If there are multiple matches, are they different to each other?
      if (length(matches) > 1) {
         differences_within_col <- find_within_match_differences( de_table.ref.marked, matches, the_test_group, the_test_dataset, the_ref_dataset, the_pval)
      }
      
   } 
   
   # Defaine the short name now.
   # <query cluster>: match1|match2
   # <query cluster>: match1(reciprocalmatchfrommatch2) 
   # <query cluster>: No similarity
   reportable_matches = matches
   shortlab_col=paste0(the_test_group,":")
   if (length(matches) > 0 ) {
      if (not_random_pval_col <= the_pval) {
         # Shortlab col won't' include matches unless the last is above random. 
         # (even if sig to next)
         shortlab_col = paste0(shortlab_col,paste(matches, collapse="|"))
      } else {
         # there are matches, but we're not putting them in shortlab
         reportable_matches = character(0) 
      }
   }
   
   extra_recip_matches = base::setdiff(reciprocal_match_list, reportable_matches)
   if (length(extra_recip_matches) > 0 ) {
      shortlab_col = paste0(shortlab_col,"(", paste(extra_recip_matches, collapse="|"), ")")
   }
   if ( (length(matches) == 0 | (length(matches) > 0 & not_random_pval_col > the_pval) ) 
        & length(extra_recip_matches) == 0 ){
      shortlab_col = paste0(shortlab_col,"No similarity")
   }
   
   
   
   # Other cols
   matches_col            <- paste(matches, collapse="|")
   reciprocal_matches_col <- paste0(reciprocal_match_list, collapse="|")
   stepped_pvals_col      <- get_stepped_pvals_str(mwtest_res_table.this)
   
   
   labels_by_similarity_table <- tibble::tibble(
      test_group               = the_test_group,
      shortlab                 = shortlab_col,
      pval                     = pval_col,
      stepped_pvals            = stepped_pvals_col,
      pval_to_random           = not_random_pval_col,
      matches                  = matches_col,
      reciprocal_matches       = reciprocal_matches_col,
      similar_non_match        = sim_outside_col,
      similar_non_match_detail = sim_outside_detail_col,
      differences_within       = differences_within_col ) 
   
   return(labels_by_similarity_table)
   
   
}



   
#' get_ranking_and_test_results
#'
#' Internal function to get reference group similarity contrasts for an individual query qroup. 
#' 
#' For use by \bold{make_ref_similarity_names_using_marked}, see that function for parameter details.
#' This function just runs this for a single query group \bold{the_test_group}
#'
#' @param de_table.ref.marked see \link[celaref]{make_ref_similarity_names_using_marked}
#' @param the_test_group The group to calculate the stats on.
#' @param the_test_dataset see \link[celaref]{make_ref_similarity_names_using_marked}
#' @param the_ref_dataset see \link[celaref]{make_ref_similarity_names_using_marked}
#' @param pval see \link[celaref]{make_ref_similarity_names_using_marked}
#' @param num_steps see \link[celaref]{make_ref_similarity_names_using_marked}
#' 
#' @return Table of similarity contrast results/assigned names e.t.c for a single group. 
#' Used internally for populating mwtest_res_table tables.
#' 
#' @seealso \code{\link[celaref]{make_ref_similarity_names_using_marked}} which calls this. 
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
#'     de_table.demo_query, de_table.demo_ref)
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
#' For use by make_ref_similarity_names_using_marked
#'
#' @param de_table.ref.marked The output of \code{\link{get_the_up_genes_for_all_possible_groups}} for the contrast of interest.
#' @param the_test_group Name of the test group in query dataset.
#' @param groupA One of the reference group names
#' @param groupB Another of the reference group names
#' @param enforceAgtB Do a one tailed test of A 'less' B (more similar)? Or two-tailed. Default = TRUE.
#'
#' @return A tibble of wilcox / man-whitneyU test results for this contrast.
#'
#' @seealso  \code{\link[celaref]{make_ref_similarity_names_using_marked}} 
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






#' get_stepped_pvals_str
#'
#' Internal function to construct the string of stepped pvalues reported by 
#' make_ref_similarity_names_using_marked
#' 
#' For use by make_ref_similarity_names_for_group
#'
#' @param mwtest_res_table.this Combined output of \code{\link{get_ranking_and_test_results}}
#'
#' @return Stepped pvalues string
#'
#' @seealso  \code{\link[celaref]{make_ref_similarity_names_for_group}} 
#'
#' @importFrom magrittr %>%
get_stepped_pvals_str  <- function(mwtest_res_table.this) {
   
   # dunno:0.000264,Weird subtype:0.226,Exciting:0.00705,Mystery cell type:NA
   mwtest_res_table.this.stepped <- mwtest_res_table.this %>% 
      dplyr::filter(.data$step == 1) %>% 
      dplyr::arrange(.data$grouprank) 
   
   stepped_pvals_str=paste0(mwtest_res_table.this.stepped$group, ":", base::signif(mwtest_res_table.this.stepped$pval, 3),
                            collapse=",")
   stepped_pvals_str=paste0(stepped_pvals_str, ",",
                            utils::tail(mwtest_res_table.this.stepped$next_group, 1),":NA")
   return(stepped_pvals_str)
}




#' get_matched_stepped_mwtest_res_table
#'
#' Internal function to grab a table of the matched group(s).
#' 
#' For use by make_ref_similarity_names_for_group
#'
#' @param mwtest_res_table.this Combined output of \code{\link{get_ranking_and_test_results}}
#' @param the_pval Pvalue threshold
#'
#' @return Stepped pvalues string
#'
#' @seealso \code{\link[celaref]{make_ref_similarity_names_for_group}} 
#'
#' @importFrom magrittr %>%
get_matched_stepped_mwtest_res_table <- function(mwtest_res_table.this, the_pval) {
   mwtest_res_table.this.stepped <- mwtest_res_table.this %>% 
      dplyr::filter(.data$step == 1) %>% 
      dplyr::arrange(.data$grouprank)
   
   mwtest_res_table.this.stepped.sig <- mwtest_res_table.this.stepped %>% dplyr::filter( .data$pval <= the_pval) 
   
   # Significance in adjacest steps defines matches
   # last in match is the last reference group that matches (usually 1)
   if (base::nrow(mwtest_res_table.this.stepped.sig) >  0) {
      last_of_match <- mwtest_res_table.this.stepped.sig$grouprank[1]
      return(mwtest_res_table.this.stepped[1:last_of_match,])
   } else {
      # Or NA for no match.
      return(NA)
   }
}



#' find_within_match_differences
#'
#' Internal function to find if there are significant difference between the
#' distribitions, when there are multiple match groups. 
#' 
#' For use by make_ref_similarity_names_for_group
#'
#' @param de_table.ref.marked see make_ref_similarity_names_for_group
#' @param matches see make_ref_similarity_names_for_group
#' @param the_test_group  see make_ref_similarity_names_for_group
#' @param the_test_dataset see make_ref_similarity_names_for_group
#' @param the_ref_dataset  see make_ref_similarity_names_for_group
#' @param the_pval  see make_ref_similarity_names_for_group
#' 
#' @return String of within match differences
#'
#' @seealso  \code{\link[celaref]{make_ref_similarity_names_for_group}} 
#'
#' @importFrom magrittr %>%
find_within_match_differences <- function(de_table.ref.marked, matches, the_test_group, the_test_dataset, the_ref_dataset, the_pval) {
   # Any significant (non-stepped) differences within match (from higher to lower)
   de_table.ref.marked.thismatch <- de_table.ref.marked %>%
      dplyr::filter(.data$test_group==the_test_group, .data$group %in% matches) 
   
   # Do an all vs all test, just within the matches.
   # Care if sig in either direction
   mwtest_res_table.inmatches <- get_ranking_and_test_results(
      the_test_group      = the_test_group,
      de_table.ref.marked = de_table.ref.marked.thismatch , 
      the_test_dataset    = the_test_dataset, 
      the_ref_dataset     = the_ref_dataset,
      num_steps=NA,
      pval=the_pval)
   
   diffs_in_match <- mwtest_res_table.inmatches %>% dplyr::filter(.data$pval <= the_pval)
   if (base::nrow(diffs_in_match) > 0 ) {
      return(paste0(diffs_in_match$group, " > ",diffs_in_match$next_group, 
                    " (p=",base::signif(diffs_in_match$pval,3),")", collapse="|"))
   } else { 
      return(NA)
   }
}




#' get_vs_random_pval
#'
#' Internal function to run a bionomial test of median test rank > 0.5 (random).
#' 
#' For use by make_ref_similarity_names_for_group
#'
#' @param de_table.ref.marked see make_ref_similarity_names_for_group
#' @param the_group Reference group name
#' @param the_test_group Test group name
#' #'
#' @return Pvalue result of a binomial test of each 'top gene' being greater than
#' the theoretical random median rank of 0.5 (halfway).
#'
#' @seealso \code{\link[celaref]{make_ref_similarity_names_for_group}} 
#'
#' @importFrom magrittr %>%
get_vs_random_pval <- function(de_table.ref.marked, the_group, the_test_group){
   # Also - last of match should be above (well, below) 0.5 - theoritical random
   # NB: power-wise there must be 6 or more genes (all true) for this to be sig.
   last_of_match_rank_dist <- de_table.ref.marked %>% 
      dplyr::filter( .data$group == the_group & .data$test_group == the_test_group) %>%  
      dplyr::pull(.data$rescaled_rank)
   not_random_pval_binom <- stats::binom.test(base::sum(last_of_match_rank_dist < 0.5), 
                                              n=base::length(last_of_match_rank_dist), 
                                              alternative = "greater")
   return(signif(not_random_pval_binom$p.value, 3))
   
}




#' get_reciprocal_matches 
#'
#' Internal function to run a bionomial test of median test rank > 0.5 (random).
#' 
#' For use by make_ref_similarity_names_using_marked
#'
#' @param mwtest_res_table.recip Combined output of \code{\link{get_ranking_and_test_results}}
#' for reciprocal test - ref vs query.
#' @param de_table.recip.marked Recriprocal ref vs query de_table.ref.marked
#' @param the_pval See make_ref_similarity_names_using_marked
#' 
#' @return List of table of reciprocal matches tested from reference to query.  
#'
#' @seealso \code{\link[celaref]{make_ref_similarity_names_using_marked}} 
#'
#' @importFrom magrittr %>%
get_reciprocal_matches <- function(mwtest_res_table.recip, de_table.recip.marked, the_pval) {
   
   the_test_groups <- unique(mwtest_res_table.recip$test_group) # test here being ref
   
   
   get_reciprocal_matches_group <- function(the_test_group, mwtest_res_table.recip, de_table.recip.marked,  the_pval) {
      # Just this group
      mwtest_res_table.recip.this <- mwtest_res_table.recip %>% dplyr::filter(.data$test_group == the_test_group)
      
      # Get a bunch of matches (stepped), or NA
      mwtest_res_table.matches <- get_matched_stepped_mwtest_res_table(mwtest_res_table.recip.this, the_pval)
      
      if (! all(is.na(mwtest_res_table.matches))) {
         # Pval to a random distribution of ranks.
         not_random_pval <- get_vs_random_pval(de_table.recip.marked, the_group=mwtest_res_table.matches$group[nrow(mwtest_res_table.matches)], the_test_group)
         
         if (not_random_pval <= the_pval) {
            # Ref + test relative to original query this is the reciprocal of
            return( dplyr::bind_cols( test_group=mwtest_res_table.matches$group,
                                      ref_group =mwtest_res_table.matches$test_group ))
         }
      }
      
      return(NULL)
   }  
   
   
   reciprocal_matches <- dplyr::bind_rows(lapply(FUN=get_reciprocal_matches_group, X=the_test_groups, 
                                          mwtest_res_table.recip=mwtest_res_table.recip, de_table.recip.marked=de_table.recip.marked, the_pval=the_pval))
   
   if(is.null(reciprocal_matches)) {reciprocal_matches <- NA}
   return(reciprocal_matches)
}



