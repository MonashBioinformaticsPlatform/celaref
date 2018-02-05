





#' contrast_each_group_to_the_rest
#'
#' Produces a table of within-experiment differential expression results (for
#' either query or reference experiment), where each group (cluster) is
#' compared to the rest of the cells.
#'
#' Note that this function is \emph{slow}, because it runs the differential
#' expression. It only needs to be run once per dataset though (unless group labels
#' change). Having package \pkg{parallel} installed is highly recomended.
#'
#' Both reference and query datasets should be processed with this
#' function.
#'
#' The tables produced by this function (usually named something like
#' \emph{de_table.datasetname}) contain summarised results of MAST results.
#' Each group is compared versus cells in the group, versus not in the group,
#' (Ie. always a 2-group contrast, other groups information is ignored). As per
#' MAST reccomendataions, the proportion of genes seen in each cell is included
#' in the model.
#'
#' @param dataset_se Summarised experiment object containing count data. Also
#' requires 'ID' and 'group' to be set within the cell information
#' (see \code{colData()})
#' @param dataset_name Short, meaningful name for this dataset/experiment.
#' @param groups2test An optional character vector specificing specific groups to check.
#' By default (set to NA), all groups will be tested.
#' @param num_cores Number of cores to use to run MAST jobs in parallel.
#' Ignored if parallel package not available. Set to 1 to avoid parallelisation.
#' Default = 4
#'
#'
#' @return A tibble the within-experiment de_table (differential expression table).
#' This is a core summary of the individual experiment/dataset, which is
#' used for the cross-dataset comparisons.
#'
#' The table feilds won't neccesarily match across datasets, as they include
#' cell annotations information. Important columns (used in downstream analysis) are:
#'
#' \describe{
#'\item{ID}{Gene identifier}
#'\item{ci_inner}{ Inner (conservative) 95\% confidence interval of log2 fold change. }
#'\item{fdr}{Multiple hypothesis corrected p-value (using BH/FDR method)}
#'\item{group}{Cells from this group were compared to everything else}
#'\item{sig_up}{Significnatly differentially expressed (fdr < 0.01), with a postivie fold change?}
#'\item{rank}{Rank position (within group), ranked by CI inner, highest to lowest. }
#'\item{rescaled_rank}{Rank scaled 0(top most overrepresented genes in group)-1(top most not-present genes)}
#'\item{dataset}{Name of dataset/experiment}
#'}
#'
#'
#' @examples
#' de_table.demo_query  <- contrast_each_group_to_the_rest(demo_query_se, "a_demo_query")
#' de_table.demo_ref    <- contrast_each_group_to_the_rest(demo_ref_se, "a_demo_ref", num_cores=2)
#'
#' @export
contrast_each_group_to_the_rest <- function(dataset_se, dataset_name, groups2test=NA, num_cores=4) {
   
   # Which groups to look at? Default all in query dataset.
   if (length(groups2test) == 1 && is.na(groups2test)) {
      groups2test = levels(dataset_se$group)
   } else { # Check its sensible before long processing steps
      if (! all(groups2test %in% levels(dataset_se$group))) {
         print(paste("can't find all of ",groups2test))
         stop("Can't find all test groups in dataset")
         return(NA)
      }
   }
   
   ## Add the 'proportion of genes covered in this factor' variable.
   # Not really a proprotion, buts in the model (and is proportional)
   # see MAST doco/paper - its used in the model for reasons.
   colData(dataset_se)$pofgenes <- scale(Matrix::colSums(as.matrix(assay(dataset_se)) > 0 ) )
   
   ## For each group, test it versus evyerthing else (paralallised, with lapply fallback )
   de_table_list <- NA
   if (requireNamespace("parallel", quietly = TRUE) & num_cores >= 1) {
      de_table_list <- parallel::mclapply(groups2test, FUN=contrast_the_group_to_the_rest, dataset_se=dataset_se, mc.cores=num_cores)
   } else {
      de_table_list <- base::lapply(groups2test, FUN=contrast_the_group_to_the_rest, dataset_se=dataset_se) 
   }
   de_table.allvsrest <- dplyr::bind_rows(de_table_list)

   # Factorise group,and add dataset name
   de_table.allvsrest$group   <- factor(de_table.allvsrest$group)
   de_table.allvsrest$dataset <- dataset_name
   
   return(de_table.allvsrest)
}











#' contrast_the_group_to_the_rest
#'
#' Internal function to calculate differential expression within an experiment
#' between a specified group and cells not in that group.
#'
#'
#' This function should only be called by \code{contrast_each_group_to_the_rest}
#' (which can be passed a single group name if desired). Else 'pofgenes' will not
#' be defined.
#'
#' MAST is supplied with log2(counts + 1.1), and zlm called with model
#' '~ TvsR + pofgenes' . The p-values reported are from the hurdle model. FDR is
#' with default fdr/BH method.
#'
#'
#' @param dataset_se Datast summarisedExperiment object.
#' @param the_group group to test
#' @param pvalue_threshold Default = 0.01
#'
#' @return A tibble, the within-experiment de_table (differential expression
#' table), for the group specified.
#'
#' @seealso \code{\link[celaref]{contrast_each_group_to_the_rest}}
#'
contrast_the_group_to_the_rest <- function( dataset_se, the_group, pvalue_threshold=0.01) {
   
   TEST = 'test'
   REST = 'rest'
   # must have pofgenes set.
   stopifnot("pofgenes" %in% colnames(colData(dataset_se)))
   
   # Define test vs rest factor for model (TvsR only applicable for a single contrast.)
   col_data_for_MAST      <- colData(dataset_se)
   col_data_for_MAST$TvsR <- factor(ifelse(col_data_for_MAST$group == the_group, TEST, REST), levels=c(REST, TEST)) #Note order for MAST
   
   
   ## Log2 transform the counts
   
   # *** NB: Adding 1.1 to each count before log! Not 1! ***
   # When there is *no* expression of a given gene within a group - coeficcients aren't generated.= (NA)
   # But they can still be sinigicant - and for scRNA these could be real+relevant! Particularly with small groups.
   # .. I think this is just a quirk of how mast is calculating them.
   # If I add 1.1 instead of 1 pre-log transformation, means log2(0+1.1) min is 0.014,
   # not 0 - and coef (log2FC) is generated
   # Adding 1.1 isn't going to have any real difference to adding 1.
   # ***
   TO_ADD <- 1.1 # Defining it here so it no get changed.
   logged_counts  <- log2( assay(dataset_se) + TO_ADD)
   
   
   # MAST uses a SE internally (inherits) but needs the parts to make it itself.
   ## Make sca object for MAST
   sca <- MAST::FromMatrix(logged_counts, cData=col_data_for_MAST, fData=rowData(dataset_se))
   
   
   ## Make Zlm model and run contrast:
   # Slow step.
   # MAST uses pofgene in the model, see doco/paper.
   # Would have liked to model group level and then to contrasts of groupA - (avg of all other groups).
   # Would allow this step to run only once per experiment, not per each group.
   # But that is complicated, and the 'summary' method (what actually does calcs) only supports
   # simple single coeffient test.
   zlm.TvsR        <- MAST::zlm(~ TvsR + pofgenes, sca)
   the_contrast    <- 'TvsRtest'
   summary.TvsR    <- MAST::summary(zlm.TvsR, doLRT = 'TvsRtest')
   summary.TvsR.df <- as.data.frame(summary.TvsR$datatable)
   
   # Grab the p-values form the 'H' hurdle model. And the logFC from the calcs
   # all for testvsrest term of formula (vs intercept) (equates to test-rest )
   mast_res <- merge(summary.TvsR.df[summary.TvsR.df$contrast==the_contrast  & summary.TvsR.df$component=='H',
                                     c("primerid", "Pr(>Chisq)")], #hurdle P-values
                     summary.TvsR.df[summary.TvsR.df$contrast==the_contrast  & summary.TvsR.df$component=='logFC',
                                     c("primerid", "coef", "ci.hi", "ci.lo")], by='primerid') #logFC coefficients
   #> head(mast_res)
   #       primerid  Pr(>Chisq)        coef       ci.hi       ci.lo
   #1 0610007P14Rik 0.785451166  0.08833968  0.33665327 -0.15997391
   #2 0610009B22Rik 0.411386651  0.15749376  0.38770946 -0.07272195
   #3 0610009O20Rik 0.869118657 -0.03563795  0.09585160 -0.16712750
   
   
   ## Reformat the results with more infomative colnames. primerid=>gene, Pr(>Chisq)=>pval
   de_table <- data.frame(mast_res)
   colnames(de_table) <- c("ID", "pval", "log2FC","ci.hi", "ci.lo")
   
   # Join on any rowData (usually other gene names if present, also pofgenes)
   de_table <- merge(x=data.frame(rowData(dataset_se)), y=de_table,  by="ID")
   
   
   # Order by Innermost/conservative CI
   # upper/lower mean numerically, actually want to use inner/outer rel to 0
   de_table$ci_inner  <- mapply(FUN=get_inner_or_outer_ci, MoreArgs = list(get_inner=TRUE),  de_table$log2FC, de_table$ci.hi, de_table$ci.lo)
   de_table$ci_outer  <- mapply(FUN=get_inner_or_outer_ci, MoreArgs = list(get_inner=FALSE), de_table$log2FC, de_table$ci.hi, de_table$ci.lo)
   de_table <- de_table[,! colnames(de_table) %in% c("ci.hi", "ci.lo")]
   de_table <- de_table[order(de_table$ci_inner, decreasing = TRUE),]
   
   
   de_table$fdr           <- p.adjust(de_table$pval, 'fdr')
   de_table$group         <- the_group
   de_table$sig           <- de_table$fdr <= pvalue_threshold
   de_table$sig_up        <- de_table$sig & de_table$log2FC > 0
   de_table$gene_count    <- nrow(de_table)
   de_table$rank          <- 1:nrow(de_table)
   de_table$rescaled_rank <- 1:nrow(de_table) / nrow(de_table)
   
   
   #       gene         pval   log2FC ci_inner ci_outer          fdr      group  sig sig_up gene_count rank rescaled_rank
   #1914  Epcam 6.926312e-37 2.210175 1.981575 2.438775 4.678723e-33 Epithelial TRUE   TRUE       6755    1  0.0001480385
   #2228  Fxyd3 4.724345e-33 2.069778 1.833179 2.306377 1.063765e-29 Epithelial TRUE   TRUE       6755    2  0.0002960770
   #3386    Mif 6.042683e-22 2.148053 1.805283 2.490823 5.831189e-19 Epithelial TRUE   TRUE       6755    3  0.0004441155
   #1902   Eno1 6.245552e-20 2.029549 1.681717 2.377380 4.687634e-17 Epithelial TRUE   TRUE       6755    4  0.0005921540
   #4256    Pkm 4.906642e-17 1.953987 1.579322 2.328652 2.117840e-14 Epithelial TRUE   TRUE       6755    5  0.0007401925
   
   return(de_table)
}





#' get_inner_or_outer_ci
#'
#' Given a fold-change, and high and low confidence interval (where lower <
#' higher), pick the innermost/most conservative one.
#'
#'
#' @param fc Fold-change
#' @param ci.hi Higher fold-change CI (numerically)
#' @param ci.lo smaller fold-change CI (numerically)
#' @param get_inner If TRUE, get the more conservative inner CI, else the bigger
#' outside one.
#'
#' @return inner or outer CI from \bold{ci.hi} or \bold{ci.low}
#'
get_inner_or_outer_ci<- function(fc, ci.hi, ci.lo, get_inner=TRUE)      {
   
   if (is.na(fc)) {
      return(0)
   }
   # if get_inner == false then get outer
   if (fc > 0) {
      # inner is lower in up, outer is high
      if (get_inner) { return(ci.lo) } else { return(ci.hi) }
   }
   # Else, this is a downwards FC, inner and outer CI are switched.
   else {
      if (get_inner) { return(ci.hi) } else { return(ci.lo) }
   }
}








#' get_the_up_genes_for_group 
#'
#' For the most overrepresented genes of the specified group in the test 
#' dataset, get their rankings in all the groups of the reference dataset. 
#' 
#' This is effectively a subset of the reference data, 'marked' with the 'top'
#' genes that represent the group of interest in the query data. The 
#' distribution of the \emph{rescaled ranks} of these marked genes in each 
#' reference data group indicate how similar they are to the query group. 
#' 
#'
#' @param the_group The group (from the test/query experiment) to examine. 
#' @param de_table.test A differential expression table of the query experiment,
#'  as generated from \code{\link[celaref]{contrast_each_group_to_the_rest}}
#' @param de_table.ref A differential expression table of the reference dataset,
#'  as generated from \code{\link[celaref]{contrast_each_group_to_the_rest}}
#' @param rankmetric Placeholder for support of different ranking methods, but only the default supported. Omit. 
#' 
#' @return \emph{de_table.marked} This will be a subset of \bold{de_table.ref}, 
#' with an added column \emph{test_group} set to \bold{the_group}. If nothing passes the rankmetric criteria, NA.
#'
#' @examples
#' de_table.marked.Group3vsRef <- get_the_up_genes_for_group(the_group="Group3",
#' de_table.test=de_table.demo_query, de_table.ref=de_table.demo_ref)
#'
#' @seealso  
#' \code{\link[celaref]{contrast_each_group_to_the_rest}} For prepraring the de_table.* tables.
#' \code{\link[celaref]{get_the_up_genes_for_all_possible_groups}} For running all query groups at once.
#
#'
#'@export
get_the_up_genes_for_group <- function(the_group, de_table.test, de_table.ref, rankmetric='TOP100_LOWER_CI_GTE1') {
   
   # Trialled a bunch. Decided on TOP100_LOWER_CI_GTE, but leave option here. (might be useful for other data types)
   the_up_genes <- c()
   if (rankmetric=='TOP100_LOWER_CI_GTE1') {
      the_up_genes <- de_table.test$ID[de_table.test$group == the_group & de_table.test$rank <= 100 & de_table.test$ci_inner >= 1]
   } else {stop("Unknown rank metric")}
   
   # Nothing selected. Or nothing selected thats in ref data. NA gets rbinded harmlesslessy
   if (length(the_up_genes) == 0 ) {
      warning(paste0("Found no mark-able genes meeting criteria when looking for 'up' genes in ",the_group,
                     " Cannot test it for similarity. Could occur if: cell groups are very similar,",
                     "a lack of statistical power (e.g. small number of cells in group),",
                     "or for a heterogenous cluster. It may only affect one/some groups, continuing with the rest.") )
      return(NA)
      }
   if (sum(de_table.ref$ID %in% the_up_genes) == 0) {return(NA)}
   
   de_table.marked <- de_table.ref[de_table.ref$ID %in% the_up_genes,]
   de_table.marked$test_group <- the_group
   
   return(de_table.marked)
   
}




#' get_the_up_genes_for_group 
#'
#' For the most overrepresented genes of each group in the test 
#' dataset, get their rankings in all the groups of the reference dataset. 
#' 
#' This is effectively a subset of the reference data, 'marked' with the 'top'
#' genes that represent the groups in the query data. The 
#' distribution of the \emph{rescaled ranks} of these marked genes in each 
#' reference data group indicate how similar they are to the query group. 
#' 
#' This function is simply a conveinent wrapper for 
#' \code{\link[celaref]{get_the_up_genes_for_group}} that merges output for 
#' each group in the query into one table.
#' 
#'
#' @param de_table.test A differential expression table of the query experiment,
#'  as generated from \code{\link[celaref]{contrast_each_group_to_the_rest}}
#' @param de_table.ref A differential expression table of the reference dataset,
#'  as generated from \code{\link[celaref]{contrast_each_group_to_the_rest}}
#' @param test_dataset_name A short meaningful name for the query experiment.   
#' @param rankmetric Placeholder for support of different ranking methods, but only the default supported. Omit. 
#' 
#' @return \emph{de_table.marked} This will alsmost be a subset of \bold{de_table.ref}, 
#' with an added column \emph{test_group} set to the query groups, and 
#' \emph{test_dataset} set to \bold{test_dataset_name}.
#' 
#' If nothing passes the rankmetric criteria, a warning is thrown and NA is 
#' returned. (This can be a genuine inability to pick out the 
#' representative 'up' genes, or due to some problem in the analysis)
#'
#' @examples
#' de_table.marked.query_vs_ref <- get_the_up_genes_for_all_possible_groups(
#' de_table.test=de_table.demo_query , 
#' de_table.ref=de_table.demo_ref ) 
#'
#' @seealso  \code{\link[celaref]{get_the_up_genes_for_group}} Function for testing a single group.
#'
#'
#'@export
get_the_up_genes_for_all_possible_groups <- function(de_table.test, de_table.ref, test_dataset_name, rankmetric='TOP100_LOWER_CI_GTE1'){
   
   # Run get_the_up_genes foreach group, 
   # But remove NAs before contsructing table
   #   e.g. no sig for myoepithelail vs rest. :. it will produce a NA result.
   de_table.marked.list <- lapply(FUN=get_the_up_genes_for_group, 
                                  X = levels(de_table.test$group), rankmetric=rankmetric,
                                  de_table.test=de_table.test, de_table.ref=de_table.ref)
   

   if (all(is.na(de_table.marked.list))) { warning("Found no mark-able genes meeting criteria when looking for 'up' genes in each group. Cannot test for similarity. Could be an error, or occur if cell groups are very similar or due to a lack of statistical power.") }
   
   # remove the NAs.
   de_table.marked.list <- de_table.marked.list[! is.na(de_table.marked.list)]
   de_table.marked      <- as.data.frame(dplyr::bind_rows(de_table.marked.list), stringsAsFactors=FALSE)
   de_table.marked$test_dataset <-  test_dataset_name
   
   
   return(de_table.marked)
}




