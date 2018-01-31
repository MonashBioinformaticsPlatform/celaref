





#' contrast_each_group_to_the_rest
#'
#' \code{contrast_each_group_to_the_rest} produces a table of differential
#' expression results, where each group (cluster) is compared to the rest of the
#' cells.
#'
#' Note that this function is \emph{slow}, because it runs the differential
#' expression. It only needs to be run once per dataset (unless group labels
#' change)
#'
#' @param dataset_se Summarised experiment object containing count data. Also
#' requires 'ID' and 'group' to be set within the cell information
#' (see \code{colData()})
#' @param dataset_name Short, meaningful name for this dataset/experiment.
#' @param groups2test An optional character vector specificing specific groups to check.
#' By default (set to NA), all groups will be tested.
#'
#' @return A tibble, the de_table (differential expression table).  This table
#' is a core interal summary of each individual dataset, which is
#' used for the cross-dataset comparisons. So save this output to pass to other
#' functions.
#'
#'
#' @examples
#' de_table.demo_query  <- contrast_each_group_to_the_rest(demo_query_se, "a_demo_query")
#'
#' @export
contrast_each_group_to_the_rest <- function(dataset_se, dataset_name, groups2test=NA) {

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

    ## For each group, test it versus evyerthing else (paralallised)
    de_table_list <- mclapply(groups2test, FUN=contrast_the_group_to_the_rest, dataset_se=dataset_se, mc.cores=NUM_CORES)
    de_table.allvsrest <- bind_rows(de_table_list)

    # Factorise group,and add dataset name
    de_table.allvsrest$group   <- factor(de_table.allvsrest$group)
    de_table.allvsrest$dataset <- dataset_name

    return(de_table.allvsrest)
}

