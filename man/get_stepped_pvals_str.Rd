% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group_labelling_functions.r
\name{get_stepped_pvals_str}
\alias{get_stepped_pvals_str}
\title{get_stepped_pvals_str}
\usage{
get_stepped_pvals_str(mwtest_res_table.this)
}
\arguments{
\item{mwtest_res_table.this}{Combined output of 
\code{\link{get_ranking_and_test_results}}}
}
\value{
Stepped pvalues string
}
\description{
Internal function to construct the string of stepped pvalues reported by 
make_ref_similarity_names_using_marked
}
\details{
For use by make_ref_similarity_names_for_group
}
\seealso{
\code{\link{make_ref_similarity_names_for_group}}
}
