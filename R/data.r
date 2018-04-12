#' Demo query de table
#' 
#' Small example dataset that is the output of
#' \link[celaref]{contrast_each_group_to_the_rest}. It contains the results
#' of each group compared to the rest of the sample (ie within sample 
#' differential expression)
#'
"de_table.demo_query"


#' Demo ref de table
#' 
#' Small example dataset that is the output of
#' \link[celaref]{contrast_each_group_to_the_rest}. It contains the results
#' of each group compared to the rest of the sample (ie within sample 
#' differential expression)
#'
"de_table.demo_ref" 


#' Demo cell info table
#' 
#' Sample sheet table listing each cell, its assignd cluster/group, and 
#' any other information that might be interesting (replicate, individual e.t.c)
#'
"demo_cell_info_table"



#' Demo count matrix
#' 
#' Counts matrix for a small, demo example datasets. Raw counts of 
#' reads per gene (row) per cell (column).
"demo_counts_matrix" 


#' Demo gene info table
#'
#' Extra table of gene-level information for the demo example dataset.
#' Can contain anything as long as theres a unique gene id.
"demo_gene_info_table" 

#' Demo microarray expression table
#' 
#' Microarray-style expression table for the demo example dataset. 
#' Rows are genes, columns are samples, as per counts matrix.
"demo_microarray_expr"


#' Demo microarray sample sheet table
#' 
#' Microarray sample sheet table for the demo example dataset. 
#' Contains array identifiers, their group and any other information that could
#' be useful.
"demo_microarray_sample_sheet" 


#' Demo query se (summarizedExperiment)
#' 
#' A summarisedExperiment object loaded from demo info tables, for a query set. 
"demo_query_se" 

#' Demo reference se (summarizedExperiment)
#' 
#' A summarisedExperiment object loaded from demo info tables, for a reference
#' set. 
"demo_ref_se"




