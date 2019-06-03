context("Loading functions")
library(celaref)


test_that("Load se from files, tables, 10X", {
  

   expect_something_in_demo_se <- function(test_se) {
      
      # any 0-length (or 1 length) dimensions are a fail, 
      # and is usual fail case.
      # but don't check what's actually there, because it could change
      # These are different sized datasets anyway.
      
      expect_gt(base::ncol(test_se), 1) # cells
      expect_gt(nrow(test_se),  1) # genes # 1 gene would be wrong too.
      expect_gt(sum(assays(test_se)[[1]]), 1) #total counts aren't all 0
      
      expect_gt(nrow(colData(test_se)) , 1 ) 
      expect_gt(ncol(colData(test_se)) , 1 ) 
      
      expect_gt(nrow(rowData(test_se)) , 1 ) 
      expect_gt(ncol(rowData(test_se)) , 1 ) 
   }
   
   
   counts_filepath    <- system.file("extdata", "sim_query_counts.tab",    package = "celaref")
   cell_info_filepath <- system.file("extdata", "sim_query_cell_info.tab", package = "celaref")
   gene_info_filepath <- system.file("extdata", "sim_query_gene_info.tab", package = "celaref")
   
   demo_se.files <- load_se_from_files(counts_filepath, 
                                 cell_info_file = cell_info_filepath, 
                                 gene_info_file = gene_info_filepath)
   expect_something_in_demo_se(demo_se.files)
   

   
   demo_se.tables <- load_se_from_tables(counts_matrix=demo_counts_matrix, 
                                  cell_info_table=demo_cell_info_table, 
                                  gene_info_table=demo_gene_info_table)
   expect_something_in_demo_se(demo_se.tables)

   
   example_10X_dir <- system.file("extdata", "sim_cr_dataset", package = "celaref")
   dataset_se.10X <- load_dataset_10Xdata(example_10X_dir, dataset_genome="GRCh38", 
        clustering_set="kmeans_4_clusters", gene_id_cols_10X=c("gene")) 
   
   expect_something_in_demo_se(dataset_se.10X)
   
})



test_that("Filtering low expression genes and groups", {
   
   demo_ref_se.trim <- trim_small_groups_and_low_expression_genes(
      dataset_se=demo_ref_se, 
      min_lib_size=1000, min_group_membership=50,
      min_reads_in_sample=1, min_detected_by_min_samples=20 
   )

   expect_equal(length(levels(colData(demo_ref_se.trim)$group)), 3)
   expect_equal(nrow(demo_ref_se.trim), 199)
   expect_equal(ncol(demo_ref_se.trim), 489)
   
})




test_that("Converting gene ids",{
   
   dataset_se <- demo_ref_se[1:10, 1:10]
   rowData(dataset_se)$dummyname <- c(rep("A",5), rep("B",5))
   rowData(dataset_se)$total_not_count <- 1:10
   
   dataset_se.2 <- convert_se_gene_ids(dataset_se, new_id='dummyname', eval_col='total_not_count')
   
   expect_equal(rowData(dataset_se.2)["A","total_not_count"], 5)
      
})