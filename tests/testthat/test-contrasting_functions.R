context("test-contrasting_functions")


test_that("MAST contrasts - dense, sparse and empty", {
   
   # Checkt the first 5 genes are the same. 
   asubset   <- seq_len(30)
   top_check <- seq_len(5)
   de_genes   <- c("Gene3",  "Gene23", "Gene10", "Gene25", "Gene30") 
   de_genes.3 <- c("Gene1",  "Gene4",  "Gene3",  "Gene30", "Gene10") #altered data.
   
   # Densse
   demo_query_se.1 <- demo_query_se[asubset,asubset]
   de_table1.demo_query  <- contrast_each_group_to_the_rest(
      demo_query_se.1, "a_demo_query", num_cores=1)
   
   expect_equal(de_table1.demo_query$ID[top_check], de_genes )
   
   
   #now sparse
   demo_query_se.2 <- demo_query_se.1
   assays(demo_query_se.2)[[1]] <-  Matrix::Matrix(assay(demo_query_se.1), sparse=TRUE)
   de_table2.demo_query  <- contrast_each_group_to_the_rest(
      demo_query_se.2, "a_demo_query", num_cores=1)
   expect_equal(de_table1.demo_query$ID[top_check], de_genes )
   
   
   # what is one gene's expression for entire group totally empty?
   # (previously caused errors and there's a workaround now.)
   demo_query_se.3 <- demo_query_se.1
   assays(demo_query_se.3)[[1]][1,demo_query_se.1$group =="Group2"] <- 0
   de_table3.demo_query  <- contrast_each_group_to_the_rest(
      demo_query_se.3, "a_demo_query", num_cores=1)
   expect_equal(de_table3.demo_query$ID[top_check], de_genes.3 )
   #ID       pval      log2FC ci_inner  ci_outer       fdr  group   sig sig_up gene_count rank rescaled_rank      dataset
   #1  Gene1 0.06701189 -1.98643911 1.734379 -5.707257 0.4505910 Group1 FALSE  FALSE         30    1    0.03333333 a_demo_query
   #2  Gene4 0.99903639 -0.04353777 1.726446 -1.813522 0.9990364 Group1 FALSE  FALSE         30    2    0.06666667 a_demo_query
   #3  Gene3 0.50595300 -2.71865608 1.394369 -6.831681 0.5837919 Group1 FALSE  FALSE         30    3    0.10000000 a_demo_query
   #4 Gene30 0.67979358 -0.94039185 1.281941 -3.162725 0.7279509 Group1 FALSE  FALSE         30    4    0.13333333 a_demo_query
   #5 Gene10 0.68075516 -0.46524523 1.147719 -2.078209 0.7279509 Group1 FALSE  FALSE         30    5    0.16666667 a_demo_query
})




test_that("MAST contrasts - hdf5-backed assays and SCE objects", {
   # Just test that runs - 
   # these things are succeptible to format / object changes.

   # dense sce 
   d.sce.den <- as(demo_query_se, "SingleCellExperiment")
   expect_equal( 10, nrow(
      contrast_each_group_to_the_rest(d.sce.den[1:10,],'test', 
      groups2test = "Group2", n.group = 20, num_cores = 1)))
   
   # sparse SCE   
   d.sce.sp        <- d.sce.den 
   assays(d.sce.sp)[[1]] <- Matrix::Matrix(assays(d.sce.sp)[[1]], sparse=TRUE)
   expect_equal( 10, nrow(
      contrast_each_group_to_the_rest(d.sce.sp[1:10,],'test',     
      groups2test = "Group2", n.group = 20, num_cores = 1)))


   # hdf5 SCE
   d.sce.hdf       <- HDF5Array::saveHDF5SummarizedExperiment( d.sce.sp , replace=TRUE)
   expect_equal( 10, nrow(
      contrast_each_group_to_the_rest(d.sce.hdf[1:10,],'test',     
      groups2test = "Group2", n.group = 20, num_cores = 1)))

})





test_that("Microarray reference", {
   
   top5 <- c("Gene100", "Gene150", "Gene57",  "Gene80",  "Gene21" )
   de_table.ma <- contrast_each_group_to_the_rest_for_norm_ma_with_limma(
      norm_expression_table=demo_microarray_expr, 
      sample_sheet_table=demo_microarray_sample_sheet,
      dataset_name="DemoSimMicroarrayRef", 
      sample_name="cell_sample", group_name="group") 
   
   
   expect_equal(de_table.ma$ID[seq_len(5)], top5)
   
})




test_that("Rankmetrics", {
   
   # Ask for just 10 genes and check them. Actually same for both mehods.
   genes.TOP100_LOWER_CI_GTE1 <- 
      c("Gene100", "Gene150", "Gene57",  "Gene80",  "Gene21",  
        "Gene30",  "Gene23",  "Gene65",  "Gene101", "Gene10")
   genes.TOP100_SIG <- genes.TOP100_LOWER_CI_GTE1 # are same

   
   de_table.marked.Group3vsRef.TOP100_LOWER_CI_GTE1 <- 
      get_the_up_genes_for_group(
                     the_group="Group3",
                     de_table.test=de_table.demo_query, 
                     de_table.ref=de_table.demo_ref,
                     rankmetric = "TOP100_LOWER_CI_GTE1",
                     n=10)
   expect_equal(de_table.marked.Group3vsRef.TOP100_LOWER_CI_GTE1$ID[
      de_table.marked.Group3vsRef.TOP100_LOWER_CI_GTE1$group == "Dunno"],
                genes.TOP100_LOWER_CI_GTE1)
   
   
   
   de_table.marked.Group3vsRef.TOP100_SIG <- 
      get_the_up_genes_for_group(
         the_group="Group3",
         de_table.test=de_table.demo_query, 
         de_table.ref=de_table.demo_ref,
         rankmetric = 'TOP100_SIG', n=10)

   expect_equal(de_table.marked.Group3vsRef.TOP100_SIG$ID[
      de_table.marked.Group3vsRef.TOP100_SIG$group == "Dunno"],
      genes.TOP100_SIG)

   
})




test_that("Subsetting ses", {
   
   dataset_se.30pergroup <- subset_cells_by_group(demo_query_se, n.group=30)
   expect_equal(sum(dataset_se.30pergroup$group == "Group3"),30) 
   expect_equal(sum(dataset_se.30pergroup$group == "Group1"),28) 
   
   demo_query_se.subset2 <- subset_se_cells_for_group_test(demo_query_se, 
                             the_group="Group3", 
                             n.group=20, 
                             n.other=30)
   expect_equal(sum(demo_query_se.subset2$group == "Group3"),20) 
   expect_equal(sum(demo_query_se.subset2$group != "Group3"),30) 
   
})



#test_that("Finding counts", {
#   
#})

