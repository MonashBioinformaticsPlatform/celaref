context("Contrasting functions")
library(celaref)

test_that("MAST contrasts - dense and sparse", {
   
   asubset   <- seq_len(30)
   top_check <- seq_len(5)
   de_genes  <- c("Gene3",  "Gene23", "Gene10", "Gene25", "Gene30") 
   
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
   
})

