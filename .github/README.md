# celaref


__***DEV VERSION : Method and doco still under active development***__


### Function  

Single cell RNA sequencing (scRNAseq) has made it possible to examine the 
cellular heterogeny within a tissue or sample, and observe changes and 
characteristics in specific cell types. To do this, we need to group the cells
into clusters and figure out what they are.

The celaref (*ce*ll *la*belling by *ref*erence) package aims to streamline the cell-type identification step, by 
suggesting cluster labels on the basis of similarity to an already-characterised
reference dataset - wheather that's from a similar experiment performed 
previously in the same lab, or from a public dataset from a similar sample. 

### Input

To look for cluster similarities celaref needs:

* The query dataset :
    - a table of read counts per cell per gene
    - a list of which cells belong in which cluster
   
* A reference dataset:
    - a table of read counts per cell per gene
    - a list of which cells belong in which *annotated* cluster
   
### Output



![](../vignettes/images/violin_plot_example.png) 


Query Group | Short Label                        | pval    |
------------|------------------------------------|---------|
cluster_1   |cluster_1:astrocytes_ependymal      |2.98e-23 |
cluster_2   |cluster_2:endothelial-mural         |8.44e-10 |
cluster_3   |cluster_3:no_similarity             |NA       |
cluster_4   |cluster_4:microglia                 |2.71e-19 |
cluster_5   |cluster_5:pyramidal SS\|interneurons|3.49e-10 |
cluster_6   |cluster_6:oligodendrocytes          |2.15e-28 |




This is a comparison of brain scRNAseq data from :

 * Zeisel, A., Manchado, A. B. M., Codeluppi, S., Lonnerberg, P., La Manno, G., Jureus, A., … Linnarsson, S. (2015). *Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq.* Science, 347(6226), 1138–42. http://doi.org/10.1126/science.aaa1934
 * Darmanis, S., Sloan, S. A., Zhang, Y., Enge, M., Caneda, C., Shuer, L. M., … Quake, S. R. (2015). *A survey of human brain transcriptome diversity at the single cell level.* Proceedings of the National Academy of Sciences, 112(23), 201507125. http://doi.org/10.1073/pnas.1507125112


### More information?

Full details in the vignette [html](http://bioinformatics.erc.monash.edu/home/sarah.williams/projects/cell_groupings/doco/celaref_doco.html) - method description, manual and example analyses.

