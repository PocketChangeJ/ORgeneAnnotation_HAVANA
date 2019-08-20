---
title: "Single-cell RNA-seq of OSNs"
date: '20 August, 2019'
output:
  html_document:
    keep_md: true
    fig_width: 5
    fig_height: 5
    code_folding: hide
    toc: true
    toc_float: 
      collapsed: false
---



Single-cell RNA-seq data from 34 single OSNs picked by hand from a GFP-OMP mouse, postnatal day 3. Cells that looked healthy and with high GFP expression were selected and processed further to produce libraries for single-cell RNA-seq.

Sequencing reads were mapped to the mouse genome `mm10` using `STAR version 2.6.0c` with gene quantification enabled, using `Ensembl version 93` annotation.

The counts from each sample were compiled into a combined count matrix with the script `readGeneCounts.R`. These can be downloaded from: ``

------

Here, we analyse the data with the purpose of finding OSNs that express split OR genes, that is, OR genes that encode a predicted functional protein across two exons.

First, we load the count matrix and mapping statistics (parsed from the STAR `Log.final.out` files). We also separate the statistics reported by `STAR` from the gene counts. We have data for 34 cells, across the whole transcriptome.


```r
ann <- read.table(paste0(dir, "geneAnnotation.Ensembl_v93.tsv"), stringsAsFactors = FALSE)

data <- read.table(paste0(dir, "geneCounts_P3OMPcells.RAW.tsv"), stringsAsFactors = FALSE)
mapping.stats <- read.table(paste0(dir, "mappingStats.tsv"), stringsAsFactors = FALSE, header = TRUE, row.names = 1)

counting.stats <- as.data.frame(t(data[1:4,-c(1:2)]))
data <- data[-c(1:4),]
dim(data[,-c(1:2)])
```

```
## [1] 54232    34
```

### Quality-control

Before analysing the data itself we need to check the quality of each library and remove any that are not good enough. Since the cells were hand picked and selected to be healthy looking cells we expect that most libraries will be of good quality.

To assess quality we check the general mapping statistics, the proportion of reads in mitochondrial reads and total number of detected genes per cell.

Samples were sequenced to a median depth of ~3 million fragments. The distribution of library sizes is quite uniform, with no obvious outliers.


```r
counting.stats <- cbind(counting.stats, colSums(data[,-c(1:2)]))
colnames(counting.stats)[5] <- "N_inExons"

ggplot(mapping.stats, aes(1, total/1e6)) + geom_violin() + geom_boxplot(width=0.1, col="grey") + geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, fill="darkgrey") + ylab("total fragments (millions)") + xlab("") + th + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
```

![](singleCell_splitORexpression_files/figure-html/libSize-1.png)<!-- -->

Overall, the mapping statistics for all libraries are quite uniform. There is one outlier each in the proportion of unmapped and multimapped fragments, and four libraries with a higher fraction of reads mapping outside annotated exons. However, the samples showing outlying behaviour in the different metrics are all different, with no individual sample showing erratic metrics in two or more categories. Bad quality samples usually fail QC thresholds in most categories.


```r
mapping.stats$unmapped <- mapping.stats$total-(mapping.stats$unique+mapping.stats$multimapped)

stopifnot(identical(row.names(mapping.stats), row.names(counting.stats)))

plots <- list()
plots[[1]] <- ggplot(mapping.stats, aes(1, unmapped/total*100)) + geom_violin(trim=FALSE) + geom_jitter(aes(1,unmapped/total*100, col=mapping.stats$total/1e6), width = 0.01) + ylab("% unmapped") + scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Blues")[-1])) + labs(colour="libsize") + xlab("") + th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + geom_hline(yintercept = 5, lty=2, col="darkgrey")

plots[[2]] <- ggplot(mapping.stats, aes(1, multimapped/total*100)) + geom_violin(trim=FALSE) + geom_jitter(aes(1,multimapped/total*100, col=mapping.stats$total/1e6), width = 0.01) + ylab("% multimapped") + scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Blues")[-1])) + labs(colour="libsize") + xlab("") + th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + geom_hline(yintercept = 4.5, lty=2, col="darkgrey")

plots[[3]] <- ggplot(counting.stats, aes(1, N_noFeature/mapping.stats$total*100)) + geom_violin(trim=FALSE) + geom_jitter(aes(1,N_noFeature/mapping.stats$total*100, col=mapping.stats$total/1e6), width = 0.01) + ylab("% outside exons") + scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Blues")[-1])) + labs(colour="libsize") + xlab("") + th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + geom_hline(yintercept = 10, lty=2, col="darkgrey")

plots[[4]] <- ggplot(counting.stats, aes(1, N_ambiguous/mapping.stats$total*100)) + geom_violin(trim=FALSE) + geom_jitter(aes(1,N_ambiguous/mapping.stats$total*100, col=mapping.stats$total/1e6), width = 0.01) + ylab("% ambiguous") + scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Blues")[-1])) + labs(colour="libsize") + xlab("") + th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggarrange(plotlist = plots, ncol = 2, nrow = 2, common.legend = TRUE, legend = "right")
```

![](singleCell_splitORexpression_files/figure-html/mapping-1.png)<!-- -->

The proportion of fragments mapped to mitochondrial genes is very low for all samples. And in average, cells have around seven thousand genes detected; there is one outlier with only 3,164 detected genes.


```r
mt <- colSums(data[data$chr=="MT",-c(1:2)])
nGenes <- apply(data[,-c(1:2)], 2, function(x) sum(x>0))

plots <- list()
plots[[1]] <- ggplot(mapping.stats, aes(1, mt/total*100)) + geom_violin(trim=FALSE) + geom_jitter(aes(1, mt/total*100, col=mapping.stats$total/1e6), width = 0.01) + ylab("% mapped to MT genes") + scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Blues")[-1])) + labs(colour="libsize") + xlab("") + th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

plots[[2]] <- ggplot(mapping.stats, aes(1, nGenes)) + geom_violin(trim=FALSE) + geom_jitter(aes(1, nGenes, col=mapping.stats$total/1e6), width = 0.01) + ylab("total detected genes") + scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Blues")[-1])) + labs(colour="libsize") + xlab("") + th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + geom_hline(yintercept = 4e3, lty=2, col="darkgrey")

ggarrange(plotlist = plots, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
```

![](singleCell_splitORexpression_files/figure-html/mt-1.png)<!-- -->

Overall, all samples are of good-quality. We will remove only the sample with low number of genes detected (which also has a small library size). We also remove all genes that are not expressed in any cell.


```r
bad.qual <- row.names(mapping.stats)[which.min(nGenes)]
data <- data[,-which(colnames(data)==bad.qual)]
data <- data[rowSums(data[,-c(1:2)])>0,]
dim(data[,-c(1:2)])
```

```
## [1] 17511    33
```

```r
mapping.stats <- mapping.stats[-which(row.names(mapping.stats)==bad.qual),]
counting.stats <- counting.stats[-which(row.names(counting.stats)==bad.qual),]
```

### Normalisation

Next we normalise the data to remove composition biases between samples. We use the method implemented in `scran`. The calculated size factors are well correlated with the library sizes.


```r
## create object
sce <- SingleCellExperiment(assays=list(counts=as.matrix(data[,-c(1:2)])))
rowData(sce)$gene <- data$gene
row.names(sce) <- uniquifyFeatureNames(row.names(assay(sce)), rowData(sce)$gene)

## size factors
sce <- computeSumFactors(sce, sizes=c(10,15,20))

plot(colSums(assay(sce))/1e6, sizeFactors(sce), pch=16, bty="l", xlab="library size (million fragments)", ylab="size factor")
abline(lm(sizeFactors(sce)~as.vector(colSums(assay(sce))/1e6)), col="red")
```

![](singleCell_splitORexpression_files/figure-html/norm-1.png)<!-- -->

```r
## normalise
sce <- normalize(sce)
```

### OR expression

OSNs express one OR gene only, and silence all others. This rule is only violated when the OR chosen for expression is a pseudogene, in which case a different OR gene is expressed until a functional one is found. All but one of the 33 OSNs express a functional OR at high levels.


```r
ors <- logcounts(sce)
ors <- ors[grep("Olfr", row.names(ors)),]

or.expr <- matrix(ncol=5, nrow=ncol(ors))
or.names <- or.expr
or.class <- or.expr

for(i in 1:ncol(ors)){
  tmp <- ors[,i][order(ors[,i], decreasing = TRUE)]
  or.expr[i,] <- tmp[1:5]
  or.names[i,] <- names(tmp[1:5])
  or.class[i,] <- ann[match(names(tmp[1:5]), ann$gene),]$biotype
}
n <- c("first", "second", "third", "fourth", "fifth")
row.names(or.expr) <- colnames(ors); colnames(or.expr) <- n; or.expr <- as.data.frame(or.expr)
row.names(or.names) <- colnames(ors); colnames(or.names) <- n;
row.names(or.class) <- colnames(ors); colnames(or.class) <- n;

ggplot(as.data.frame(or.expr), aes(1:nrow(or.expr), first, col=mapping.stats$total/1e6)) + geom_point(size=2.5) + xlab("single OSNs") + ylab("log2 normalised expression") + ggtitle("Highest expressed OR gene") + th + labs(col="libSize") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
```

![](singleCell_splitORexpression_files/figure-html/topORs-1.png)<!-- -->

```r
# median(or.expr$first) # 15.58
# sd(or.expr$first) # 2.48
```

The one cell without abundant OR expression expresses instead a TAAR gene, a different class of receptor protein, also found in a small proportion of OSNs.


```r
other <- logcounts(sce)
other <- other[grep("Taar|Vmn|Gucy2d", row.names(other)),]

tmp <- other[,which(colnames(other) == row.names(or.expr)[which.min(or.expr$first)])]
tmp[tmp>0]
```

```
##    Taar4 
## 13.27016
```

So we add this gene as the receptor with highest expression for that cell.


```r
## add Taar4 to the or.expr data
idx <- which.min(or.expr$first)
or.expr[idx,][-1] <- or.expr[idx,][-5]
or.expr[idx,][1] <- tmp[tmp>0]

or.names[idx,][-1] <- or.names[idx,][-5]
or.names[idx,][1] <- names(tmp[tmp>0])

ggplot(as.data.frame(or.expr), aes(1:nrow(or.expr), first, col=mapping.stats$total/1e6, shape=ifelse(grepl("Taar", or.names[,1]), "TAAR", "OR"))) + geom_point(size=2.5) + xlab("single OSNs") + ylab("log2 normalised expression") + ggtitle("Highest expressed OR|TAAR gene") + th + labs(col="libSize", shape="") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylim(c(3,18))
```

![](singleCell_splitORexpression_files/figure-html/addTAAR-1.png)<!-- -->

All 33 cells express a different receptor gene.


```r
or.names[,1][order(or.names[,1])]
```

```
##        P3OMP31         P3OMP1        P3OMP39        P3OMP33        P3OMP38 
##     "Olfr1014"     "Olfr1048"     "Olfr1167"     "Olfr1195"     "Olfr1263" 
##        P3OMP32        P3OMP35        P3OMP36         P3OMP7        P3OMP34 
##     "Olfr1284"     "Olfr1299"     "Olfr1340"     "Olfr1348"     "Olfr1366" 
##        P3OMP13         P3OMP5         P3OMP8        P3OMP17        P3OMP29 
## "Olfr1369-ps1"     "Olfr1382"     "Olfr1444"     "Olfr1491"       "Olfr17" 
##         P3OMP4        P3OMP19        P3OMP27        P3OMP41         P3OMP6 
##       "Olfr19"       "Olfr43"      "Olfr481"      "Olfr536"      "Olfr576" 
##        P3OMP11        P3OMP18         P3OMP2        P3OMP14        P3OMP20 
##      "Olfr591"      "Olfr670"      "Olfr691"  "Olfr718-ps1"       "Olfr73" 
##        P3OMP12        P3OMP22        P3OMP21         P3OMP9        P3OMP15 
##      "Olfr740"      "Olfr743"  "Olfr766-ps1"      "Olfr799"        "Olfr9" 
##        P3OMP16        P3OMP28        P3OMP23 
##      "Olfr975"      "Olfr987"        "Taar4"
```

And the OR gene tends to be one of the most abundant genes expressed in the cell, usually within the top five most highly expressed genes (29 of the 33 cells).


```r
dataNorm <- logcounts(sce)

ranks <- sapply(1:nrow(or.expr), function(x) which(dataNorm[order(dataNorm[,x], decreasing = TRUE),x] == or.expr[x,1]))

ggplot(as.data.frame(ranks), aes(1:length(ranks), ranks)) + geom_point() + xlab("OSNs") + ylab("OR rank in transcriptome") + th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + geom_hline(yintercept = 5, lty=2, col="grey")
```

![](singleCell_splitORexpression_files/figure-html/ORrank-1.png)<!-- -->


And all but one show clear monogenic expression, with the second highest OR gene expressed at least a hundred times less than the top OR gene.


```r
or.expr <- 2^or.expr-1

ggplot(as.data.frame(or.expr), aes(1, log10(or.expr$first/or.expr$second))) + geom_boxplot() + geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) + xlab("single OSNs") + ylab("first / second OR") + th + labs(col="exprFirstOR", shape="biotype") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + scale_y_continuous(breaks = 1:5, labels = c(10,100,1000,1e4,1e5))
```

![](singleCell_splitORexpression_files/figure-html/secondOR-1.png)<!-- -->

The exception is P3OMP33:


```r
rbind(or.names[which.max(or.expr$second),], or.expr[which.max(or.expr$second),])
```

```
##                    first           second            third
## 1               Olfr1195         Olfr1193     Olfr1372-ps1
## P3OMP33 12515.8431917587 3288.98904946998 627.618307803594
##                   fourth            fifth
## 1               Olfr1284           Olfr17
## P3OMP33 1.92226127964347 1.92226127964347
```

*Olfr1195* and *Olfr1193* are adjacent genes. *Olfr1195* shows robust expression. Additional reads map uniquely to *Olfr1193*, including the coding sequence, but these do no support the annotated gene model. Instead, the data shows an additional isoform that spans both genes. Thus, it is likely that this cell expresses *Olfr1195* monogenically, and that some of its reads have been incorrectly assigned to *Olfr1195* due to lack of annotation of additional isoforms. Since the two genes are in opposite strands, strand-specific RNA-seq could resolve this issue.

Most other cells have only tens of counts on the second most abundant OR gene. There are six cells with more than 45 counts from protein-coding genes, generally a few hundred. In four of the six, the second OR is adjacent to the top expressed gene. In most cases, the second receptor is a close paralogue of the top receptor. In these cases, the second OR gets a certain amount of reads multimapped to both receptors; these are not counted and thus the expression of the top receptor is underestimated. Often a small number of reads also map uniquely, but with several mismatches that are consistent across all reads, a sign of mismapping due to sequencing errors. Additionally, some of the UTRs are too short and thus reads from the top receptor spill over to the second receptor, just because of proximity. In all these cases, these estimates of hundreds of counts are spurious. Thus, after accounting for mapping and annotation inaccuracies, mongenic expression is strongly observed in all OSNs.

Seven cells show expression of annotated pseudogenes, and in five of these, the pseudogene is the second most abundant OR locus expressed. Many are expressed at low levels (around 20 normalised counts) but three have several hundred normalised counts. Such expression levels, however, are still very low compared to the expression of the functional ORs, which have at least ~10,000 counts.

### Split OR expression

Within the 33 cells, two of them express abundantly OR genes that we identified as split OR genes, that is, intact open reading frames coded across two exons. Both of these OR genes -*Olfr718-ps1* and *Olfr766-ps1*- are expressed at similar levels to OR genes contained within a single exon.


```r
split <- unique(c("Olfr104-ps","Olfr104-ps","Olfr105-ps","Olfr106-ps","Olfr1116","Olfr1117-ps1","Olfr1118","Olfr1123","Olfr1174-ps","Olfr1175-ps","Olfr1177-ps","Olfr1183","Olfr1289","Olfr1291-ps1","Olfr1293-ps","Olfr1331","Olfr1333","Olfr1358","Olfr1480","Olfr239","OR5BS1P","Olfr286","Olfr287","Olfr288","Olfr324","Olfr324","Olfr55","Olfr560","Olfr592","Olfr607","Olfr680-ps1","Olfr682-ps1","Olfr718-ps1","Olfr735","Olfr745","Olfr764-ps1","Olfr766-ps1","Olfr844","Olfr857","Olfr869","Olfr872","Olfr873","Olfr873","Olfr18","Olfr94"))

is.split <- rep("singleExon", nrow(or.expr))
is.split[which(or.names[,1] %in% split)] <- "multiExon"

ggplot(or.expr, aes(1, first/1e3)) + geom_boxplot() + geom_jitter(width = 0.03, aes(col=is.split)) + xlab("") + ylab("normalised counts x 1000") + labs(col="") + th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
```

![](singleCell_splitORexpression_files/figure-html/split-1.png)<!-- -->

And both are able to induce monogenic expression, since the second highest OR is expressed several hundred to several thousand times lower.


```r
tmp <- data.frame(rank=rep(colnames(or.expr), each=nrow(or.expr)), expr=c(or.expr$first, or.expr$second, or.expr$third, or.expr$fourth, or.expr$fifth))
tmp$rank <- factor(tmp$rank, levels=colnames(or.expr)) # relevel

is.pseudo <- c(or.class[,1], or.class[,2], or.class[,3], or.class[,4], or.class[,5])
is.pseudo <- ifelse(is.pseudo == "protein_coding", "protein_coding", "pseudogene")

ggplot(tmp, aes(rank, log10(expr))) + geom_boxplot() + geom_jitter(width=0.05, aes(shape=is.pseudo, col=rep(is.split,5))) + labs(col="", shape="") + xlab("") + ylab("log10 normalised counts") + th
```

![](singleCell_splitORexpression_files/figure-html/monogenic-1.png)<!-- -->

These data suggest that the split OR genes encode functional receptors that behave in the same way as those encoded within a single exon.



```r
sessionInfo()
```

```
## R version 3.5.3 (2019-03-11)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: OS X El Capitan 10.11.6
## 
## Matrix products: default
## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] RColorBrewer_1.1-2          ggpubr_0.2                 
##  [3] magrittr_1.5                scater_1.10.1              
##  [5] ggplot2_3.1.0               scran_1.10.2               
##  [7] SingleCellExperiment_1.4.1  SummarizedExperiment_1.12.0
##  [9] DelayedArray_0.8.0          matrixStats_0.54.0         
## [11] Biobase_2.42.0              GenomicRanges_1.34.0       
## [13] GenomeInfoDb_1.18.2         IRanges_2.16.0             
## [15] S4Vectors_0.20.1            BiocGenerics_0.28.0        
## [17] BiocParallel_1.16.6        
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.0               locfit_1.5-9.1          
##  [3] lattice_0.20-38          assertthat_0.2.0        
##  [5] digest_0.6.18            R6_2.4.0                
##  [7] plyr_1.8.4               dynamicTreeCut_1.63-1   
##  [9] evaluate_0.13            pillar_1.3.1            
## [11] zlibbioc_1.28.0          rlang_0.3.1             
## [13] lazyeval_0.2.1           Matrix_1.2-15           
## [15] rmarkdown_1.12           labeling_0.3            
## [17] BiocNeighbors_1.0.0      statmod_1.4.30          
## [19] stringr_1.4.0            igraph_1.2.4            
## [21] RCurl_1.95-4.12          munsell_0.5.0           
## [23] HDF5Array_1.10.1         vipor_0.4.5             
## [25] compiler_3.5.3           xfun_0.5                
## [27] pkgconfig_2.0.2          ggbeeswarm_0.6.0        
## [29] htmltools_0.3.6          tidyselect_0.2.5        
## [31] tibble_2.0.1             gridExtra_2.3           
## [33] GenomeInfoDbData_1.2.0   edgeR_3.24.3            
## [35] viridisLite_0.3.0        withr_2.1.2             
## [37] crayon_1.3.4             dplyr_0.8.0.1           
## [39] bitops_1.0-6             grid_3.5.3              
## [41] gtable_0.2.0             scales_1.0.0            
## [43] stringi_1.4.3            reshape2_1.4.3          
## [45] XVector_0.22.0           viridis_0.5.1           
## [47] limma_3.38.3             DelayedMatrixStats_1.4.0
## [49] cowplot_0.9.4            Rhdf5lib_1.4.2          
## [51] tools_3.5.3              glue_1.3.0              
## [53] beeswarm_0.2.3           purrr_0.3.1             
## [55] yaml_2.2.0               colorspace_1.4-0        
## [57] rhdf5_2.26.2             knitr_1.22
```

