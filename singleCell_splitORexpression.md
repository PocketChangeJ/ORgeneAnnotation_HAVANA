---
title: "Single-cell RNA-seq of OSNs"
date: '26 August, 2019'
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
ann <- read.table(paste0(dir, "data/geneAnnotation.Ensembl_v93.tsv"), stringsAsFactors = FALSE)

data <- read.table(paste0(dir, "data/geneCounts_P3OMPcells.RAW.tsv"), stringsAsFactors = FALSE)
mapping.stats <- read.table(paste0(dir, "data/mappingStats.tsv"), stringsAsFactors = FALSE, header = TRUE, row.names = 1)

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

The one cell without abundant OR expression (`row.names(or.expr)[which.min(or.expr$first)]`) expresses instead a TAAR gene, a different class of receptor protein, also found in a small proportion of OSNs.


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

#### Index switching

These cells were sequenced on the HiSeq X Ten platform, which suffers from index switching. This leads to some reads from a given library being incorrectly assigned to other libraries that share one of the indices used for multiplexing. Below is a heatmap of all the OR genes detected, and the coloured bar at the top indicates one of the indices used for multiplexing. It is clear that among the libraries sharing the same index there is index switching, i.e., a small number of reads from the very highly expressed OR gene from each cell bleeding into the others cells with the same index.


```r
multiplexing <- read.table(paste0(dir, "data/indices.tsv"), header = TRUE, stringsAsFactors = FALSE)
multiplexing <- multiplexing[match(colnames(ors), multiplexing$cell),]

heatmap.2(as.matrix(log10(ors+1)), trace="none", col=rev(brewer.pal(n=10, "RdYlBu")), Colv = FALSE, ColSideColors = as.character(factor(multiplexing$index2, labels = rainbow(n=3))), dendrogram = "row", key.title = "", key.xlab = "log10 conts")
```

![](singleCell_splitORexpression_files/figure-html/multiplexing-1.png)<!-- -->

Thus, for all the most highly expressed OR genes per cell, we set the counts in other cells to zero.


```r
## remove the counts of top expressed genes in other cells, since they are bleed through
for(i in or.names[,1][-21]){
  ors[i,-which.max(ors[i,])] <- 0
}

## recompute ranking
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

## add Taar4 to the or.expr data
tmp <- other[,which(colnames(other) == row.names(or.expr)[which.min(or.expr$first)])]
or.expr[idx,][-1] <- or.expr[idx,][-5]
or.expr[idx,][1] <- tmp[tmp>0]

or.names[idx,][-1] <- or.names[idx,][-5]
or.names[idx,][1] <- names(tmp[tmp>0])
```

#### Monogenic OR expression

All but one cell show clear monogenic expression, with the second highest OR gene expressed at least a hundred times less than the top OR gene.


```r
or.expr <- 2^or.expr-1

ggplot(as.data.frame(or.expr), aes(1, log10(or.expr$first/(or.expr$second+1)))) + geom_boxplot() + geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) + xlab("single OSNs") + ylab("first / second OR") + th + labs(col="exprFirstOR", shape="biotype") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + scale_y_continuous(breaks = 1:5, labels = c(10,100,1000,1e4,1e5))
```

![](singleCell_splitORexpression_files/figure-html/secondOR-1.png)<!-- -->

The exception is P3OMP33:


```r
rbind(or.names[which.max(or.expr$second),], or.expr[which.max(or.expr$second),])
```

```
##                    first           second            third          fourth
## 1               Olfr1195         Olfr1193     Olfr1372-ps1        Olfr1383
## P3OMP33 12515.8433055288 3288.98907936719 627.618313508701 0.9611306485585
##                   fifth
## 1              Olfr1370
## P3OMP33 0.9611306485585
```

*Olfr1195* and *Olfr1193* are adjacent genes, on opposite strands. *Olfr1195* shows robust expression. Additional reads map uniquely to *Olfr1193*, including the coding sequence, but these do no support the annotated gene model. Instead, the data shows an additional isoform that spans both genes. Thus, it is likely that this cell expresses *Olfr1195* monogenically, and that some of its reads have been incorrectly assigned to *Olfr1195* due to lack of annotation of additional isoforms. If indeed a new isoform covers both genes, and is transcribed from the reverse strand, then only OLFR1195 would be translated; the splice data supports this. 

![](pics/Olfr1195.png)

Nonetheless, since the missing isoform does cover *Olfr1193* coding sequence, we cannot rule out that this gene is being expressed from the forward strand, and a second OR protein generated. 

Most other cells have only tens of counts on the second most abundant OR gene. 


```r
tmp <- data.frame(cell=rep(row.names(or.expr),2), rank=rep(colnames(or.expr)[1:2], each=nrow(or.expr)), expr=c(or.expr$first, or.expr$second))
tmp$rank <- factor(tmp$rank, levels=colnames(or.expr)) # relevel

is.pseudo <- c(or.class[,1], or.class[,2])
is.pseudo <- ifelse(is.pseudo == "protein_coding", "protein_coding", "pseudogene")

ggplot(tmp, aes(rank, log10(expr+1))) + geom_boxplot() + geom_jitter(width=0.05, aes(col=is.pseudo)) + labs(col="") + xlab("") + ylab("log10 normalised counts") + th #+ geom_text_repel(data = tmp[tmp$rank=="second" & tmp$expr > 30,], label=tmp[tmp$rank=="second" & tmp$expr > 30,]$cell, cex=2.5)
```

![](singleCell_splitORexpression_files/figure-html/top2-1.png)<!-- -->

There are seven cells with more than 30 counts from protein-coding genes, generally a few hundred, except for P3OMP33 which we have already discussed.


```r
tmp <- or.expr[or.expr$second>45,]
tmp[or.class[or.expr$second>45,2]=="protein_coding",1:2]
```

```
##            first     second
## P3OMP2  81182.42  225.98529
## P3OMP7  67784.54   79.44556
## P3OMP11 49401.91   47.98269
## P3OMP19 80479.42  319.48424
## P3OMP21 59303.92  153.90643
## P3OMP31 78982.71  584.09473
## P3OMP33 12515.84 3288.98908
```

In five of the seven, the second OR is adjacent to the top expressed gene, and often these are paralogues with high identity. This is suggestive of mismapping from the top to the second OR. 


```r
tmp <- or.names[or.expr$second>45,]
tmp <- as.data.frame(tmp[or.class[or.expr$second>45,2]=="protein_coding",1:2])
tmp$adjacent <- c("yes", "yes", "skip1", "yes", "no", "yes", "yes")
tmp$identity <- c(71.11,92.63,89.81,98.72,NA,88.85,57.14) # from Ensembl's paalogue data
tmp
```

```
##               first   second adjacent identity
## P3OMP2      Olfr691  Olfr690      yes    71.11
## P3OMP7     Olfr1348 Olfr1347      yes    92.63
## P3OMP11     Olfr591  Olfr593    skip1    89.81
## P3OMP19      Olfr43  Olfr403      yes    98.72
## P3OMP21 Olfr766-ps1 Olfr1420       no       NA
## P3OMP31    Olfr1014 Olfr1013      yes    88.85
## P3OMP33    Olfr1195 Olfr1193      yes    57.14
```

This is the case for P3OMP7, P3OMP11, P3OMP19 and P3OMP31, where the second OR gene has mostly multimapped reads (unique reads are coloured blue; multimapped in any other colour) and a few uniquely mapped ones, but with several mismatches that are supported by several reads.

*Olfr1348*-*Olfr1347* in P3OMP7:

![](pics/Olfr1348.png)

*Olfr591*-*Olfr593* in P3OMP11:
![](pics/Olfr591.png)

*Olfr43*-*Olfr403* in P3OMP19:
![](pics/Olfr43.png)

*Olfr1014*-*Olfr1013* in P3OMP31:
![](pics/Olfr1014.png)

In the case of P3OMP2, it looks like the UTR from *Olfr691* is not fully annotated and the reads spill over to the first exon of *Olfr690*, but correspond still to *Olfr691*.

![](pics/Olfr691.png)

Finally, for P3OMP21, the second OR is not a neighbour of the top OR and they are not paralogous. The reads mapping to the second OR, *Olfr1420*, look like high-quality alignments. Thus, it seems that this second OR is transcribed, but at much lower levels than the top OR: 154 versus 59,304 normalised counts.

![](pics/Olfr1420.png)

Thus, in six out of the seven cases, we can dismiss the counts from the second most highly expressed OR.

In the case of the third OR gene, only one cell has more than 30 normalised counts for protein-coding OR genes.


```r
rbind(or.names[or.expr$third>30 & or.class[,3]=="protein_coding",], or.expr[or.expr$third>30 & or.class[,3]=="protein_coding",])
```

```
##                    first           second            third
## 1                 Olfr43          Olfr403          Olfr855
## P3OMP19 80479.4204080011 319.484244749838 113.828218937816
##                   fourth             fifth
## 1            Olfr408-ps1          Olfr1508
## P3OMP19 1.91307930987927 0.956539654939634
```

The second OR in this cell is product of mismapping, but the third OR shows good quality alignments. *Olfr855* is a distant paralogue of *Olfr43*, with only 42.49% identity between the two. Thus, these counts seem genuine.

![](pics/Olfr855.png)

Thus, we remove the counts of the second highest ORs where we've shown mismapping or spill-over. This leaves only two cells with a few hundred fragments mapped to a second OR gene, and P3OMP33 with several thousand, that are most likely from the top OR (but we cannot be sure, so we don't remove them).


```r
## for these six cells, we shift the counts from the third, to second and so on
or.expr[row.names(tmp)[-c(5,7)],2:4] <- or.expr[row.names(tmp)[-c(5,7)],3:5]
or.names[row.names(tmp)[-c(5,7)],2:4] <- or.names[row.names(tmp)[-c(5,7)],3:5]
or.class[row.names(tmp)[-c(5,7)],2:4] <- or.class[row.names(tmp)[-c(5,7)],3:5]
## recover the sixth OR
for(i in row.names(tmp)[-c(5,7)]){
  tmp <- ors[,i][order(ors[,i], decreasing = TRUE)]
  or.expr[i,5] <- tmp[6]
  or.names[i,5] <- names(tmp[6])
  or.class[i,5] <- ann[match(names(tmp[6]), ann$gene),]$biotype
}

## plot
tmp <- data.frame(cell=rep(row.names(or.expr),2), rank=rep(colnames(or.expr)[1:2], each=nrow(or.expr)), expr=c(or.expr$first, or.expr$second))
tmp$rank <- factor(tmp$rank, levels=colnames(or.expr)) # relevel

is.pseudo <- c(or.class[,1], or.class[,2])
is.pseudo <- ifelse(is.pseudo == "protein_coding", "protein_coding", "pseudogene")

ggplot(tmp, aes(rank, log10(expr+1))) + geom_boxplot() + geom_jitter(width=0.05, aes(col=is.pseudo)) + labs(col="") + xlab("") + ylab("log10 normalised counts") + th #+ geom_text_repel(data = tmp[tmp$rank=="second" & tmp$expr > 30,], label=tmp[tmp$rank=="second" & tmp$expr > 30,]$cell, cex=2.5)
```

![](singleCell_splitORexpression_files/figure-html/top2_mismapping-1.png)<!-- -->

Ten cells show expression of annotated pseudogenes, and in seven the pseudogene is the second most abundant OR locus expressed. Many are expressed at low levels (around 20 normalised counts) but three have several hundred normalised counts. Such expression levels, however, are still very low compared to the expression of the functional ORs, which have at least ~10,000 counts.


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

is.split <- factor(is.split, levels = c("singleExon", "multiExon"))

ggplot(tmp, aes(rank, log10(expr+1))) + geom_boxplot() + geom_jitter(width=0.05, aes(col=is.pseudo, shape=rep(is.split,5))) + labs(col="", shape="") + xlab("") + ylab("log10 normalised counts") + th
```

![](singleCell_splitORexpression_files/figure-html/monogenic-1.png)<!-- -->

These data suggest that the split OR genes encode functional receptors that behave in the same way as those encoded within a single exon.


```r
sessionInfo()
```

```
## R version 3.6.1 (2019-07-05)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS High Sierra 10.13.6
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] gplots_3.0.1.1              RColorBrewer_1.1-2         
##  [3] ggrepel_0.8.1               ggpubr_0.2.2               
##  [5] magrittr_1.5                scater_1.12.2              
##  [7] ggplot2_3.2.1               scran_1.12.1               
##  [9] SingleCellExperiment_1.6.0  SummarizedExperiment_1.14.1
## [11] DelayedArray_0.10.0         BiocParallel_1.18.1        
## [13] matrixStats_0.54.0          Biobase_2.44.0             
## [15] GenomicRanges_1.36.0        GenomeInfoDb_1.20.0        
## [17] IRanges_2.18.1              S4Vectors_0.22.0           
## [19] BiocGenerics_0.30.0        
## 
## loaded via a namespace (and not attached):
##  [1] viridis_0.5.1            dynamicTreeCut_1.63-1   
##  [3] edgeR_3.26.7             BiocSingular_1.0.0      
##  [5] viridisLite_0.3.0        DelayedMatrixStats_1.6.0
##  [7] gtools_3.8.1             assertthat_0.2.1        
##  [9] statmod_1.4.32           dqrng_0.2.1             
## [11] GenomeInfoDbData_1.2.1   vipor_0.4.5             
## [13] yaml_2.2.0               pillar_1.4.2            
## [15] lattice_0.20-38          glue_1.3.1              
## [17] limma_3.40.6             digest_0.6.20           
## [19] XVector_0.24.0           ggsignif_0.6.0          
## [21] colorspace_1.4-1         cowplot_1.0.0           
## [23] htmltools_0.3.6          Matrix_1.2-17           
## [25] pkgconfig_2.0.2          zlibbioc_1.30.0         
## [27] purrr_0.3.2              scales_1.0.0            
## [29] gdata_2.18.0             tibble_2.1.3            
## [31] withr_2.1.2              lazyeval_0.2.2          
## [33] crayon_1.3.4             evaluate_0.14           
## [35] beeswarm_0.2.3           tools_3.6.1             
## [37] stringr_1.4.0            munsell_0.5.0           
## [39] locfit_1.5-9.1           irlba_2.3.3             
## [41] compiler_3.6.1           rsvd_1.0.2              
## [43] caTools_1.17.1.2         rlang_0.4.0             
## [45] grid_3.6.1               RCurl_1.95-4.12         
## [47] BiocNeighbors_1.2.0      igraph_1.2.4.1          
## [49] labeling_0.3             bitops_1.0-6            
## [51] rmarkdown_1.15           gtable_0.3.0            
## [53] R6_2.4.0                 gridExtra_2.3           
## [55] knitr_1.24               dplyr_0.8.3             
## [57] KernSmooth_2.23-15       stringi_1.4.3           
## [59] ggbeeswarm_0.6.0         Rcpp_1.0.2              
## [61] tidyselect_0.2.5         xfun_0.8
```

