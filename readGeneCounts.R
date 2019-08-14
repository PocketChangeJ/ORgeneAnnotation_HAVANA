### script to collate the counts per gene from all cells

dir <- "/user01/group_folders/Personal/Ximena/ORs/singleCell/Khan/"

## gene annotation, parsed from Ensembl's GTF
ann <- read.table("/user01/group_folders/Personal/Ximena/REFERENCE/Mus_musculus.GRCm38.93.ann", stringsAsFactors=FALSE, row.names=1)
colnames(ann) <- c("gene", "chr", "start", "end", "strand", "biotype")

## cells processed, each in a separate directory. Get the folder names.
files <- list.files(paste0(dir, "STAR/"))
files <- files[c(1,12,28,30:34,2:3,5:11,13:27,29,4)] # reorder in numeric order

## read the data for the first cell and create the count matrix
## output is from STAR read counts per gene
data <- read.table(paste0(dir, "STAR/", files[1], "/ReadsPerGene.out.tab"), stringsAsFactors=FALSE)[,1:2]
stopifnot(identical(row.names(ann), data$V1[-c(1:4)]))
row.names(data) <- data$V1
data$V1 <- NULL
data$gene <- ann[match(row.names(data), row.names(ann)),1]
data <- data[,2:1]
colnames(data)[2] <- files[1]

## read files for all other cells and append to count matrix
for(f in files[-1]){
  data <- cbind(data, read.table(paste0(dir, "STAR/", f, "/ReadsPerGene.out.tab"), stringsAsFactors=FALSE)[,2])
}
colnames(data)[-1] <- files

## save
write.table(data, paste0(dir, "geneCounts_P3OMPcells.RAW.tsv"), quote=FALSE, sep="\t")
write.table(ann, paste0(dir, "geneAnnotation.Ensembl_v93.tsv"), quote=FALSE, sep="\t")
