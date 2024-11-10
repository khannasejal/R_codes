library(readxl)
library(dplyr)
library(DESeq2)

#tnbc cell lines:: BT549, 76NF2V, SUM149, MB231, MB436, HCC1937
#non_tnbc cell lines:: MCF7, UACC812, T47D, MB361, ZR751

folder <- "D:\\sejal\\bam_bai_folder\\tnbc\\counts"
prefixes = c("BT549", "76NF2V", "SUM149", "MB231", "MB436", "HCC1937")
for (prefix in prefixes) {file_pattern <- paste0("^", prefix, ".*\\.bed")

#List all files in the folder with the given prefix and ".bed" extension
bed_files <- list.files(path = folder, pattern = file_pattern, full.names = TRUE)
cts <- read.table(bed_files, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
cts_matrix = as.matrix(cts)
View(cts_matrix)
cts_mx_unique = cts_matrix[!duplicated(cts_matrix[,4]),]    #remove the duplicate enhancer regions
rownames(cts_mx_unique) = cts_mx_unique[,4]  #put enhancer names as row names
View(cts_mx_unique)

#read the metadata file with condition information
metadata <- read.csv("metadata.csv",row.names = 1)
all(rownames(metadata) %in% colnames(cts_mx_unique))
cts_mx_unique <- cts_mx_unique[, rownames(metadata)]
View(cts_mx_unique)
all(rownames(metadata) == colnames(cts_mx_unique))
View(metadata)
mode(cts_mx_unique) <- "integer"
metadata$condition <- factor(metadata$condition)
metadata$cell_line <- factor(metadata$cell_line)

dds <- DESeqDataSetFromMatrix(countData = cts_mx_unique, colData = metadata, design = ~condition)
View(counts(dds))
dds_deseq <- DESeq(dds)
print(dds_deseq)
resultsNames(dds_deseq)
res<- results(dds_deseq)
View(res)
diffse <- paste0(prefix,"_deseq.txt")
write.table(res, file = diffse, sep ="\t", quote = FALSE)}