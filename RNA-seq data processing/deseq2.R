#####DESeq2 Analysis#####

library("DESeq2")

mainDir <- "/home/lailab/bxin/data"   #working directory
method <- "deseq2"

setwd(mainDir)

if (file.exists(method)){
    out_folder <- file.path(mainDir, method)
} else {
    dir.create(file.path(mainDir, method))
    out_folder <- file.path(mainDir, method)
}

threshold <- 0.05

### define functions
preparationD2 <- function (table){
  countData <- read.table(table,row.names=1,header=TRUE)
  colData <- data.frame(row.names= paste0("maize",1:ncol(countData)), condition= c(rep("high",33), rep("low",34)))
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~ condition)
  colData(dds)$condition <- factor(colData(dds)$condition,
                                   levels=c("high","low"))
  library("BiocParallel")
  register(MulticoreParam(20))
  dds <- DESeq(dds, parallel=TRUE)
  res <- as.data.frame(results(dds))
  return(res)
}

processingD2thre <- function(res,thre){
  test <- subset(res, res$padj< thre)
  return(test)
}

processingD2ranked <- function(res){
  test <- res[order(res$padj,res$pvalue),]
  return(test)
}


table_file <- "combined_htseq.txt"
res <- preparationD2(table_file)
sorted_res <- processingD2ranked(res)
DEresult <- processingD2thre(res,threshold)

save(sorted_res, file=paste0(out_folder, "/sorted_res.RData"))

write.table(DEresult,paste0(out_folder, "/", method, "_result.txt"),sep="\t",row.names=TRUE)
