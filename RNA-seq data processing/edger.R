#####edgeR Analysis#####

library(limma)
library(edgeR)

method <- "edger"
mainDir <- "/home/lailab/bxin/data"   #working directory

setwd(mainDir)

if (file.exists(method)){
    out_folder <- file.path(mainDir, method)
} else {
    dir.create(file.path(mainDir, method))
    out_folder <- file.path(mainDir, method)
}

threshold <- 0.05

#####Classic exactTest:
preparationE <- function (table){
  countTable <- read.table(table,header=TRUE, row.names=1)
  group <- factor(c(rep("high", 33), rep("low",34)))
  dge <- DGEList (counts = countTable , group = group )

  dge <- calcNormFactors(dge,method="TMM")
  dge <- estimateCommonDisp(dge,verbose=TRUE)
  dge <- estimateTagwiseDisp(dge)
  res <- exactTest(dge)
}


processingEthre <- function(res, thre){
  tops <- topTags(res , n=nrow(res))
  #topTags() default:
  #topTags(object, n=10, adjust.method="BH", sort.by="PValue")
  tops <- tops[tops$table$FDR<thre,]
  return(tops)
}

processingEranked <-function(res){
  tops <- topTags(res , n=nrow(res))
  return(tops)
}

table_file <- "combined_htseq.txt"
res <- preparationE(table_file)
sorted_res <- processingEranked(res)
DEresult <- processingEthre(res,threshold)

save(sorted_res, file=paste0(out_folder, "/sorted_res.RData"))

write.table(DEresult,paste0(out_folder, "/",method, "_result.txt"),sep="\t",row.names=TRUE)
