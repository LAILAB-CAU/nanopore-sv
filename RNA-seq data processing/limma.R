library(limma)
library(edgeR)

method <- "limma"
mainDir <- "/home/lailab/bxin/data"   #working directory

setwd(mainDir)

if (file.exists(method)){
    out_folder <- file.path(mainDir, method)
} else {
    dir.create(file.path(mainDir, method))
    out_folder <- file.path(mainDir, method)
}

threshold <- 0.05

###Vector indicating number of samples in each subset:
preparationL <- function(table){
  counts <- read.table(table, header = TRUE, row.names = 1)
  mType <- factor(c(rep("high", 33), rep("low",34)))
  nf <- calcNormFactors(counts)

  design <- model.matrix(~mType)
  y <- voom (counts, design, lib.size = colSums(counts)*nf)
  fit <- lmFit(y,design)
  fit <- eBayes(fit)
  return(fit)
}


processingLthre <- function(res,coef,thre){
  test <- topTable(res,coef=coef,n=nrow(res))
  #topTable() default:
  #topTable(fit, coef=NULL, number=10, genelist=fit$genes, adjust.method="BH",
  #sort.by="B", resort.by=NULL, p.value=1, lfc=0, confint=FALSE)
  #coef is the column of interest,here :2 Male or Female
  test <- subset(test,test$adj.P.Val <thre)
  #thre: FDR threshold
  return(test)
}


processingLranked <- function(res,coef){
  test <- topTable(res,coef=coef,n=nrow(res))
  return(test)
}


table_file <- "combined_htseq.txt"
res <- preparationL(table_file)
sorted_res <- processingLranked(res,2)
DEresult <- processingLthre(res,2,threshold)

save(sorted_res, file=paste0(out_folder, "/sorted_res.RData"))

write.table(DEresult,paste0(out_folder, "/",method, "_result.txt"),sep="\t",row.names=TRUE)
