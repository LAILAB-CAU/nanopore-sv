#####SAMseq Analysis#####

library(samr)

method <- "samseq"
mainDir <- "/home/lailab/bxin/data"   #working directory

setwd(mainDir)

if (file.exists(method)){
    out_folder <- file.path(mainDir, method)
} else {
    dir.create(file.path(mainDir, method))
    out_folder <- file.path(mainDir, method)
}

threshold <- 0.05

preparationSthre <- function(table){
  counts <- read.table(table, header = TRUE, row.names= 1)
  group <- c(rep(1, 33), rep(2, 34))
  samfit <- SAMseq(counts, group ,geneid=as.character(rownames(counts)),genenames=as.character(rownames(counts)), resp.type= "Two class unpaired", nperms= 1000, random.seed=1234567, fdr.output=0.05)
  #SAMseq() default:
  #SAMseq <- function(x, y, censoring.status = NULL, 
  #resp.type = c("Quantitative", "Two class unpaired", "Survival", "Multiclass", "Two class paired"), 
  #geneid = NULL, genenames = NULL, nperms = 100, random.seed = NULL, nresamp = 20, fdr.output = 0.20)
  return(samfit)
}

namedThre <- function(samseq)
{
  temporary <- samseq[[4]]
  temporary <- rbind(temporary$genes.up,temporary$genes.lo)
  temporary <- temporary[order(as.numeric(temporary[,5])),,drop=FALSE]
  temporary <- temporary[,-1,drop=FALSE]
  #Double check for FDR less than 0.05 
  truth <- NULL
  for (i in 1:nrow(temporary)){
    if (as.numeric(temporary[i,4]) <= 5){
      truth <- rbind(truth, temporary[i,])
    }}
  return(truth)
}

preparationSranked <- function(table){
  counts <- read.table(table, header = TRUE, row.names= 1)
  group <- c(rep(1, 33), rep(2, 34))
  samfit <- SAMseq(counts, group ,geneid=as.character(rownames(counts)),genenames=as.character(rownames(counts)), resp.type= "Two class unpaired", nperms= 1000, random.seed=1234567,fdr.output=1)
  return(samfit)
}

namedRank <- function(samseq){
  temporary <- samseq[[4]]
  temporary <- rbind(temporary$genes.up, temporary$genes.lo)
  temporary <- temporary[order(as.numeric(temporary[, 5])),,drop=FALSE ]
  temporary <- temporary[,-1,drop=FALSE]
  return(temporary)
}


table_file <- "combined_htseq.txt"
res_thre <- preparationSthre(table_file)
named_res_thre <- namedThre(res_thre)

sorted_res <- preparationSranked(table_file)
named_sorted_res <- namedRank(sorted_res)

save(named_sorted_res, file=paste0(out_folder, "/sorted_res.RData"))

write.table(named_res_thre,paste0(out_folder, "/",method, "_result.txt"),sep="\t",row.names=TRUE)