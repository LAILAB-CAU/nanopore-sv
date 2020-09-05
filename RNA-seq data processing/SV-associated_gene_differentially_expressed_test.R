require(ggplot2)
require(scales)

rm(list = ls())

## make sure you add column names to the genotype files
ins = read.table("INS.genotype_gene_sv_pair.txt", header = TRUE)
del = read.table("DEL.genotype_gene_sv_pair.txt", header = TRUE)

data = read.table("DTA73MX-mBsnpREF-200425-expressed.csv",sep=",", head=T,row.names=1)
dd <- matrix(as.numeric(as.matrix(data)), nrow=dim(data)[1],ncol=dim(data)[2])
rownames(dd) <- rownames(data)
colnames(dd) <- as.character(lapply(colnames(data), function(x) substr(x,1,4)))
rm(data)

del_p_value_greater = c()
del_p_value_less = c()
del_sv_chrs = c()
del_sv_start = c()
del_gene_name = c()
del_gene_region = c()

support_threshold = 5

for (i in 1:dim(del)[1]){
  gene_name = as.character(del[i,7])
  if (gene_name %in% rownames(dd)){
    genotype = del[i,12:117]
    having_samples= colnames(genotype[which(genotype=="1/1")])
    not_having_samples= colnames(genotype[which(genotype=="0/0")])
    # 如果106个个体中support rate 不足5
    if (length(having_samples)< support_threshold || length(not_having_samples)< support_threshold){
      next
    }
    else{
      having_express = c()
      # Not all 106 accessions have gene expression profiles, only 73 accessions have RNA-seq data
      for (g in having_samples){
        if (g %in% colnames(dd)){having_express=c(having_express,dd[gene_name, g])}
      }
      not_having_express = c()
      for (g in not_having_samples){
        if (g %in% colnames(dd)){not_having_express=c(not_having_express,dd[gene_name, g])}
      }
      # if a group have less than 5 accessions that have RNA-seq data
      if (length(having_express)< support_threshold || length(not_having_express)< support_threshold){
        next
      }
      del_p_value_greater <- c(del_p_value_greater, wilcox.test(having_express,not_having_express, alternative = "greater")$p.value)
      del_p_value_less <- c(del_p_value_less, wilcox.test(having_express,not_having_express, alternative = "less")$p.value)
      del_sv_chrs <- c(del_sv_chrs,as.character(del[i,1]))
      del_sv_start <- c(del_sv_start, del[i,2])
      del_gene_name <- c(del_gene_name, as.character(del[i,7]))
      del_gene_region <- c(del_gene_region, as.character(del[i,8]))
    }
  }
}

del_df <- data.frame(cbind(del_sv_chrs, del_sv_start, del_gene_name, del_gene_region))
del_df$del_p_value_greater <- del_p_value_greater
del_df$del_p_value_less <- del_p_value_less
save(del_df, file=paste0(folder, "/del_df.RData"))

del_df = na.omit(del_df)
## adjust p_values
del_df$BH_greater = signif(p.adjust(del_df$del_p_value_greater,method = "BH"),4)
del_df$BH_less = signif(p.adjust(del_df$del_p_value_less,method = "BH"),4)
significant_del_df_bh = del_df[del_df$BH_greater<0.05|del_df$BH_less<0.05,]

write.table(significant_del_df_bh, "del_gene_pair_result_less0.05_BH.txt",sep="\t",row.names=FALSE, quote = FALSE)

ins_p_value_greater = c()
ins_p_value_less = c()
ins_sv_chrs = c()
ins_sv_start = c()
ins_gene_name = c()
ins_gene_region = c()

support_threshold = 5

for (i in 1:dim(ins)[1]){
  gene_name = as.character(ins[i,7])
  if (gene_name %in% rownames(dd)){
    genotype = ins[i,12:117]
    having_samples= colnames(genotype[which(genotype=="1/1")])
    not_having_samples= colnames(genotype[which(genotype=="0/0")])
    if (length(having_samples)< support_threshold || length(not_having_samples)< support_threshold){
      next
    }
    else{
      having_express = c()
      for (g in having_samples){
        if (g %in% colnames(dd)){having_express=c(having_express,dd[gene_name, g])}
      }
      not_having_express = c()
      for (g in not_having_samples){
        if (g %in% colnames(dd)){not_having_express=c(not_having_express,dd[gene_name, g])}
      }
      if (length(having_express)< support_threshold || length(not_having_express)< support_threshold){
        next
      }
      ins_p_value_greater <- c(ins_p_value_greater, wilcox.test(having_express,not_having_express, alternative = "greater")$p.value)
      ins_p_value_less <- c(ins_p_value_less, wilcox.test(having_express,not_having_express, alternative = "less")$p.value)
      ins_sv_chrs <- c(ins_sv_chrs,as.character(ins[i,1]))
      ins_sv_start <- c(ins_sv_start, ins[i,2])
      ins_gene_name <- c(ins_gene_name, as.character(ins[i,7]))
      ins_gene_region <- c(ins_gene_region,as.character(ins[i,8]))
    }
  }
}


ins_df <- data.frame(cbind(ins_sv_chrs, ins_sv_start, ins_gene_name, ins_gene_region))
ins_df$ins_p_value_greater <- ins_p_value_greater
ins_df$ins_p_value_less <- ins_p_value_less
ins_df = na.omit(ins_df)

## adjust p_values
ins_df$BH_greater = signif(p.adjust(ins_df$ins_p_value_greater,method = "BH"),4)
ins_df$BH_less = signif(p.adjust(ins_df$ins_p_value_less,method = "BH"),4)
significant_ins_df_bh = ins_df[ins_df$BH_greater<0.05|ins_df$BH_less<0.05,]
write.table(significant_ins_df_bh,"ins_gene_pair_result_less0.05_BH.txt",sep="\t",row.names=FALSE, quote = FALSE)
