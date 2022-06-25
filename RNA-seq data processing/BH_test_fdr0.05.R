fdr = 0.05

d2 = read.table("permutations_covariates.all.chunks.txt.gz", header=F, stringsAsFactors=F)
colnames(d2) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope", "ppval", "bpval")
plot(d2$ppval, d2$bpval, xlab="Direct method", ylab="Beta approximation", main="Check plot")
abline(0, 1, col="red")

# BH correction for multiple phenotypes tested. 
d2$bh = p.adjust(d2$bpval, method="fdr")

#Then, you can extract all significant MP-QTL pairs with a 10% false discovery rate (FDR) for instance and write them to a file permutations.all.chunks.benjamini.txt using:
write.table(d2[which(d2$bh <= fdr),], "permutations_covariates.all.chunks.benjamini.txt", quote=F, row.names=F, col.names=T)