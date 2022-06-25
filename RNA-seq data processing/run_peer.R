library(peer, quietly=TRUE)  # https://github.com/PMBio/peer

WriteTable <- function(data, filename, index.name) {
	datafile <- file(filename, open = "wt")
	on.exit(close(datafile))
	header <- c(index.name, colnames(data))
	writeLines(paste0(header, collapse="\t"), con=datafile, sep="\n")
	write.table(data, datafile, sep="\t", col.names=FALSE, quote=FALSE)
}

expr.file = "seedling.normalized.expression.txt"
prefix = "seedling"
n = 30
alphaprior_a=0.001
alphaprior_b=0.01
epsprior_a=0.1
epsprior_b=10
max_iter=1000
output_dir="./"


cat("PEER: loading expression data ... ")
if (grepl('.gz$', expr.file)) {
    nrows <- as.integer(system(paste0("zcat ", expr.file, " | wc -l | cut -d' ' -f1 "), intern=TRUE, wait=TRUE))
} else {
    nrows <- as.integer(system(paste0("wc -l ", expr.file, " | cut -d' ' -f1 "), intern=TRUE, wait=TRUE))
}
if (grepl('.bed$', expr.file) || grepl('.bed.gz$', expr.file)) {
    df <- read.table(expr.file, sep="\t", nrows=nrows, header=TRUE, check.names=FALSE, comment.char="")
    row.names(df) <- df[, 4]
    df <- df[, 5:ncol(df)]
} else {
    df <- read.table(expr.file, sep="\t", nrows=nrows, header=TRUE, check.names=FALSE, comment.char="", row.names=1)
}
M <- t(as.matrix(df))
cat("done.\n")

# run PEER
cat(paste0("PEER: estimating hidden confounders (", n, ")\n"))
model <- PEER()
invisible(PEER_setNk(model, n))
invisible(PEER_setPhenoMean(model, M))
invisible(PEER_setPriorAlpha(model, alphaprior_a, alphaprior_b))
invisible(PEER_setPriorEps(model, epsprior_a, epsprior_b))
invisible(PEER_setNmax_iterations(model, max_iter))
# if(!is.null(covs)) {
#   invisible(PEER_setCovariates(model, covs))
# }
time <- system.time(PEER_update(model))

X <- PEER_getX(model)  # samples x PEER factors
A <- PEER_getAlpha(model)  # PEER factors x 1
R <- t(PEER_getResiduals(model))  # genes x samples

# add relevant row/column names
c <- paste0("InferredCov",1:ncol(X))
rownames(X) <- rownames(M)
colnames(X) <- c
rownames(A) <- c
colnames(A) <- "Alpha"
A <- as.data.frame(A)
A$Relevance <- 1.0 / A$Alpha
rownames(R) <- colnames(M)
colnames(R) <- rownames(M)

# write results
cat("PEER: writing results ... ")
WriteTable(t(X), file.path(output_dir, paste0(prefix, "_PEER_covariates.txt")), "ID")  # format(X, digits=6)
WriteTable(A, file.path(output_dir, paste0(prefix, "_PEER_alpha.txt")), "ID")
WriteTable(R, file.path(output_dir, paste0(prefix, "_PEER_residuals.txt")), "ID")
cat("done.\n")
