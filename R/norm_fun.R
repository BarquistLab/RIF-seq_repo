#' logcpm values, required for normalization,
#' otherwise, the natural logarithm is used
#' throughout the analysis.
logcounts <- function(Y,N){
  t(log2(t(Y+0.5)/(N+1)*1e+06))
}

wgk <- function(Ygk,Nk){
  (1-Ygk/Nk)/Ygk
}

#' normalize sequencing libraries with trimmed mean of M values (TMM)
#' similar to implementaion in edgeR
#' ERCC spike-ins are not expected to have large log-fold changes, therefore
#' trimM and trimA are set to zero
#' @export
normTMM <- function(ercc_reads, N, trimM=0, trimA=0){

  trimmed_ercc <- subset(ercc_reads, apply(ercc_reads,1,function(x) all(x!=0)))

  # calculate logcpm values + mean as reference
  logcpm <- logcounts(trimmed_ercc, N)
  Nref <- mean(N)

  logref <- rowMeans(logcpm)
  ref_reads <- 2^logref*Nref * 1e-06
  
  # M (logFC) values
  # remove genes with largest M values
  logFC <- logcpm - logref
  mean_logFC <- rowMeans(logFC)
  n_row <- nrow(logFC)
  n_keep_M <- ceiling(n_row*(1-trimM))
  trim_logFC <- order(abs(mean_logFC))[1:n_keep_M]
  trimmed_ercc <- trimmed_ercc[trim_logFC, ]
  
  # re-calculate logcpm values
  logcpm <- logcounts(trimmed_ercc, N)
  Nref <- mean(N)
  
  logref <- rowMeans(logcpm)
  ref_reads <- 2^logref*Nref * 1e-06
  
  # A (logA) values
  # remove genes with largest A values (intensity) to avoid outliers
  logA <- (logcpm + logref)/2
  mean_logA <- rowMeans(logA)
  n_row <- nrow(logA)
  n_keep_A <- ceiling(n_row*(1-trimA))
  trim_logA <- order(abs(mean_logA))[1:n_keep_A]
  trimmed_ercc <- trimmed_ercc[trim_logA, ]
  
  # re-calculate logcpm values
  logcpm <- logcounts(trimmed_ercc, N)
  Nref <- mean(N)
  
  logref <- rowMeans(logcpm)
  ref_reads <- 2^logref*Nref * 1e-06

  # calculate weights and logFC
  logFC <- logcpm - logref
  weights <- matrix(nrow=nrow(trimmed_ercc), ncol=ncol(trimmed_ercc))
  for(i in 1:ncol(trimmed_ercc)){
    weights[,i] <- wgk(trimmed_ercc[,i],N[i]) + wgk(ref_reads,Nref)
  }

  # calculate and return TMM values
  logTMM <- colSums(logFC/weights)/colSums(1/weights)
  2^logTMM/2^mean(logTMM)
}

#' remove genes/ercc spike-ins with very low/no counts, return TMM normalization 
#' with ERCC spike-ins keeping the original library size
#' @export
get_norm_facs <- function(rif_reads, ercc_reads, count_min=10, sample_min=3) {
  N <- colSums(rif_reads)
  keep <- which(rowSums(cpm(rif_reads, N) > count_min) >= sample_min)
  rif_reads <- rif_reads[keep,]
  keep <- which(rowSums(cpm(ercc_reads, N) > count_min) >= sample_min)
  ercc_reads <- ercc_reads[keep,]
  N <- colSums(rif_reads)
  normTMM(ercc_reads, N)
}
