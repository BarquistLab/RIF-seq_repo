#' log cpm values, only for normalization,
#' otherwise, the natural logarithm is used
#' throughout the analysis.
logcounts <- function(Y,N){
  t(log2(t(Y+0.5)/(N+1)*1e+06))
}

wgk <- function(Ygk,Nk){
  (1-Ygk/Nk)/Ygk
}

#' This normalization method is similar to normTMM, but it uses a reference sample
#' instead of using the geometric mean as reference.
normTMMr <- function(ercc_reads, N, nref){

  trimmed_ercc <- subset(ercc_reads, apply(ercc_reads,1,function(x) all(x!=0)))

  logcpm <- logcounts(trimmed_ercc, N)

  weights <- matrix(nrow=nrow(trimmed_ercc), ncol=ncol(trimmed_ercc))
  for(i in 1:ncol(trimmed_ercc)){
    weights[,i] <- wgk(trimmed_ercc[,i],N[i]) + wgk(trimmed_ercc[,nref],N[nref])
  }

  logFC <- logcpm - logcpm[,nref]

  logTMM <- colSums(logFC/weights)/colSums(1/weights)
  2^logTMM/2^mean(logTMM)
}

#' @export
normTMM <- function(ercc_reads, N, trimM=0.3, trimA=0.05){

  trimmed_ercc <- subset(ercc_reads, apply(ercc_reads,1,function(x) all(x!=0)))

  # calculate logcpm values
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

#' @export
normTMM_nonzero <- function(ercc_reads, N, trimM=0.3, trimA=0.05){
  
  trimmed_ercc <- ercc_reads
  
  # calculate logcpm values
  logcpm <- logcounts(trimmed_ercc, N)
  logcpm[trimmed_ercc<=0] <- NA
  Nref <- mean(N)
  
  logref <- rowMeans(logcpm, na.rm=T)
  ref_reads <- 2^logref*Nref * 1e-06
  
  # M (logFC) values
  # remove genes with largest M values
  logFC <- logcpm - logref
  logFC[ercc_reads<=0] <- NA
  mean_logFC <- rowMeans(logFC, na.rm=T)
  n_row <- nrow(logFC)
  n_keep_M <- ceiling(n_row*(1-trimM))
  trim_logFC <- order(abs(mean_logFC))[1:n_keep_M]
  trimmed_ercc <- trimmed_ercc[trim_logFC, ]
  
  # re-calculate logcpm values
  logcpm <- logcounts(trimmed_ercc, N)
  logcpm[trimmed_ercc<=0] <- NA
  Nref <- mean(N)
  
  logref <- rowMeans(logcpm, na.rm=T)
  ref_reads <- 2^logref*Nref * 1e-06
  
  # A (logA) values
  # remove genes with largest A values (intensity) to avoid outliers
  logA <- (logcpm + logref)/2
  logA[trimmed_ercc<=0] <- NA
  mean_logA <- rowMeans(logA, na.rm=T)
  n_row <- nrow(logA)
  n_keep_A <- ceiling(n_row*(1-trimA))
  trim_logA <- order(abs(mean_logA))[1:n_keep_A]
  trimmed_ercc <- trimmed_ercc[trim_logA, ]
  
  # re-calculate logcpm values
  logcpm <- logcounts(trimmed_ercc, N)
  logcpm[trimmed_ercc<=0] <- NA
  Nref <- mean(N, na.rm=T)
  
  logref <- rowMeans(logcpm, na.rm=T)
  ref_reads <- 2^logref*Nref * 1e-06
  
  # calculate weights and logFC
  logFC <- logcpm - logref
  logFC[trimmed_ercc<=0] <- NA
  weights <- matrix(nrow=nrow(trimmed_ercc), ncol=ncol(trimmed_ercc))
  for(i in 1:ncol(trimmed_ercc)){
    weights[,i] <- wgk(trimmed_ercc[,i],N[i]) + wgk(ref_reads,Nref)
  }
  weights[trimmed_ercc<=0] <- NA
  # calculate and return TMM values
  logTMM <- colSums(logFC/weights, na.rm=T)/colSums(1/weights, na.rm=T)
  2^logTMM/2^mean(logTMM, na.rm=T)
}

#' This function is similar to normTMM, but it uses the natural logarithm instead of log2.
#' In practice, the results for the normalization factors hardly depend on the choice of basis.
normEXP <- function(ercc_reads, N){

  trimmed_ercc <- subset(ercc_reads, apply(ercc_reads,1,function(x) all(x!=0)))

  logcpm <- log(cpm(trimmed_ercc, N))
  Nref <- mean(N)

  logref <- rowMeans(logcpm)
  ref_reads <- exp(logref)*Nref * 1e-06

  weights <- matrix(nrow=nrow(trimmed_ercc), ncol=ncol(trimmed_ercc))
  for(i in 1:ncol(trimmed_ercc)){
    weights[,i] <- wgk(trimmed_ercc[,i],N[i]) + wgk(ref_reads,Nref)
  }
  print(1/weights[3,])

  logFC <- logcpm - logref

  logTMM <- colSums(logFC/weights)/colSums(1/weights)
  print(wgk(ref_reads,Nref))
  exp(logTMM)/exp(mean(logTMM))
}

#' remove genes/ercc spike-ins with very low/no counts, return normTMM
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

#' remove genes/ercc spike-ins with very low/no counts, return normTMM not using zero counts
#' @export
get_norm_facs_nonzero <- function(rif_reads, ercc_reads, count_min=10, sample_min=3) {
  N <- colSums(rif_reads)
  keep <- which(rowSums(cpm(rif_reads, N) > count_min) >= sample_min)
  rif_reads <- rif_reads[keep,]
  keep <- which(rowSums(cpm(ercc_reads, N) > count_min) >= sample_min)
  ercc_reads <- ercc_reads[keep,]
  N <- colSums(rif_reads)
  normTMM_nonzero(ercc_reads, N)
}

#' plot library size
#' @export
plot_lib_size <- function(plot_name, dfit, adj=1, labels=c("WT", expression(paste(Delta, "proQ")), "proQ+")) {
  dplot <- dfit %>% group_by(batch, condition, time) %>% summarise(N=sum(raw), n_f=mean(nf), n_f2=mean(nf2))
  dplot$sample <- 1:nrow(dplot)
  if (adj==1) dplot$N <- dplot$N*dplot$n_f
  else if (adj==2) dplot$N <- dplot$N*dplot$n_f2
  p <- ggplot(data=dplot, aes(x=sample, y=N, fill=condition)) +
    geom_bar(stat="identity") +
    scale_fill_manual("", values=cb_palette, labels=labels) +
    labs(x="batch", y = "N") +
    scale_y_continuous("",label=fancy_scientific)  +
    scale_x_discrete() +
    annotation_custom(grid::textGrob("1", vjust=1.5, gp = grid::gpar(fontsize = 20)), xmin=7.5, xmax=7.5, ymin=0, ymax=0) +
    annotation_custom(grid::textGrob("2", vjust=1.5, gp = grid::gpar(fontsize = 20)), xmin=22.5, xmax=22.5, ymin=0, ymax=0) +
    annotation_custom(grid::textGrob("3", vjust=1.5, gp = grid::gpar(fontsize = 20)), xmin=37.5, xmax=37.5, ymin=0, ymax=0) +
    annotation_custom(grid::textGrob("4", vjust=1.5, gp = grid::gpar(fontsize = 20)), xmin=52.5, xmax=52.5, ymin=0, ymax=0) +
    annotation_custom(grid::textGrob("5", vjust=1.5, gp = grid::gpar(fontsize = 20)), xmin=67.5, xmax=67.5, ymin=0, ymax=0) +
    annotation_custom(grid::textGrob("6", vjust=1.5, gp = grid::gpar(fontsize = 20)), xmin=82.5, xmax=82.5, ymin=0, ymax=0) +
    theme_classic() +
    theme(axis.text=element_text(size=25),
          axis.title=element_text(size=29),
          plot.title=element_text(size=32,face="bold"),
          legend.text=element_text(size=25),
          legend.title=element_text(size=29,face="bold"),
          legend.position = c(0.8,0.9))
  print(p)
  cairo_pdf(plot_name, width=12, height=6.5)
  print(p)
  dev.off()
}

#' principle component analysis
#' @export
plot_pca <- function(plot_name, dfit, cond, adj=T) {
  if (!adj) dfit$logcpm <- log((dfit$raw + 0.5) / (dfit$N + 1) * 1e6)
  variation <- dfit %>% group_by(gene) %>% summarise(avg=mean(logcpm), var=sd(logcpm))
  dplot <- reshape(dfit[c("logcpm", "sample", "gene")],
                   v.names="logcpm", timevar="sample", idvar="gene", direction="wide")
  row.names(dplot) <- dplot$gene
  dplot <- dplot[, 2:ncol(dplot)]
  colnames(dplot) <- sub("logcpm.", "", colnames(dplot))
  dplot <- dplot[order( - variation$var ), ]
  
  rif_pca <- prcomp(dplot[1:1000,],
                    center = TRUE,
                    scale. = TRUE)
  
  plot(rif_pca, type="l")
  
  cond <- cbind(cond, data.frame(rif_pca$rotation)[match(row.names(cond), row.names(rif_pca$rotation))])
  
  p <- ggplot(cond, aes(x=PC1, y=PC2, col=factor(time))) +
    geom_point(aes(pch=condition), size=2) +
    scale_color_manual("time", values=cb_palette, labels=c(paste(unique(cond$time), " min"))) +
    scale_shape_discrete("condition", labels=c("WT", expression(paste(Delta, "proQ")), "proQ+")) +
    theme_classic() +
    theme(axis.text=element_text(size=25),
          axis.title=element_text(size=29),
          plot.title=element_text(size=32,face="bold"),
          legend.text=element_text(size=25),
          legend.title=element_text(size=29,face="bold"))
  print(p)
  cairo_pdf(plot_name, width = 9, height = 6.5)
  print(p)
  dev.off()
}

#' @export
plot_MA <- function(dfit, dfit_ercc, sample_ids, adj=T) {
  if (!adj) {
    dfit$logcpm <- log((dfit$raw + 0.5) / (dfit$N + 1) * 1e6)
    dfit_ercc$logcpm <- log((dfit_ercc$raw + 0.5) / (dfit_ercc$N + 1) * 1e6)
  }
  
  dfit$type <- "gene"
  dfit_ercc$type <- "ercc"
  sample1 <- rbind(subset(dfit, sample==sample_ids[1]),
                   subset(dfit_ercc, sample==sample_ids[1]))
  sample2 <- rbind(subset(dfit, sample==sample_ids[2]),
                   subset(dfit_ercc, sample==sample_ids[2]))
  
  dplot <- sample2
  dplot$M <- dplot$logcpm - sample1$logcpm
  dplot$A <- (dplot$logcpm + sample1$logcpm)/2
  dplot$type <- factor(dplot$type, levels=c("gene", "ercc"))
  
  ercc_avg <- mean(subset(dplot, type=="ercc")$M)
  
  p <- ggplot(dplot) +
    rasterise(geom_point(aes(x=A, y=M, col=type, size=type)), dpi=250) +
    scale_color_manual("", values=cb_palette[c(2,6)]) +
    scale_size_manual("", values=c(2, 2.5)) +
    geom_hline(yintercept=ercc_avg, size=1, col=cb_palette[6]) +
    theme_classic() +
    theme(axis.text=element_text(size=25),
          axis.title=element_text(size=29),
          plot.title=element_text(size=32,face="bold"),
          legend.text=element_text(size=25),
          legend.title=element_text(size=29,face="bold"))
  return(p)
}

#' @export
calc_batch_effect <- function(dfit) {
  dfit_stats <- dfit %>% group_by(gene, time, condition) %>% summarise(logcpm_avg=mean(logcpm))
  dfit$logrel2 <- dfit$logcpm - dfit_stats$logcpm_avg[match(
    paste(dfit$gene, dfit$time, dfit$condition, sep=":"),
    paste(dfit_stats$gene, dfit_stats$time, dfit_stats$condition, sep=":"))]
  dfit_stats2 <- dfit %>% group_by(sample) %>% summarise(nf2=exp(mean(logrel2)), nf1=mean(nf))
  
  dfit$nf2 <- dfit$nf*dfit_stats2$nf2[match(dfit$sample, dfit_stats2$sample)]
  dfit$logrel3 <- dfit$logrel - log(dfit_stats2$nf2[match(dfit$sample, dfit_stats2$sample)])
  
  return(dfit)
}
