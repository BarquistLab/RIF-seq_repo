#' Gene-wise decay plots: fitted curve + data points
#' @export
plot_gene <- function(data, results, lab, res_dir, addNorm=F, addDect=F, rm.cond=NULL, show.int=T, col=viridis(4)){
  
  if(any(grepl("dproQ", data$condition))) con_lab <- c("WT", expression(paste(Delta,"proQ")), "proQ+", expression(paste(Delta,"cspCE")))
  else if(any(grepl("DeltaCE", data$condition))) con_lab <-  c("WT", expression(paste(Delta,"cspCE")))
  else if(any(grepl("drppH", data$condition))) con_lab <-  c("WT", expression(paste(Delta,"drppH")))
  else con_lab <- unique(data$it_c)
  
  if(addNorm) data$logrel <- data$logrel3
  
  if (!is.null(rm.cond)) {
    data <- subset(data, it_c!=rm.cond)
    results <- subset(results, it_c!=rm.cond)
    con_lab <- con_lab[-rm.cond]
  }
  
  # zero counts are set to 0.5, i.e. we should mark the detection limit of log(0.5/10) - logcpm0
  data_stats <- subset(data, time==0) %>% group_by(it_c) %>% summarise(logcpm0=mean(logcpm))
  data_stats$detect_lim <- log(0.05) - data_stats$logcpm0
    
  p <- ggplot(data) +
    geom_point(aes(x=time, y=logrel, col=as.factor(it_c), shape=as.factor(it_c)), size=3) +
    geom_line(aes(x=time, y=mu, col=as.factor(it_c), group=as.factor(it_c), linetype=as.factor(it_c)),
              inherit.aes=F, data=results, size=1.5) +
    scale_color_manual(name="condition", labels=con_lab, values=col) +
    scale_fill_manual(name="condition", labels=con_lab, values=col) +
    scale_linetype_discrete(name="condition", labels=con_lab) +
    scale_shape_discrete(name="condition", labels=con_lab) +
    xlab("t [min]") +
    ylab("rel. log-abundance") +
    theme_classic() +
    theme(axis.text=element_text(size=18),
          axis.title=element_text(size=21),
          plot.title=element_text(size=21, margin = margin(t = 10, b = -20)),
          legend.text=element_text(size=18),
          legend.title=element_blank(),
          legend.position="top",
          legend.key.width=unit(0.05, "cm")) +
    annotate("text", x=Inf, y=Inf, label=unique(data$com_id),
             hjust="inward", vjust="inward", size=9)
  
  if (addDect) {
    p <- p +
      geom_hline(aes(yintercept=detect_lim, col=as.factor(it_c)), data=data_stats, linetype="dashed")
  }
  
  if (show.int) {
    p <- p +
      geom_ribbon(aes(x=time, ymin=mu_m, ymax=mu_p, group=as.factor(it_c), fill=as.factor(it_c)),
                  inherit.aes=F, alpha=0.2, data=results)
  }
  print(p)
  cairo_pdf(paste(res_dir,"decay_rate_",lab,".pdf",sep = ""), width = 4, height = 3.5)
  print(p)
  dev.off()
}

#' calculate decay curve for a set of parameters
#' @export
calc_mu <- function(pib, beta, gamma, time, model_id) {
  
  mu <- 0
  if (grepl("LNM", model_id)) {
    mu <- log( pib + (1 - pib) * exp(-beta*(time - gamma)) )
  } else if (model_id=="NB") {
    mu <- log( pib + (1 - pib) * exp(-beta*(time - gamma)) )
  } else if (model_id=="PLM") {
    mu <- -beta*(time - gamma)
  } else {
    mu <- -beta*time
  }
  return(mu)
}

#' gene-wise plots of (exponential) decay
#' @export
plot_genes <- function(res_dir, dfit, param, id_RNA=NULL, model_id="", addNorm=F, addDect=F, rm.cond=NULL, show.int=T, col=viridis(4)){
  
  if(is.null(id_RNA)) {
    id_RNA <- unique(dfit$locus_tag)
  }
  else {
    id_RNA <- id_RNA[which(id_RNA %in% unique(dfit$locus_tag))]
  }
  
  dfit$it_c = dfit$condition %>% as.numeric
  dfit$it_b = dfit$replicate
  dfit$it_s = dfit$sample %>% as.numeric
  dfit$it_g = dfit$locus_tag %>% as.numeric
  dfit$it_t = dfit$time %>% as.numeric
  
  N_g <- nrow(param)
  N_con <- length(unique(dfit$it_c))
  t_max <- max(dfit$time)
  
  beta_list <- vector("list", length=N_con)
  beta_list[[1]] <- data.frame(beta05=param$beta05_WT, beta=param$beta_WT, beta95=param$beta95_WT)
  
  pib_list <- vector("list", length=N_con)
  pib_list[[1]] <- data.frame(pib05=param$pib05_WT, pib=param$pib_WT, pib95=param$pib95_WT)
  
  ic <- 1
  while (ic < N_con) {
    betas <- paste(c("beta05_", "beta_", "beta95_"), ic, sep="")
    pibs <- paste(c("pib05_", "pib_", "pib95_"), ic, sep="")
    
    beta_list[[ic+1]] <- param %>% select(betas)
    pib_list[[ic+1]] <- param %>% select(pibs)
    
    ic <- ic + 1
  }
  all_results <- vector("list", length=length(id_RNA))
  for (i in 1:length(id_RNA)) {
    id <- which(param$locus_tag==id_RNA[i])
    results <- data.frame(time = numeric(),
                          it_c = numeric(),
                          beta = double(),
                          beta_m = double(),
                          beta_p = double(),
                          gamma = double(),
                          pib = double(),
                          mu=double(), mu_m=double(), mu_p=double())
    for (j in 1:N_con) {
      results_con <- data.frame(time = (1:(t_max*10)/10),
                                it_c = j,
                                beta = beta_list[[j]][id, 2] %>% unlist(),
                                beta_m = beta_list[[j]][id, 1] %>% unlist(),
                                beta_p = beta_list[[j]][id, 3] %>% unlist(),
                                gamma = param$gamma[id],
                                pib = pib_list[[j]][id, 2] %>% unlist(),
                                mu=0, mu_m=0, mu_p=0)
      if(model_id=="LM") results_con$gamma <- 0
      decay <- which(results_con$time>=results_con$gamma)
      results_con$mu[decay] <-
        do.call(function(pib, beta, gamma, time, ...) calc_mu(pib, beta, gamma, time, model_id=model_id), subset(results_con, time>gamma))
      results_con$mu_m[decay] <-
        do.call(function(pib, beta_m, gamma, time, ...) calc_mu(pib, beta_m, gamma, time, model_id=model_id), subset(results_con, time>gamma))
      results_con$mu_p[decay] <-
        do.call(function(pib, beta_p, gamma, time, ...) calc_mu(pib, beta_p, gamma, time, model_id=model_id), subset(results_con, time>gamma))
      results <- rbind(results, results_con)
    }
    all_results[[i]] <- results
  }
  
  for (i in 1:length(id_RNA)) {
    results <- all_results[[i]]
    dplot <- subset(dfit, locus_tag==id_RNA[i])
    if (is.null(dplot$com_id)) {
      lab <- paste(model_id, unique(dplot$locus_tag), sep="_")
    } else {
      lab <- paste(model_id, unique(dplot$com_id), sep="_")
    }
    plot_gene(dplot, results, lab, res_dir, addNorm=addNorm, addDect=addDect, rm.cond=rm.cond, show.int=show.int, col=col)
  }
}
