#' Extract results from cmdstanr (LNM-1.0) and return data.frame with the parameters which are relevant for downstream analysis
#' @export
extract_params <- function(param_df, sig=5e-2, nThread=1, delta_thr=2e-2, pi_thr=1e-4, hl_thr=5e-2, hl_thr2=0.1) {

  sigm1 <- 1 - sig
  	
  beta_WT <- param_df[,grep("beta_WT", colnames(param_df))]
  delta_beta <- param_df[,grep("delta_beta[\\.|\\[]", colnames(param_df))]
	gamma_fit <- param_df[,grep("^gamma", colnames(param_df))] %>% apply(2,function(x) quantile(x, probs=c(sig,0.5,sigm1)))
	# required for 68% interval
	gamma_fit_sd <- param_df[,grep("^gamma", colnames(param_df))] %>% apply(2,function(x) quantile(x, probs=c(0.16,0.5,0.84)))
	sigma_fit <- param_df[,grep("^sigma_g", colnames(param_df))] %>% apply(2,function(x) quantile(x, probs=c(sig,0.5,sigm1)))
	pib_b <- param_df[,grep("^pib_b", colnames(param_df))]

	N_g <- ncol(gamma_fit)
	N_con <- ncol(delta_beta)/N_g + 1
	N_t <- ncol(sigma_fit)/N_g

	betaWT_fit <- beta_WT %>% apply(2,function(x) quantile(x, probs=c(sig,0.5,sigm1)))
	
	hlWTs <- log(2)/beta_WT
	hlWT_fit <- hlWTs %>% apply(2,function(x) quantile(x, probs=c(sig, 0.16, 0.5, 0.84, sigm1)))
	pibWT_fit <- pib_b[,1:N_g] %>% apply(2,function(x) quantile(x, probs=c(sig,0.5,sigm1)))
	if (N_con>1) {
	  delta_fit <- delta_beta %>% apply(2,function(x) quantile(x, probs=c(sig,0.5,sigm1)))
	  pib_b_fit <- pib_b %>% apply(2,function(x) quantile(x, probs=c(sig,0.5,sigm1)))
	}

	param_fit <- data.frame(beta_WT=betaWT_fit[2,],
	                        beta05_WT=betaWT_fit[1,],
	                        beta95_WT=betaWT_fit[3,],
	                        hl_WT=hlWT_fit[3,],
	                        hl05_WT=hlWT_fit[1,],
	                        hl95_WT=hlWT_fit[5,],
	                        hl16_WT=hlWT_fit[2,],
	                        hl84_WT=hlWT_fit[4,],
	                        gamma=gamma_fit[2,],
	                        gamma05=gamma_fit[1,],
	                        gamma95=gamma_fit[3,],
	                        gamma16=gamma_fit_sd[1,],
	                        gamma84=gamma_fit_sd[3,],
	                        pib_WT=pibWT_fit[2,],
	                        pib05_WT=pibWT_fit[1,],
	                        pib95_WT=pibWT_fit[3,], row.names=1:N_g)
	ic <- 1
	while (ic < N_con) {
	  delta_cond <- delta_fit[, ((ic-1)*N_g+1):(ic*N_g)]
	  beta_cond <- apply( beta_WT + delta_beta[,((ic-1)*N_g+1):(ic*N_g)], 2, function(x) quantile(x, probs=c(sig,0.5,sigm1)))
	  
	  hls <- log(2)/(beta_WT + delta_beta[,((ic-1)*N_g+1):(ic*N_g)])
	  hl_cond <- hls %>% apply(2, function(x) quantile(x, probs=c(sig,0.5,sigm1)))
	  delta_hl_cond <- (hls - hlWTs) %>% apply(2, function(x) quantile(x, probs=c(sig, 0.16, 0.5, 0.84, sigm1)))
	  pib_cond <- pib_b_fit[, ((ic)*N_g+1):((ic+1)*N_g)]
	  
	  # calculate p values
	  p_cond <- delta_beta[,((ic-1)*N_g+1):(ic*N_g)] %>% apply(2, function(x) p_value(x, int0=delta_thr))
	  p_hl <- (hls - hlWTs) %>% apply(2, function(x) p_value(x, int0=hl_thr))
	  p_hl2 <- ((hls - hlWTs)/hlWTs) %>% apply(2, function(x) p_value(x, int0=hl_thr2))
	  
	  ratio_hl_cond <- sapply(1:N_g, function(ig) {
	    ratio <- 0
	    if(delta_hl_cond[3,ig]>0) {
	      ratio <- delta_hl_cond[3,ig] / (delta_hl_cond[5,ig] - delta_hl_cond[3,ig])
	    } else {
	      ratio <- delta_hl_cond[3,ig] / (delta_hl_cond[3,ig] - delta_hl_cond[1,ig])
	    }
	    return(ratio)
	  })

	  data_cond <- data.frame(delta_cond[2,], delta_cond[1,], delta_cond[3,], p_cond,
	                          beta_cond[2,], beta_cond[1,], beta_cond[3,],
	                          hl_cond[2,], hl_cond[1,], hl_cond[3,],
	                          delta_hl_cond[3,], delta_hl_cond[1,], delta_hl_cond[5,], delta_hl_cond[2,], delta_hl_cond[4,], p_hl, p_hl2, ratio_hl_cond,
	                          pib_cond[2,], pib_cond[1,], pib_cond[3,], row.names=1:N_g)
	  colnames(data_cond) <- c(paste("delta_", ic, sep=""), paste("delta05_", ic, sep=""), paste("delta95_", ic, sep=""), paste("p_value_", ic, sep=""),
	                           paste("beta_", ic, sep=""), paste("beta05_", ic, sep=""), paste("beta95_", ic, sep=""),
	                           paste("hl_", ic, sep=""), paste("hl05_", ic, sep=""), paste("hl95_", ic, sep=""), 
	                           paste("delta_hl_", ic, sep=""), paste("delta_hl05_", ic, sep=""), paste("delta_hl95_", ic, sep=""), paste("delta_hl16_", ic, sep=""), paste("delta_hl84_", ic, sep=""), paste0("p_hl_", ic), paste0("p_hl2_", ic), paste0("ratio_hl_", ic),
	                           paste("pib_", ic, sep=""), paste("pib05_", ic, sep=""), paste("pib95_", ic, sep=""))
	  param_fit <- cbind(param_fit, data_cond)
	  ic <- ic+1
	}
	data_sigma <- data.frame(sigma_fit[2, 1:N_g], row.names=1:N_g)
	N_t <- (ncol(sigma_fit)/N_g)
	if(N_t > 1) {
	  for (it in 2:N_t) {
	    data_sigma <- cbind(data_sigma, data.frame(sigma_fit[2, ((it-1)*N_g+1):(it*N_g)]))
	  }
	}
	colnames(data_sigma) <- paste("sigma_g", 1:N_t, sep="")
	param_fit <- cbind(param_fit, data_sigma)
  return(param_fit)
}

#' Extract results from cmdstanr (LNM-1.0 with gamma=0) and return data.frame with the parameters which are relevant for downstream analysis
#' Adapted routine for a statistical model where gamma is 0, the corresponding statistical model is not part of this repository.
#' @export
extract_params_g0 <- function(param_df, sig=5e-2, nThread=1, delta_thr=2e-2, pi_thr=1e-4, hl_thr=5e-2, hl_thr2=0.1) {
  
  sigm1 <- 1 - sig
  
  beta_WT <- param_df[,grep("beta_WT", colnames(param_df))]
  delta_beta <- param_df[,grep("delta_beta", colnames(param_df))]
  # required for 68% interval
  sigma_fit <- param_df[,grep("^sigma_g", colnames(param_df))] %>% apply(2,function(x) quantile(x, probs=c(sig,0.5,sigm1)))
  pib_b <- param_df[,grep("^pib_b", colnames(param_df))]
  
  N_g <- ncol(beta_WT)
  N_con <- ncol(delta_beta)/N_g + 1
  N_t <- ncol(sigma_fit)/N_g
  
  betaWT_fit <- beta_WT %>% apply(2,function(x) quantile(x, probs=c(sig,0.5,sigm1)))
  
  hlWTs <- log(2)/beta_WT
  hlWT_fit <- hlWTs %>% apply(2,function(x) quantile(x, probs=c(sig, 0.16, 0.5, 0.84, sigm1)))
  pibWT_fit <- pib_b[,1:N_g] %>% apply(2,function(x) quantile(x, probs=c(sig,0.5,sigm1)))
  if (N_con>1) {
    delta_fit <- delta_beta %>% apply(2,function(x) quantile(x, probs=c(sig,0.5,sigm1)))
    pib_b_fit <- pib_b %>% apply(2,function(x) quantile(x, probs=c(sig,0.5,sigm1)))
  }
  
  param_fit <- data.frame(beta_WT=betaWT_fit[2,],
                          beta05_WT=betaWT_fit[1,],
                          beta95_WT=betaWT_fit[3,],
                          hl_WT=hlWT_fit[3,],
                          hl05_WT=hlWT_fit[1,],
                          hl95_WT=hlWT_fit[5,],
                          hl16_WT=hlWT_fit[2,],
                          hl84_WT=hlWT_fit[4,],
                          gamma=0,
                          gamma05=0,
                          gamma95=0,
                          gamma16=0,
                          gamma84=0,
                          pib_WT=pibWT_fit[2,],
                          pib05_WT=pibWT_fit[1,],
                          pib95_WT=pibWT_fit[3,], row.names=1:N_g)
  ic <- 1
  while (ic < N_con) {
    delta_cond <- delta_fit[, ((ic-1)*N_g+1):(ic*N_g)]
    beta_cond <- apply( beta_WT + delta_beta[,((ic-1)*N_g+1):(ic*N_g)], 2, function(x) quantile(x, probs=c(sig,0.5,sigm1)))
    
    hls <- log(2)/(beta_WT + delta_beta[,((ic-1)*N_g+1):(ic*N_g)])
    hl_cond <- hls %>% apply(2, function(x) quantile(x, probs=c(sig,0.5,sigm1)))
    delta_hl_cond <- (hls - hlWTs) %>% apply(2, function(x) quantile(x, probs=c(sig, 0.16, 0.5, 0.84, sigm1)))
    pib_cond <- pib_b_fit[, ((ic)*N_g+1):((ic+1)*N_g)]
    
    # calculate p values
    p_cond <- delta_beta[,((ic-1)*N_g+1):(ic*N_g)] %>% apply(2, function(x) p_value(x, int0=delta_thr))
    p_hl <- (hls - hlWTs) %>% apply(2, function(x) p_value(x, int0=hl_thr))
    p_hl2 <- ((hls - hlWTs)/hlWTs) %>% apply(2, function(x) p_value(x, int0=hl_thr2))
    
    ratio_hl_cond <- sapply(1:N_g, function(ig) {
      ratio <- 0
      if(delta_hl_cond[3,ig]>0) {
        ratio <- delta_hl_cond[3,ig] / (delta_hl_cond[5,ig] - delta_hl_cond[3,ig])
      } else {
        ratio <- delta_hl_cond[3,ig] / (delta_hl_cond[3,ig] - delta_hl_cond[1,ig])
      }
      return(ratio)
    })
    
    data_cond <- data.frame(delta_cond[2,], delta_cond[1,], delta_cond[3,], p_cond,
                            beta_cond[2,], beta_cond[1,], beta_cond[3,],
                            hl_cond[2,], hl_cond[1,], hl_cond[3,],
                            delta_hl_cond[3,], delta_hl_cond[1,], delta_hl_cond[5,], delta_hl_cond[2,], delta_hl_cond[4,], p_hl, p_hl2, ratio_hl_cond,
                            pib_cond[2,], pib_cond[1,], pib_cond[3,], row.names=1:N_g)
    colnames(data_cond) <- c(paste("delta_", ic, sep=""), paste("delta05_", ic, sep=""), paste("delta95_", ic, sep=""), paste("p_value_", ic, sep=""),
                             paste("beta_", ic, sep=""), paste("beta05_", ic, sep=""), paste("beta95_", ic, sep=""),
                             paste("hl_", ic, sep=""), paste("hl05_", ic, sep=""), paste("hl95_", ic, sep=""), 
                             paste("delta_hl_", ic, sep=""), paste("delta_hl05_", ic, sep=""), paste("delta_hl95_", ic, sep=""), paste("delta_hl16_", ic, sep=""), paste("delta_hl84_", ic, sep=""), paste0("p_hl_", ic), paste0("p_hl2_", ic), paste0("ratio_hl_", ic),
                             paste("pib_", ic, sep=""), paste("pib05_", ic, sep=""), paste("pib95_", ic, sep=""))
    param_fit <- cbind(param_fit, data_cond)
    ic <- ic+1
  }
  data_sigma <- data.frame(sigma_fit[2, 1:N_g], row.names=1:N_g)
  for (it in 2:N_t) {
    data_sigma <- cbind(data_sigma, data.frame(sigma_fit[2, ((it-1)*N_g+1):(it*N_g)]))
  }
  colnames(data_sigma) <- paste("sigma_g", 1:N_t, sep="")
  param_fit <- cbind(param_fit, data_sigma)
  return(param_fit)
}

#' Extract results from cmdstanr (NB-1.0) and return data.frame with the parameters which are relevant for downstream analysis
#' NB-1.0 is not part of the github repository
#' @export
extract_params_nb <- function(param_df, sig=5e-2, nThread=1, delta_thr=2e-2, pi_thr=1e-4, hl_thr=5e-2, hl_thr2=0.1) {
  
  sigm1 <- 1 - sig
  
  beta_WT <- param_df[,grep("beta_WT", colnames(param_df))]
  delta_beta <- param_df[,grep("delta_beta", colnames(param_df))]
  gamma_fit <- param_df[,grep("^gamma", colnames(param_df))] %>% apply(2,function(x) quantile(x, probs=c(sig,0.5,sigm1)))
  # required for 68% interval
  gamma_fit_sd <- param_df[,grep("^gamma", colnames(param_df))] %>% apply(2,function(x) quantile(x, probs=c(0.16,0.5,0.84)))
  sigma_fit <- param_df[,grep("^inv_rho_g", colnames(param_df))] %>% apply(2,function(x) quantile(x, probs=c(sig,0.5,sigm1)))
  pib_b <- param_df[,grep("^pib_b", colnames(param_df))]
  
  N_g <- ncol(gamma_fit)
  N_con <- ncol(delta_beta)/N_g + 1
  N_t <- ncol(sigma_fit)/N_g
  
  betaWT_fit <- beta_WT %>% apply(2,function(x) quantile(x, probs=c(sig,0.5,sigm1)))
  
  hlWTs <- log(2)/beta_WT
  hlWT_fit <- hlWTs %>% apply(2,function(x) quantile(x, probs=c(sig, 0.16, 0.5, 0.84, sigm1)))
  pibWT_fit <- pib_b[,1:N_g] %>% apply(2,function(x) quantile(x, probs=c(sig,0.5,sigm1)))
  if (N_con>1) {
    delta_fit <- delta_beta %>% apply(2,function(x) quantile(x, probs=c(sig,0.5,sigm1)))
    pib_b_fit <- pib_b %>% apply(2,function(x) quantile(x, probs=c(sig,0.5,sigm1)))
  }
  
  param_fit <- data.frame(beta_WT=betaWT_fit[2,],
                          beta05_WT=betaWT_fit[1,],
                          beta95_WT=betaWT_fit[3,],
                          hl_WT=hlWT_fit[3,],
                          hl05_WT=hlWT_fit[1,],
                          hl95_WT=hlWT_fit[5,],
                          hl16_WT=hlWT_fit[2,],
                          hl84_WT=hlWT_fit[4,],
                          gamma=gamma_fit[2,],
                          gamma05=gamma_fit[1,],
                          gamma95=gamma_fit[3,],
                          gamma16=gamma_fit_sd[1,],
                          gamma84=gamma_fit_sd[3,],
                          pib_WT=pibWT_fit[2,],
                          pib05_WT=pibWT_fit[1,],
                          pib95_WT=pibWT_fit[3,], row.names=1:N_g)
  ic <- 1
  while (ic < N_con) {
    delta_cond <- delta_fit[, ((ic-1)*N_g+1):(ic*N_g)]
    beta_cond <- apply( beta_WT + delta_beta[,((ic-1)*N_g+1):(ic*N_g)], 2, function(x) quantile(x, probs=c(sig,0.5,sigm1)))
    
    hls <- log(2)/(beta_WT + delta_beta[,((ic-1)*N_g+1):(ic*N_g)])
    hl_cond <- hls %>% apply(2, function(x) quantile(x, probs=c(sig,0.5,sigm1)))
    delta_hl_cond <- (hls - hlWTs) %>% apply(2, function(x) quantile(x, probs=c(sig, 0.16, 0.5, 0.84, sigm1)))
    pib_cond <- pib_b_fit[, ((ic)*N_g+1):((ic+1)*N_g)]
    
    # calculate p values
    p_cond <- delta_beta[,((ic-1)*N_g+1):(ic*N_g)] %>% apply(2, function(x) p_value(x, int0=delta_thr))
    p_hl <- (hls - hlWTs) %>% apply(2, function(x) p_value(x, int0=hl_thr))
    p_hl2 <- ((hls - hlWTs)/hlWTs) %>% apply(2, function(x) p_value(x, int0=hl_thr2))
    
    ratio_hl_cond <- sapply(1:N_g, function(ig) {
      ratio <- 0
      if(delta_hl_cond[3,ig]>0) {
        ratio <- delta_hl_cond[3,ig] / (delta_hl_cond[5,ig] - delta_hl_cond[3,ig])
      } else {
        ratio <- delta_hl_cond[3,ig] / (delta_hl_cond[3,ig] - delta_hl_cond[1,ig])
      }
      return(ratio)
    })
    
    data_cond <- data.frame(delta_cond[2,], delta_cond[1,], delta_cond[3,], p_cond,
                            beta_cond[2,], beta_cond[1,], beta_cond[3,],
                            hl_cond[2,], hl_cond[1,], hl_cond[3,],
                            delta_hl_cond[3,], delta_hl_cond[1,], delta_hl_cond[5,], delta_hl_cond[2,], delta_hl_cond[4,], p_hl, p_hl2, ratio_hl_cond,
                            pib_cond[2,], pib_cond[1,], pib_cond[3,], row.names=1:N_g)
    colnames(data_cond) <- c(paste("delta_", ic, sep=""), paste("delta05_", ic, sep=""), paste("delta95_", ic, sep=""), paste("p_value_", ic, sep=""),
                             paste("beta_", ic, sep=""), paste("beta05_", ic, sep=""), paste("beta95_", ic, sep=""),
                             paste("hl_", ic, sep=""), paste("hl05_", ic, sep=""), paste("hl95_", ic, sep=""), 
                             paste("delta_hl_", ic, sep=""), paste("delta_hl05_", ic, sep=""), paste("delta_hl95_", ic, sep=""), paste("delta_hl16_", ic, sep=""), paste("delta_hl84_", ic, sep=""), paste0("p_hl_", ic), paste0("p_hl2_", ic), paste0("ratio_hl_", ic),
                             paste("pib_", ic, sep=""), paste("pib05_", ic, sep=""), paste("pib95_", ic, sep=""))
    param_fit <- cbind(param_fit, data_cond)
    ic <- ic+1
  }
  data_sigma <- data.frame(sigma_fit[2, 1:N_g], row.names=1:N_g)
  for (it in 2:N_t) {
    data_sigma <- cbind(data_sigma, data.frame(sigma_fit[2, ((it-1)*N_g+1):(it*N_g)]))
  }
  colnames(data_sigma) <- paste("inv_rho_g", 1:N_t, sep="")
  param_fit <- cbind(param_fit, data_sigma)
  return(param_fit)
}
