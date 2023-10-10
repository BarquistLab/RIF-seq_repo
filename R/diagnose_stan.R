bayesplot_theme_update(text = element_text(size = 10, family = "sans"))

#'
#' @export
extract_np <- function(mfit_diag, n_iter=1000) {
  
  mfit_diag <- mfit_diag[, grep("__|Chain", colnames(mfit_diag))]
  n_chains <- mfit_diag$Chain %>% unique %>% length
  mfit_diag$Iteration <- rep(1:n_iter, n_chains)
  np <- pivot_longer(mfit_diag, !c("Iteration", "Chain"), names_to="Parameter", values_to="Value")
  np <- np[,c("Chain", "Iteration", "Parameter", "Value")]
  
  np <- subset(np, Parameter != "lp__")
  np$Parameter <- factor(np$Parameter, levels=c("accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__"))
  
  return(np)
}

#'
#' @export
plot_pairs <- function(plot_name, mfit, np, pars, width=8, height=8) {
  
  color_scheme_set("darkgray")
  p <- mcmc_pairs(mfit, np = np, pars=pars,
                  np_style = pairs_style_np(div_color="green", div_alpha=0.8))
  print(p)
  ggsave(plot_name, p, width=width, height=height)
}

#'
#' @export
plot_trace <- function(plot_name, mfit, np, pars, width=8, height=4) {
  
  color_scheme_set("viridis")
  p <- mcmc_trace(mfit, np = np, pars=pars,
                  facet_args = list(ncol = 1, strip.position = "left"))
  print(p)
  ggsave(plot_name, p, width=width, height=height)
}

#'
#' @export
plot_energy <- function(plot_name, np, width=6, height=3) {
  p <- mcmc_nuts_energy(np) + theme_classic()
  print(p)
  cairo_pdf(plot_name, width=width, height=height)
  print(p)
  dev.off()
}
