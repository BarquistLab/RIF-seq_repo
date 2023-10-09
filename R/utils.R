#' custum colors
cb_palette_8 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cb_palette <- c("#193cbcff", "#1473afff", "#589acfff", "#89c3efff", "#ea594eff", "#e5b039ff", "#ede65aff")
cb_palette_4 <- c(viridis(4)[1:3], cb_palette[6])

#' calculate counts per million (cpm)
#' @export
cpm <- function(Y,N){
  t(t(Y+0.5)/(N+1)*1e+06)
}

#' convert decay rate to half-lifes
#' @export
dr_to_hl <- function(beta) {
  log(2)/beta
}

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

#' @export
nzmean <- function(x) {
  zvals <- x==0
  if (all(zvals)) 0 else mean(x[!zvals])
}

#' @export
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

#' Fraction of samples that fall into the interval [-int0, int0], which corresponds to the null hypothesis H0.
#' The p value corresponds to the probability of the null hypothesis.
#' @export
p_value <- function(x, int0=0, h0=0) {
  x_med <- median(x)
  ifelse(x_med>h0, length(which(x<int0+h0))/length(x), length(which(x>-int0+h0))/length(x))
}

#' @export
p_value_sample <- function(x, int0=0, h0=0, n_x=2^18, n_sam=1e6) {
  x_med <- median(x)
  p <- 1
  if(x_med>h0){
    den <- density(x, n=n_x)
    samples <- sample(den$x, n_sam, replace=T, prob=den$y)
    p <- length(which(samples<int0+h0))/n_sam
  } else {
    den <- density(x, n=n_x)
    samples <- sample(den$x, n_sam, replace=T, prob=den$y)
    p <- length(which(samples>-int0+h0))/n_sam
  }
  return(p)
}

#' @export
add_FDR <- function(param){
  
  N_con <- ncol(param %>% select(starts_with("p_value")))
  
  for (i in 1:N_con) {
    param <- param %>% arrange(param %>% select(paste0("p_value_", i)))
    FDR <- cumsum(param %>% pull(paste0("p_value_", i)))/(1:nrow(param))
    param <- cbind(param, FDR=FDR)
    names(param)[which(names(param)=="FDR")] <- paste0("FDR_",i)
  }
  
  param <- param %>% arrange(FDR_1)
  return(param)
}

#' @export
n_c <- function(x) {
  if(x>1) return(x+2)
  else return(x)
}

#' @export
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}
