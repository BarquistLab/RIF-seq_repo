# RIF-seq_repo

This repository contains Stan models for the analysis of RNA decay in bacteria which have been treated with the transcription initiation inhibitor rifampicin followed by sequencing (RIF-seq) which were used to extract RNA decay parameters in [^f1]. It was shown that models based on log-counts can be used to describe RNA-seq data when the mean-variance relationship is modeled [^f2]. Because of computational efficiency, all models in this repository are based on log-counts rather than raw counts. The Stan models are in the directory *stan_models*. The directory *R* contains some R functions which may be helpful for the downstream analysis, but which are not required for running the Stan models. An example on how to run the Stan model with cmdstanr can be found the *examples* directory. The required data files are provided in *data*. Furthermore, a log-normal model (LNM) with and without modeling the mean-variance relationship are compared in *model_comparison.Rmd* including using Baysian model comparison techniques like leave-one-out cross validation (LOO-CV)[^f3] and posterior predictive checks [^f4]. Alternatively, the data can be exported using ``stan_rdump``` as described below and run with cmdstan, the command line version of Stan (or with PyStan).

## Stan models

Stan is a probabilistic programming language. Extensive documentation on how to develop and run Stan models can be found on the Stan website: https://mc-stan.org/ . The implementation of Stan models is very efficient, because Stan separates the model implementation from the sampling process.

The directory *stan_models* contains Stan model files to analyze RIF-seq data:

<div class="columns-2">
  
  - **LNM-1.0.stan** This Stan file contains a log-normal model (LNM) which fits RNA decay curves of rifampicin treated bacteria. Its parameters have been described in [^f1]. It is compatible with Stan (>=2.26). Since it is based on normalized log-counts, it is suitable for sequencing experiments with similar library sizes. If the library sizes are very different, we recommend using LNMcdv-1.0.stan, which models the variance-mean relationship.
  - **LNMcdv-1.0.stan** In this Stan model, the variance depends on the raw sequencing counts. It is therefore suitable for sequencing experiments with large differences in library size, while it is still computationally more efficient than a count-based (negative binomial) model. In *model_comparison.Rmd*, the LNM with and without mean-variance modeling are compared.
  - **LNM_loo-1.0.stan** Same as LNM-1.0.stan, but the quantities requried for LOO-CV and posterior predictive checks are computed in the generated quantities block. This makes the model computationally less efficient. Therefore, it should only be used for model comparison.
  - **LNMcdv_loo-1.0.stan** Same as LNMcdv-1.0.stan, but the quantities requried for LOO-CV and posterior predictive checks are computed in the generated quantities block. This makes the model computationally less efficient. Therefore, it should only be used for model comparison.
  - **LNM.stan** This is the original model used to extract RNA decay rates in [^f1]. It is compatible and tested with Stan 2.28 - 2.31. Starting from Stan 2.32, the old syntax of array definitions will no longer work https://mc-stan.org/docs/2_29/reference-manual/brackets-array-syntax.html . We recommend using LNM-1.0.stan instead, where the array definitions have been updated to the new syntax (Stan >=2.26).
  - **LNMsim.stan** simulates RNA decay curves for multiple experimental conditions or bacterial strains (Stan <=2.31).

</div>
If you have questions on how to adapt the statistical models to your sequencing data or if you are interested in using count-based models, please contact the authors directly.

Here, we provide useful links regarding installing and running Stan models, in addition to instructions on how to prepare the data with **R** before running LNM.stan with cmdstan. For an example with cmdstanr, we have included *example_cmdstanr.Rmd* in the *examples* directory.

## Installing and running the Stan model LNM.stan

The Stan models can be run using various interfaces. Depending on the size of the analyzed dataset, fitting the parameters can take several days even with multithreading. It is therefore recommended to use **cmdstan** (example workflow in *example_cmdstan.Rmd*): https://mc-stan.org/docs/2_31/cmdstan-guide/cmdstan-installation.html

Alternatively, smaller (test) data sets can be analyzed using **cmdstanr**. This requires the installation of the the R package cmdstanr **and** the installation of cmdstan. An installation guide can be found here (example workflow in *example_cmdstanr.Rmd*):
https://mc-stan.org/cmdstanr/
and explanations on how to install cmdstan and get started with cmdstanr are provided here:
https://mc-stan.org/cmdstanr/articles/cmdstanr.html

> **Note**
> 
> The Stan models use the function ```reduce_sum```. To reduce computation time, multithreading should be enabled following the guide https://mc-stan.org/docs/2_31/cmdstan-guide/parallelization.html .

After installing cmdstan, a directory for the Stan model can be created within the cmdstan directory and the Stan model can be installed following the cmdstan instructions: https://mc-stan.org/docs/2_31/cmdstan-guide/compiling-a-stan-program.html .

After creating a RIF-seq_data.R file, the decay rates can be fitted with ``` ./LNM sample data file=RIF-seq_data.R``` following the cmdstan instructions: https://mc-stan.org/docs/2_31/cmdstan-guide/mcmc-intro.html. The model was tested using 1000 warmup iterations and 1000 sampling iterations.

## Required R packages for workflow with cmdstan/cmdstanr

<div class="columns-2">
  
  - **CRAN** knitr, data.table, tidyverse, ggplot2, bayesplot
  - **github** cmdstanr (only required for workflow with cmdstanr): https://mc-stan.org/cmdstanr/
</div>

Running the model comparison script requires installing additional packages (see the file "DESCRIPTION" for a list).

## Exporting the read counts and metadata to cmdstan format with R

The read counts and metadata can be exported for cmdstan using the rstan function ```stan_rdump```. The following quantities are required:

<div class="columns-2">
  
  - ```N_con``` Number of experimental conditions or bacterial strains.
  - ```N_s``` Number of samples.
  - ```N_b``` An array containing the number of replicates for every condition.
  - ```N_t``` Number of time points.
  - ```N_tot``` Total number of data points (number of rows of the long data frame)
  - ```N_g``` Number of genetic features.
</div>

Then, the read counts and metadata should be organized in a long data frame ```stan_df``` with the following columns:

<div class="columns-2">
  
  - ```raw``` raw counts, numeric.
  - ```time``` time (in minutes) which corresponds to the sample index ```it_s```, numeric.
  - ```sample``` unique sample ID, character.
  - ```condition``` label of condition, i.e. WT, dRBP, character.
  - ```locus_tag``` locus tag of genetic feature, character.
  - ```it_b``` replicate index, runs from 1,...,```N_b[it_c]```, numeric.
</div>

Furthermore, ```N_s``` normalization factors ```n_f``` (e.g. calculated from ERCC spike-ins) and RNA-seq library sizes ```N_lib``` (total number of reads per sample) are required.

The data is then exported via ```stan_rdump``` (see *example_cmdstan.Rmd* for a working example):

```
# Define iterators required by the Stan models
# condition, sample and locus_tag are expected to be factors, replicate and time numeric.
it_c = stan_df$condition %>% as.numeric
it_b = stan_df$replicate
it_s = stan_df$sample %>% as.numeric
it_g = stan_df$locus_tag %>% as.numeric
it_t = stan_df$time %>% factor(levels=c("0", "3", "6", "12", "24")) %>% as.numeric

# number of genetic features
N_g = unique(stan_df$locus_tag) %>% length
# combined iterators
it_gc = (it_c - 1)*N_g + it_g
N_t = unique(stan_df$time) %>% length

# raw read counts, will be converted to log-counts in the Stan model
raw = stan_df$raw
t = stan_df$time

it_gct = (it_gc - 1)*N_t + it_t
it_gt = (it_g - 1)*N_t + it_t

# number of conditions/strains, including the reference condition/strain
N_con = max(it_c)
# number of samples/sequencing libraries
N_s = max(it_s)
# number of replicates per condition, same order as in it_c
N_b = stan_df %>% group_by(condition) %>% summarize(N_b=unique(it_b) %>% length) %>% pull(N_b)
N_t = N_t
N_g = N_g
# total number of data points
N_tot = nrow(stan_df)

# TMM normalization factors and library sizes, same order as in it_s
n_f = nf_proQ
N_lib = colSums(proQ_reads)
```

Two switches have to be set before exporting the data via ```stan_rdump```:

```
# Remove batch effects by using center mean normalization (yes=1, no=0)
batch_effects <- 1
# Choose model:
# 1: LNM
# 2: LM
# 3: PLM
model_id <- 1

rstan::stan_rdump(c("N_con", "N_s", "N_b", "N_t", "N_g", "N_tot"), file="RIF-seq_data.R")
rstan::stan_rdump(c("raw", "t", "N_lib", "n_f", "batch_effects", "model_id"), file="RIF-seq_data.R", append=T)
rstan::stan_rdump(c("it_g", "it_gc", "it_gct", "it_gt", "it_t", "it_s", "it_b", "it_c"), file="RIF-seq_data.R", append=T)
```

As explained in [^f1], we recommend using the log-normal model for the analysis of RIF-seq data, because it accounts for the delay because of ongoing transcription as well as for a finite baseline concentration of stable RNA which can be observed in RIF-seq data.

## Importing the cmdstan results with R

cmdstan saves the results from the HMC sampler in a csv file, e.g. RIF-seq_results.csv. Since the csv files can be quite large, it is recommended to import them using ```fread```, e.g.

```
stan_samples <- fread(cmd = 'grep -v "#" RIF-seq_results.csv', data.table = F)
```

The WT (control) decay rates of gene ```it_g``` are saved in ```betas.1.it_g``` (LNM) or ```beta_WT.it_g``` (LNM-1.0, LNMcdv-1.0), the differences in decay rate for condition/strain ```it_c``` are saved in ```betas.it_c.it_g``` (LNM) or ```delta_beta.(it_c-1).it_g``` (LNM-1.0, LNMcdv-1.0). Then, the medians and quantiles can be obtained from the data frame ```stan_samples``` as follows:

```
beta_WT <- stan_samples[,grep("^betas\\.1\\.", colnames(stan_samples)] %>% apply(2, function(x) quantile(x, probs=c(0.05, 0.5, 0.95))
```

> **Note**
> 
> If the files are too large to import them with fread, the columns containing the relevant parameters can be extracted from the csv file using command line tools.

## Simulating decay curves

The RIF-seq_data.R file can also be used to simulate decay curves with the same number of replicates/samples/conditions/genes. Alternatively, one can create a RIF-seq_sim_data.R file containing only the required information following the same steps as for the RIF-seq_data.R file, omitting the variables ```raw, N_lib, n_f, batch_effects, model_id```.

```
rstan::stan_rdump(c("N_con", "N_s", "N_b", "N_t", "N_g", "N_tot"), file=RIF-seq_data.R)
rstan::stan_rdump(c("t"), file=RIF-seq_data.R, append=T)
rstan::stan_rdump(c("it_g", "it_gc", "it_gct", "it_gt", "it_t", "it_s", "it_b", "it_c"), file=RIF-seq_data.R, append=T)
```

The model simulates normalized relative log-counts (i.e. the mean log-count at t=0 min is set to zero). These can be converted to raw counts using library sizes and normalization factors.

[^f1]: L. Jenniches, C. Michaux, S. Reichardt, J. Vogel, A. J. Westermann, L. Barquist. Improved RNA stability estimation through Bayesian modeling reveals most bacterial transcripts have sub-minute half-lives. bioRxiv 2023.06.15.545072; doi: https://doi.org/10.1101/2023.06.15.545072
[^f2]: Law, C.W., Chen, Y., Shi, W. et al. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol 15, R29 (2014). https://doi.org/10.1186/gb-2014-15-2-r29
[^f3]: A. Vehtari, A. Gelman, J. Gabry, Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. Stat.
1146 Comput. 27, 1413–1432 (2017).
[^f4]: http://avehtari.github.io/BDA_R_demos/demos_rstan/ppc/poisson-ppc.html
