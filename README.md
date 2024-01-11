# RNA Decay Analysis with Stan

Welcome to the RNA Decay Analysis repository! This collection of Stan models is designed for studying RNA decay in bacteria treated with the transcription initiation inhibitor rifampicin, followed by sequencing (RIF-seq). These models, discussed in detail in [^f1], aim to extract essential RNA decay parameters. Our approach is based on log-counts rather than raw counts, as it has been demonstrated that log-count models effectively describe RNA-seq data when accounting for the mean-variance relationship [^f2].

## Directory Structure:

<div class="columns-2">
  
  - **stan_models:** This directory contains all Stan models.
  - **R:** The R directory houses helpful R functions. While not mandatory for running Stan models, these functions can enhance downstream analysis.
  - **examples:** Discover examples illustrating how to run the Stan models using cmdstanr/cmdstan. Clear instructions and sample data are provided.
  - **data:** Essential data files required for running the example workflows are conveniently stored here.
  - **annotations:** This directory contains the annotations used for the analysis of RNA decay rates in *Salmonella* [^f1].

</div>

# Getting Started

## Installing and running the Stan models

To run the Stan models, various interfaces are available. The choice of interface depends on the dataset size, as fitting parameters might take several days, even with multithreading. For optimal performance, we recommend using the command-line interface of Stan, **cmdstan** (detailed example workflow provided in *example_cmdstan.Rmd*). Alternatives are the R interface cmdstanr and the Python interace PyStan.

### Using cmdstan

1. **Install cmdstan**

Follow the comprehensive guide https://mc-stan.org/docs/2_31/cmdstan-guide/cmdstan-installation.html . The models have been tested with Stan v2.31.0.

2. **Enable multithreading**

The Stan models use the function ```reduce_sum```. To enhance computation speed, enable multithreading following the guide https://mc-stan.org/docs/2_31/cmdstan-guide/parallelization.html .

3. **Compile the Stan model**

After installing cmdstan, create a directory for the Stan model within the cmdstan directory. To compile the model, follow the instructions here: https://mc-stan.org/docs/2_31/cmdstan-guide/compiling-a-stan-program.html .

4. **Fitting the Stan model**

Below, we explain how to export the sequencing data to the R dump file RIF-seq_data.R. Alternatively, the JSON format can be used. After formatting the data, fit the decay rates using the following command

```
./LNM sample data file=RIF-seq_data.R
```

Refer to the cmdstan instructions: https://mc-stan.org/docs/2_31/cmdstan-guide/mcmc-intro.html. The models have been tested with 1000 warmup iterations and 1000 sampling iterations.

### Using cmdstanr

For smaller (test) datasets, consider using **cmdstanr**. Install the R package cmdstanr **and** cmdstan following the installation guide (example workflow in *example_cmdstanr.Rmd*):
https://mc-stan.org/cmdstanr/
Explanations on how to install cmdstan and get started with cmdstanr are provided here:
https://mc-stan.org/cmdstanr/articles/cmdstanr.html



## Example workflows

Get started with Bayesian modeling with our three comprehensive example workflows, each providing a step-by-step example of a Bayesian analysis workflow of RNA decay in rifampicin-treated bacteria.

### Prerequisites

Before running the examples, ensure you have the necessary R packges installed.

<div class="columns-2">
  
  - **CRAN** knitr, data.table, tidyverse, ggplot2, bayesplot
  - **github** cmdstanr (only required for workflow with cmdstanr): https://mc-stan.org/cmdstanr/
</div>

For additional packages required for the model comparison script, consult the "DESCRIPTION" file in the repository.

### List of workflows

<div class="columns-2">
  
  - **example_cmdstan.Rmd** Normalize and prepare RIF-seq read counts for 50 genes, 3 replicates of WT *Salmonella* and a ProQ deletion strain. Export the data to the format required by the command-line version of Stan, **cmdstan**. Fit the Stan model separately using cmdstan. Copy the cmdstan output file to the *data* directory before running Bayesian diagnostics and plotting the fit results.
   - **example_cmdstanr.Rmd** Similar to *example_cmdstanr.Rmd*, this workflow employs the R package cmdstanr. Find instructions on installing cmdstanr [here](https://mc-stan.org/cmdstanr/) and in the example script. If cmdstanr has not been used previously, run ```check_cmdstan_toolchain``` and ```install_cmdstan``` before fitting the Stan model.
   - **model_comparison.Rmd** This example extends the cmdstanr workflow, fitting both the log-normal model (LNM) as well as the LNM with count-dependent variance (LNMcdv). Use the Stan model files which calculate the pointwise log-likelihood ```log_lik``` and saves a random draw from the log-normal distribution ```y_rep```. These quantities are essential for the Bayesian model comparison techniques leave-one-out cross validation (LOO-CV)[^f3] and posterior predictive checks[^f4]. After fitting the models, explore model comparisons with LOO-CV, posterior predictive checking, and correlation coefficients. Observe how the width of the log-normal distribution (the variation unexplained by the model) behaves across different models. The initial run requires fitting the models and saving them as RDS files for efficient reuse. If you are interested in comparing cpu times, comment in and run the respective code.
  
 </div>

## Stan models

Stan is a probabilistic programming language. Find extensive documentation on how to develop and run Stan models on the Stan website: https://mc-stan.org/ . The implementation of Stan models is very efficient, because Stan separates model implementation from the sampling process.

The directory *stan_models* contains Stan model files to analyze RIF-seq data:

<div class="columns-2">
  
  - **LNM-1.0.stan** A log-normal model (LNM) which fits RNA decay curves of rifampicin treated bacteria. Its parameters have been described in [^f1]. It is compatible with Stan (>=2.26). Since it is based on normalized log-counts, it is suitable for sequencing experiments with similar library sizes. If the library sizes are very different, we recommend using LNMcdv-1.0.stan, which models the variance-mean relationship.
  - **LNMcdv-1.0.stan** Similar to the LNM, the variance depends on the raw sequencing counts in this model. Its use is recommanded for sequencing experiments with large differences in library size. Computationally, it is more efficient than a count-based (negative binomial) model. Run the example *model_comparison.Rmd*, to compare the LNM with and without mean-variance modeling.
  - **LNM_loo-1.0.stan** Same as LNM-1.0.stan, but the quantities requried for LOO-CV and posterior predictive checks are computed in the generated quantities block. This makes the model computationally less efficient. Therefore, it should only be used for model comparison purposes.
  - **LNMcdv_loo-1.0.stan** Same as LNMcdv-1.0.stan, but the quantities requried for LOO-CV and posterior predictive checks are computed in the generated quantities block. This makes the model computationally less efficient. Therefore, it should only be used for model comparison purposes.
  - **LNM.stan** This is the original model used to extract RNA decay rates in [^f1]. It is compatible and tested with Stan 2.28 - 2.31. Starting from Stan 2.32, the old syntax of array definitions is no longer suported https://mc-stan.org/docs/2_29/reference-manual/brackets-array-syntax.html . We recommend using LNM-1.0.stan instead, where the array definitions have been updated to the new syntax (Stan >=2.26).
  - **LNMsim.stan** simulates RNA decay curves for multiple experimental conditions or bacterial strains (Stan <=2.31).

</div>
If you have questions on how to adapt the statistical models to your sequencing data or if you are interested in using count-based models, please contact the authors directly.

## Exporting the read counts and metadata to cmdstan format with R

Export read counts and metadata for cmdstan using the rstan function ```stan_rdump```. The following quantities are required:

<div class="columns-2">
  
  - ```N_con``` Number of experimental conditions or bacterial strains.
  - ```N_s``` Number of samples.
  - ```N_b``` An array containing the number of replicates for every condition.
  - ```N_t``` Number of time points.
  - ```N_tot``` Total number of data points (number of rows of the long data frame)
  - ```N_g``` Number of genetic features.
</div>

Organize the read counts and metadata in a long data frame ```stan_df``` with the following columns:

<div class="columns-2">
  
  - ```raw``` raw counts, numeric.
  - ```time``` time (in minutes) which corresponds to the sample index ```it_s```, numeric.
  - ```sample``` unique sample ID, character.
  - ```condition``` label of condition, i.e. WT, dRBP, character.
  - ```locus_tag``` locus tag of genetic feature, character.
  - ```it_b``` replicate index, runs from 1,...,```N_b[it_c]```, numeric.
</div>

Furthermore, ```N_s``` normalization factors ```n_f``` (e.g. calculated from ERCC spike-ins) and RNA-seq library sizes ```N_lib``` (total number of reads per sample) are required.

Export the data via ```stan_rdump``` (see *example_cmdstan.Rmd* for a working example):

```
# Define iterators required by the Stan models
# condition, sample and locus_tag have to be factors, replicate and time numeric.
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

As explained in [^f1], we recommend using the log-normal model for the analysis of RIF-seq data, because it accounts for common confounding factors of rifampicin treated bacteria. It includes the delay because of ongoing transcription as well as a finite baseline concentration of stable RNA observable in RIF-seq data.

## Importing the cmdstan results with R

cmdstan saves the results from the HMC sampler in a csv file, e.g. RIF-seq_results.csv. Since the csv files can be quite large, import them using ```fread```, e.g.

```
stan_samples <- fread(cmd = 'grep -v "#" RIF-seq_results.csv', data.table = F)
```

The WT (control) decay rates of gene ```it_g``` are saved in ```betas.1.it_g``` (LNM) or ```beta_WT.it_g``` (LNM-1.0, LNMcdv-1.0), the differences in decay rate for condition/strain ```it_c``` are saved in ```betas.it_c.it_g``` (LNM) or ```delta_beta.(it_c-1).it_g``` (LNM-1.0, LNMcdv-1.0). Calculate medians and quantiles from the data frame ```stan_samples``` as follows:

```
beta_WT <- stan_samples[,grep("^betas\\.1\\.", colnames(stan_samples)] %>% apply(2, function(x) quantile(x, probs=c(0.05, 0.5, 0.95))
```

> **Note**
> 
> If the files are too large to import them with fread, extract the columns containing the relevant parameters using command line tools.

## Simulating decay curves

Use the RIF-seq_data.R file to simulate decay curves with the same number of replicates/samples/conditions/genes as your real data. Alternatively, create a RIF-seq_sim_data.R file containing only the required information following the same steps as for the RIF-seq_data.R file, omitting the variables ```raw, N_lib, n_f, batch_effects, model_id```.

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
