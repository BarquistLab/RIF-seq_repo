# RIF-seq_repo

The directory *stan_models* contains Stan model files to analyze RIF-seq data:


<div class="columns-2">
  
  - **LNM.stan** extracts RNA decay rates from RNA sequencing data over a time course following treatment with the transcription initiation inhibitor rifampicin. Original model used in
  - **LNMsim.stan** simulates RNA decay curves for multiple experimental conditions or bacterial strains.
</div>

Extensive documentation on how to develop and run Stan models can be found on the Stan website: https://mc-stan.org/ . Here, useful links regarding installing and running Stan models can be found in addition to instructions on how to prepare the data with **R** before running LNM.stan with cmdstan.

## Installing and running the Stan model LNM.stan

The Stan models are compatible and tested with Stan 2.28 - 2.31. Starting from Stan 2.32, the syntax of array definition has changed https://mc-stan.org/docs/2_29/reference-manual/brackets-array-syntax.html . To run the models with Stan 2.32 and higher, the array definitions have to be adjusted.

The Stan models can be run using various interfaces. Depeding on the size of the analyzed dataset, fitting the parameters can take several days even with multithreading. It is therefore recommended to use cmdstan: https://mc-stan.org/docs/2_31/cmdstan-guide/cmdstan-installation.html .

> **Note**
> 
> LNM.stan uses the Stan function ```reduce_sum```. To reduce computation time, multithreading should be used following the guide https://mc-stan.org/docs/2_31/cmdstan-guide/parallelization.html .

After installing cmdstan, a directory for the Stan model can be created within the cmdstan directory and the Stan model can be installed following the cmdstan instructions: https://mc-stan.org/docs/2_31/cmdstan-guide/compiling-a-stan-program.html .

After creating a RIF-seq_data.R file, the decay rates can be fitted with ``` ./LNM sample data file=RIF-seq_data.R``` following the cmdstan instructions: https://mc-stan.org/docs/2_31/cmdstan-guide/mcmc-intro.html. The model was tested using 1000 warmup iterations and 1000 sampling iterations.

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

The data is then exported via ```stan_rdump```:

```
raw <- stan_df$raw
t <- stan_df$time
it_g <- stan_df$locus_tag %>% factor %>% as.numeric

# relevel stan_df$condition to the control condition.
# Differences in decay rate will be given relative to the control.
stan_df$condition <- factor(stan_df$condition) %>% relevel(ref="WT")
it_c <- stan_df$condition %>% as.numeric
it_t <- stan_df$time %>% factor %>% as.numeric
it_b <- stan_df$it_b
it_s <- stan_df$sample %>% factor %>% as.numeric

it_gc = (it_c - 1)*N_g + it_g
it_gct = (it_gc - 1)*N_t + it_t
it_gt = (it_g - 1)*N_t + it_t

N_g <- unique(stan_df$locus_tag) %>% length
N_con <- max(stan_df$it_c)
N_s <- max(stan_df$it_s)
N_b <- stan_df %>% group_by(condition) %>% summarize(N_b=unique(it_b) %>% length) %>% pull(N_b)
N_t <- unique(stan_df$time) %>% length
N_tot <- nrow(stan_df)
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

rstan::stan_rdump(c("N_con", "N_s", "N_b", "N_t", "N_g", "N_tot"), file=RIF-seq_data.R)
rstan::stan_rdump(c("raw", "t", "N_lib", "n_f", "batch_effects", "model_id"), file=RIF-seq_data.R, append=T)
rstan::stan_rdump(c("it_g", "it_gc", "it_gct", "it_gt", "it_t", "it_s", "it_b", "it_c"), file=RIF-seq_data.R, append=T)
```

## Importing the cmdstan results with R

cmdstan saves the results from the HMC sampler in a csv file, e.g. RIF-seq_results.csv. Since the csv files can be quite large, it is recommended to import them using ```fread```, e.g.

```
stan_samples <- fread(cmd = 'grep -v "#" RIF-seq_results.csv', data.table = F)
```

The WT (control) decay rates are saved in ```betas.1.it_g```, the differences in decay rate for condition/strain ```it_c``` are saved in ```betas.it_c.it_g```. Then, the medians and quantiles can be obtained from the data frame ```stan_samples``` as follows:

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
