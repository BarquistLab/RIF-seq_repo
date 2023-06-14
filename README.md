# RIF-seq_repo

The directory *stan_models* contains to Stan model files:

**LNM.stan** extracts RNA decay rates from RNA sequencing data over a time course following treatment with the transcription initiation inhibitor rifampicin.
**LNMsim.stan** simulates RNA decay curves for multiple experimental conditions or bacterial strains.

Extensive documentation on how to develop and run Stan models can be found on the Stan website: https://mc-stan.org/ . Here, useful links regarding installing and running Stan models can be found in addition to instructions on how to prepare the data with **R** before running LNM.stan with cmdstan.

## Installing and running the Stan model

The Stan models are compatible and tested with Stan 2.28 - 2.31. Starting from Stan 2.32, the syntax of array definition has changed https://mc-stan.org/docs/2_29/reference-manual/brackets-array-syntax.html . To run the models with Stan 2.32 and higher, the array definitions have to be adjusted.

The Stan models can be run using various interfaces. Depeding on the size of the analyzed dataset, fitting the parameters can take several days even with multithreading. It is therefore recommended to use cmdstan: https://mc-stan.org/docs/2_31/cmdstan-guide/cmdstan-installation.html . Furthermore, LNM.stan uses the Stan function ```reduce_sum```. To reduce computation time, multithreading should be used following the guide https://mc-stan.org/docs/2_31/cmdstan-guide/parallelization.html .

After installing cmdstan, a directory for the Stan model can be created within the cmdstan directory and the Stan model can be installed following the cmdstan instructions: https://mc-stan.org/docs/2_31/cmdstan-guide/compiling-a-stan-program.html .

After creating a RIF-seq_data.R file, the decay rates can be fitted following the cmdstan instructions: https://mc-stan.org/docs/2_31/cmdstan-guide/mcmc-intro.html.

## Exporting the read counts and metadata to cmdstan format with R

The read counts and metadata can be exported for cmdstan using the rstan function ```stan_rdump```.

## Importing the cmdstan results with R
