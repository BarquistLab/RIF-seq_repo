# RIF-seq_repo

The directory *stan_models* contains to Stan model files:

**LNM.stan** extracts RNA decay rates from RNA sequencing data over a time course following treatment with the transcription initiation inhibitor rifampicin.
**LNMsim.stan** simulates RNA decay curves for multiple experimental conditions or bacterial strains.

## Installing and running the Stan model

The Stan models are compatible and tested with Stan 2.28 - 2.31. Starting from Stan 2.32, the syntax of array definition has changed [https://mc-stan.org/docs/2_29/reference-manual/brackets-array-syntax.html https://mc-stan.org/docs/2_29/reference-manual/brackets-array-syntax.html]
