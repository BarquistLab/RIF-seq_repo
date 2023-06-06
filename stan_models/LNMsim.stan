data{
// number of conditions, samples, time points, genes, data points
    int<lower=1> N_con;
    int<lower=1> N_s;
    int<lower=1> N_b[N_con];
    int<lower=1> N_t;
    int<lower=1> N_g;
    int<lower=1> N_tot;
// iterator indicating gene
    int<lower=1> it_g[N_tot];
    int<lower=1> it_c[N_tot];
// iterator indicating gene + condition
    int<lower=1> it_gc[N_tot];
// iterator indicating gene + condition + time
    int<lower=1> it_gct[N_tot];
// iterator indicating gene + time
    int<lower=1> it_gt[N_tot];
// iterator indicating time
    int<lower=1> it_t[N_tot];
// sample id's
    int<lower=1> it_s[N_tot];
// time
    vector<lower=0>[N_tot] t;
}

transformed data{
    // number of samples
    vector[N_con-1] f[N_g];
    for (ic in 2:N_con) {
      f[, ic-1] = rep_array(1, N_g);
    }
}

parameters{
    // chose parameters such that all priors are centered at zero
    // to increase efficiency
    vector<lower=0>[N_g] betaWT; // WT decay rate
    vector[N_con-1] delta_beta[N_g]; // difference in decay rate compared to WT
    vector<lower=0, upper=12>[N_g] gamma; // elongation time (min)
    vector<lower=0>[N_g] sigma_g; // variation (unexplained by the model)
    vector<lower=0, upper=1>[N_con] pib_b[N_g]; // fraction of baseline concentration (as compared to t=0 min)
}

transformed parameters{
}

model{
    // priors
    betaWT ~ normal(0.75, 0.3);
    gamma ~ cauchy(0,2);
    for (i in 1:N_con) pib_b[,i] ~ normal(0, 0.05);
    
    for (i in 1:N_con-1) delta_beta[,i] ~ normal(0, 0.08);
    
    sigma_g ~ normal(0.35, 0.25);

}

generated quantities{
    // The generated quantity block generates the decay curves with the predicted relative log-counts for every genetic feature and genotype.
    vector[N_tot] logrel_pred;
    {
      vector[N_tot] x1; // indicator variable t>gamma
  
      for (i in 1:N_tot) {
        if (t[i] < gamma[it_g[i]]) {
          x1[i] = 0;
        } else {
          x1[i] = 1;
        }
      }
      
      for (i in 1:N_tot) {
        real pi_g;
        if (it_gc[i] == it_g[i]) {
          pi_g =  pib_b[it_g[i], it_c[i]] +  (1 - pib_b[it_g[i], it_c[i]]) * exp( - x1[i] * betaWT[it_g[i]] * ( t[i] - gamma[it_g[i]] ) );
        }
        else {
          pi_g =  pib_b[it_g[i], it_c[i]] +  (1 - pib_b[it_g[i], it_c[i]]) * exp( - x1[i] * (betaWT[it_g[i]] + delta_beta[it_g[i],it_c[i]-1]) * ( t[i] - gamma[it_g[i]] ) );
        }
        logrel_pred[i] = normal_rng(log(pi_g), sigma_g[it_g[i]]);
      }
    }
}
