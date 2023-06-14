functions{
    // statistical model, running in parallel for individual genes
    vector ln_model(vector param, vector theta, real[] X, int[] y){
      int N_con = y[2];  // number of conditions
      int N_b [N_con] = y[4:(N_con+3)];  // number of batches/replicates
      int N_t = y[1];  // number of different time points
      int N_s = y[3];  // number of samples
      int model_id = y[2*N_s + N_con + 4]; // model switch, see definition below
      
      vector[N_s] t = to_vector(X[N_s+1:2*N_s]);  // list of time values
      vector[N_s] logrel = to_vector(X[1:N_s]);  // corresponding normalized log-cpm values (log-cpm at t=0 subtracted)
      int it_t[N_s] = y[(4 + N_con):(N_s + N_con + 3)];  // iterator indicating time
      int con[N_s] = y[(N_s + N_con + 4):(2*N_s + N_con + 3)];  // condition identification number (WT: con=1)
      
      vector[N_s] mu_g;  // logrel is distributed around mu_g
      vector[N_s] pi_g;  // mu_g=log(pi_g)
      
      vector[N_con] betas = to_vector(theta[1:N_con]);  // WT decay rate and differences to WT decay rate
      real gamma = theta[N_con+1];  // elongation time
      vector[N_con] pib = theta[(N_con+2):(2*N_con+1)]; // baseline value of exponential decay
      vector[N_t] sigma_g = theta[(2*N_con+2):(2*N_con + 1 + N_t)];  // fitted standard deviation
      
      real ll;  //log likelihood of one gene
      
      {
        vector[N_s] x1; // indicator variable t>gamma
        
        for (i in 1:N_s) {
          if (t[i] < gamma) {
            x1[i] = 0;
          } else {
            x1[i] = 1;
          }
        }
        
        for (i in 1:N_s) {
          if (con[i]==1) {  // WT
              if (model_id == 1) { // LNM
                pi_g[i] =  pib[con[i]] +  (1 - pib[con[i]]) * exp( - x1[i] * betas[1] * ( t[i] - gamma ) );
              } else if(model_id == 2) { // LM
                pi_g[i] = exp( - betas[1] * t[i] );
              } else { // PLM
                pi_g[i] =  exp( - x1[i] * betas[1] * ( t[i] - gamma ) );
              }
          }
          else {  // other conditions
              if (model_id == 1) { // LNM
                pi_g[i] =  pib[con[i]] +  (1 - pib[con[i]]) * exp( - x1[i] * (betas[1] + betas[con[i]]) * ( t[i] - gamma ) );
              } else if (model_id == 2) { // LM
                pi_g[i] =  exp( - (betas[1] + betas[con[i]]) * t[i] );
              } else { // PLM
                pi_g[i] =  exp( - x1[i] * (betas[1] + betas[con[i]]) * ( t[i] - gamma ) );
              }
          }
        }
      mu_g = log(pi_g);
      }

      // log-likelihood contribution of one gene
      ll = normal_lpdf(logrel | mu_g, sigma_g[it_t]);
      return [ll]';
    }
}

data{
// number of conditions, samples, replicates (per condition/strain), time points, genes, data points
    int<lower=1> N_con;
    int<lower=1> N_s;
    int<lower=1> N_b[N_con];
    int<lower=1> N_t;
    int<lower=1> N_g;
    int<lower=1> N_tot;
// count data
    vector[N_tot] raw; // raw counts
    vector[N_s] N_lib; // library sizes
    vector[N_s] n_f; // normalization factors
// iterator indicating gene
    int<lower=1> it_g[N_tot];
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
// batch id's
    int<lower=1> it_b[N_tot];
// condition id's
    int<lower=1> it_c[N_tot];
// time
    vector<lower=0>[N_tot] t;
// switch use center mean normalization yes=1/no=0
    int<lower=0, upper=1> batch_effects;
// choose model:
// 1: LNM
// 2: LM
// 3: PLM
    int<lower=1, upper=3> model_id;
}

transformed data{
    // number of time courses and corresponding iterator
    int<lower=1> N_b_max = max(N_b);
    vector[N_tot] logcpm;
    // logrel: normalized logcpm values, logcpm0 is subtracted
    vector[N_tot] logrel;
    // logcpm value at t=0, depends on condition and gene
    vector[N_con*N_g] logcpm_avg[N_t];
    vector[N_s] logcpm_sample;
    // matrices required for multithreading (to store "shards" of data)
    // order data by gene because most parameters are fitted gene-wise
    // real values are stored in X, integer values in y
    real X[N_g, 2*N_s];
    int y[N_g, 2*N_s + N_con + 4];
    
    logcpm = log( (raw + 0.5) ./ (N_lib[it_s] .*  n_f[it_s] + 1) * 1e06);

    // calculate gene-wise average logcpm at t=0
    for (j in 1:N_g) {
      for (con in 1:N_con) {
        for (it in 1:N_t) {
          real tmp = 0.;
          int counter = 0;
          for (i in 1:N_tot) {
            if (it_gc[i] == N_g*(con-1)+j && it_t[i]==it) {
              tmp += logcpm[i];
              counter += 1;
            }
          }
          logcpm_avg[it, (con-1)*N_g+j] = tmp/counter;
        }

      }
    }
    logrel = logcpm - logcpm_avg[1, it_gc];
    
    // center mean normalization
    for (is in 1:N_s) {
      real tmp=0;
      for (i in 1:N_tot) {
        if (it_s[i]==is) {
          tmp += logcpm[i] - logcpm_avg[it_t[i], it_gc[i]];
        }
      }
      logcpm_sample[is] = tmp/N_g;
    }
    // subtract center mean normalization constant if batch_effects==1
    logrel -= logcpm_sample[it_s]*batch_effects;
    
    // fill logrel and iterators into the matrices X and y
    // for parallel computing
    for (j in 1:N_g) {
      vector[N_s] glogcpm;  // logrel values for gene it_g
      real gt [N_s];  // corresponding time points
      int git_t [N_s]; // corresponding time iterator
      int git_c [N_s];
      int counter = 1;
      
      for (i in 1:N_tot) {
        if (it_g[i] == j) {
          glogcpm[counter] = logrel[i];
          gt[counter] = t[i];
          git_t[counter] = it_t[i];
          git_c[counter] = it_c[i];
          counter += 1;
        }
      }
      X[j, 1:N_s] = to_array_1d(glogcpm);
      X[j, N_s+1:2*N_s] = gt;
      y[j, 1] = N_t;
      y[j, 2] = N_con;
      y[j, 3] = N_s;
      y[j, 4:(N_con+3)] = N_b;
      y[j, (N_con + 4):(N_s + N_con + 3)] = git_t;
      y[j, (N_s + N_con + 4):(2*N_s + N_con + 3)] = git_c;
      y[j, 2*N_s + N_con + 4] = model_id;
    }
}

parameters{
    // chose parameters such that all priors are centered at zero
    // to increase efficiency
    vector[N_con] betas[N_g]; // decay rates (WT) and differences in decay rate
    vector<lower=0, upper=12>[N_g] gamma; // elongation time (min)
    vector<lower=0>[N_t] sigma_g[N_g]; // variation (unexplained by the model)
    vector<lower=0, upper=0.2>[N_con] pib_b[N_g]; // fraction of baseline concentration (as compared to t=0 min)
    
    // hyperparameters
    real<lower=-0.5> mu_beta;
    real<lower=0> sigma_beta;
    vector<lower=0>[N_t] mu_sigma;
    vector<lower=0>[N_t] sigma_sigma;
}

model{
    // priors of hyperparameters
    mu_sigma ~ cauchy(0, 1);
    sigma_sigma ~ normal(0.3, 0.3);
    mu_beta ~ cauchy(0, 1);
    sigma_beta ~ cauchy(0, 1);
    
    // priors of parameters of
    if (model_id == 1 || model_id == 3) gamma ~ cauchy(0,2);
    if (model_id == 1) for (i in 1:N_con) pib_b[,i] ~ normal(0, 0.25);
    
    betas[, 1] ~ normal(mu_beta + 0.5, sigma_beta); // WT decay rates
    for (i in 2:N_con) betas[,i] ~ normal(0, 0.2); // Differences in decay rate (compared to WT)
    
    for (i in 1:N_t) sigma_g[,i] ~ normal(mu_sigma[i], sigma_sigma[i]);

    
    {
      // store parameters in vectors for parallelization
      vector[2*N_t + 2] param; // parameters shared by all genes
      vector[2*N_con + 1 + N_t] theta[N_g]; // gene-specific parameters
      
      param[1:N_t] = mu_sigma;
      param[(N_t + 1):(2*N_t)] = sigma_sigma;
      param[(2*N_t + 1):(2*N_t + 2)] = [mu_beta, sigma_beta]';
      
      theta[,1:N_con] = betas;
      theta[,(N_con+1)] = to_array_1d(gamma);
      theta[,(N_con+2):(2*N_con + 1)] = pib_b;
      theta[,(2*N_con+2):(2*N_con + 1 + N_t)] = sigma_g;
  
      // sum over gene-wise contributions to the log-likelihood
      target += sum(map_rect(ln_model, param, theta, X, y));
    }
}

generated quantities{
  // comment in for model comparison
  //   vector[N_tot] log_lik;
  // {
  //   vector[N_tot] x1; // indicator variable t>gamma
  // 
  //   for (i in 1:N_tot) {
  //     if (t[i] < gamma[it_g[i]]) {
  //       x1[i] = 0;
  //     } else {
  //       x1[i] = 1;
  //     }
  //   }
  //   
  //   for (i in 1:N_tot) {
  //     real pi_g;
  //     if (it_gc[i] == it_g[i]) {
  //       if (model_id == 1) {
  //         pi_g =  pib_b[it_g[i], it_c[i]] +  (1 - pib_b[it_g[i], it_c[i]]) * exp( - x1[i] * betas[it_g[i],1] * ( t[i] - gamma[it_g[i]] ) );
  //       } else if(model_id == 2) {
  //         pi_g = exp( - betas[it_g[i],1] * t[i] );
  //       } else {
  //         pi_g =  exp( - x1[i] * betas[it_g[i],1] * ( t[i] - gamma[it_g[i]] ) );
  //       }
  //     }
  //     else {
  //       if (model_id == 1) {
  //         pi_g =  pib_b[it_g[i], it_c[i]] +  (1 - pib_b[it_g[i], it_c[i]]) * exp( - x1[i] * (betas[it_g[i],1] + betas[it_g[i],it_c[i]]) * ( t[i] - gamma[it_g[i]] ) );
  //       } else if (model_id == 2) {
  //         pi_g =  exp( - (betas[it_g[i],1] + betas[it_g[i],it_c[i]]) * t[i] );
  //       } else {
  //         pi_g =  exp( - x1[i] * (betas[it_g[i],1] + betas[it_g[i],it_c[i]]) * ( t[i] - gamma[it_g[i]] ) );
  //       }
  //     }
  //     log_lik[i] = normal_lpdf(logrel[i] | log(pi_g), sigma_g[it_g[i],it_t[i]]);
  //   }
  // }
}


