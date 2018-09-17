functions {

  vector luminosity(real gamma, vector ep, real norm){
    // relation between "true" L and "true" Epeak (for each timebin, ofc)
    // -2 is for pivot of 100 keV (beacuse it is in log space)
    return norm + gamma * (ep - 2);
  }



}


data {		

  int<lower=0> N;	// number of all timebins of the whole sample
  int<lower=0> N_grbs;

  //  int bin_id[N]; // array of timebin id-s for all grbs (all sample)	  
  int grb_length[N_grbs];
  
  //input data from csv file		
  vector[N] ep_obs; // Epeak
  vector[N] lum_obs; // Luminosity
  vector<lower=0>[N] ep_err; // E_peak errors
  vector<lower=0>[N] lum_err; // Luminosity errors


  int<lower=0> N_model;
  row_vector[N_model] ep_model;
  
}

transformed data {


}

parameters {	
  // definitions of the parameters of the model
  
  // definitions of parameters used to perform non-centered
  // parametrization of "original" parameters gamma i norm		
 
  
  real gamma_mu;		
  real norm_mu;
  real<lower=0> gamma_sigma;
  real<lower=0> norm_sigma;
  

  vector[N_grbs] norm_raw;  
  vector[N_grbs] gamma_raw; 


  
  vector<lower=0>[N] Ep_latent; // not sure why this is here and L_latent there 

  // definitions of hyper-parameters
  
}

transformed parameters {
  // definitions of "original" parameters
  vector[N_grbs] norm;
  vector[N_grbs] gamma;
  vector<lower=0>[N] L_latent;
  real norm_shift;
  real gamma_shift;
  // shift the norm so that it is
  // sampling on a proper scale
  norm_shift = norm_mu +52;
  gamma_shift = gamma_mu +2;

  // non-centered  parameterization
  gamma = gamma_shift + gamma_sigma * gamma_raw;
  norm = norm_shift + norm_sigma * norm_raw;

  // loop through the data and set the latent luminosity
  {  
    int k = 1;
    for (n in 1:N_grbs) {
      L_latent[k:k+grb_length[n] -1] = luminosity(gamma[n], Ep_latent[k:k+grb_length[n] -1], norm[n]);

      // increment the length selection
      k += grb_length[n];
    }
  }
}


model {

  // hyper parameters
  gamma_mu ~ normal(0,1);
  gamma_sigma ~ normal(0,1);
  
  norm_mu ~ normal(0,1);
  norm_sigma ~ normal(0,1);

  // non-centered  parameterization
  gamma_raw ~ normal(0,1);
  norm_raw ~ normal(0,1);
  
  // distribution of "true", unknown value of Epeak (for all grbs and all timebins)		
  Ep_latent ~ normal(2,1); 	

 
  // distribution of input data -> Epeaks
  ep_obs ~ normal(Ep_latent, ep_err);  

  // distribution of input data -> Luminosities
  lum_obs ~ normal(L_latent, lum_err);
}


generated quantities {

  real gamma_true;
  matrix[N_grbs, N_model] fitted_lines;
  //  vector[N_model] group_line;

  
  for (n in 1:N_grbs) {
    
    row_vector[N_model] line = norm[n] + gamma[n] * (ep_model-2);
    
    fitted_lines[n] = line;

    
  }  

}

