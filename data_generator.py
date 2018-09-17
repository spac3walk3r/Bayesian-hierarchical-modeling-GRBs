import numpy as np
import scipy.stats as stats
import pystan
import matplotlib.pyplot as plt




class DataGenerator(object):

    def __init__(self, N_grbs, norm_mu=52, norm_sigma=1, gamma_mu=1.5, gamma_sigma=0.2, min_length=3, max_length=15):



        self._N_grbs = N_grbs
        self._norm_mu = norm_mu
        self._norm_sigma = norm_sigma

        self._gamma_mu = gamma_mu
        self._gamma_sigma = gamma_sigma

        self._min_length = min_length
        self._max_length = max_length


    def generate_sample(self, seed=123):

        # set the seed for reproducability
        np.random.seed(seed)

        self._grb_lengths = [int(x) for x in np.random.uniform(self._min_length, self._max_length, size = self._N_grbs)]

        gamma_rng = stats.norm(self._gamma_mu, self._gamma_sigma)

        norm_rng = stats.norm(self._norm_mu, self._norm_sigma)


        self._gammas = []
        self._norms = []
        self._epeak = []
        self._luminosity = []
        self._ep_obs = []
        self._luminosity_obs = []

        self._ep_err = []
        self._luminosity_err = []
        
        for n in range(self._N_grbs):

            gamma = gamma_rng.rvs()
            norm = norm_rng.rvs()


            epeaks = np.random.uniform(1,3, size=self._grb_lengths[n] )
            luminosities = norm + gamma * (epeaks - 2)

            ep_err = np.random.exponential(scale=0.05, size=self._grb_lengths[n])

            ep_obs = epeaks + stats.norm.rvs(loc=0, scale=ep_err, size=self._grb_lengths[n]) 

            luminoisty_err = np.random.exponential(scale=0.1, size=self._grb_lengths[n])

            luminosity_obs = luminosities + stats.norm.rvs(loc=0, scale=luminoisty_err, size=self._grb_lengths[n]) 

            # add on the values

            self._luminosity.extend(luminosities)
            self._epeak.extend(epeaks)
            self._ep_obs.extend(ep_obs)
            self._ep_err.extend(ep_err)
            self._luminosity_obs.extend(luminosity_obs)
            self._luminosity_err.extend(luminoisty_err)

            self._gammas.append(gamma)
            self._norms.append(norm)

    def write_data(self, filename='sim.data.R', n_model=100):


        data = dict(
            N = len(self._ep_obs),
            N_grbs = self._N_grbs,
            grb_length = self._grb_lengths,
            gamma_mu = self._gamma_mu,
            gamma_sigma = self._gamma_sigma,
            norm_mu = self._norm_mu,
            norm_sigma = self._norm_sigma,
            gamma = self._gammas,
            norm = self._norms,
            ep_obs = self._ep_obs,
            ep_err = self._ep_err,
            lum_obs = self._luminosity_obs,
            lum_err = self._luminosity_err,
            N_model = int(n_model),
            ep_model=np.linspace(.5,3.5,n_model)
        )

        pystan.stan_rdump(data,filename)

    def plot(self):

        fig, ax = plt.subplots()


        i = 0
        for n in range(self._N_grbs):

            length = self._grb_lengths[n]
            
            ax.errorbar(self._ep_obs[i:i+length],
                        self._luminosity_obs[i:i+length],
                        xerr = self._ep_err[i:i+length],
                        yerr = self._luminosity_err[i:i+length],
                        fmt='.'
            )

            i+=length


            

    
