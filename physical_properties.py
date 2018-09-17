import warnings
warnings.simplefilter('ignore')


# Scientific libraries
import numpy as np
import scipy.integrate as integrate

# Import Pandas
import pandas as pd

# Astro
import astropy.io.fits as fits
from astropy import units as u
from astropy.cosmology import WMAP9 as cosmo

# Graphic libraries

#%matplotlib notebook
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#optional 3ML imports

from threeML import *
from threeML.io.file_utils import *
from threeML.random_variates import RandomVariates

import os
import re
import sys
import traceback
import fnmatch

from ipywidgets import *

from models_params import Variates, Params


def my_band(x,K,alpha,xc,beta):
    
    piv = 100.0
    
    if (alpha < beta):
        raise ModelAssertionViolation("Alpha cannot be less than beta")

    if x <= (alpha - beta) * xc :
        out = x * K * np.power(x / piv, alpha) * np.exp(-x / xc)
    else:
        out = x * K * np.power((alpha - beta) * xc / piv, alpha - beta) * np.exp(beta - alpha) * \
                np.power(x / piv, beta)

    return out


def my_cutoffpl(x,K,index,xc):
    
    piv = 100.0
    # Compute it in logarithm to avoid roundoff errors, then raise it
    log_v = index * np.log(x / piv) - x / xc

    return x * K * np.exp(log_v)


def get_betas(results,GRB_name):
    
    betas = []
    
    for timebin in results:
        
        timebin_index = results.index(timebin) +1    
        # get the name of the model
        model_name = timebin.optimized_model.source.spectrum.to_dict()['main'].keys()[0]

        if model_name == 'Band_grbm':      
            beta = timebin.get_variates('source.spectrum.main.Band_grbm.beta')
            beta_med = beta.median
        else:
            beta_med = -100.0
        
        betas.append(beta_med)
    
    betas = np.array(betas)
    
    return betas


def Epeak_rest_frame(results,GRB_name,z):
    
    Epeaks = []
    Epeak_errs = []

    for timebin in results:
        
        # to get the name of the model, fancy way! heheheh B-)
        model_name = timebin.optimized_model.source.spectrum.to_dict()['main'].keys()[0]

        if model_name == 'Band_grbm':
            Ecut = timebin.get_variates('source.spectrum.main.Band_grbm.xc')
            alpha = timebin.get_variates('source.spectrum.main.Band_grbm.alpha')
        else:
            Ecut = timebin.get_variates('source.spectrum.main.Cutoff_powerlaw.xc')
            alpha = timebin.get_variates('source.spectrum.main.Cutoff_powerlaw.index')

        Epeak_obs = (alpha + 2)*Ecut

        #transpose to rest frame!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Epeak = Epeak_obs *(1+z)

        log_Epeak = np.log10(Epeak)

        hpd = log_Epeak.highest_posterior_density_interval()
        err_temp = hpd[1]-log_Epeak.median, log_Epeak.median-hpd[0]

        log_err = (np.array(err_temp)).max()


        timebin_index = results.index(timebin) +1

        #plot!
        fig, ax = plt.subplots()
        fig.suptitle('log_Epeak(rest) distribution - %s_timebin_%d'%(GRB_name,timebin_index), fontsize=12)
        ax.hist(log_Epeak, 50, alpha = 0.5)
        ax.scatter(hpd,[1,1], label='hpd', marker='v', s=100)
        ax.scatter(log_Epeak.median, 1, label='median', marker='v', s=100)
        plt.legend()
        fig.savefig('physical_properties_distribution_plots/%s_timebin_%d_distribution_logEpeak.png'%(GRB_name,timebin_index), bbox_inches="tight", frameon=True, overwrite=True)
        plt.close(fig)

        Epeaks.append(log_Epeak.median)
        Epeak_errs.append(log_err)
    
    Epeaks = np.array(Epeaks)
    Epeak_errs = np.array(Epeak_errs)
        
    Epeak_w_errs = Epeaks, Epeak_errs
    return Epeak_w_errs



def Flux_and_lumin_rest_frame(results,GRB_name,z,Epeak_rest):
    
    #epeak is log value, so we need to recalculate it back
    Ep = np.power(10,Epeak_rest)
    
    
    Fluxes = []
    Flux_errs = []
    
    Lumins = []
    Lumin_errs = []

    for timebin in results:
        
        ep_index = results.index(timebin)
        
        e1 = 0.5*Ep[ep_index]
        e2 = Ep[ep_index]
        
        variate = Variates(timebin)
        loop_len = variate.length
        
        if variate.name == 'Band_grbm':
            flux = []
            for id in range(loop_len):
                param = Params(id,variate)
                params = param.K,param.alpha,param.xc,param.beta
                temp_flux = integrate.quad(my_band,e1,e2, args=params)
                temp_flux = temp_flux[0] * 1.60218e-9 # from keV to ergs! (/cm2/s)
                flux.append(temp_flux)
            
        if variate.name == 'Cutoff_powerlaw':
            flux = []
            for id in range(loop_len):
                param = Params(id,variate)
                params = param.K,param.index,param.xc
                temp_flux = integrate.quad(my_cutoffpl,e1,e2, args=params)
                temp_flux = temp_flux[0] * 1.60218e-9 # from keV to ergs! (/cm2/s)
                flux.append(temp_flux)
        
        flux = RandomVariates(flux)
        
        log_flux = np.log10(flux)
        
        hpd = log_flux.highest_posterior_density_interval()
        err_temp = hpd[1]-log_flux.median, log_flux.median-hpd[0]
        
        log_flux_err = (np.array(err_temp)).max()
        
        timebin_index = results.index(timebin) +1
        
        #plot!
        fig, ax = plt.subplots()
        fig.suptitle('log_Flux distribution - %s_timebin_%d'%(GRB_name,timebin_index), fontsize=12)
        ax.hist(log_flux, 50, alpha = 0.5)
        ax.scatter(hpd,[1,1], label='hpd', marker='v', s=100)
        ax.scatter(log_flux.median, 1, label='median', marker='v', s=100)
        plt.legend()
        fig.savefig('physical_properties_distribution_plots/%s_timebin_%d_distribution_logFlux.png'%(GRB_name,timebin_index), bbox_inches="tight", frameon=True, overwrite=True)
        plt.close(fig)
        
        
        # now calculate luminosities:
        
        dL = cosmo.luminosity_distance(z)
        dL = dL.to(u.cm)
        dL = dL.value

        lumin = flux * 4 * math.pi * dL * dL
        
        log_lumin = np.log10(lumin)
        
        hpd = log_lumin.highest_posterior_density_interval()
        err_temp = hpd[1]-log_lumin.median, log_lumin.median-hpd[0]
        
        log_lumin_err = (np.array(err_temp)).max()

        #plot!
        fig, ax = plt.subplots()
        fig.suptitle('log_Luminosity distribution - %s_timebin_%d'%(GRB_name,timebin_index), fontsize=12)
        ax.hist(log_lumin, 50, alpha = 0.5)
        ax.scatter(hpd,[1,1], label='hpd', marker='v', s=100)
        ax.scatter(log_lumin.median, 1, label='median', marker='v', s=100)
        plt.legend()
        fig.savefig('physical_properties_distribution_plots/%s_timebin_%d_distribution_logLumin.png'%(GRB_name,timebin_index),bbox_inches="tight", frameon=True, overwrite=True)
        plt.close(fig)
        
        Fluxes.append(log_flux.median)
        Flux_errs.append(log_flux_err)
    
        Lumins.append(log_lumin.median)
        Lumin_errs.append(log_lumin_err)
        
    Fluxes = np.array(Fluxes)
    Flux_errs = np.array(Flux_errs)
    Lumins = np.array(Lumins)
    Lumin_errs = np.array(Lumin_errs)
    
    F_and_L_w_errs = Fluxes, Flux_errs, Lumins, Lumin_errs
    return F_and_L_w_errs