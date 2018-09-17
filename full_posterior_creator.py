import warnings
warnings.simplefilter('ignore')

import logging
logging.basicConfig(filename='GRESKE.log')

# Scientific libraries
import numpy as np

# Import Pandas
import pandas as pd

# Astro
import astropy.io.fits as fits
from astropy import units as u
from astropy.cosmology import WMAP9 as cosmo

# Graphic libraries
# if script:
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#if notebook:
#%matplotlib notebook
#import matplotlib.pyplot as plt

import scipy
from  scipy.stats import gaussian_kde
import cmath

#optional 3ML imports

from threeML import *
from threeML.io.file_utils import *
from threeML.random_variates import RandomVariates

import os
import re

# my new imports
from grb_list import besties_list, multies_list, pulses_list, names_list
from fit_results_getter import ResultGetter, ResultGetter_single, ResultGetter_pulse
from models_params import Variates, Params

import vapeplot



colors_short = []
for color in vapeplot.palette('vaporwave'):
    colors_short.append(color)

colors = []
for color in vapeplot.palette('sunset'):
    colors.append(color)
for color in vapeplot.palette('cool'):
    colors.append(color)
for color in vapeplot.palette('jazzcup'):
    colors.append(color)
for color in vapeplot.palette('seapunk'):
    colors.append(color)
for color in vapeplot.palette('macplus'):
    colors.append(color)
colors.append('gray')


besties = besties_list()
multies = multies_list()
pulsies = pulses_list()
grb_names = names_list()


decays_list = pd.read_csv('rises_and_decays.txt', delimiter='\t')

names = decays_list.GRB_name
lengths = decays_list.all_bins
peak_bins_id = decays_list.peak-1 # position is bin number - 1



# get all timebins results, new way; each result is list of timebins of that grb

fit_results = []

for trig, GRB_name, z in besties:
    pulse = 0    
    fit_result = ResultGetter_pulse(GRB_name, pulse)
    fit_results.append(fit_result)
    
    
    
rootdir = '/data29s/fermi/abacelj/GRB_correlations/GRB_builder/'

for trig, GRB_name, z in multies:
    
    with within_directory(rootdir):

        with within_directory(GRB_name):
            number_of_pulses = 0
            for i in range(1,10):
                pulse_dir = 'prepared_pha_files_%d'%i
                if(path_exists_and_is_directory(pulse_dir)==True):
                    number_of_pulses += 1

        for pulse in range(1,number_of_pulses+1):
            
            fit_result = ResultGetter_pulse(GRB_name, pulse)
            fit_results.append(fit_result)
            
            
# get decay phase of fit_results

decay_results = []

for fit_result in fit_results:
    
    index = fit_results.index(fit_result)
    
    decay_result = fit_result[peak_bins_id[index]:]
    decay_results.append(decay_result)
    
    
# get POSTERIOR model for each grb (grb == pulse) and each of it's decay timebins

### timebin_model_sample is actually list of all samples of that timebin


models_all = []
models_all_best = []
Epeaks_all = []


# this loop goues through 38 of grbs
for decay_result in decay_results:
    
    models_grb = []
    models_grb_best = []
    Epeaks_grb = []
    
    # this loop goes through timebins of current grb
    for timebin in decay_result:
        
        model_timebin_samples = [] 
        
        variate = Variates(timebin)
        loop_len = variate.length

        if variate.name == 'Band_grbm':

            for id in range(loop_len):
                param = Params(id,variate)
                band = Band_grbm(alpha=param.alpha, beta=param.beta, K=param.K, xc=param.xc, piv=100.0)
                model_timebin_samples.append(band)
            
            band_best = Band_grbm()
            band_best.alpha = np.median(variate.alpha)
            band_best.beta = np.median(variate.beta)
            band_best.K = np.median(variate.K)
            band_best.xc = np.median(variate.xc)
            band_best.piv = 100.0
            models_grb_best.append(band_best)
            
            Ecut_dist = variate.xc
            alpha_dist = variate.alpha

            
        if variate.name == 'Cutoff_powerlaw':
            
            for id in range(loop_len):
                param = Params(id,variate)
                cutoffpl = Cutoff_powerlaw(index = param.index, K = param.K, xc = param.xc, piv = 100.0)
                model_timebin_samples.append(cutoffpl)
                
            cutoffpl_best = Cutoff_powerlaw()
            cutoffpl_best.index = np.median(variate.index)
            cutoffpl_best.K = np.median(variate.K)
            cutoffpl_best.xc = np.median(variate.xc)
            cutoffpl_best.piv = 100.0
            models_grb_best.append(cutoffpl_best)
                
            Ecut_dist = variate.xc
            alpha_dist = variate.index
        
        
        Epeak = (alpha_dist + 2.0) * Ecut_dist
        Epeak = np.median(Epeak)
        
        models_grb.append(model_timebin_samples)
        Epeaks_grb.append(Epeak)

    models_all.append(models_grb)
    models_all_best.append(models_grb_best)
    Epeaks_all.append(Epeaks_grb)


    
# plot each grb on its own - full posterior

from matplotlib.patches import Patch
from matplotlib.lines import Line2D


dat = np.linspace(10,30000,100000)
x_range = np.linspace(0,1e5,100000)


for models_grb_best in models_all_best:
    
    fig, ax = plt.subplots()
    fig.set_size_inches(10,7)
    
    legend_elements = [Line2D([0], [0], c='black', alpha=0.6, lw=2, label='fit spectra'),
                   Line2D([0], [0], marker='o', c='black', alpha=0.8, label='peak energies',linewidth=0 , markersize=8),
                   Patch(facecolor='gray', alpha=0.3 ,label='data range\n[10keV-30MeV]')]
    
    index = models_all_best.index(models_grb_best)
    Epeaks_grb = Epeaks_all[index]
    models_grb = models_all[index]
 
    for i in range(len(models_grb_best)):
        
        model_timebin_best = models_grb_best[i]
        model_timebin_samples = models_grb[i]
        Epeak = Epeaks_grb[i]
        
        #color sheme
        if (len(models_grb_best))>14:
            color = colors[i]
        else:
            color = colors_short[i]  
        
            
        for model_sample in model_timebin_samples:

            # plot model sample over wide range
            ax.plot(x_range, x_range*x_range*model_sample(x_range), c=color, alpha=0.01, linewidth=0.6)


        # plot best (median) model over wide range
        ax.plot(x_range, x_range*x_range*model_timebin_best(x_range), c='black', alpha=0.8, linewidth=0.8, zorder=3)

        #plot peaks of spectrums
        ax.scatter(Epeak, np.max(dat*dat*model_timebin_best(dat)),color='black', s=30, alpha=0.6, zorder = 4)
        
        marker = Line2D([0], [0], marker='o', c=color, alpha=0.8, label=('timebin %d'%(i+1)), linewidth=0.8, markersize=6)
        legend_elements.append(marker)

    
    #highlight area where i have data
    ax.fill_between(dat,y1=1e0,y2=1e4, color='gray', alpha=0.2)

    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.set_ylim(1e0,1e4)
    ax.set_xlim(1e0,1e5)

    ax.set_xlabel('Energy [keV]')
    ax.set_ylabel(r'$E^2 N_E \; [erg \cdot erg\;cm^{-2}\;s^{-1}\;keV^{-}]$')
    
    
    #ax.legend(handles=legend_elements, loc='upper right');
    
    #Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax.legend(handles=legend_elements,loc='center left', bbox_to_anchor=(1, 0.5))
    

    #ax.legend(handles=legend_elements, loc='upper left');
    plt.title('%s'%(names[index]), fontsize=14)

    fig.savefig('plots_from_posteriors/%s_full.png'%(names[index]), bbox_inches="tight", frameon=True, overwrite=True)
    plt.close(fig)



