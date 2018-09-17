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


write_dir = '/data29s/fermi/abacelj/GRB_correlations/K_corrections_new_way/'


def models_writer(GRB_name,results):

    file_path = '%s/timebin_models.txt'%write_dir
    if file_existing_and_readable(file_path)==False:
        writer = open(file_path,'w')
        writer.write('GRB_name\ttimebin\tmodel\n')
    else:
        writer = open(file_path,'a')
        
    for timebin in results:
        
        timebin_index = results.index(timebin) +1
        # to get the name of the model, fancy way! heheheh B-)
        model_name = timebin.optimized_model.source.spectrum.to_dict()['main'].keys()[0]
        
        writer.write('%s\t%d\t%s\n'%(GRB_name,timebin_index,model_name))
  
    writer.close()
    
    
def write_grb_data(GRB_name, Epeak_rest, Epeak_rest_err, Lumin_rest, Lumin_rest_err, betas):

    file_path = '%s/correlation_data_Epeak_max.txt'%write_dir
    
    if file_existing_and_readable(file_path)==False:
        writer = open(file_path,'w')
        writer.write('GRB_name\ttimebin\tE_peak\tEpeak_err\tLuminosity\tLuminosity_err\tbeta\n')
    else:
        writer = open(file_path,'a')
        
    for id in range(len(Epeak_rest)):
    
        Epeak = Epeak_rest[id]
        Lumin = Lumin_rest[id]
        Epeak_err = Epeak_rest_err[id]
        Lumin_err = Lumin_rest_err[id]
        beta = betas[id]
        
        writer.write('%s\t%d\t%f\t%f\t%e\t%e\t%f\n'%(GRB_name,id,Epeak,Epeak_err,Lumin,Lumin_err,beta))
    
    writer.close()