# Scientific libraries
import numpy as np

# Import Pandas
import pandas as pd

# Astro
import astropy.io.fits as fits

#optional 3ML imports

from threeML import *
from threeML.io.file_utils import *

import os
import re

def ResultGetter(besties,multies):
    
    rootdir = '/data29s/fermi/abacelj/GRB_correlations/GRB_builder/'
    
    with within_directory(rootdir):
        
        results = []
    
        for trig, GRB_name, z in besties:

            file_list = []
            GRB_directory = "%s"%(GRB_name)     

            with within_directory(GRB_directory):

                if path_exists_and_is_directory('fit_results_BAYESIAN')==True:
                    for subdir, dirs, files in os.walk('fit_results_BAYESIAN'):
                        for file in sorted(files):
                            file_list.append(file)
                            file_list.sort(key=lambda var:[int(x) if x.isdigit() 
                                                           else x for x in re.findall(r'[^0-9]|[0-9]+', var)])

                    with within_directory('fit_results_BAYESIAN'):
                        for file in file_list:
                            result = load_analysis_results(file)
                            results.append(result)
                                    
                                    
                                    
        for trig, GRB_name, z in multies:
    
            GRB_directory = "%s"%(GRB_name)     
    
            with within_directory(GRB_directory):
                number_of_pulses = 0
                for i in range(1,10):
                    pulse_dir = 'prepared_pha_files_%d'%i
                    if(path_exists_and_is_directory(pulse_dir)==True):
                        number_of_pulses += 1

            for pulse in range(1,number_of_pulses+1):
             
                file_list = []
                with within_directory(GRB_directory):
                    if path_exists_and_is_directory('fit_results_BAYESIAN_%d'%pulse)==True:
                        for subdir, dirs, files in os.walk('fit_results_BAYESIAN_%d'%pulse):
                            for file in sorted(files):
                                file_list.append(file)
                                file_list.sort(key=lambda var:[int(x) if x.isdigit() 
                                                               else x for x in re.findall(r'[^0-9]|[0-9]+', var)])

                        with within_directory('fit_results_BAYESIAN_%d'%pulse):
                            for file in file_list:
                                result = load_analysis_results(file)
                                results.append(result)
                                
    return results


def ResultGetter_single(GRB_name):
        
    rootdir = '/data29s/fermi/abacelj/GRB_correlations/GRB_builder/'
    
    with within_directory(rootdir):
            
            results = []
            file_list = []
            GRB_directory = "%s"%(GRB_name)     

            with within_directory(GRB_directory):

                if path_exists_and_is_directory('fit_results_BAYESIAN')==True:
                    for subdir, dirs, files in os.walk('fit_results_BAYESIAN'):
                        for file in sorted(files):
                            file_list.append(file)
                            file_list.sort(key=lambda var:[int(x) if x.isdigit() 
                                                           else x for x in re.findall(r'[^0-9]|[0-9]+', var)])

                    with within_directory('fit_results_BAYESIAN'):
                        for file in file_list:
                            result = load_analysis_results(file)
                            results.append(result)
                    
    return results





def ResultGetter_pulse(GRB_name, pulse):
        
    rootdir = '/data29s/fermi/abacelj/GRB_correlations/GRB_builder/'
    
    results = []
    file_list = []
    
    GRB_directory = "%s"%(GRB_name)
    
    with within_directory(rootdir):
        

            with within_directory(GRB_directory):

                if pulse == 0:
   
                    if path_exists_and_is_directory('fit_results_BAYESIAN')==True:
                        for subdir, dirs, files in os.walk('fit_results_BAYESIAN'):
                            for file in sorted(files):
                                file_list.append(file)
                                file_list.sort(key=lambda var:[int(x) if x.isdigit() 
                                                               else x for x in re.findall(r'[^0-9]|[0-9]+', var)])

                        with within_directory('fit_results_BAYESIAN'):
                            for file in file_list:
                                result = load_analysis_results(file)
                                results.append(result)
                                        
                else:

                    if path_exists_and_is_directory('fit_results_BAYESIAN_%d'%pulse)==True:
                        for subdir, dirs, files in os.walk('fit_results_BAYESIAN_%d'%pulse):
                            for file in sorted(files):
                                file_list.append(file)
                                file_list.sort(key=lambda var:[int(x) if x.isdigit() 
                                                               else x for x in re.findall(r'[^0-9]|[0-9]+', var)])

                        with within_directory('fit_results_BAYESIAN_%d'%pulse):
                            for file in file_list:
                                result = load_analysis_results(file)
                                results.append(result)

    return results

