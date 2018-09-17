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

#optional 3ML imports

from threeML import *
from threeML.io.file_utils import *
from threeML.random_variates import RandomVariates

import os
import re
import sys
import traceback
import fnmatch


class Variates:
    
    def __init__(self, timebin):
        self._timebin = timebin
        # to get the name of the model, fancy way! heheheh B-)
        self.model_name = timebin.optimized_model.source.spectrum.to_dict()['main'].keys()[0]
    @property
    def length(self):
        return self._timebin.samples.shape[1]
    @property
    def name(self):
        return self._timebin.optimized_model.source.spectrum.to_dict()['main'].keys()[0]
    @property
    def K(self):
        return self._timebin.get_variates('source.spectrum.main.%s.K'%self.model_name)
    @property
    def xc(self):
        return self._timebin.get_variates('source.spectrum.main.%s.xc'%self.model_name)
    @property
    def alpha(self):
        return self._timebin.get_variates('source.spectrum.main.Band_grbm.alpha')
    @property
    def beta(self):
        return self._timebin.get_variates('source.spectrum.main.Band_grbm.beta')
    @property
    def index(self):
        return self._timebin.get_variates('source.spectrum.main.Cutoff_powerlaw.index')
    
    
class Params:
    
    def __init__(self,id,variate):
        self._variate = variate
        self.id = id
    @property   
    def name(self):
        return self._variate.name()
    @property
    def K(self):
        return self._variate.K[self.id]
    @property
    def xc(self):
        return self._variate.xc[self.id]
    @property
    def alpha(self):
        return self._variate.alpha[self.id]
    @property
    def beta(self):
        return self._variate.beta[self.id]
    @property
    def index(self):
        return self._variate.index[self.id]