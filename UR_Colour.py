#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 00:37:37 2025

@author: dingyuancao
"""

import os
import sep
import csv
import glob
from mpdaf.obj import WCS
import aplpy
import fitsio
import pandas as pd
import numpy as np
import warnings


import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.image as mpimg
from matplotlib.patches import Ellipse
from matplotlib.ticker import PercentFormatter

from PyAstronomy import pyasl
from scipy.interpolate import interp1d
from scipy.integrate import trapezoid
from scipy.stats import ks_2samp
from scipy.stats import norm
from photutils.detection import DAOStarFinder
from photutils.centroids import centroid_2dg

from astropy.utils.exceptions import AstropyWarning
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from astropy.stats import mad_std
import astropy.units as u
import astropy.constants as conts
import astropy.coordinates as coord
from astropy.table import Table
from astropy.io import fits
import astropy.wcs as AsWCS
#from astropy import wcs
import scipy.integrate
from astropy.cosmology import FlatLambdaCDM
import requests
from PIL import Image
from io import BytesIO

import hostphot

from hostphot._constants import workdir
from hostphot._constants import font_family
from hostphot.cutouts import download_images
from hostphot.cutouts import check_existing_images
from hostphot.coadd import coadd_images
from hostphot.image_masking import create_mask
from hostphot.image_cleaning import remove_nan
from hostphot.utils import plot_fits, get_survey_filters
from hostphot.utils import suppress_stdout
from hostphot.sed_plotting import plot_sed
import hostphot.local_photometry as lp
import hostphot.global_photometry as gp

##

def computeColour(SED, filter1, filter2, system="AB", Vega=None):
    c=conts.c.value
    h=(conts.h).value
    # We compute the AB magnitudes by default
    # One can add Vega magnitudes if needed
    # Interpolate the SED and filter curves over steps in wavelength
    
    # Filter 1
    wmin=np.min(filter1['wavelength'])
    wmax=np.max(filter1['wavelength'])
    wstep=1.0
    
    wavelength1=np.arange(wmin,wmax+wstep,wstep)
    
    T_filter1=interp1d(filter1['wavelength'],filter1['transmission'],
                       bounds_error=False, fill_value=0.0)(wavelength1)
 
    SED_interp=interp1d(SED['wavelength'],SED['luminosity'],
                       bounds_error=False, fill_value=0.0)(wavelength1)
    
    
    # Compute the upper intergral
    upper=trapezoid(SED_interp * T_filter1 / c / h * wavelength1,wavelength1)
    
    # Compute the lower intergral
    if system == "AB":
        lower=trapezoid(3631e-23 * T_filter1 / c / h * wavelength1,wavelength1)
    else:
        Vega_interp=interp1d(Vega['wavelength'],Vega['luminosity'],
                       bounds_error=False, fill_value=0.0)(wavelength1)
    
        lower=trapezoid( Vega_interp * T_filter1 / c / h * wavelength1,wavelength1)   
 
    mag1=-2.5*np.log10(upper/lower)
    

    # Filter 2
    wmin=np.min(filter2['wavelength'])
    wmax=np.max(filter2['wavelength'])
    wstep=1.0

    wavelength2=np.arange(wmin,wmax+wstep,wstep)
    
    T_filter2=interp1d(filter2['wavelength'],filter2['transmission'],
                       bounds_error=False, fill_value=0.0)(wavelength2)

    SED_interp=interp1d(SED['wavelength'],SED['luminosity'],
                       bounds_error=False, fill_value=0.0)(wavelength2)


    # Compute the intergrals
    upper=trapezoid(SED_interp * T_filter2 / c / h * wavelength2, wavelength2)

    
    # Compute the lower intergral
    if system == "AB":
        lower=trapezoid(3631e-23 * T_filter2 / c / h * wavelength2,wavelength2)
    else:
        Vega_interp=interp1d(Vega['wavelength'],Vega['luminosity'],
                       bounds_error=False, fill_value=0.0)(wavelength2)
    
        lower=trapezoid(Vega_interp  * T_filter2 / c / h * wavelength2,wavelength2)
    
    mag2=-2.5*np.log10(upper/lower)
 
    colour=mag1-mag2
    
    return colour


def ApplycomputeColour(folderPath):
    # the folder should contain bes-fit SED data for each host.
    # Used Bessell filter curves used in SALT2
    u_filter="sux.dat" 
    r_filter="sr.dat"
    Vega="alpha_lyr_stis_005.ascii"

    Vega=np.genfromtxt(Vega, names = ['wavelength','luminosity'])
    u_filterCurve=np.genfromtxt(u_filter, names = ['wavelength','transmission'])
    r_filterCurve=np.genfromtxt(r_filter, names = ['wavelength','transmission'])


    # Contains the best-fitting SED files for each host galaxy
    results = []
    for filename in os.listdir(folderPath):
        if filename.endswith('.txt'):
            SED=np.genfromtxt(os.path.join(folderPath, filename), names = ['wavelength','luminosity'])
            colour=computeColour(SED, u_filterCurve, r_filterCurve, "Vega", Vega=Vega)
            results.append((filename[:-4], colour))

    Host_color_Vega = pd.DataFrame(results, columns=['snid', 'U-R Color'])
    return Host_color_Vega

