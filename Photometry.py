#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 15:48:32 2025

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

# Define function for Global Photometry
riz_lst = ['r','i','z']
surveys=['PS1','DES','SDSS']
coadd_dict={'PS1':'riz','DES':'riz','SDSS':'riz'}

##
def mask_survey(surveyy, name, host_ra, host_dec, exist_survey):
    filters = get_survey_filters(surveyy)
    coadd_mask_params = create_mask(name, host_ra, host_dec,
                                    filt='riz', survey= surveyy,
                                    extract_params=True, crossmatch=True)
    
    sigma_dict = {survey:8 if survey!='GALEX' else 4 for survey in surveys}
    for survey in exist_survey:
        filters = get_survey_filters(survey)
        for filt in filters:
            create_mask(name, host_ra, host_dec,
                        filt, survey=survey,
                        common_params=coadd_mask_params,
                        sigma=sigma_dict[survey])


def global_photometry(name,z,ra,dec,host_ra,host_dec):
    # Download Image cutouts
    for survey in surveys:
        download_images(name,host_ra,host_dec, size = 10, survey=survey, overwrite=True)
    exist_survey = []
    for survey in surveys:
        filters_tocheck = get_survey_filters(survey)
        non_exist_filter = check_existing_images(name, filters=filters_tocheck , survey=survey)
        #print("The survey is ", survey, "which has images of filters ", non_exist_filter)
        compress_lst = []
        for single_filter in riz_lst:
            if single_filter not in non_exist_filter:
                #print("THe appended survey is ", single_filter)
                compress_lst.append(single_filter)
        if compress_lst:
            exist_survey.append(survey)
            filter_names = ""
            for filt in compress_lst:
                filter_names += filt
            #print("THe name is ", filter_names)
            coadd_images(name, filters=filter_names, survey=survey)
    if not exist_survey:
        return

    print("Coadding completed")
    # Image Masking
    # Global Photometry
    eps = 0.0005
    if 'PS1' in exist_survey:
        mask_survey('PS1', name, host_ra, host_dec, exist_survey)
        aperture_params = gp.extract_kronparams(name, host_ra, host_dec,
                                            filt = 'riz', survey = 'PS1',
                                            ra=ra, dec=dec, use_mask=True,
                                            optimize_kronrad=True, eps=eps,
                                            save_plots=True,
                                            save_aperture_params=True,threshold=15)
    else:
        mask_survey(exist_survey[0], name, host_ra, host_dec, exist_survey)
        aperture_params = gp.extract_kronparams(name, host_ra, host_dec,
                                            filt = 'riz', survey = exist_survey[0],
                                            ra=ra, dec=dec, use_mask=True,
                                            optimize_kronrad=True, eps=eps,
                                            save_plots=True,
                                            save_aperture_params=True,threshold=15) 
        
    print("Masking completed")
    
    for survey in exist_survey:
        gp.multi_band_phot(name, host_ra, host_dec,
                            survey=survey, ra=ra, dec=dec,
                            use_mask=True, correct_extinction=True,
                            aperture_params=aperture_params,
                            save_plots=True, save_results=True,
                            raise_exception=True,threshold=15)
    # SED Plotting
    plot_sed(name, 'global', z=z)

    return aperture_params

