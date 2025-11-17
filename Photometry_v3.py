#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 28 15:37:15 2025

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
from pathlib import Path

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
from hostphot.cutouts import download_images
from hostphot.photometry import global_photometry as gp
from hostphot.photometry import local_photometry as lp
from hostphot.photometry.sed_plotting import plot_sed
from hostphot.processing import coadd_images
from hostphot.processing import masking
# from hostphot.processing import objects_detection
from hostphot.utils import plot_fits, plot_image
from hostphot.surveys_utils import get_survey_filters, check_filters_validity, check_survey_validity
from hostphot.utils import check_work_dir
from astropy.utils.data import clear_download_cache

phot_type = 'global'
riz_list = ['r','i','z']
#surveys = ['PanSTARRS','SDSS','GALEX']
coadd_dict = {'PanSTARRS':'riz', 'SDSS':'riz', 'DES':'riz'} 
# surveys = ['PS1']
# coadd_dict = {'PS1':'riz'} 

#FN_Host_properties = pd.read_csv('FN_Host_properties.csv')

def check_existing_images(name, filters, survey):
    check_survey_validity(survey)

    check_work_dir(workdir)
    obj_dir = os.path.join(workdir, name + '/' + survey) 

    filters_without_image = []
    for filt in filters:
        filt_image = os.path.join(obj_dir, f'{survey}_{filt}.fits')
        if os.path.isfile(filt_image) is False:
            filters_without_image.append(filt)

    return filters_without_image

##
def mask_and_photometry(name,host_ra,host_dec,ra,dec,exist_survey,ref_filt,ref_survey):
    #filters = get_survey_filters(surveyy)
    masking.create_mask(name, host_ra, host_dec,
                            filt=ref_filt, survey=ref_survey,threshold=1.5,r=2,
                            ra=ra, dec=dec,  # to plot the SN position
                            save_plots=True, save_mask_params=True,
                            save_input=False)
    
    aperture_params = gp.extract_aperture(name, host_ra, host_dec,
                                          filt=ref_filt, survey=ref_survey,
                                          ra=ra, dec=dec,
                                          threshold=1.5, use_mask=True,
                                          optimize_kronrad=True, eps=0.001,
                                          gal_dist_thresh=-1,  
                                          save_aperture_params=True,
                                          save_plots=True)
    
    for survey in exist_survey:
        filters = get_survey_filters(survey)
        for filt in filters:
            try:
                masking.create_mask(name, host_ra, host_dec,
                                    filt=filt, survey=survey,
                                    ra=ra, dec=dec,  # to plot the SN position   
                                    threshold=1.5,
                                    ref_filt=ref_filt, # to get the mask from
                                    ref_survey=ref_survey,
                                    save_plots=True, save_mask_params=False,
                                    save_input=False) 
            except Exception as exc:
                print(f"{survey}-{filt}: {exc}")
    for survey in exist_survey: 
        gp.multi_band_phot(name,host_ra, host_dec, survey=survey,
                           ra=ra, dec=dec, use_mask=True,
                           correct_extinction=False,common_aperture=True,
                           ref_filt=ref_filt, ref_survey=ref_survey,
                           save_plots=True, save_results=True,
                           save_aperture_params=False,
                           raise_exception=False)
    return aperture_params

def global_photometry(name,z,ra,dec,host_ra,host_dec,surveys):
    # Download Image cutouts
    for survey in surveys:
        download_images(name,host_ra,host_dec,survey=survey, overwrite=False,
                                 save_input=False)
    # Image Coadding
    exist_survey = []
    exist_riz_survey = []
    ps1_riz_exist = False
    for survey in surveys:
        if survey == "PanSTARRS":
            ps1_riz_exist = False
        filters_tocheck = get_survey_filters(survey)
        non_exist_filter = check_existing_images(name, 
                                                 filters=filters_tocheck, 
                                                 survey=survey)
        # print("The survey is ", survey, "which has images of filters ", non_exist_filter)
        
        exist_filter = []
        compress_riz_lst = []
        
        for single_filter in filters_tocheck:
            if single_filter not in non_exist_filter:
                # print("THe appended survey is ", single_filter)
                exist_filter.append(single_filter)
        for single_filter in riz_list:
            if single_filter not in non_exist_filter and single_filter in filters_tocheck:
                # print("THe appended survey is ", single_filter)
                compress_riz_lst.append(single_filter)
        if exist_filter:
            exist_survey.append(survey)
        if compress_riz_lst == riz_list:
            exist_riz_survey.append(survey)
            if survey == "PanSTARRS":
                ps1_riz_exist = True
                coadd_images(name, filters='riz', survey=survey)
            # if survey == PS1, 继续下面的步骤: “if survey != "GALEX": ...”
            # else:
            #     if ps1_riz_exist == True, 回到 “for survey in surveys:”
            #     else: 继续下面的步骤: “if survey != "GALEX": ...”
            if ps1_riz_exist:
                continue # True => continue
            
            coadd_images(name, filters='riz', survey=survey)
    if not exist_survey:
        return
    print("Coadding completed")
    
    # Image Masking
    # Global Photometry
    # eps = 0.0005
    if 'PanSTARRS' in exist_survey:
        extracted_params = mask_and_photometry(name,host_ra,host_dec,ra,dec,exist_survey,'riz','PanSTARRS')

    else:
        extracted_params = mask_and_photometry(name,host_ra,host_dec,ra,dec,exist_survey,'riz',exist_riz_survey[0])
        
    print("Grabing shape params completed")
    

    # SED Plotting
    out_dir = Path(workdir,name)
    plot_sed(name, 'global', z=z, outfile = out_dir/"sed_global_v3.jpg")

    return extracted_params


    

    
   
    
   
    