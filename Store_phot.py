#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 23:19:35 2025

@author: dingyuancao
"""
import os
import sep
import csv
import glob
import queue
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


surveys=['PS1','DES','SDSS']

# Check existing survey and filter
def check_survey_filter(name):
    exits_survey=""
    for survey in surveys: # 取决于该星系在哪个/哪几个survey中有数据
        filters_tocheck = get_survey_filters(survey)
        non_exist_filter = check_existing_images(name, filters=filters_tocheck , survey=survey)
        if not non_exist_filter:
            exits_survey+=survey
    return exits_survey


# Extract and filter magnitudes and errors for host galaxies
def extract_phot(name): 
    phot_data = []
    obj_dir = os.path.join(workdir,name)
    exist_survey =[]
    for survey in surveys: # 取决于该星系在哪个/哪几个survey中有数据
        filters_tocheck = get_survey_filters(survey)
        non_exist_filter = check_existing_images(name, filters=filters_tocheck , survey=survey)
        if not non_exist_filter:
                exist_survey.append(survey)
                #print(exist_survey)
                exist_filters = get_survey_filters(survey)
                exist_csv = survey + '_global.csv' 
                exist_data = pd.read_csv(os.path.join(obj_dir,exist_csv))
                for filter_name in exist_filters:
                    phot_data.append(exist_data.at[0,filter_name])
                    phot_data.append(exist_data.at[0,filter_name + '_err'])
    return phot_data


def filterred_phot(name , z):
    filterred_mag_err = [name, z, -9,  -9, -9]
    phot_data = extract_phot(name)
    for m in range(0, len(phot_data),2):
        magnitude = phot_data[m]
        error = phot_data[m+1] 
        if magnitude >= 0 and error <= 1:
            filterred_mag_err.append(round(magnitude, 3))
            filterred_mag_err.append(round(error, 3))  
    return filterred_mag_err


# Generate all combinations of surveys
def generate_all_combinations(n):
    output_lst = []
    q = queue.Queue()
    q.put([0])
    q.put([1])
    while not q.empty():
        a = q.get()
        if len(a) < n:
            q.put(a + [0])
            q.put(a + [1])
        else:
            output_lst.append(a)
    # Remove [0, 0, 0] from the output list
    if [0] * n in output_lst:
        output_lst.remove([0] * n)                   
    return output_lst


# Create a dictionary with the names of the galaxies for each type of combination of surveys
def create_dict(num_surveys, data_file):
    # Define dictionary
    txt_dict = {}
    survey_possibilities = generate_all_combinations(num_surveys)

    for p in survey_possibilities:
        survey_combination = ''
        for idx in range(len(p)):
            if p[idx] == 1:
                survey_combination += surveys[idx]
        txt_dict[survey_combination] = []

    for name in data_file[301:321]['snid']:
        exist_survey = check_survey_filter(name)
        for key in txt_dict.keys():
            if key == exist_survey:
                txt_dict[key].append(name)
    return txt_dict
                

# Create header for each key in the dictionary
def create_header(survey_combination):
    survey_in_key = [s for s in surveys if s in survey_combination]
    PS1_FILTERS = ['PS1g', 'PS1r', 'PS1i', 'PS1z', 'PS1y']
    DES_FILTERS = ['DESg', 'DESr', 'DESi', 'DESz', 'DESy']
    SDSS_FILTERS = ['sdss_dr7_u','sdss_dr7_g','sdss_dr7_r','sdss_dr7_i','sdss_dr7_z']
    CALIBTYPES = ['#', 'CALIBTYPES '] +['AB,']*(5*len(survey_in_key)-1) + ['AB']
    for survey in survey_in_key:
        header = '# FILTERS '
        if 'PS1' in survey_in_key:
            for filt in PS1_FILTERS:
                header += filt + ', '
        if 'DES' in survey_in_key:
            for filt in DES_FILTERS:
                header += filt + ', '
        if 'SDSS' in survey_in_key: 
            for filt in SDSS_FILTERS:
                header += filt + ', '
        header = header[:-2]
        return header, CALIBTYPES


# Write the photometry data to a text file
def write_txt(num_surveys, data_file):
    txt_dict = create_dict(num_surveys, data_file)
    for key in txt_dict.keys():
        if txt_dict[key]:
            with open(f'{key}.txt', 'w') as file:
                # write the header
                FILTERS, CALIBTYPES = create_header(key)
                file.write(FILTERS)
                file.write('\n')
                file.write(' '.join(str(item) for item in CALIBTYPES))
                file.write('\n')

                # write the data
                for name in txt_dict[key]:
                    phot_data = filterred_phot(name, data_file[data_file['snid'] == name]['Redshift'].values[0])
                    file.write(' '.join(str(item) for item in phot_data))
                    file.write('\n')
                    file.write('\n')







