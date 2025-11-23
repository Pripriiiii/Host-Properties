#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  1 15:14:39 2025

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


# debass = pd.read_csv('/Users/dingyuancao/Desktop/DEBASS_allHost/DEBASS_allHost.csv')
# abs_path = '/Users/dingyuancao/Desktop/RSAA-Summer/Data-Analysing/images/'
# file_name = 'global_photometry.csv'

def mean_band_error(df, bands):
    errs = []
    for b in bands:
        col = f"{b}_err"
        if df is not None and col in df.columns and len(df) > 0:
            v = df[col].iloc[0]
            if pd.notna(v):
                errs.append(float(v))
    # If no valid error values, treat as "inf" so it loses comparisons
    return np.mean(errs) if errs else float("inf")


def append_placeholders(new_line, bands):
    """Append -99 for value and error per band."""
    new_line.extend(["-99 "] * (2 * len(bands)))

def append_cleaned_data(new_line, df, bands):
    """
    Append magnitude + error for each band.
    If magnitude is NaN OR error is NaN OR error > 1 → both replaced with -99.
    """
    for b in bands:
        mag = None
        err = None
        if df is not None and len(df) > 0:
            if b in df.columns:
                mag = df[b].iloc[0]
            if f"{b}_err" in df.columns:
                err = df[f"{b}_err"].iloc[0]

        # condition: missing magnitude or invalid error
        if pd.isna(mag) or pd.isna(err) or (err > 1):
            new_line.extend(["-99 ", "-99 "])
        else:
            new_line.extend([str(round(mag,3))+ ' ', str(round(err,3))+ ' '])
            
            
survey_filters = {
    'GALEX': ['NUV'], 
    'SDSS': ['u', 'g', 'r', 'i', 'z'], 
    'PanSTARRS': ['g', 'r', 'i', 'z', 'y'], 
    'DES': ['g', 'r', 'i', 'z', 'Y'], 
    'LegacySurvey': ['g', 'r', 'z']}

SDSS_remain_filters= ['g', 'r', 'i', 'z']
surveys = ['SDSS','PanSTARRS','DES','GALEX','LegacySurvey']
header = ['# FILTERS '] 
CALIBTYPES = ['# ', 'CALIBTYPES ']
for survey in survey_filters.keys():
    for filter in survey_filters[survey]:
        if survey == 'SDSS':
            header.append(f'sdss_dr7_{filter}, ')
            CALIBTYPES.append('AB, ')
        elif survey == 'PanSTARRS':
            header.append(f'PS1{filter},')
            CALIBTYPES.append('AB, ')
        else:
            header.append(f'{survey}{filter}, ')
            CALIBTYPES.append('AB, ')

def derive_PS1DES(new_line):
    GALEX_bands = survey_filters['GALEX']
    append_cleaned_data(new_line, GALEX_df, GALEX_bands)
    new_line.extend(["-99 "]*10)
    ps1_bands  = survey_filters["PanSTARRS"]
    des_bands  = survey_filters["DES"]
    ps1_mean = mean_band_error(PS1_df, ps1_bands)
    des_mean = mean_band_error(DES_df, des_bands)
    if ps1_mean < des_mean:
        # Choose PanSTARRS; DES becomes -99
        append_cleaned_data(new_line, PS1_df, ps1_bands)
        append_placeholders(new_line, des_bands)
    else:
        # Choose DES; PanSTARRS becomes -99
        append_placeholders(new_line, ps1_bands)
        append_cleaned_data(new_line, DES_df, des_bands)

def write_txt(data_file, txt_name, abs_path, file_name):
    with open(f'{txt_name}.txt', 'w') as f:
        f.writelines(header)
        f.write('\n')
        f.writelines(CALIBTYPES)
        f.write('\n')
        for name in list(data_file['snid']):
            new_line = [f'{name} ']
            new_line.append(str(round(data_file[data_file['snid']==name]['Redshift'].values[0],3))+ ' ')
            new_line.extend(['-9 ']*3)
            SDSS_file_path = os.path.join(abs_path, name, 'SDSS', file_name)
            PS1_file_path = os.path.join(abs_path, name, 'PanSTARRS', file_name)
            DES_file_path = os.path.join(abs_path, name, 'DES', file_name)
            GALEX_file_path = os.path.join(abs_path, name, 'GALEX', file_name)
            Legacy_file_path = os.path.join(abs_path, name, 'LegacySurvey', file_name)
            # Initialize to None so they exist even if the file is missing
            SDSS_df = PS1_df = DES_df = GALEX_df = LEGACY_df = None
    
            if os.path.exists(SDSS_file_path):
                SDSS_df = pd.read_csv(SDSS_file_path, nrows=1)  # only need first row
            if os.path.exists(PS1_file_path):
                PS1_df = pd.read_csv(PS1_file_path, nrows=1)
            if os.path.exists(DES_file_path):
                DES_df = pd.read_csv(DES_file_path, nrows=1)
            if os.path.exists(GALEX_file_path):
                GALEX_df = pd.read_csv(GALEX_file_path, nrows=1)
            if os.path.exists(Legacy_file_path):
                LEGACY_df = pd.read_csv(Legacy_file_path, nrows=1)
    
            if SDSS_df is not None and not SDSS_df.empty:
                all_fail = False
    
                # Helper to read a scalar safely
                def get0(df, col, default=-99):
                    if col in df.columns and len(df) > 0:
                        val = df[col].iloc[0]
                        return default if pd.isna(val) else val
                    else:
                        return default
    
                u_err = get0(SDSS_df, 'u_err')
    
                if pd.notna(u_err): 
                    if u_err > 0.5:
                        # If u_err > 0.5, check each band’s err and mark fail if any > 0.5
                        for band in SDSS_remain_filters:
                            berr = get0(SDSS_df, f'{band}_err')
                            if pd.isna(berr) or berr > 0.5:
                                # you likely want to add 38 placeholders as *38 columns*
                                derive_PS1DES(new_line)
                                Legacy_bands = survey_filters['LegacySurvey']
                                append_cleaned_data(new_line, LEGACY_df, Legacy_bands)
                                all_fail = True
                                break
                        if all_fail:
                            f.writelines(new_line)
                            f.write('\n')
                            f.write('\n')
                            continue
                        else:
                            # Add two placeholders as *two columns*
                            GALEX_bands = survey_filters['GALEX']
                            append_cleaned_data(new_line, GALEX_df, GALEX_bands)
                            new_line.extend(['-99 ', '-99 '])
                            for band in SDSS_remain_filters:
                                new_line.append(str(round(get0(SDSS_df, band), 3)) + ' ')
                                new_line.append(str(round(get0(SDSS_df, f'{band}_err'), 3))+ ' ')
                            new_line.extend(['-99 '] * 20)
                    elif  u_err <= 0.5:
                        # u_err <= 0.5 
                        new_line.extend(['-99 ', '-99 '])
                        for band in survey_filters['SDSS']:
                            new_line.append(str(round(get0(SDSS_df, band),3))+ ' ')
                            new_line.append(str(round(get0(SDSS_df, f'{band}_err'),3))+ ' ')
                        # Add 28 placeholders as *28 columns*
                        new_line.extend(['-99 '] * 20)  # ← adjust 28 to your real column count
                        # Do NOT break; continue lets you process other surveys or write the line
                else:
                    for band in SDSS_remain_filters:
                            berr = get0(SDSS_df, f'{band}_err')
                            if pd.isna(berr) or berr > 0.5:
                                # you likely want to add 38 placeholders as *38 columns*
                                new_line.extend(['-99 '] * 38)   # ← adjust 38 to your real column count
                                all_fail = True
                                break
                    GALEX_bands = survey_filters['GALEX']
                    append_cleaned_data(new_line, GALEX_df, GALEX_bands)
                    new_line.extend(['-99 ', '-99 '])
                    if not all_fail:
                        for band in SDSS_remain_filters:
                            new_line.append(str(round(get0(SDSS_df, band), 3)) + ' ')
                            new_line.append(str(round(get0(SDSS_df, f'{band}_err'), 3))+ ' ')
                    else:
                        print(name, "SDSS u band and some of other bands fail")
                        continue
            else:
                derive_PS1DES(new_line)
    
            Legacy_bands = survey_filters['LegacySurvey']
            append_cleaned_data(new_line, LEGACY_df, Legacy_bands)
            f.writelines(new_line)
            f.write('\n')
            f.write('\n')

    

