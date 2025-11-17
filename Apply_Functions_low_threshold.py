#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 14:37:30 2025

@author: dingyuancao
"""

import os
import time
import glob
import math
import pandas as pd
import numpy as np
from matplotlib.colors import LogNorm
import astropy.io.fits as fits
import matplotlib as plt
from Photometry import global_photometry
from astropy import units as u
from astropy.coordinates import SkyCoord, Distance
from hostphot._constants import workdir


# Define variables and load dataset
surveys=['PanSTARRS','SDSS','DES','GALEX']
coadd_dict={'PanSTARRS':'riz', 'SDSS':'riz', 'DES':'riz'} 
riz_lst = ['r','i','z']

# surveys=['LegacySurvey']
# coadd_dict={'LegacySurvey':'grz'}
# grz_lst = ['g','r','z']

# surveys=['DES']
# coadd_dict={'DES':'riz'}
# riz_lst = ['r','i','z']

debass = pd.read_csv('/Users/dingyuancao/Desktop/DEBASS_allHost/DEBASS_allHost.csv')
names = list(debass['snid'])
redshifts = list(debass['Redshift'])
ras = list(debass['SN_RA'])
decs = list(debass['SN_DEC'])
host_ras = list(debass['HOST_RA'])
host_decs = list(debass['HOST_DEC'])

# desirt = pd.read_csv('DESIRT copy.csv')
# CID = list(desirt['CID'])
# names = [str(i) for i in CID]
# redshifts = list(desirt['REDSHIFT_CMB'])
# ras = list(desirt['RA'])
# decs = list(desirt['DEC'])
# host_ras = list(desirt['Host Legacy Survey RA'])
# host_decs = list(desirt['Host Legacy Survey Dec'])




# shape_para = []
# # Compute Global Photometry len(debass)
# for n in range(len(debass)): 
#     if np.isnan(host_ras[n]): 
#         n += 1
#         shape_para.append([names[n],'','',''])
#     else:
#         parameters = global_photometry(name=names[n],z=redshifts[n], 
#                                         ra = ras[n], dec = decs[n], 
#                                         host_ra=host_ras[n],host_dec=host_decs[n],
#                                         thresh1 = 15, thresh2 = 1.5)
#         try:
#             a = parameters[0]['a'][0]
#             b = parameters[0]['b'][0]
#             theta = parameters[0]['theta'][0]
#             shape_para.append([names[n],a,b,theta])
#             Shape_Params = pd.DataFrame(shape_para, columns=['snid','a','b','theta'])
#         except:
#             print(names[n],"has no images and params")
#         n += 1



index_error_names = []
shape_para = []
n = 0
surveys = ['PanSTARRS','SDSS','DES','GALEX']
repeat_time = 0
while n < len(names):
    try:
        print("The current n is ", n)
        if np.isnan(host_ras[n]):
            n += 1
            repeat_time = 0
            shape_para.append([names[n],'','',''])
        else:
            parameters = global_photometry(names[n],redshifts[n],ras[n],decs[n],host_ras[n],host_decs[n],surveys)
            try:
                a = parameters[0]['a'][0]
                b = parameters[0]['b'][0]
                theta = parameters[0]['theta'][0]
                shape_para.append([names[n],a,b,theta])
                Shape_Params = pd.DataFrame(shape_para, columns=['snid','a','b','theta'])
            except:
                print(names[n],"has no images and params")
            n += 1
            repeat_time = 0
    except IndexError as IE:
        print(f"IndexError at n = {n}: {IE}, Skipping to n = {n+1}")
        index_error_names.append(names[n])
        n += 1
        repeat_time = 0
    except Exception as e:
        time.sleep(2)
        repeat_time += 1
        if repeat_time >= 10:
            n += 1
            repeat_time = 0


