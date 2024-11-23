import os
import glob
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

import astropy.units as u
import astropy.constants as c

import hostphot
from hostphot.cutouts import download_images
from hostphot.cutouts import check_existing_images
from hostphot.coadd import coadd_images
from hostphot.image_masking import create_mask
from hostphot.utils import plot_fits, get_survey_filters
from hostphot.sed_plotting import plot_sed
import hostphot.local_photometry as lp
import hostphot.global_photometry as gp

Only_PS1 = pd.read_csv('Only_PS1.csv') # Includes targets that are only available in PS1 survey
In_Both = pd.read_csv('In_Both.csv') # Includes targets that are available in both PS1 and SDSS surveys
DEBASS_Host = pd.read_csv('DEBASS_Host.csv') # Includes filtered targets in DEBASS, with host ra&dec

all_names = DEBASS_Host['snid']
PS1_names = list(Only_PS1['snid'])
both_names = list(In_Both['snid'])

PS1_z = list(Only_PS1['Redshift'])
both_z = list(In_Both['Redshift'])
redshifts = list(DEBASS_Host['Redshift'])

ras = list(DEBASS_Host['SN_RA'])
decs = list(DEBASS_Host['SN_DEC'])
host_ras = list(DEBASS_Host['HOST_RA'])
host_decs = list(DEBASS_Host['HOST_DEC'])

def global_photometry(name,z,ra,dec,host_ra,host_dec, surveys, coadd_dict):
    #Download Image cutouts
    for survey in surveys:
        download_images(name,host_ra,host_dec, size = 6, survey=survey, overwrite=True)
    
    # Image Coadding
    for survey, coadd_filters in coadd_dict.items():
        coadd_images(name, filters=coadd_filters, survey=survey)
    
    # Image Masking
    surveyy = 'PS1'
    filters = get_survey_filters(surveyy)
    coadd_filters = coadd_dict[surveyy]
    coadd_mask_params = create_mask(name, host_ra, host_dec,
                                    filt='riz', survey= 'PS1',
                                    extract_params=True, crossmatch=True)
    sigma_dict = {survey:8 if survey!='GALEX' else 4 for survey in surveys}
    for survey, coadd_filters in coadd_dict.items():
        filters = get_survey_filters(survey)
        for filt in filters:
            create_mask(name, host_ra, host_dec,
                        filt, survey=survey,
                        common_params=coadd_mask_params,
                        sigma=sigma_dict[survey])

    
    # Global Photometry
    eps = 0.0005
    aperture_params = gp.extract_kronparams(name, host_ra, host_dec,
                                            filt = 'riz', survey = 'PS1',
                                            ra=ra, dec=dec, use_mask=True,
                                            optimize_kronrad=True, eps=eps,
                                            save_plots=True,
                                            save_aperture_params=True,threshold=1.5)
    for survey in surveys:
        gp.multi_band_phot(name, host_ra, host_dec,
                            survey=survey, ra=ra, dec=dec,
                            use_mask=True, correct_extinction=True,
                            aperture_params=aperture_params,
                            save_plots=True, save_results=True,
                            raise_exception=True)
        
    # SED Plotting
    plot_sed(name, 'global', z=z)
    return aperture_params

# NOTE that to obtain the best shape parameters, manually adjustments on thresholds in functions `gp.extract_kronparams` and `create_mask`
n = 0
shape_para = []
while n < len(all_names):
    if all_names[n] in PS1_names:
        parameters = global_photometry(name = all_names[n], z=redshifts[n],ra=ras[n], 
                          dec=decs[n],host_ra=host_ras[n], host_dec=host_decs[n],
                          surveys=['PS1'], coadd_dict={'PS1':'riz'})
        a = parameters[0]['a'][0]
        b = parameters[0]['b'][0]
        theta = parameters[0]['theta'][0]
        shape_para.append([all_names[n],a,b,theta])
        n += 1
    elif all_names[n] in both_names:
        parameters = global_photometry(name = all_names[n], z=redshifts[n],ra=ras[n], 
                          dec=decs[n],host_ra=host_ras[n], host_dec=host_decs[n],
                          surveys=['PS1','SDSS'], 
                          coadd_dict={'PS1':'riz','SDSS':'riz'})  
        a = parameters[0]['a'][0]
        b = parameters[0]['b'][0]
        theta = parameters[0]['theta'][0]
        shape_para.append([all_names[n],a,b,theta])
        n += 1
    else:
        n+=1
shape_para = pd.DataFrame(shape_para,columns=["snid","a","b","theta"])































    
    
    
    
    
