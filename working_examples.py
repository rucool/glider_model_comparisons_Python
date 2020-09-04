#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 10:21:51 2019

@author: aristizabal
"""

#%% Cell #1: Read glider data from the IOOS thredds server and plots a    
#   scatter plot of a glider transect for the entire length of the deployment

from read_glider_data import read_glider_data_thredds_server

#url_thredds = 'https://data.ioos.us/thredds/dodsC/deployments/rutgers/ru33-20180801T1323/ru33-20180801T1323.nc3.nc'
url_thredds = 'http://gliders.ioos.us/thredds/dodsC/deployments/aoml/SG668-20190819T1217/SG668-20190819T1217.nc3.nc'
var_name = 'temperature'
#var = 'salinity'
scatter_plot = 'yes'

varg, timeg, latg, long, depthg, dataset_id = \
             read_glider_data_thredds_server(url_thredds,var_name,scatter_plot)             
             

#%% Cell #2: Read glider data from the IOOS thredds server and plots a    
#   scatter plot of a glider transect for a specific time window

from read_glider_data import read_glider_data_thredds_server

#url_glider = 'https://data.ioos.us/thredds/dodsC/deployments/rutgers/ru33-20180801T1323/ru33-20180801T1323.nc3.nc'
#date_ini = '2018/09/01/00' # year/month/day/hour
#date_end = '2018/09/10/00' # year/month/day/hour

url_thredds = 'http://gliders.ioos.us/thredds/dodsC/deployments/aoml/SG668-20190819T1217/SG668-20190819T1217.nc3.nc'
var_name = 'temperature'
#var = 'salinity'
date_ini = '2019/09/01/00' # year/month/day/hour
date_end = '2019/09/10/00' # year/month/day/hour
scatter_plot = 'yes'

kwargs = dict(date_ini=date_ini,date_end=date_end)

varg, timeg, latg, long, depthg, dataset_id = \
             read_glider_data_thredds_server(url_thredds,var_name,scatter_plot,**kwargs)

#%% Cell #3: Same as cell #2, in addition to interpolating the variable 
# of interest (temperature, salinity or density) to regular depth levels
# (gridding variables in the vertical)                      
    
from read_glider_data import read_glider_data_thredds_server
from process_glider_data import grid_glider_data

#url_thredds = 'https://data.ioos.us/thredds/dodsC/deployments/rutgers/ru33-20180801T1323/ru33-20180801T1323.nc3.nc'
#date_ini = '2018/09/01/00' # year/month/day/hour
#date_end = '2018/09/10/00' # year/month/day/hour

url_thredds = 'http://gliders.ioos.us/thredds/dodsC/deployments/aoml/SG668-20190819T1217/SG668-20190819T1217.nc3.nc'
var_name = 'temperature'
#var = 'salinity'
date_ini = '2019/09/01/00' # year/month/day/hour
date_end = '2019/09/10/00' # year/month/day/hour
scatter_plot = 'no'
kwargs = dict(date_ini=date_ini,date_end=date_end)

tempg, timeg, latg, long, depthg, dataset_id = \
             read_glider_data_thredds_server(url_thredds,var_name,scatter_plot,**kwargs)    

delta_z = 0.4 # bin size in the vertical when gridding the variable vertical profile 
              # default value is 0.3   
contour_plot = 'yes'  # default value is 'yes'   
tempg_gridded, timegg, depthg_gridded = \
                    grid_glider_data(var_name,dataset_id,tempg,timeg,latg,long,depthg,delta_z,contour_plot)   

#%% Cell #4: Search for glider data sets given 
#   a latitude and longitude box and time window 
          
from read_glider_data import retrieve_dataset_id_erddap_server

# Server location
url_erddap = 'https://data.ioos.us/gliders/erddap'

'''
# MAB
lon_lim = [-75.0,-72.0]
lat_lim = [38.0,40.0]

# date limits
date_ini = '2018/09/01/00'
date_end = '2018/09/10/00'
'''

# Caribbean
lon_lim = [-80,-60.0]
lat_lim = [10.0,30.0]

# date limits
date_ini = '2019/09/01/00'
date_end = '2019/09/10/00'

gliders = retrieve_dataset_id_erddap_server(url_erddap,lat_lim,lon_lim,date_ini,date_end)
print('The gliders found are ')
print(gliders)
             
#%% Cell #5: Search for glider data sets given a 
#   latitude and longitude box and time window, choose one those data sets 
#   (glider_id), plot a scatter plot of the chosen glider transect, grid 
#   and plot a contour plot of the chosen glider transect 
         
from read_glider_data import retrieve_dataset_id_erddap_server
from read_glider_data import read_glider_data_erddap_server
from process_glider_data import grid_glider_data

# Server location
url_erddap = 'https://data.ioos.us/gliders/erddap'

'''
# MAB
lon_lim = [-75.0,-72.0]
lat_lim = [38.0,40.0]

# date limits
date_ini = '2018/09/01/00'
date_end = '2018/09/10/00'
'''

# Caribbean
lon_lim = [-80,-60.0]
lat_lim = [10.0,30.0]

# date limits
date_ini = '2019/09/01/00'
date_end = '2019/09/10/00'

gliders = retrieve_dataset_id_erddap_server(url_erddap,lat_lim,lon_lim,date_ini,date_end)

dataset_id = gliders[10]

# variable to retrieve
var_name = 'temperature'
#var_name = 'salinity'

kwargs = dict(date_ini=date_ini,date_end=date_end)
scatter_plot = 'yes'

tempg, saltg, timeg, latg, long, depthg = read_glider_data_erddap_server(url_erddap,dataset_id,\
                                   lat_lim,lon_lim,scatter_plot,**kwargs)
    
#tempg, saltg, timeg, latg, long, depthg = read_glider_data_erddap_server(url_erddap,dataset_id,\
#                                   lat_lim,lon_lim,scatter_plot)  

contour_plot = 'yes' # default value is 'yes'
delta_z = 0.4     # default value is 0.3

tempg_gridded, timegg, depthg_gridded = \
                    grid_glider_data(var_name,dataset_id,tempg,timeg,latg,long,depthg,delta_z,contour_plot)
                                                   
#%% cell #6: Search for glider data sets given a 
#    latitude and longitude box and time window, choose one those data sets 
#    (dataset_id), grid in the vertical the glider transect, get the glider
#    transect in the GOFS 3.1 grid, and plot both the transect from the glider
#    deployment and GOFS 3.1 output

from read_glider_data import read_glider_data_erddap_server
from read_glider_data import retrieve_dataset_id_erddap_server
from process_glider_data import grid_glider_data
from glider_transect_model_com import get_glider_transect_from_GOFS

# Servers location
url_erddap = 'https://data.ioos.us/gliders/erddap'
url_GOFS = 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_93.0/ts3z'

'''
# MAB
lon_lim = [-75.0,-72.0]
lat_lim = [38.0,40.0]
# date limits
date_ini = '2018/09/01/00'
date_end = '2018/09/10/00'
'''

# Caribbean
lon_lim = [-80,-60.0]
lat_lim = [10.0,30.0]

# date limits
date_ini = '2019/09/01/00'
date_end = '2019/09/10/00'
kwargs = dict(date_ini=date_ini,date_end=date_end)

scatter_plot = 'yes'
contour_plot = 'yes' # default value is 'yes'
delta_z = 0.4     # default value is 0.3

# model variable name
model_name = 'GOFS 3.1'
var_name_model = 'water_temp'
var_name_glider = 'temperature'

gliders = retrieve_dataset_id_erddap_server(url_erddap,lat_lim,lon_lim,date_ini,date_end)
dataset_id = gliders[10]

tempg, saltg, timeg, latg, long, depthg = read_glider_data_erddap_server(url_erddap,dataset_id,\
                                   lat_lim,lon_lim,scatter_plot,**kwargs)
    
tempg_gridded, timegg, depthg_gridded = \
                    grid_glider_data(var_name_glider,dataset_id,tempg,timeg,latg,long,depthg,delta_z,contour_plot)

# Get temperature transect from model    
temp_GOFS, time_GOFS, depth_GOFS, lat_GOFS, lon_GOFS = \
              get_glider_transect_from_GOFS(url_GOFS,var_name_model,model_name,\
                                        tempg,timeg,latg,long,depthg,contour_plot)
                  
#%% cell #7: Search for glider data sets given a 
#    latitude and longitude box and time window, choose one those data sets 
#    (dataset_id), grid in the vertical the glider transect, get the glider
#    transect in the AmSeas grid, and plot both the transect from the glider
#    deployment and the AmSeas output

from read_glider_data import read_glider_data_erddap_server
from read_glider_data import retrieve_dataset_id_erddap_server
from process_glider_data import grid_glider_data
from glider_transect_model_com import get_glider_transect_from_Amseas

# Servers location
url_erddap = 'https://data.ioos.us/gliders/erddap'
url_amseas = 'https://www.ncei.noaa.gov/thredds-coastal/dodsC/amseas/amseas_20130405_to_current/' #'20190901/ncom_relo_amseas_u_2019090100_t003.nc'

# Caribbean
lon_lim = [-80,-60.0]
lat_lim = [10.0,30.0]

# date limits
date_ini = '2019/09/10/00'
date_end = '2019/09/15/00'
kwargs = dict(date_ini=date_ini,date_end=date_end)

scatter_plot = 'yes'
contour_plot = 'yes' # default value is 'yes'
delta_z = 0.4     # default value is 0.3

# model variable name
model_name = 'Amseas'
var_name_model = 'water_temp'
var_name_glider = 'temperature'

gliders = retrieve_dataset_id_erddap_server(url_erddap,lat_lim,lon_lim,date_ini,date_end)
dataset_id = gliders[10]

tempg, saltg, timeg, latg, long, depthg = read_glider_data_erddap_server(url_erddap,dataset_id,\
                                   lat_lim,lon_lim,scatter_plot,**kwargs)
    
tempg_gridded, timegg, depthg_gridded = \
                    grid_glider_data(var_name_glider,dataset_id,tempg,timeg,latg,long,depthg,delta_z,contour_plot)

temp_amseas, time_amseas, depth_amseas, lat_amseas, lon_amseas = \
              get_glider_transect_from_Amseas(url_amseas,var_name_model,model_name,\
                                        tempg,timeg,latg,long,depthg,contour_plot='yes') 
