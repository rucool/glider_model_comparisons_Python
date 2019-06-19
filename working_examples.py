#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 10:21:51 2019

@author: aristizabal
"""

#%% read_glider_data_thredds_server

from read_glider_data import read_glider_data_thredds_server

url_glider = 'https://data.ioos.us/thredds/dodsC/deployments/rutgers/ru33-20180801T1323/ru33-20180801T1323.nc3.nc'
var = 'temperature'
scatter_plot = 'yes'

varg, latg, long, depthg, timeg, inst_id = \
             read_glider_data_thredds_server(url_glider,var,scatter_plot)             
             

#%% read_glider_data_thredds_server

from read_glider_data import read_glider_data_thredds_server

url_glider = 'https://data.ioos.us/thredds/dodsC/deployments/rutgers/ru33-20180801T1323/ru33-20180801T1323.nc3.nc'
var = 'temperature'
#var = 'salinity'
date_ini = '2018/09/01/00' # year/month/day/hour
date_end = '2018/09/10/00' # year/month/day/hour
scatter_plot = 'yes'
kwargs = dict(date_ini=date_ini,date_end=date_end)

varg, latg, long, depthg, timeg, inst_id = \
             read_glider_data_thredds_server(url_glider,var,scatter_plot,**kwargs)
             
#%% retrieve_glider_id_erddap_server 
#read_glider_data_erddap_server
# grid_glider_data_erddap             

from read_glider_data import retrieve_glider_id_erddap_server
from read_glider_data import read_glider_data_erddap_server
from process_glider_data import grid_glider_data_erddap

# Server location
url_server = 'https://data.ioos.us/gliders/erddap'

# MAB
lon_lim = [-75.0,-72.0]
lat_lim = [38.0,40.0]

# date limits
date_ini = '2018-09-01T00:00:00Z'
date_end = '2018-09-10T00:00:00Z'

gliders = retrieve_glider_id_erddap_server(url_server,lat_lim,lon_lim,date_ini,date_end)

dataset_id = gliders[0]

# variable to retrieve
var = 'temperature'
#var = 'salinity'

scatter_plot = 'yes'

df = read_glider_data_erddap_server(url_server,dataset_id,var,\
                                   lat_lim,lon_lim,date_ini,date_end,\
                                   scatter_plot)

contour_plot = 'yes'
delta_z = 0.3

depthg_gridded, tempg_gridded, timeg, latg, long = \
                          grid_glider_data_erddap(df,var,delta_z,contour_plot)
                          
#%% read_glider_data_thredds_server
# grid_glider_data_thredd                        
    
from read_glider_data import read_glider_data_thredds_server
from process_glider_data import grid_glider_data_thredd

url_glider = 'https://data.ioos.us/thredds/dodsC/deployments/rutgers/ru33-20180801T1323/ru33-20180801T1323.nc3.nc'
var_name = 'temperature'
date_ini = '2018/09/01/00' # year/month/day/hour
date_end = '2018/09/10/00' # year/month/day/hour
scatter_plot = 'yes'
kwargs = dict(date_ini=date_ini,date_end=date_end)

varg, latg, long, depthg, timeg, inst_id = \
             read_glider_data_thredds_server(url_glider,var,scatter_plot,**kwargs)    
    
contour_plot='yes'    
depthg_gridded, varg_gridded, timegg = \
                    grid_glider_data_thredd(timeg,latg,long,depthg,varg,var_name,inst_id)                          

#%% glider_transect_model_comp

from read_glider_data import retrieve_glider_id_erddap_server
from glider_transect_model_com import glider_transect_model_com_erddap_server

# Servers location
url_glider = 'https://data.ioos.us/gliders/erddap'
url_model = 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_93.0/ts3z'

# MAB
lon_lim = [-75.0,-72.0]
lat_lim = [38.0,40.0]

# date limits
date_ini = '2018-09-01T00:00:00Z'
date_end = '2018-09-05T00:00:00Z'

# glider variable to retrieve
var_glider = 'temperature'

# model variable name
model_name = 'GOFS 3.1'
var_model = 'water_temp'

gliders = retrieve_glider_id_erddap_server(url_glider,lat_lim,lon_lim,date_ini,date_end)

for glid in gliders:
    dataset_id = glid
    print(glid)
    glider_transect_model_com_erddap_server(url_glider,dataset_id,url_model,lat_lim,lon_lim,\
                              date_ini,date_end,var_glider,var_model,model_name)


                    