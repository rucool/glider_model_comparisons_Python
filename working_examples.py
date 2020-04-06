#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 10:21:51 2019

@author: aristizabal
"""

#%% Cell #1: reads glider data from the IOOS thredds server and plots a    
#   scatter plot of a glider transect for the entire length of the deployment

from read_glider_data import read_glider_data_thredds_server

#url_thredds = 'https://data.ioos.us/thredds/dodsC/deployments/rutgers/ru33-20180801T1323/ru33-20180801T1323.nc3.nc'
url_thredds = 'http://gliders.ioos.us/thredds/dodsC/deployments/aoml/SG668-20190819T1217/SG668-20190819T1217.nc3.nc'
var_name = 'temperature'
#var = 'salinity'
scatter_plot = 'yes'

varg, latg, long, depthg, timeg, inst_id = \
             read_glider_data_thredds_server(url_thredds,var_name,scatter_plot)             
             

#%% Cell #2: reads glider data from the IOOS thredds server and plots a    
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

varg, latg, long, depthg, timeg, inst_id = \
             read_glider_data_thredds_server(url_thredds,var_name,scatter_plot,**kwargs)

#%% Cell #3: same as cell #2, in addition to interpolating the variable 
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

varg, latg, long, depthg, timeg, inst_id = \
             read_glider_data_thredds_server(url_thredds,var_name,scatter_plot,**kwargs)    

delta_z = 0.4 # bin size in the vertical when gridding the variable vertical profile 
              # default value is 0.3   
contour_plot = 'yes'  # default value is 'yes'   
depthg_gridded, varg_gridded, timegg = \
                    grid_glider_data(timeg,latg,long,depthg,varg,var_name,inst_id,delta_z,contour_plot)   

#%% Cell #4: search for glider data sets given 
#   a latitude and longitude box and time window 
          
from read_glider_data import retrieve_glider_id_erddap_server

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

gliders = retrieve_glider_id_erddap_server(url_erddap,lat_lim,lon_lim,date_ini,date_end)
print('The gliders found are ')
print(gliders)
             
#%% Cell #5: search for glider data sets given a 
#   latitude and longitude box and time window, choose one those data sets 
#   (glider_id), plot a scatter plot of the chosen glider transect, grid 
#   and plot a contour plot of the chosen glider transect 
         
from read_glider_data import retrieve_glider_id_erddap_server
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

gliders = retrieve_glider_id_erddap_server(url_erddap,lat_lim,lon_lim,date_ini,date_end)

dataset_id = gliders[10]

# variable to retrieve
var = 'temperature'
#var = 'salinity'

kwargs = dict(date_ini=date_ini,date_end=date_end)
scatter_plot = 'yes'

varg, latg, long, depthg, timeg = read_glider_data_erddap_server(url_erddap,dataset_id,var,\
                                   lat_lim,lon_lim,scatter_plot,**kwargs)
    
#varg, latg, long, depthg, timeg = read_glider_data_erddap_server(url_erddap,dataset_id,var,\
#                                   lat_lim,lon_lim,scatter_plot)    

contour_plot = 'yes' # default value is 'yes'
delta_z = 0.4     # default value is 0.3

depthg_gridded, varg_gridded, timegg = \
                    grid_glider_data(timeg,latg,long,depthg,varg,var_name,inst_id,delta_z,contour_plot)
                                                   
                       
#%% cell #6: search for glider data sets given a 
#    latitude and longitude box and time window, choose one those data sets 
#    (glider_id), grid in the vertical the glider transect, get the glider
#    transect in the GOFS 3.1 grid, and plot both the transect from the glider
#    deployment and GOFS 3.1 output

from read_glider_data import retrieve_glider_id_erddap_server
from glider_transect_model_com import glider_transect_erddap_server_vs_model

# Servers location
url_erddap = 'https://data.ioos.us/gliders/erddap'
url_model = 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_93.0/ts3z'

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


# glider variable to retrieve
var_name_glider = 'temperature'
#var_glider = 'salinity'
delta_z = 0.4 # bin size in the vertical when gridding the variable vertical profile 
              # default value is 0.3  

# model variable name
model_name = 'GOFS 3.1'
var_name_model = 'water_temp'
#var_model = 'salinity'

gliders = retrieve_glider_id_erddap_server(url_erddap,lat_lim,lon_lim,date_ini,date_end)

'''
for glid in gliders:
    dataset_id = glid
    print(glid)
    timeg,depthg_gridded,varg_gridded,timem,depthm,target_varm = \
    glider_transect_model_com_erddap_server(url_glider,dataset_id,url_model,\
                              lat_lim,lon_lim,\
                              date_ini,date_end,var_glider,var_model,model_name,delta_z=0.4)
''' 
                     
dataset_id = gliders[10]
kwargs = dict(date_ini=date_ini,date_end=date_end)
print(dataset_id)
timeg,depthg_gridded,varg_gridded,timem,depthm,target_varm = \
glider_transect_erddap_server_vs_model(url_erddap,dataset_id,url_model,\
            lat_lim,lon_lim,var_name_glider,var_name_model,model_name,delta_z=0.4,**kwargs) 
 
'''    
# Run "glider_transect_erddap_server_vs_model" without the "kwargs" argument
# if you want to get the entire deployment    
dataset_id = gliders[3]    
timeg,depthg_gridded,varg_gridded,timem,depthm,target_varm = \
glider_transect_erddap_server_vs_model(url_erddap,dataset_id,url_model,\
            lat_lim,lon_lim,var_name_glider,var_name_model,model_name,delta_z=0.4) 
'''
    
#%% cell #7: Read glider data from a specific glider from the thredds server, reads the same 
# glider transect in the GOFS 3.1 grid, and plot both the transect from the glider
# deployment and GOFS 3.1 output
    
from glider_transect_model_com import glider_transect_thredds_server_vs_model

url_thredds = 'http://gliders.ioos.us/thredds/dodsC/deployments/aoml/SG668-20190819T1217/SG668-20190819T1217.nc3.nc'
url_model = 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_93.0/ts3z'

# Caribbean
lon_lim = [-80,-60.0]
lat_lim = [10.0,30.0]

# date limits
date_ini = '2019/09/01/00'
date_end = '2019/09/10/00'

# glider variable to retrieve
var_name_glider = 'temperature'
#var_name_glider = 'salinity'
delta_z = 0.4 # bin size in the vertical when gridding the variable vertical profile 
              # default value is 0.3  

# model variable name
model_name = 'GOFS 3.1'
var_name_model = 'water_temp'
#var_name_model = 'salinity'
    
kwargs = dict(date_ini=date_ini,date_end=date_end)

timeg,depthg_gridded,varg_gridded,timem,depthm,target_varm = \
    glider_transect_thredds_server_vs_model(url_thredds,var_name_glider,\
                                url_model,var_name_model,model_name,delta_z=0.4,**kwargs)

'''
# Run "glider_transect_erddap_server_vs_model" without the "kwargs" argument
# if you want to get the entire deployment          
timeg,depthg_gridded,varg_gridded,timem,depthm,target_varm = \
    glider_transect_thredds_server_vs_model(url_thredds,var_name_glider,\
                                url_model,var_name_model,model_name,delta_z=0.4)        
'''       