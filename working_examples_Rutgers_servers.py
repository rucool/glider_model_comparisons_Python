#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 11:07:45 2020

@author: aristizabal
"""
#%% Cell #1: Search for glider data sets given
#   a latitude and longitude box and time window

from read_glider_data import retrieve_dataset_id_erddap_server

url_erddap = "http://slocum-data.marine.rutgers.edu/erddap"

# Caribbean
lon_lim = [-80,-60.0]
lat_lim = [10.0,30.0]

# date limits
date_ini = '2020/10/17 00:00:00'
date_end = '2020/10/18 16:00:00'

gliders = retrieve_dataset_id_erddap_server(url_erddap,lat_lim,lon_lim,date_ini,date_end)

#%% Cell #2: Given a glider data set, Search for available
# variables within that data set

from read_glider_data import retrieve_variable_names_erddap_server

url_erddap = "http://slocum-data.marine.rutgers.edu/erddap"
#dataset_id = [idg for idg in gliders if idg.split('-')[3] == 'sci'][0]
dataset_id = 'ru29-20200908T1623-profile-sci-rt'

variable_list = retrieve_variable_names_erddap_server(url_erddap,dataset_id)

#%% Cell #3: Read glider data given a dataset_id, a list of variables,
# latitude and logitude limits, a time window (optional).

from read_glider_data import read_glider_variables_erddap_server

url_erddap = "http://slocum-data.marine.rutgers.edu/erddap"

# Caribbean
lon_lim = [-80,-60.0]
lat_lim = [10.0,30.0]

# date limits
date_ini = '2020/10/17 00:00:00'
date_end = '2020/10/18 16:00:00'

#dataset_id = [idg for idg in gliders if idg.split('-')[3] == 'sci'][0]
dataset_id = 'ru29-20200908T1623-profile-sci-rt'

kwargs = dict(date_ini=date_ini,date_end=date_end)

variable_names = [
            'depth',
            'latitude',
            'longitude',
            'time',
            'temperature',
            'u',
            'v'
            ]

df = read_glider_variables_erddap_server(url_erddap,dataset_id,lat_lim,lon_lim,\
                                   variable_names,**kwargs)

#%% Cell #4: Given the dataframe df from cell #3, convert the requested the variables into arrays
# according to the upcast and downcast recorded in the depth variable

from process_glider_data import find_profiles
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime

# depth
indp = [ind for ind,var in enumerate(variable_names) if var == 'depth'][0]
ddg = df[df.columns[indp]].values
okdepth = np.isfinite(ddg)
dg = ddg[okdepth]
depthg = find_profiles(dg,dg)

# time
indt = [ind for ind,var in enumerate(variable_names) if var == 'time'][0]
ttg = df[df.columns[indt]].values[okdepth]
tg = np.asarray([mdates.date2num(datetime.strptime(tt,'%Y-%m-%dT%H:%M:%SZ')) for tt in ttg])
timeg = find_profiles(dg,tg)

# latitude
indlt = [ind for ind,var in enumerate(variable_names) if var == 'latitude'][0]
ltg = df[df.columns[indlt]].values[okdepth]
latg = find_profiles(dg,ltg)

# longitude
indln = [ind for ind,var in enumerate(variable_names) if var == 'longitude'][0]
lng = df[df.columns[indln]].values[okdepth]
long = find_profiles(dg,lng)

# u velocity
indu = [ind for ind,var in enumerate(variable_names) if var == 'u'][0]
uug = df[df.columns[indu]].values[okdepth]
ug = find_profiles(dg,uug)

# v velocity
indv = [ind for ind,var in enumerate(variable_names) if var == 'v'][0]
vvg = df[df.columns[indv]].values[okdepth]
vg = find_profiles(dg,vvg)

# temperature
indt = [ind for ind,var in enumerate(variable_names) if var == 'temperature'][0]
temg = df[df.columns[indt]].values[okdepth]
tempg = find_profiles(dg,temg)

#%% cell #5: Extract a velocity transect from GOFS 3.1 given
# times, latitudes and longitudes of the transect

from glider_transect_model_com import get_transect_from_GOFS
import numpy as np
import matplotlib.dates as mdates

# Server location
url_model = 'http://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/uv3z'

# model variable name
model_name = 'GOFS 3.1'
var_name_model = 'water_u'

contour_plot = 'yes' # default value is 'yes'

tt = np.nanmean(timeg,axis=0)
okt = np.isfinite(tt)
time = mdates.num2date(tt[okt])
lat = np.nanmean(latg,axis=0)[okt]
lon = np.nanmean(long,axis=0)[okt]

var_model, ttmodel, depth_model, lat_model, lon_model = \
    get_transect_from_GOFS(url_model,var_name_model,model_name,time,lat,lon,contour_plot='yes')

#%% cell #6: Extract a temperature transect from GOFS 3.1 given
# times, latitudes and longitudes of the transect

from glider_transect_model_com import get_transect_from_GOFS
import numpy as np
import matplotlib.dates as mdates

# Server location
url_model = 'http://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/ts3z'

# model variable name
model_name = 'GOFS 3.1'
var_name_model = 'water_temp'

contour_plot = 'yes' # default value is 'yes'

tt = np.nanmean(timeg,axis=0)
okt = np.isfinite(tt)
time = mdates.num2date(tt[okt])
lat = np.nanmean(latg,axis=0)[okt]
lon = np.nanmean(long,axis=0)[okt]

var_model, ttmodel, depth_model, lat_model, lon_model = \
    get_transect_from_GOFS(url_model,var_name_model,model_name,time,lat,lon,contour_plot='yes')
