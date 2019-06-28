#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 14:08:27 2019

@author: aristizabal
"""

def glider_transect_model_com_erddap_server(url_glider,dataset_id,url_model,lat_lim,lon_lim,\
                              date_ini,date_end,var_glider,var_model,model_name):
 
    """
    Created on Wed Feb  6 11:49:24 2019

    @author: aristizabal
    
    This function reads a glider transect and find the corresponding glider transect in
    an ocean model (example GOFS 3.1)
    
    Inputs:
    url_glider: url address of glider erddap server
                Example: 'https://data.ioos.us/gliders/erddap'
    dataset_id: id of glider dataset from "retrieve_glider_id_erddap_server"
    url_model: url address of model output.
              example: 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_93.0/ts3z'
    lat_lim: latitude limits for the search. 
            Example, lat_lim = [38.0,40.0]
    lon_lim: longitude limits for the search. 
            Example, lon_lim = [-75.0,-72.0]
    date_ini: initial date of time window. 
        This function uses the data format '%Y-%m-%d T %H:%M:%S Z'. 
        Examaple: date_ini = '2018-08-02T00:00:00Z' 
    date_end: initial date of time window. 
        This function uses the data format '%Y-%m-%d T %H:%M:%S Z'. 
        Examaple: date_ini = '2018-08-10T00:00:00Z'
    var_glider: name of variable in glider dataset. Example "temperature"
    var_model: name of variable in model output. Example "water_temp"
    model_name: short name of model. Just use to title plots
                    
    Outputs:
    timeg: glider time vector. it represents the average time of each 
          down or upcast
    depthg_gridded: gridded glider depth vector
    varg_gridded: gridded glider variable matrix
    timem: model time vector
    depthm: model depth vector
    target_varm: model variable matrix
    
    """
    
    import xarray as xr
    import netCDF4
    import numpy as np
    import datetime
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    import cmocean

    from read_glider_data import read_glider_data_erddap_server
    from process_glider_data import grid_glider_data_erddap

    # Read GOFS 3.1 output
    print('Retrieving coordinates from model')
    model = xr.open_dataset(url_model,decode_times=False)
    
    latm = model.lat[:]
    lonm = model.lon[:]
    depthm = model.depth[:]
    ttm = model.time
    tm = netCDF4.num2date(ttm[:],ttm.units) 

    tmin = datetime.datetime.strptime(date_ini,'%Y-%m-%dT%H:%M:%SZ')
    tmax = datetime.datetime.strptime(date_end,'%Y-%m-%dT%H:%M:%SZ')

    oktimem = np.where(np.logical_and(tm >= tmin, tm <= tmax))
    
    timem = tm[oktimem]
    
    # Read and process glider data
    print('Reading glider data')
    df = read_glider_data_erddap_server(url_glider,dataset_id,var_glider,\
                                        lat_lim,lon_lim,date_ini,date_end,\
                                        scatter_plot='no')

    depthg_gridded, varg_gridded, timeg, latg, long = \
                       grid_glider_data_erddap(df,var_glider,contour_plot='no')

    # Conversion from glider longitude and latitude to GOFS convention
    target_lon = np.empty((len(long),))
    target_lon[:] = np.nan
    for i,ii in enumerate(long):
        if ii < 0: 
            target_lon[i] = 360 + ii
        else:
            target_lon[i] = ii
    target_lat = latg

    # Changing times to timestamp
    tstamp_glider = [mdates.date2num(timeg[i]) for i in np.arange(len(timeg))]
    tstamp_model = [mdates.date2num(timem[i]) for i in np.arange(len(timem))]

    # interpolating glider lon and lat to lat and lon on model time
    sublonm=np.interp(tstamp_model,tstamp_glider,target_lon)
    sublatm=np.interp(tstamp_model,tstamp_glider,target_lat)

    # getting the model grid positions for sublonm and sublatm
    oklonm=np.round(np.interp(sublonm,lonm,np.arange(len(lonm)))).astype(int)
    oklatm=np.round(np.interp(sublatm,latm,np.arange(len(latm)))).astype(int)
    
    # Getting glider transect from model
    print('Getting glider transect from model. If it breaks is because GOFS 3.1 server is not responding')
    target_varm = np.empty((len(depthm),len(oktimem[0])))
    target_varm[:] = np.nan
    for i in range(len(oktimem[0])):
        print(len(oktimem[0]),' ',i)
        target_varm[:,i] = model.variables[var_model][oktimem[0][i],:,oklatm[i],oklonm[i]]

    target_varm[target_varm < -100] = np.nan

    # plot
    if var_glider == 'temperature':
        color_map = cmocean.cm.thermal
    else:
        if var_glider == 'salinity':
            color_map = cmocean.cm.haline
        else:
            color_map = 'RdBu_r'
        
    okg = depthg_gridded <= np.max(depthg_gridded) 
    okm = target_varm <= np.max(depthg_gridded) 
    min_val = np.round(np.min([np.nanmin(varg_gridded[okg]),np.nanmin(target_varm[okm])]))
    max_val = np.round(np.max([np.nanmax(varg_gridded[okg]),np.nanmax(target_varm[okm])]))
    
    if var_glider == 'salinity':
        kw = dict(levels = np.arange(min_val,max_val,0.25))
    else:
        nlevels = max_val - min_val + 1
        kw = dict(levels = np.linspace(min_val,max_val,nlevels))
    

    # plot
    fig, ax = plt.subplots(figsize=(12, 6))

    ax = plt.subplot(211)        
    #plt.contour(timeg,-depthg_gridded,varg_gridded,colors = 'lightgrey',**kw)
    cs = plt.contourf(timeg,-depthg_gridded,varg_gridded,cmap=color_map,**kw)
    plt.contour(timeg,-depthg_gridded,varg_gridded,[26],colors='k')

    cs = fig.colorbar(cs, orientation='vertical') 
    cs.ax.set_ylabel(var_glider[0].upper()+var_glider[1:],fontsize=14,labelpad=15)
        
    ax.set_xlim(df.index[0], df.index[-1])
    ax.set_ylim(-np.max(depthg_gridded), 0)
    ax.set_ylabel('Depth (m)',fontsize=14)
    ax.set_xticklabels(' ')
    
    
    plt.title('Along Track ' + var_glider[0].upper() + var_glider[1:] + ' Profile ' + dataset_id.split('-')[0])
    
    ax = plt.subplot(212)        
    #plt.contour(timeg,-depthg_gridded,varg_gridded,colors = 'lightgrey',**kw)
    cs = plt.contourf(timem,-depthm,target_varm,cmap=color_map,**kw)

    cs = fig.colorbar(cs, orientation='vertical') 
    cs.ax.set_ylabel(var_glider[0].upper()+var_glider[1:],fontsize=14,labelpad=15)

    ax.set_xlim(df.index[0], df.index[-1])
    ax.set_ylim(-np.max(depthg_gridded), 0)
    ax.set_ylabel('Depth (m)',fontsize=14)
    xfmt = mdates.DateFormatter('%H:%Mh\n%d-%b')
    ax.xaxis.set_major_formatter(xfmt)

    plt.title('Along Track ' + var_glider[0].upper() + var_glider[1:] + ' Profile ' + model_name)     

    return timeg,depthg_gridded,varg_gridded,timem,depthm,target_varm