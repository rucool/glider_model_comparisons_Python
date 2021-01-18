#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 12:59:28 2020

@author: aristizabal
"""

#%%

def read_glider_data_erddap_Rutgers_server(url_erddap,dataset_id,\
                                   lat_lim,lon_lim,scatter_plot,**kwargs):    
    
    from erddapy import ERDDAP
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    import cmocean
    import numpy as np
    
    date_ini = kwargs.get('date_ini', None)
    date_end = kwargs.get('date_end', None)
    
    # Find time window of interest    
    if np.logical_or(date_ini==None,date_end==None):
        constraints = {
            'latitude>=': lat_lim[0],
            'latitude<=': lat_lim[1],
            'longitude>=': lon_lim[0],
            'longitude<=': lon_lim[1],
            }
    else:
        constraints = {
            'time>=': date_ini,
            'time<=': date_end,
            'latitude>=': lat_lim[0],
            'latitude<=': lat_lim[1],
            'longitude>=': lon_lim[0],
            'longitude<=': lon_lim[1],
            }
     
    variables = [
            'depth',
            'latitude',
            'longitude',
            'time',
            'temperature',
            'salinity'
            ]
    
    e = ERDDAP(
            server=url_erddap,
            protocol='tabledap',
            response='nc'
            )
    
    e.dataset_id = dataset_id
    e.constraints = constraints
    e.variables = variables
        
    # Converting glider data to data frame
    # Cheching that data frame has data
    df = e.to_pandas()
    if len(df) != 0: 
    
        df = e.to_pandas(
            index_col='time (UTC)',
            parse_dates=True,
            skiprows=(1,)  # units information can be dropped.
            ).dropna()
        
    dg = df['depth (m)'].values
    tg = df.index.values
    vg1 = df[df.columns[3]].values
    vg2 = df[df.columns[4]].values
    
    upcast = np.where(np.diff(dg) < 0)[0]
    oku = np.where(np.diff(upcast)>1)[0]
    end_upcast = upcast[oku]
    
    downcast = np.where(np.diff(dg) > 0)[0]
    okd = np.where(np.diff(downcast)>1)[0]
    end_downcast = downcast[okd]
    
    ind = np.hstack([0,np.unique(np.hstack([end_upcast,end_downcast])),len(dg)])
    zn = np.max(np.diff(ind))
    
    depthg = np.empty((zn,len(ind)))
    depthg[:] = np.nan
    timeg = np.empty((zn,len(ind)))
    timeg[:] = np.nan
    tempg = np.empty((zn,len(ind)))
    tempg[:] = np.nan
    saltg = np.empty((zn,len(ind)))
    saltg[:] = np.nan
    
    for i in np.arange(len(ind)):
        if i == 0:
            indd = np.argsort(dg[ind[i]:ind[i+1]+2])
            depthg[0:len(dg[ind[i]:ind[i+1]+2]),i] = dg[ind[i]:ind[i+1]+2][indd]
            timeg[0:len(dg[ind[i]:ind[i+1]+2]),i] = mdates.date2num(tg[ind[i]:ind[i+1]+2][indd])
            tempg[0:len(vg1[ind[i]:ind[i+1]+2]),i] = vg1[ind[i]:ind[i+1]+2][indd]
            saltg[0:len(vg2[ind[i]:ind[i+1]+2]),i] = vg2[ind[i]:ind[i+1]+2][indd]
        if i < len(ind)-1:
            indd = np.argsort(dg[ind[i]+1:ind[i+1]+2])
            depthg[0:len(dg[ind[i]+1:ind[i+1]+2]),i] = dg[ind[i]+1:ind[i+1]+2][indd]
            timeg[0:len(dg[ind[i]+1:ind[i+1]+2]),i] = mdates.date2num(tg[ind[i]+1:ind[i+1]+2][indd])
            tempg[0:len(vg1[ind[i]+1:ind[i+1]+2]),i] = vg1[ind[i]+1:ind[i+1]+2][indd]
            saltg[0:len(vg2[ind[i]+1:ind[i+1]+2]),i] = vg2[ind[i]+1:ind[i+1]+2][indd]
        else:
            indd = np.argsort(dg[ind[i]+1:len(dg)])
            depthg[0:len(dg[ind[i]+1:len(dg)]),i] = dg[ind[i]+1:len(dg)][indd]
            timeg[0:len(dg[ind[i]+1:len(dg)]),i] = mdates.date2num(tg[ind[i]+1:len(dg)][indd])
            tempg[0:len(vg1[ind[i]+1:len(vg1)]),i] = vg1[ind[i]+1:len(vg1)][indd]
            saltg[0:len(vg2[ind[i]+1:len(vg2)]),i] = vg2[ind[i]+1:len(vg2)][indd]   
    
    # Scatter plot
    if scatter_plot == 'yes':
        
        color_map = cmocean.cm.thermal
        varg = tempg
        #timeg_matrix = np.tile(timeg.T,(depthg.shape[0],1))
        ttg = np.ravel(timeg)
        dg = np.ravel(depthg)
        teg = np.ravel(varg)
    
        kw = dict(c=teg, marker='*', edgecolor='none')
    
        fig, ax = plt.subplots(figsize=(10, 3))
        cs = ax.scatter(ttg,-dg,cmap=color_map,**kw)
        #fig.colorbar(cs)
        ax.set_xlim(np.nanmin(ttg), np.nanmax(ttg))
    
        ax.set_ylabel('Depth (m)',fontsize=14)
        cbar = plt.colorbar(cs)
        cbar.ax.set_ylabel('Temperature ($^oC$)',fontsize=14)
        ax.set_title(dataset_id,fontsize=16)
        xfmt = mdates.DateFormatter('%H:%Mh\n%d-%b')
        ax.xaxis.set_major_formatter(xfmt)
        plt.ylim([-np.nanmax(dg),0])
        
        color_map = cmocean.cm.haline
        varg = saltg
        #timeg_matrix = np.tile(timeg.T,(depthg.shape[0],1))
        ttg = np.ravel(timeg)
        dg = np.ravel(depthg)
        teg = np.ravel(varg)
    
        kw = dict(c=teg, marker='*', edgecolor='none')
    
        fig, ax = plt.subplots(figsize=(10, 3))
        cs = ax.scatter(ttg,-dg,cmap=color_map,**kw)
        #fig.colorbar(cs)
        ax.set_xlim(np.nanmin(ttg), np.nanmax(ttg))
    
        ax.set_ylabel('Depth (m)',fontsize=14)
        cbar = plt.colorbar(cs)
        cbar.ax.set_ylabel('Salinity',fontsize=14)
        ax.set_title(dataset_id,fontsize=16)
        xfmt = mdates.DateFormatter('%H:%Mh\n%d-%b')
        ax.xaxis.set_major_formatter(xfmt)
        plt.ylim([-np.nanmax(dg),0])
        
    return tempg, saltg, timeg, latg, long, depthg
        

#%% Cell #4: Search for glider data sets given 
#   a latitude and longitude box and time window 
          
from read_glider_data import retrieve_dataset_id_erddap_server
#from read_glider_data import read_glider_data_erddap_server

# Server location
url_erddap = 'http://slocum-data.marine.rutgers.edu/erddap'

# Caribbean
lon_lim = [-80,-60.0]
lat_lim = [10.0,30.0]

# date limits
date_ini = '2020/10/01/00'
date_end = '2020/10/06/00'

gliders = retrieve_dataset_id_erddap_server(url_erddap,lat_lim,lon_lim,date_ini,date_end)
print('The gliders found are ')
print(gliders)

dataset_id = gliders[2]

kwargs = dict(date_ini=date_ini,date_end=date_end)
scatter_plot = 'yes'

tempg, saltg, timeg, latg, long, depthg = read_glider_data_erddap_Rutgers_server(url_erddap,dataset_id,\
                                   lat_lim,lon_lim,scatter_plot,**kwargs)

#%% 
import numpy as np    
from process_glider_data import grid_glider_data    
    
# variable to retrieve
var_name = 'temperature'
    
contour_plot = 'yes' # default value is 'yes'
delta_z = 0.4     # default value is 0.3

okt = np.isfinite(np.nanmean(timeg,axis=0))
timeg_vec = np.nanmean(timeg,axis=0)[okt]
tempgg = tempg[:,okt]
depthgg = depthg[:,okt]

tempg_gridded, timegg, depthg_gridded = \
                    grid_glider_data(var_name,dataset_id,tempgg,timeg_vec,latg,long,depthgg,delta_z,contour_plot)
                    