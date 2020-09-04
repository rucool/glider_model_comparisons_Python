#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 14:08:27 2019

@author: aristizabal
"""

#%%
def get_glider_transect_from_GOFS(url_model,var_name_model,model_name,varg,timeg,latg,long,depthg,contour_plot='yes'):

    """
    Created on April 23 2020

    @author: aristizabal

    This function finds the corresponding glider transect in
    the Global Ocean Forecasting System (GOFS 3.1)

    Inputs:
    url_model: url address of model output.
              example: 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_93.0/ts3z'
    var_name_model: name of variable in model output. Example "water_temp"
    model_name: short name of model. Just use to title plots
    varg: glider profiles matrix for the variable chosen
    latg: glider latitude vector
    long: glider longitude vector
    timeg: glider time vector
    depthg: depth matrix for all profiles

    Outputs:
    var_model
    time_model
    depth_model
    lat_model
    lon_model

    """

    import xarray as xr
    import netCDF4
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from datetime import datetime
    import cmocean

    print('Retrieving coordinates from model')
    model = xr.open_dataset(url_model,decode_times=False)

    lat_model = np.asarray(model.lat[:])
    lon_model = np.asarray(model.lon[:])
    depth_model = np.asarray(model.depth[:])
    ttm = model.time
    tm = netCDF4.num2date(ttm[:],ttm.units)

    tini = mdates.num2date(mdates.date2num(timeg[0]))
    tend = mdates.num2date(mdates.date2num(timeg[-1]))
    #tini = timeg[0]
    #tend = timeg[-1]

    oktimem = np.where(np.logical_and(tm >= tini,tm <= tend))
    
    time_model = tm[oktimem]

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
    ttmodel = np.asarray([datetime(time_model[i].year,time_model[i].month,time_model[i].day,\
                    time_model[i].hour) for i in np.arange(len(time_model))])

    tstamp_model = [mdates.date2num(ttmodel[i]) for i in np.arange(len(ttmodel))]

    # interpolating glider lon and lat to lat and lon on model time
    sublonm=np.interp(tstamp_model,tstamp_glider,target_lon)
    sublatm=np.interp(tstamp_model,tstamp_glider,target_lat)

    # getting the model grid positions for sublonm and sublatm
    oklonm=np.round(np.interp(sublonm,lon_model,np.arange(len(lon_model)))).astype(int)
    oklatm=np.round(np.interp(sublatm,lat_model,np.arange(len(lat_model)))).astype(int)

    # Getting glider transect from model
    print('Getting glider transect from '+ model_name)
    var_model = np.empty((len(depth_model),len(oktimem[0])))
    var_model[:] = np.nan
    for i in range(len(oktimem[0])):
        print(len(oktimem[0]),' ',i)
        var_model[:,i] = model.variables[var_name_model][oktimem[0][i],:,oklatm[i],oklonm[i]]

    # Countour plot
    if contour_plot == 'yes':

        if var_name_model == 'water_temp':
            color_map = cmocean.cm.thermal
            clabel = var_name_model[0].upper()+var_name_model[1:] + ' ($^oC$)'
        else:
            if var_name_model == 'salinity':
                color_map = cmocean.cm.haline
                clabel = var_name_model[0].upper()+var_name_model[1:]
            else:
                color_map = 'RdBu_r'

        okg = depthg <= np.nanmax(depthg)
        okm = depth_model <= np.nanmax(depthg)
        min_val = np.int(np.floor(np.min([np.nanmin(varg[okg]),np.nanmin(var_model[okm])])))
        max_val = np.int(np.ceil(np.max([np.nanmax(varg[okg]),np.nanmax(var_model[okm])])))
    
        if var_name_model == 'salinity':
            kw = dict(levels = np.arange(min_val,max_val+0.25,0.25))
        else:
            nlevels = max_val - min_val + 1
            kw = dict(levels = np.linspace(min_val,max_val,nlevels))
        
        fig, ax = plt.subplots(figsize=(10, 3))
        cs = plt.contourf(ttmodel,-depth_model,var_model,cmap=color_map,**kw)
        plt.contour(ttmodel,-depth_model,var_model,[26],colors='k')
        cs = fig.colorbar(cs, orientation='vertical')
        cs.ax.set_ylabel(clabel,fontsize=14,labelpad=15)
    
        ax.set_xlim(timeg[0], timeg[-1])
        ax.set_ylim(-np.nanmax(depthg), 0)
        ax.set_ylabel('Depth (m)',fontsize=14)
        xfmt = mdates.DateFormatter('%H:%Mh\n%d-%b')
        ax.xaxis.set_major_formatter(xfmt)
        plt.title('Along Track ' + var_name_model[0].upper() + var_name_model[1:] +\
                  ' Profile ' + model_name,fontsize=16)

    return var_model, ttmodel, depth_model, lat_model, lon_model

#%%

def get_glider_transect_from_Amseas(url_amseas,var_name_model,model_name,\
                                        varg,timeg,latg,long,depthg,contour_plot='yes'):

    from datetime import datetime, timedelta
    import numpy as np
    import xarray as xr
    import netCDF4
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    import cmocean

    date_ini = str(timeg[0])
    date_end = str(timeg[-1])

    yini = int(date_ini.split('-')[0])
    mini = int(date_ini.split('-')[1])
    dini = int(date_ini.split('-')[2].split('T')[0])
    #hini = int(date_ini.split('-')[2].split('T')[1][0:2])

    yend = int(date_end.split('-')[0])
    mend = int(date_end.split('-')[1])
    dend = int(date_end.split('-')[2].split('T')[0])
    #hend = int(date_end.split('-')[2].split('T')[1][0:2])

    tini = datetime(yini,mini,dini,0)
    tend = datetime(yend,mend,dend+1,0)

    yini = str(tini).split('-')[0]
    mini = str(tini).split('-')[1]
    dini = str(tini).split('-')[2].split(' ')[0]
    hini = str(tini).split('-')[2].split(' ')[1][0:2]

    n_files = (tend-tini).days *8
    time_vec = [tini + n*timedelta(hours=3) for n in np.arange(n_files)]

    file_amseas0 = url_amseas + yini + mini + dini + \
            '/ncom_relo_amseas_u_' + yini + mini + dini + '00' + '_t0' + hini + '.nc'

    amseas = xr.open_dataset(file_amseas0,decode_times=False)

    print('Getting depth, lat and lon from Amseas')
    lat_amseas = np.asarray(amseas.lat[:])
    lon_amseas = np.asarray(amseas.lon[:])
    depth_amseas = np.asarray(amseas.depth[:])

    print('Getting time from Amseas')
    tm = []
    for t,tt in enumerate(time_vec):
        print(tt)
        y = str(tt).split('-')[0]
        m = str(tt).split('-')[1]
        d = str(tt).split('-')[2].split(' ')[0]
        h = str(tt).split('-')[2].split(' ')[1][0:2]

        file_amseas = url_amseas + y + m + d + \
            '/ncom_relo_amseas_u_' + y + m + d + '00' + '_t0' + h + '.nc'

        amseas = xr.open_dataset(file_amseas,decode_times=False)

        ttm = amseas.time
        tm.append(datetime.strptime(str(netCDF4.num2date(ttm[:],ttm.units)[0]),'%Y-%m-%d %H:%M:%S'))
    time_amseas = np.asarray(tm)

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
    tstamp_model = [mdates.date2num(time_amseas[i]) for i in np.arange(len(time_amseas))]

    # interpolating glider lon and lat to lat and lon on model time
    sublon_amseas = np.interp(tstamp_model,tstamp_glider,target_lon)
    sublat_amseas = np.interp(tstamp_model,tstamp_glider,target_lat)

    # getting the model grid positions for sublonm and sublatm
    oklon_amseas = np.round(np.interp(sublon_amseas,lon_amseas,np.arange(len(lon_amseas)))).astype(int)
    oklat_amseas = np.round(np.interp(sublat_amseas,lat_amseas,np.arange(len(lat_amseas)))).astype(int)

    var_amseas = np.empty((len(depth_amseas),len(time_vec)))
    var_amseas[:] = np.nan
    for t,tt in enumerate(time_vec):
        print(tt)
        y = str(tt).split('-')[0]
        m = str(tt).split('-')[1]
        d = str(tt).split('-')[2].split(' ')[0]
        h = str(tt).split('-')[2].split(' ')[1][0:2]

        file_amseas = url_amseas + y + m + d + \
            '/ncom_relo_amseas_u_' + y + m + d + '00' + '_t0' + h + '.nc'

        amseas = xr.open_dataset(file_amseas,decode_times=False)

        var_amseas[:,t] = amseas.variables[var_name_model][0,:,oklat_amseas[t],oklon_amseas[t]]

    # Countour plot
    if contour_plot == 'yes':

        if var_name_model == 'water_temp':
            color_map = cmocean.cm.thermal
            clabel = var_name_model[0].upper()+var_name_model[1:] + ' ($^oC$)'
        else:
            if var_name_model == 'salinity':
                color_map = cmocean.cm.haline
                clabel = var_name_model[0].upper()+var_name_model[1:]
            else:
                color_map = 'RdBu_r'

        okg = depthg <= np.nanmax(depthg)
        okm = depth_amseas <= np.nanmax(depthg)
        min_val = np.int(np.floor(np.min([np.nanmin(varg[okg]),np.nanmin(var_amseas[okm])])))
        max_val = np.int(np.ceil(np.max([np.nanmax(varg[okg]),np.nanmax(var_amseas[okm])])))

        if var_name_model == 'salinity':
            kw = dict(levels = np.arange(min_val,max_val+0.25,0.25))
        else:
            nlevels = max_val - min_val + 1
            kw = dict(levels = np.linspace(min_val,max_val,nlevels))

        fig, ax = plt.subplots(figsize=(10, 3))
        cs = plt.contourf(mdates.date2num(time_amseas),-depth_amseas,var_amseas,cmap=color_map,**kw)
        plt.contour(mdates.date2num(time_amseas),-depth_amseas,var_amseas,[26],colors='k')
        cs = fig.colorbar(cs, orientation='vertical')
        cs.ax.set_ylabel(clabel,fontsize=14,labelpad=15)
    
        ax.set_xlim(timeg[0], timeg[-1])
        ax.set_ylim(-np.nanmax(depthg), 0)
        ax.set_ylabel('Depth (m)',fontsize=14)
        xfmt = mdates.DateFormatter('%H:%Mh\n%d-%b')
        ax.xaxis.set_major_formatter(xfmt)
        plt.title('Along Track ' + var_name_model[0].upper() + var_name_model[1:] +\
                  ' Profile ' + model_name,fontsize=16)

    return var_amseas, time_amseas, depth_amseas, lat_amseas, lon_amseas
