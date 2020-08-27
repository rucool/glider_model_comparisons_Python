#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%

def read_glider_data_thredds_server(url_thredds,var_name,scatter_plot,**kwargs):
    
    """
    Created on Tue Feb  5 10:05:37 2019

    @author: aristizabal

    This function reads glider data from the IOOS Data Assembly Center (DAC).
    
    Inputs:
    url_thredds: url address or directory on local computer where the netcdf 
                file with the glider data resides. Example:
                url_glider = 'https://data.ioos.us/thredds/dodsC/deployments/rutgers/ru33-20180801T1323/ru33-20180801T1323.nc3.nc'
    var_name: variable to plot. Ex: 'temperature', 'salinity'. Make sure
        to use the same name as defined in the netcdf file
    scatter_plot: if equal to 'yes' then a scatter plot 
            of the glider transect is plotted
    kwargs = dict(date_ini='2018/09/01/00',date_end='2018/09/10/00')      
    date_ini: initial date the user wish to visualize the data. 
        This function uses the data format '%Y/%m/%d/%H'. 
        Examaple: date_ini = '2018/09/01/00'
        if it is not passed date_ini is the initial time of deployment
    date_end: final date the user wish to visualize the data. 
        This function uses the data format '%Y/%m/%d/%H'. 
        Examaple: date_ini = '2018/09/10/00' 
        if it is not passed date_ini is the initial time of deployment
                    
    Outputs:
    varg: all the glider profiles for the variable chosen within the user defined time window
    latg: latitude within the user defined time window
    long: longitude within the user defined time window
    timeg: user defined time window
    depthg: depth vector for all profiles
    """

    import xarray as xr
    import netCDF4
    import datetime
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    import cmocean
    
    date_ini = kwargs.get('date_ini', None)
    date_end = kwargs.get('date_end', None)

    gdata = xr.open_dataset(url_thredds,decode_times=False)
    
    dataset_id = gdata.id.split('_')[0]

    variable = np.asarray(gdata.variables[var_name][0][:])
    latitude = np.asarray(gdata.latitude[0])
    longitude = np.asarray(gdata.longitude[0])
    depth = np.asarray(gdata.depth[0])
    
    time = gdata.time[0]
    time = netCDF4.num2date(time,time.units)
    
    # Find time window of interest    
    if date_ini==None:
        tti = time[0]
    else:
        tti = datetime.datetime.strptime(date_ini,'%Y/%m/%d/%H')
        
    if date_end==None:
        tte = time[-1]
    else:
        tte = datetime.datetime.strptime(date_end,'%Y/%m/%d/%H')
        
    oktimeg = np.logical_and(time >= tti,time <= tte)
        
    # Fiels within time window
    varg =  variable[oktimeg,:].T
    latg = latitude[oktimeg]
    long = longitude[oktimeg]
    depthg = depth[oktimeg,:].T
    timeg = time[oktimeg]
    
    # Scatter plot
    if scatter_plot == 'yes':
        
        if var_name == 'temperature':
            color_map = cmocean.cm.thermal
            clabel = var_name[0].upper()+var_name[1:] + ' ($^oC$)'
        else:
            if var_name == 'salinity':
                color_map = cmocean.cm.haline
                clabel = var_name[0].upper()+var_name[1:]
            else:
                color_map = 'RdBu_r'
        
        timeg_matrix = np.tile(timeg.T,(depthg.shape[0],1))
        ttg = np.ravel(timeg_matrix)
        dg = np.ravel(depthg)
        teg = np.ravel(varg)

        kw = dict(c=teg, marker='*', edgecolor='none')

        fig, ax = plt.subplots(figsize=(10, 3))
        cs = ax.scatter(ttg,-dg,cmap=color_map,**kw)
        #fig.colorbar(cs)
        ax.set_xlim(timeg[0], timeg[-1])

        ax.set_ylabel('Depth (m)',fontsize=14)
        cbar = plt.colorbar(cs)
        cbar.ax.set_ylabel(clabel,fontsize=14)
        ax.set_title(dataset_id,fontsize=16)
        xfmt = mdates.DateFormatter('%H:%Mh\n%d-%b')
        ax.xaxis.set_major_formatter(xfmt)
        plt.ylim([-np.nanmax(dg),0])
    
    return varg, timeg, latg, long, depthg, dataset_id

#%%
    
def retrieve_dataset_id_erddap_server(url_erddap,lat_lim,lon_lim,date_ini,date_end):
    
    """
    Created on Tue Feb  5 10:05:37 2019

    @author: aristizabal

    This function retrieves glider ids from the IOOS 
    Data Assembly Center (DAC).
    
    Inputs:
    url_erddap: url address of erddap server
                Example: 'https://data.ioos.us/gliders/erddap'
    lat_lim: latitude limits for the search. 
            Example, lat_lim = [38.0,40.0]
    lon_lim: longitude limits for the search. 
            Example, lon_lim = [-75.0,-72.0]
    date_ini: initial date of time window. 
        This function accepts the data formats '%Y-%m-%d T %H:%M:%S Z' and '%Y/%m/%d/%H'. 
        Examaple: date_ini = '2018-08-02T00:00:00Z' 
    date_end: initial date of time window. 
        This function accepts the data formats '%Y-%m-%d T %H:%M:%S Z' and '%Y/%m/%d/%H'. 
        Examaple: date_ini = '2018-08-10T00:00:00Z'
                    
    Outputs:
    gliders: list of gliders ids that fall within the lat, lon and
             time constraints
    """
    
    from erddapy import ERDDAP
    import pandas as pd

    e = ERDDAP(server = url_erddap)

    # Search constraints
    kw = {
            'min_lon': lon_lim[0],
            'max_lon': lon_lim[1],
            'min_lat': lat_lim[0],
            'max_lat': lat_lim[1],
            'min_time': date_ini,
            'max_time': date_end,
            }

    search_url = e.get_search_url(response='csv', **kw)
    search = pd.read_csv(search_url)
    
    # Extract the IDs
    gliders = search['Dataset ID'].values
    
    return gliders

#%%
    
def read_glider_data_erddap_server(url_erddap,dataset_id,\
                                   lat_lim,lon_lim,scatter_plot,**kwargs):
    
    """
    Created on Tue Feb  5 10:05:37 2019

    @author: aristizabal

    This function reads glider data from the IOOS 
    Data Assembly Center (DAC).
    
    Inputs:
    url_erddap: url address of thredds server
                Example: 'https://data.ioos.us/gliders/erddap'
    dataset_id: this id is retrieved from the glider DAC using the
               function "retrieve_glider_id_erddap_server".
               Example: 'ru30-20180705T1825'
    lat_lim: latitude limits for the search. 
            Example, lat_lim = [38.0,40.0]
    lon_lim: longitude limits for the search. 
            Example, lon_lim = [-75.0,-72.0]
    date_ini: initial date of time window. 
        This function accepts the data formats '%Y-%m-%d T %H:%M:%S Z' and '%Y/%m/%d/%H'. 
        Examaple: date_ini = '2018-08-02T00:00:00Z' or '2018/08/02/00'
    date_end: initial date of time window. 
        This function uses the data format '%Y-%m-%d T %H:%M:%S Z'. 
        Examaple: date_ini = '2018-08-10T00:00:00Z' and '2018/08/10/00'
    scatter_plot: if equal to 'yes' then a scatter plot 
            of the glider transect is plotted
                    
    Outputs:
    tempg: all the glider profiles of temperature within the user defined time window
    saltg: all the glider profiles of salinity within the user defined time window
    latg: latitude within the user defined time window
    long: longitude within the user defined time window
    timeg: user defined time window
    depthg: depth vector for all profiles
    """
   
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
        
    # Coverting glider vectors into arrays
    timeg, ind = np.unique(df.index.values,return_index=True)
    latg = df['latitude (degrees_north)'].values[ind]
    long = df['longitude (degrees_east)'].values[ind]

    dg = df['depth (m)'].values
    vg1 = df[df.columns[3]].values
    vg2 = df[df.columns[4]].values

    zn = np.int(np.max(np.diff(np.hstack([ind,len(dg)]))))

    depthg = np.empty((zn,len(timeg)))
    depthg[:] = np.nan
    tempg = np.empty((zn,len(timeg)))
    tempg[:] = np.nan
    saltg = np.empty((zn,len(timeg)))
    saltg[:] = np.nan

    for i,ii in enumerate(ind):
        if i < len(timeg)-1:
            depthg[0:len(dg[ind[i]:ind[i+1]]),i] = dg[ind[i]:ind[i+1]]
            tempg[0:len(vg1[ind[i]:ind[i+1]]),i] = vg1[ind[i]:ind[i+1]]
            saltg[0:len(vg2[ind[i]:ind[i+1]]),i] = vg2[ind[i]:ind[i+1]]
        else:
            depthg[0:len(dg[ind[i]:len(dg)]),i] = dg[ind[i]:len(dg)]
            tempg[0:len(vg1[ind[i]:len(vg1)]),i] = vg1[ind[i]:len(vg1)]
            saltg[0:len(vg2[ind[i]:len(vg2)]),i] = vg2[ind[i]:len(vg2)]
                      
    # Scatter plot
    if scatter_plot == 'yes':
        
        color_map = cmocean.cm.thermal
        varg = tempg
        timeg_matrix = np.tile(timeg.T,(depthg.shape[0],1))
        ttg = np.ravel(timeg_matrix)
        dg = np.ravel(depthg)
        teg = np.ravel(varg)

        kw = dict(c=teg, marker='*', edgecolor='none')

        fig, ax = plt.subplots(figsize=(10, 3))
        cs = ax.scatter(ttg,-dg,cmap=color_map,**kw)
        #fig.colorbar(cs)
        ax.set_xlim(timeg[0], timeg[-1])

        ax.set_ylabel('Depth (m)',fontsize=14)
        cbar = plt.colorbar(cs)
        cbar.ax.set_ylabel('Temperature ($^oC$)',fontsize=14)
        ax.set_title(dataset_id,fontsize=16)
        xfmt = mdates.DateFormatter('%H:%Mh\n%d-%b')
        ax.xaxis.set_major_formatter(xfmt)
        plt.ylim([-np.nanmax(dg),0])
        
        color_map = cmocean.cm.haline
        varg = saltg
        timeg_matrix = np.tile(timeg.T,(depthg.shape[0],1))
        ttg = np.ravel(timeg_matrix)
        dg = np.ravel(depthg)
        teg = np.ravel(varg)

        kw = dict(c=teg, marker='*', edgecolor='none')

        fig, ax = plt.subplots(figsize=(10, 3))
        cs = ax.scatter(ttg,-dg,cmap=color_map,**kw)
        #fig.colorbar(cs)
        ax.set_xlim(timeg[0], timeg[-1])

        ax.set_ylabel('Depth (m)',fontsize=14)
        cbar = plt.colorbar(cs)
        cbar.ax.set_ylabel('Salinity',fontsize=14)
        ax.set_title(dataset_id,fontsize=16)
        xfmt = mdates.DateFormatter('%H:%Mh\n%d-%b')
        ax.xaxis.set_major_formatter(xfmt)
        plt.ylim([-np.nanmax(dg),0])
   
    return tempg, saltg, timeg, latg, long, depthg
    
     