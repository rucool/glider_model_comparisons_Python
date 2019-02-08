#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%

def read_glider_data_thredds_server(url_glider,var,scatter_plot,**kwargs):
    
    """
    Created on Tue Feb  5 10:05:37 2019

    @author: aristizabal

    This function reads glider data from the IOOS Data Assembly Center (DAC).
    
    Inputs:
    url_glider: url address or directory on local computer where the netcdf 
                file with the glider data resides. Example:
                url_glider = 'https://data.ioos.us/thredds/dodsC/deployments/rutgers/ru33-20180801T1323/ru33-20180801T1323.nc3.nc'
    var: variable to plot. Ex: 'temperature', 'salinity'. Make sure
        to use the same name as defined in the netcdf file
    scatter_plot: if equal to 'yes' then a scatter plot 
            of the glider transect is plotted
    kwargs:        
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
    
    date_ini = kwargs.get('date_ini', None)
    date_end = kwargs.get('date_end', None)

    gdata = xr.open_dataset(url_glider,decode_times=False)
    
    inst_id = gdata.id.split('_')[0]

    variable =gdata.variables[var][0]
    latitude = gdata.latitude[0]
    longitude = gdata.longitude[0]
    depth = gdata.depth[0]
    
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
    varg =  variable[oktimeg,:]
    latg = latitude[oktimeg]
    long = longitude[oktimeg]
    depthg = depth[oktimeg,:]
    timeg = time[oktimeg]
    
    # Scatter plot
    if scatter_plot == 'yes':
        timeg_matrix = np.transpose(np.tile(timeg.T,(depthg.shape[1],1)))
        ttg = np.ravel(timeg_matrix)
        dg = np.ravel(depthg)
        teg = np.ravel(varg)

        kw = dict(c=teg, marker='*', edgecolor='none')

        fig, ax = plt.subplots(figsize=(10, 6))
        cs = ax.scatter(ttg,-dg,cmap='RdYlBu_r',**kw)
        #fig.colorbar(cs)
        ax.set_xlim(timeg[0], timeg[-1])

        ax.set_ylabel('Depth (m)',fontsize=16)
        cbar = plt.colorbar(cs)
        cbar.ax.set_ylabel(var,fontsize=16)
        ax.set_title(inst_id.split('-')[0],fontsize=20)
        xfmt = mdates.DateFormatter('%H:%Mh\n%d-%b')
        ax.xaxis.set_major_formatter(xfmt)
    
    return varg, latg, long, depthg, timeg, inst_id

#%%
    
def retrieve_glider_id_erddap_server(url_server,lat_lim,lon_lim,date_ini,date_end):
    
    """
    Created on Tue Feb  5 10:05:37 2019

    @author: aristizabal

    This function retrieves glider ids from the IOOS 
    Data Assembly Center (DAC).
    
    Inputs:
    url_server: url address of erddap server
                Example: 'https://data.ioos.us/gliders/erddap'
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
                    
    Outputs:
    gliders: list of gliders ids that fall within the lat, lon and
             time constraints
    """
    
    from erddapy import ERDDAP
    import pandas as pd

    e = ERDDAP(server = url_server)

    # Search constraints
    kw2018 = {
            'min_lon': lon_lim[0],
            'max_lon': lon_lim[1],
            'min_lat': lat_lim[0],
            'max_lat': lat_lim[1],
            'min_time': date_ini,
            'max_time': date_end,
            }

    search_url = e.get_search_url(response='csv', **kw2018)
    search = pd.read_csv(search_url)
    
    # Extract the IDs
    gliders = search['Dataset ID'].values
    
    return gliders

#%% 

def read_glider_data_erddap_server(url_server,dataset_id,var,\
                                   lat_lim,lon_lim,date_ini,date_end,\
                                   scatter_plot):
    
    """
    Created on Tue Feb  5 10:05:37 2019

    @author: aristizabal

    This function reads glider data from the IOOS 
    Data Assembly Center (DAC).
    
    Inputs:
    url_server: url address of thredds server
                Example: 'https://data.ioos.us/gliders/erddap'
    dataset_id: this id is retrieved from the glider DAC using the
               function "retrieve_glider_id_erddap_server".
               Example: 'ru30-20180705T1825'
    var: variable that needs to be read. Example, var = 'temperature'
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
    scatter_plot: if equal to 'yes' then a scatter plot 
            of the glider transect is plotted
                    
    Outputs:
    df: dataframe that contains the glider data
    """
    
    from erddapy import ERDDAP
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    # Readind data set
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
            var,
            ]

    e = ERDDAP(
            server=url_server,
            protocol='tabledap',
            response='nc'
            )

    e.dataset_id = dataset_id
    e.constraints = constraints
    e.variables = variables
    
    # Concerting glider data to data frame
    df = e.to_pandas(
            index_col='time',
            parse_dates=True,
            skiprows=(1,)  # units information can be dropped.
            ).dropna()

    # Scatter plot
    if scatter_plot == 'yes':

        fig, ax=plt.subplots(figsize=(10, 6), facecolor='w', edgecolor='w')

        kw = dict(s=30, c=df[var], marker='*', edgecolor='none')
        cs = ax.scatter(df.index, -df['depth'], **kw, cmap='RdYlBu_r')
        
        ax.set_xlim(df.index[0], df.index[-1])
        xfmt = mdates.DateFormatter('%H:%Mh\n%d-%b')
        ax.xaxis.set_major_formatter(xfmt)
        
        cbar = fig.colorbar(cs, orientation='vertical')
        cbar.ax.set_ylabel(var)
        ax.set_ylabel('Depth (m)');

    return df
     