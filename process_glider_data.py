#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def grid_glider_data_erddap(df,var,delta_z=0.3,contour_plot='yes'):

    """
    Created on Wed Feb  6 11:49:24 2019

    @author: aristizabal

    This function grids glider data so it can be plotted as a countour plot.

    Inputs:
    df: data frame that contains glider data. It should contain at least
        time, latitude, longitude and a variable such as temperature or
        salinity. df is obtained running "read_glider_data_erddap_server"
    var: variable to be gridded. example: 'temperature', 'salinity'
    delta_z: desired spacing in meters of the vertical levels of output
             variable var_gridded. example, delta_z=0.5. Default value is 0.3
    contour_plot: if equal to 'yes' then a contour plot
            of the glider transect is plotted. Default value is 'yes'


    Outputs:
    depthg_gridded: gridded depth vector
    varg_gridded: gridded variable matrix
    timeg: time vector
    """

    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    # Coverting glider vectors into arrays
    timeg, ind = np.unique(df.index.values,return_index=True)
    latg = df['latitude (degrees_north)'].values[ind]
    long = df['longitude (degrees_east)'].values[ind]

    dg = df['depth (m)'].values
    vg = df[df.columns[3]].values

    zn = np.int(np.round(np.max(dg)/delta_z))

    depthg = np.empty((zn,len(timeg)))
    depthg[:] = np.nan
    varg = np.empty((zn,len(timeg)))
    varg[:] = np.nan

    # Grid variables
    depthg_gridded = np.arange(0,np.nanmax(dg),delta_z)
    varg_gridded = np.empty((len(depthg_gridded),len(timeg)))
    varg_gridded[:] = np.nan

    for i,ii in enumerate(ind):
        if i < len(timeg)-1:
            depthg[0:len(dg[ind[i]:ind[i+1]]),i] = dg[ind[i]:ind[i+1]]
            varg[0:len(vg[ind[i]:ind[i+1]]),i] = vg[ind[i]:ind[i+1]]
        else:
            depthg[0:len(dg[ind[i]:len(dg)]),i] = dg[ind[i]:len(dg)]
            varg[0:len(vg[ind[i]:len(vg)]),i] = vg[ind[i]:len(vg)]

    for t,tt in enumerate(timeg):
        depthu,oku = np.unique(depthg[:,t],return_index=True)
        varu = varg[oku,t]
        okdd = np.isfinite(depthu)
        depthf = depthu[okdd]
        varf = varu[okdd]
        ok = np.isfinite(varf)
        if np.sum(ok) < 3:
            varg_gridded[:,t] = np.nan
        else:
            okd = depthg_gridded < np.max(depthf[ok])
            varg_gridded[okd,t] = np.interp(depthg_gridded[okd],depthf[ok],varf[ok])


    # Countour plot
    if contour_plot == 'yes':

        fig, ax=plt.subplots(figsize=(10, 6), facecolor='w', edgecolor='w')

        nlevels = np.round(np.nanmax(varg_gridded)) - np.round(np.nanmin(varg_gridded)) + 1
        kw = dict(levels = np.linspace(np.round(np.nanmin(varg_gridded)),\
                                       np.round(np.nanmax(varg_gridded)),nlevels))
        #plt.contour(timeg,-depthg_gridded,varg_gridded,colors = 'lightgrey',**kw)
        cs = plt.contourf(timeg,-depthg_gridded,varg_gridded,cmap='RdYlBu_r',**kw)

        ax.set_xlim(df.index[0], df.index[-1])
        xfmt = mdates.DateFormatter('%H:%Mh\n%d-%b')
        ax.xaxis.set_major_formatter(xfmt)

        cbar = fig.colorbar(cs, orientation='vertical')
        cbar.ax.set_ylabel(var)
        #cbar.ax.set_ylabel('Temperature ($^\circ$C)')
        ax.set_ylabel('Depth (m)');

    return depthg_gridded, varg_gridded, timeg, latg, long

#%%

def grid_glider_data_thredd(timeg,latg,long,depthg,varg,var,inst_id,delta_z=0.3,contour_plot='yes'):

    """
    Created on Wed Feb 25 2019

    @author: aristizabal

    This function grids glider data so it can be plotted as a countour plot.

    Inputs:
    timeg: time vector
    latg: latitude within the user defined time window. length is the same as
          timeg
    long: longitude within the user defined time window. length is the same as
          timeg
    depthg: depth array for all profiles. It must be a 2D array with
          Dimensions timeg x depth
    varg: variable to be gridded. example: 'temperature', 'salinity'
          it must be a 2D array with dimension timeg x depth
    var: name of variable to be gridded
    delta_z: desired spacing in meters of the vertical levels of output
             variable var_gridded. example, delta_z=0.5. Default value is 0.3
    contour_plot: if equal to 'yes' then a contour plot
            of the glider transect is plotted. Default value is 'yes'


    Outputs:
    depthg_gridded: gridded depth vector
    varg_gridded: gridded variable matrix
    timegg: sorted time vector

    """

    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    # sort time variable
    okt = np.argsort(timeg)
    timegg = timeg[okt]
    #latgg = latg[okt]
    #longg = long[okt]
    depthgg = depthg[okt,:]
    vargg = varg[okt,:]

    # Grid variables
    depthg_gridded = np.arange(0,np.nanmax(depthgg),delta_z)
    varg_gridded = np.empty((len(timegg),len(depthg_gridded)))
    varg_gridded[:] = np.nan

    for t,tt in enumerate(timegg):
        depthu,oku = np.unique(depthgg[t,:],return_index=True)
        varu = vargg[t,oku]
        okdd = np.isfinite(depthu)
        depthf = depthu[okdd]
        varf = varu[okdd]
        ok = np.asarray(np.isfinite(varf))
        if np.sum(ok) < 3:
            varg_gridded[t,:] = np.nan
        else:
            okd = depthg_gridded < np.max(depthf[ok])
            varg_gridded[t,okd] = np.interp(depthg_gridded[okd],depthf[ok],varf[ok])

    # Countour plot
    if contour_plot == 'yes':

        fig, ax=plt.subplots(figsize=(10, 6), facecolor='w', edgecolor='w')

        nlevels = np.round(np.nanmax(varg_gridded)) - np.round(np.nanmin(varg_gridded)) + 1
        kw = dict(levels = np.linspace(np.round(np.nanmin(varg_gridded)),\
                                       np.round(np.nanmax(varg_gridded)),nlevels))
        #plt.contour(timegg,-depthg_gridded,varg_gridded.T,levels=26,colors = 'k')
        cs = plt.contourf(timegg,-depthg_gridded,varg_gridded.T,cmap='RdYlBu_r',**kw)
        plt.title(inst_id.split('-')[0],fontsize=20)

        ax.set_xlim(timeg[0], timeg[-1])
        xfmt = mdates.DateFormatter('%H:%Mh\n%d-%b')
        ax.xaxis.set_major_formatter(xfmt)

        cbar = fig.colorbar(cs, orientation='vertical')
        cbar.ax.set_ylabel(var,fontsize=16)
        ax.set_ylabel('Depth (m)',fontsize=16);

    return depthg_gridded, varg_gridded, timegg
