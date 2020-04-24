#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%

def grid_glider_data(var_name,dataset_id,varg,timeg,latg,long,depthg,delta_z=0.3,contour_plot='yes'):

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
    var_name: name of variable to be gridded
    dataset_id: obtained from running function "retrieve_dataset_id_erddap_server"
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
    import cmocean

    # sort time variable
    okt = np.argsort(timeg)
    timegg = timeg[okt]
    #latgg = latg[okt]
    #longg = long[okt]
    depthgg = depthg[:,okt]
    vargg = varg[:,okt]

    # Grid variables
    depthg_gridded = np.arange(0,np.nanmax(depthgg),delta_z)
    varg_gridded = np.empty((len(depthg_gridded),len(timegg)))
    varg_gridded[:] = np.nan

    for t,tt in enumerate(timegg):
        depthu,oku = np.unique(depthgg[:,t],return_index=True)
        varu = vargg[oku,t]
        okdd = np.isfinite(depthu)
        depthf = depthu[okdd]
        varf = varu[okdd]
        ok = np.asarray(np.isfinite(varf))
        if np.sum(ok) < 3:
            varg_gridded[:,t] = np.nan
        else:
            #okd = depthg_gridded < np.max(depthf[ok])
            okd = np.logical_and(depthg_gridded >= np.min(depthf[ok]),\
                                 depthg_gridded < np.max(depthf[ok]))
            varg_gridded[okd,t] = np.interp(depthg_gridded[okd],depthf[ok],varf[ok])

    # Countour plot
    if contour_plot == 'yes':
        
        if var_name == 'temperature':
            color_map = cmocean.cm.thermal
            clabel = var_name[0].upper()+var_name[1:] + ' ($^oC$)'
        else:
            if var_name == 'salinity':
                color_map = cmocean.cm.haline
                clabel = var_name[0].upper()+var_name[1:]
            else:
                color_map = 'RdBu_r'
                
        okg = depthg_gridded <= np.max(depthg_gridded) 
        min_val = np.int(np.floor(np.min(np.nanmin(varg_gridded[okg]))))
        max_val = np.int(np.ceil(np.max(np.nanmax(varg_gridded[okg]))))
    
        if var_name == 'salinity':
            kw = dict(levels = np.arange(min_val,max_val+0.25,0.25))
        else:
            nlevels = max_val - min_val + 1
            kw = dict(levels = np.linspace(min_val,max_val,nlevels)) 
            
        #nlevels = np.round(np.nanmax(varg_gridded)) - np.round(np.nanmin(varg_gridded)) + 1
        #kw = dict(levels = np.linspace(np.round(np.nanmin(varg_gridded)),\
        #                               np.round(np.nanmax(varg_gridded)),nlevels))

        fig, ax=plt.subplots(figsize=(10, 3))
        
        plt.contour(timegg,-depthg_gridded,varg_gridded,levels=[26],colors = 'k')
        cs = plt.contourf(timegg,-depthg_gridded,varg_gridded,cmap=color_map,**kw)
        plt.title(dataset_id,fontsize=16)

        ax.set_xlim(timeg[0], timeg[-1])
        xfmt = mdates.DateFormatter('%H:%Mh\n%d-%b')
        ax.xaxis.set_major_formatter(xfmt)

        cbar = fig.colorbar(cs, orientation='vertical')
        cbar.ax.set_ylabel(clabel,fontsize=14)
        ax.set_ylabel('Depth (m)',fontsize=14)

    return varg_gridded, timegg, depthg_gridded,