glider_model_comparisons

Author: Maria Aristizabal Vargars

This repository is a set of Python functions that retrieve and plot glider data 
from the IOOS glider DAC. There are two ways to access the data: 
from a thredds server at 'https://data.ioos.us/thredds/dodsC/deployments/', 
and from an erddap server 'https://data.ioos.us/gliders/erddap'. 

The function "glider_transect_model_com_erddap_server" also retrieves model 
output from a global operational ocean models. So far it is set up to read from 
the Global Ocean Forecasting System 
(GOFS 3.1 at 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_93.0/ts3z' 
(GOFS 3.1).

Take a look at the script 'working_examples.py' to see examples of how to use 
the different functions.
