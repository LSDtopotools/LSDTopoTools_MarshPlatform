"""
Example_Plots.py

This script produces nice plots to visualise your results

"""


#------------------------------------------------------------------
#0. Set up display environment in putty if you are working on a terminal with no graphical interface.
import matplotlib
matplotlib.use('Agg')

#----------------------------------------------------------------
#1. Load useful Python packages
import os
import sys

import numpy as np
import functools
import math as mt
import cmath
import scipy as sp
import scipy.stats as stats
from datetime import datetime
import cPickle
from pylab import *
import functools
import itertools as itt
from osgeo import gdal, osr
from osgeo import gdal, gdalconst
from osgeo.gdalconst import *
from copy import copy
from matplotlib import cm
import matplotlib.colors as colors
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import *
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap, cm
from matplotlib.patches import Rectangle
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import timeit
#------------------------------------------------------------------
# Import the marsh-finding functions
from Example_functions import ENVI_raster_binary_to_2d_array
from Example_functions import ENVI_raster_binary_from_2d_array
from Example_functions import kernel



# This is a function that makes an outline out of a raster where each object has a single given value

def Outline (Raster, Outline_value, Nodata_value):

    P1 = np.where(Raster[:,1:] != Raster[:,:-1])
    Raster[P1] = Outline_value           

    P2 = np.where(Raster[1:,:] != Raster[:-1,:])
    Raster[P2] = Outline_value
    
    for i in range(len(Raster)):
        for j in range(len(Raster[0,:])):
            if Raster[i,j] == Outline_value:
                K = kernel (Raster, 3, i, j)
                if np.mean(K) < 0:
                    Raster[i,j] = 0
    
    return Raster






#------------------------------------------------------------------
#2. Set up the important variables

#Select site names. Simply add a site in the list to analyse multiple sites simultaneously.
Sites = ["FEL"]

# Set the value for empty DEM cells
Nodata_value = -9999


# Use this if you want to time the execution
Start = timeit.default_timer()


#------------------------------------------------------------------
#3. Start Plotting

#Plot 1: Draw the platform on a DEM, superimposed on a hillshade
for site in Sites:
    fig=plt.figure(1, facecolor='White',figsize=[10,10])
    ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1, axisbg='white')

    # Name the axes
    ax1.set_xlabel('x (m)', fontsize = 12)
    ax1.set_ylabel('y (m)', fontsize = 12)

    # Load the relevant data
    HS, post_HS, envidata_HS = ENVI_raster_binary_to_2d_array ("Input/%s_hs_clip.bil" % (site), site)
    DEM, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array ("Input/%s_DEM_clip.bil" % (site), site)
    Platform, post_Platform, envidata_Platform = ENVI_raster_binary_to_2d_array ("Output/%s_Marsh.bil" % (site), site)

    
    # Make a mask!
    Platform_mask = np.ma.masked_where(Platform <=0, Platform)
    Platform_mask[Platform_mask>0] = DEM[Platform_mask>0]
    
    # Make a map!
    Map_HS = ax1.imshow(HS, interpolation='None', cmap=plt.cm.gist_gray, vmin = 100, vmax = 200)
    Map_DEM = ax1.imshow(DEM, interpolation='None', cmap=plt.cm.gist_gray, vmin = np.amin(DEM[DEM!=Nodata_value]), vmax = np.amax(DEM), alpha = 0.5)
    Map_Marsh = ax1.imshow(Platform_mask, interpolation='None', cmap=plt.cm.gist_earth, vmin=np.amin(DEM[DEM!=Nodata_value]), vmax=np.amax(DEM), alpha = 0.5)
    


plt.savefig('Output/Figure_1.png')




#Plot 2: Draw the marsh outline on a DEM, superimposed on a hillshade
for site in Sites:
    fig=plt.figure(2, facecolor='White',figsize=[10,10])
    ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1, axisbg='white')

    # Name the axes
    ax1.set_xlabel('x (m)', fontsize = 12)
    ax1.set_ylabel('y (m)', fontsize = 12)

    # Load the relevant data
    HS, post_HS, envidata_HS = ENVI_raster_binary_to_2d_array ("Input/%s_hs_clip.bil" % (site), site)
    DEM, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array ("Input/%s_DEM_clip.bil" % (site), site)
    Platform, post_Platform, envidata_Platform = ENVI_raster_binary_to_2d_array ("Output/%s_Marsh.bil" % (site), site)

    # Outline the marsh
    Platform[Platform > 0] = 1
    Marsh_outline = Outline (Platform, 2, Nodata_value)

    
    # Make a mask!
    Outline_mask = np.ma.masked_where(Marsh_outline <=1, Marsh_outline)

    # Make a map!
    Map_HS = ax1.imshow(HS, interpolation='None', cmap=plt.cm.gist_gray, vmin = 100, vmax = 200)
    Map_DEM = ax1.imshow(DEM, interpolation='None', cmap=plt.cm.gist_gray, vmin = np.amin(DEM[DEM!=Nodata_value]), vmax = np.amax(DEM), alpha = 0.5)
    
    Map_Marsh = ax1.imshow(Outline_mask, interpolation='None', cmap=plt.cm.Reds, vmin = 0, vmax = 2, alpha = 1)
    

plt.savefig('Output/Figure_2.png')



#Plot 3: Draw the marsh map and reference outline, superimposed on a hillshade
for site in Sites:
    fig=plt.figure(3, facecolor='White',figsize=[10,10])
    ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1, axisbg='white')


    # Name the axes
    ax1.set_xlabel('x (m)', fontsize = 12)
    ax1.set_ylabel('y (m)', fontsize = 12)

    # Load the relevant data
    HS, post_HS, envidata_HS = ENVI_raster_binary_to_2d_array ("Input/%s_hs_clip.bil" % (site), site)
    DEM, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array ("Input/%s_DEM_clip.bil" % (site), site)
    Platform, post_Platform, envidata_Platform = ENVI_raster_binary_to_2d_array ("Output/%s_Marsh.bil" % (site), site)
    Reference, post_Reference, envidata_Reference = ENVI_raster_binary_to_2d_array ("Input/%s_ref_DEM_clip.bil" % (site), site)

    # Outline the reference
    Reference[Reference > 0] = 1
    Ref_outline = Outline (Reference, 2, Nodata_value)

    
    # Make a mask!
    Outline_mask = np.ma.masked_where(Ref_outline <=1, Ref_outline)
    
    
    # Make a map!
    Platform_mask = np.ma.masked_where(Platform <=0, Platform)
    Platform_mask[Platform_mask>0] = DEM[Platform_mask>0]

    Map_HS = ax1.imshow(HS, interpolation='None', cmap=plt.cm.gist_gray, vmin = 100, vmax = 200)
    Map_DEM = ax1.imshow(DEM, interpolation='None', cmap=plt.cm.gist_gray, vmin = np.amin(DEM[DEM!=Nodata_value]), vmax = np.amax(DEM), alpha = 0.5)
    Map_Marsh = ax1.imshow(Platform_mask, interpolation='None', cmap=plt.cm.gist_earth, vmin=np.amin(DEM[DEM!=Nodata_value]), vmax=np.amax(DEM), alpha = 0.5)
    
    Map_Marsh = ax1.imshow(Outline_mask, interpolation='None', cmap=plt.cm.Reds, vmin = 0, vmax = 2, alpha = 1)

plt.savefig('Output/Figure_3.png')




#Plot 4: Draw the confusion map, superimposed on a hillshade
for site in Sites:
    fig=plt.figure(4, facecolor='White',figsize=[10,10])
    ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1, axisbg='white')


    # Name the axes
    ax1.set_xlabel('x (m)', fontsize = 12)
    ax1.set_ylabel('y (m)', fontsize = 12)

    # Load the relevant data
    HS, post_HS, envidata_HS = ENVI_raster_binary_to_2d_array ("Input/%s_hs_clip.bil" % (site), site)
    Confusion, post_Confusion, envidata_Confusion = ENVI_raster_binary_to_2d_array ("Output/%s_Confusion_DEM.bil" % (site), site)

    
    # Make a mask!
    Confusion_mask = np.ma.masked_where(Confusion == Nodata_value, Confusion)
    
    
    # Make a map!
    Map_HS = ax1.imshow(HS, interpolation='None', cmap=plt.cm.gist_gray, vmin = 100, vmax = 200)
    Map_DEM = ax1.imshow(Confusion_mask, interpolation='None', cmap=plt.cm.Spectral, vmin = -2, vmax = 2, alpha = 0.5)
    

plt.savefig('Output/Figure_4.png')



Stop = timeit.default_timer()
print 'Runtime = ', Stop - Start , 's'

