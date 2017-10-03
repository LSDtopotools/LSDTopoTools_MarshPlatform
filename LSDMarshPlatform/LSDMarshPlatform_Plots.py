"""
Example_Plots.py


    This contains wrapper functions that simplify plotting raster
    and vector data for publication-ready figures.
    The documentation of the examples can be found here:
    https://lsdtopotools.github.io/LSDTopoTools_ChiMudd2014/
    Simon Mudd and Fiona Clubb, June 2017
    Released under GPL3

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
from LSDMarshPlatform_functions import ENVI_raster_binary_to_2d_array
from LSDMarshPlatform_functions import ENVI_raster_binary_from_2d_array
from LSDMarshPlatform_functions import kernel
from LSDMarshPlatform_functions import Outline


#------------------------------------------------------------------
#2. Set up the important variables

# Name your data input directory
#Input_dir = "//csce.datastore.ed.ac.uk/csce/geos/users/s1563094/Software/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/"
# Name your results output directory
#Output_dir = "//csce.datastore.ed.ac.uk/csce/geos/users/s1563094/Software/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/"
# NB: When naming your work directories, make sure the syntax is compatible with the OS you are using. Above is an example in the Windows command shell



#Select site names. Simply add a site in the list to analyse multiple sites simultaneously.
Sites = ["FEL"]

# Set the value for empty DEM cells
Nodata_value = -9999

def Plot_platform_on_hillshade(Input_dir =  "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/", 
            Output_dir = "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Sites = ["FEL"]):
    """
    This plots the extracted marsh platform on a hillshade
    
    Args:
        Input_dir (str): Name your data input directory
        Output_dir (str): Name your results output directory
        Sites (str list): A list of strings. The file names are modified based on these sites
        
    Author: GCHG
    """
    #Plot 1: Draw the platform on a DEM, superimposed on a hillshade
    for site in Sites:
        fig=plt.figure(1, facecolor='White',figsize=[10,10])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1, facecolor='white')

        # Name the axes
        ax1.set_xlabel('x (m)', fontsize = 12)
        ax1.set_ylabel('y (m)', fontsize = 12)

        # Load the relevant data
        HS, post_HS, envidata_HS = ENVI_raster_binary_to_2d_array (Input_dir+"%s_hs.bil" % (site), site)
        DEM, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array (Input_dir+"%s.bil" % (site), site)
        Platform, post_Platform, envidata_Platform = ENVI_raster_binary_to_2d_array (Output_dir+"%s_Marsh.bil" % (site), site)


        # Make a mask!
        Platform_mask = np.ma.masked_where(Platform <=0, Platform)
        Platform_mask[Platform_mask>0] = DEM[Platform_mask>0]

        # Make a map!
        Map_HS = ax1.imshow(HS, interpolation='None', cmap=plt.cm.gist_gray, vmin = 100, vmax = 200)
        Map_DEM = ax1.imshow(DEM, interpolation='None', cmap=plt.cm.gist_gray, vmin = np.amin(DEM[DEM!=Nodata_value]), vmax = np.amax(DEM), alpha = 0.5)
        Map_Marsh = ax1.imshow(Platform_mask, interpolation='None', cmap=plt.cm.gist_earth, vmin=np.amin(DEM[DEM!=Nodata_value]), vmax=np.amax(DEM), alpha = 0.5)

        plt.savefig(Output_dir+'Figure_1.png')



def Plot_marsh_outline_on_hillshade(Input_dir =  "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/", 
            Output_dir = "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Sites = ["FEL"]):
    """
    This draws the marsh outline on a DEM, superimposed on a hillshade.
    
    Args:
        Input_dir (str): Name your data input directory.
        Output_dir (str): Name your results output directory.
        Sites (str list): A list of strings. The file names are modified based on these sites.
        
    Author: GCHG
    """
    
    for site in Sites:
        fig=plt.figure(2, facecolor='White',figsize=[10,10])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1, facecolor='white')

        # Name the axes
        ax1.set_xlabel('x (m)', fontsize = 12)
        ax1.set_ylabel('y (m)', fontsize = 12)

        # Load the relevant data
        HS, post_HS, envidata_HS = ENVI_raster_binary_to_2d_array (Input_dir+"%s_hs.bil" % (site), site)
        DEM, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array (Input_dir+"%s.bil" % (site), site)
        Platform, post_Platform, envidata_Platform = ENVI_raster_binary_to_2d_array (Output_dir+"%s_Marsh.bil" % (site), site)

        # Outline the marsh
        Platform[Platform > 0] = 1
        Marsh_outline = Outline (Platform, 2, Nodata_value)


        # Make a mask!
        Outline_mask = np.ma.masked_where(Marsh_outline <=1, Marsh_outline)

        # Make a map!
        Map_HS = ax1.imshow(HS, interpolation='None', cmap=plt.cm.gist_gray, vmin = 100, vmax = 200)
        Map_DEM = ax1.imshow(DEM, interpolation='None', cmap=plt.cm.gist_gray, vmin = np.amin(DEM[DEM!=Nodata_value]), vmax = np.amax(DEM), alpha = 0.5)

        Map_Marsh = ax1.imshow(Outline_mask, interpolation='None', cmap=plt.cm.Reds, vmin = 0, vmax = 2, alpha = 1)


    plt.savefig(Output_dir+'Figure_2.png')


def Plot_marsh_reference_on_hillshade(Input_dir =  "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/", 
            Output_dir = "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Sites = ["FEL"]):
    """
    This draws the marsh reference on a hillshade
    
    Args:
        Input_dir (str): Name your data input directory.
        Output_dir (str): Name your results output directory.
        Sites (str list): A list of strings. The file names are modified based on these sites.
        
    Author: GCHG
    """
    #Plot 3: Draw the marsh map and reference outline, superimposed on a hillshade
    for site in Sites:
        fig=plt.figure(3, facecolor='White',figsize=[10,10])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1, facecolor='white')


        # Name the axes
        ax1.set_xlabel('x (m)', fontsize = 12)
        ax1.set_ylabel('y (m)', fontsize = 12)

        # Load the relevant data
        HS, post_HS, envidata_HS = ENVI_raster_binary_to_2d_array (Input_dir+"%s_hs.bil" % (site), site)
        DEM, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array (Input_dir+"%s.bil" % (site), site)
        Platform, post_Platform, envidata_Platform = ENVI_raster_binary_to_2d_array (Output_dir+"%s_Marsh.bil" % (site), site)
        Reference, post_Reference, envidata_Reference = ENVI_raster_binary_to_2d_array (Input_dir+"%s_ref.bil" % (site), site)

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

    plt.savefig(Output_dir+'Figure_3.png')



def Plot_confusion_map_on_hillshade(Input_dir =  "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/", 
            Output_dir = "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Sites = ["FEL"]):
    """
    This draws the marsh reference on a hillshade
    
    Args:
        Input_dir (str): Name your data input directory.
        Output_dir (str): Name your results output directory.
        Sites (str list): A list of strings. The file names are modified based on these sites.
        
    Author: GCHG
    """
    #Plot 4: Draw the confusion map, superimposed on a hillshade
    for site in Sites:
        fig=plt.figure(4, facecolor='White',figsize=[10,10])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1, facecolor='white')


        # Name the axes
        ax1.set_xlabel('x (m)', fontsize = 12)
        ax1.set_ylabel('y (m)', fontsize = 12)

        # Load the relevant data
        HS, post_HS, envidata_HS = ENVI_raster_binary_to_2d_array (Input_dir+"%s_hs.bil" % (site), site)
        Confusion, post_Confusion, envidata_Confusion = ENVI_raster_binary_to_2d_array (Output_dir+"%s_Confusion.bil" % (site), site)


        # Make a mask!
        Confusion_mask = np.ma.masked_where(Confusion == Nodata_value, Confusion)


        # Make a map!
        Map_HS = ax1.imshow(HS, interpolation='None', cmap=plt.cm.gist_gray, vmin = 100, vmax = 200)
        Map_DEM = ax1.imshow(Confusion_mask, interpolation='None', cmap=plt.cm.Spectral, vmin = -2, vmax = 2, alpha = 0.5)


    plt.savefig(Output_dir+'Figure_4.png')



