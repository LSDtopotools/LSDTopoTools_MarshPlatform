"""
LSDMarshPlatform_Marsh_ID.py

This is your driver file to run the marsh platform extraction.

Please read the README and the instructions in this script before you run it.

Authors: Guillaume GH Goodwin and Simon Marius Mudd

"""

#------------------------------------------------------------------
#0. Set up display environment if you are working on a terminal with no GUI.
import matplotlib
matplotlib.use('Agg')

#------------------------------------------------------------------

# Useful Python packages
import numpy as np
import cPickle
import timeit
import os

# A very useful package
from LSDMarshPlatform_functions import ENVI_raster_binary_to_2d_array
from LSDMarshPlatform_functions import ENVI_raster_binary_from_2d_array

# The main functions for the marsh identification
from LSDMarshPlatform_functions import MARSH_ID
from LSDMarshPlatform_functions import Confusion

# Retained directories from Guillaume
# "//csce.datastore.ed.ac.uk/csce/geos/users/s1563094/Software/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/"


def MarshID(Input_dir =  "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Output_dir = "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Sites = ["FEL_DEM_clip"], opt1 = -2.0, opt2 = 0.85, opt3 = 8.0, compare_with_digitised_marsh = False):
    """
    This function wraps all the marsh ID scripts in one location

    Args:
        Input_dir (str): Name your data input directory
        Output_dir (str): Name your results output directory
        Sites (str list): A list of strings. The file names are modified based on these sites
        opt1 (flt): first optimisation
        opt2 (flt): 2nd optimisation
        opt3 (flt): 3rd optimisation
        compare_with_digitised_marsh (bool): If true, this will compare the data with a digitised marsh platform

    Author:
        GCHG, Modified by SMM 02/10/2017
    """
    #------------------------------------------------------------------
    print("Welcome to the marsh ID program!")
    print("I am opening the file: "+Input_dir)

    # Set the value for empty DEM cells
    Nodata_value = -9999

    # Timing
    Start = timeit.default_timer()


    for site in Sites:
        print("Loading input data from site: "+site)

        # NB: When loading input data, please make sure the naming convention shown here is respected.
        print(" Loading DEM")
        DEM, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array (Input_dir+"%s.bil" % (site), site)
        print " Loading Slopes"
        # check to get the correct slope raster
        slope_fname = site+"_slope.bil"
        if not os.path.isfile(Input_dir+slope_fname):
            slope_fname = site+"_SLOPE.bil"
        Slope, post_Slope, envidata_Slope =  ENVI_raster_binary_to_2d_array (Input_dir+slope_fname, site)


        # Here begins the detection process
        print "Identifying the platform and scarps"
        DEM_work = np.copy(DEM)
        Search_space, Scarps, Platform = MARSH_ID(DEM, Slope, Nodata_value, opt1, opt2, opt3)
        Platform_work = np.copy(Platform)
        Scarps[Scarps == 0] = Nodata_value

        # Here is where you save your output files for use in a GIS software
        print "Saving marsh features"
        new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, Output_dir+"%s_Search_space.bil" % (site), post_DEM, Search_space)
        new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, Output_dir+"%s_Scarps.bil" % (site), post_DEM, Scarps)
        new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, Output_dir+"%s_Marsh.bil" % (site), post_DEM, Platform)




        # Disable the following section if you do not wish to compare your results to a reference marsh
        if compare_with_digitised_marsh:
            # NB When loading input data, please make sure the naming convention shown here is respected.
            print " Loading detected Marsh"
            Platform_work, post_Platform, envidata_Platform =  ENVI_raster_binary_to_2d_array (Output_dir+"%s_Marsh.bil" % (site), site)
            print "Loading reference marsh"
            Reference, post_Reference, envidata_Reference =  ENVI_raster_binary_to_2d_array (Input_dir+"%s_ref.bil" % (site), site)
            print "Evaluating the performance of the detection"
            Confusion_matrix, Performance, Metrix = Confusion (Platform_work, Reference, Nodata_value)
            new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_Platform, Output_dir+"%s_Confusion.bil" % (site),                                                                                post_Platform, Confusion_matrix)
            cPickle.dump(Performance,open(Output_dir+"%s_Performance.pkl" % (site), "wb"))
            cPickle.dump(Metrix,open(Output_dir+"%s_Metrix.pkl" % (site), "wb"))




    # Comment these 2 lines if you don't want to know how long the script run for.
    Stop = timeit.default_timer()
    print 'Runtime = ', Stop - Start , 's'
