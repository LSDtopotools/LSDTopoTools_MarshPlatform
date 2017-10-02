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
#1. Load functions

# Useful Python packages
import numpy as np
import cPickle
import timeit

# A very useful package
from LSDMarshPlatform_functions import ENVI_raster_binary_to_2d_array
from LSDMarshPlatform_functions import ENVI_raster_binary_from_2d_array

# The main functions for the marsh identification
from LSDMarshPlatform_functions import MARSH_ID
from LSDMarshPlatform_functions import Confusion



#------------------------------------------------------------------
#2. Set up your work environment


# Name your data input directory
Input_dir = "//csce.datastore.ed.ac.uk/csce/geos/users/s1563094/Software/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/"
# Name your results output directory
Output_dir = "//csce.datastore.ed.ac.uk/csce/geos/users/s1563094/Software/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/"
# NB: When naming your work directories, make sure the syntax is compatible with the OS you are using. Above is an example in the Windows command shell


# Select site names. Simply add a site in the list to analyse multiple sites simultaneously.
Sites = ["FEL"]


# Set the value for detection parameters. Values are currently set to the published optimum.
opt1 = -2.0
opt2 = 0.85
opt3 = 8.0


# Set the value for empty DEM cells
Nodata_value = -9999




#------------------------------------------------------------------
#3. Run the Marsh identification script

# Comment this line if you don't want to know how long the script run for.
Start = timeit.default_timer()


for site in Sites:
    print "Loading input data"
    
    # NB: When loading input data, please make sure the naming convention shown here is respected.
    print " Loading DEM"
    DEM, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array (Input_dir+"%s_DEM_clip.bil" % (site), site)
    print " Loading Slopes"
    Slope, post_Slope, envidata_Slope =  ENVI_raster_binary_to_2d_array (Input_dir+"%s_slope_clip.bil" % (site), site)


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
    
    # NB When loading input data, please make sure the naming convention shown here is respected.
    print " Loading detected Marsh"
    Platform_work, post_Platform, envidata_Platform =  ENVI_raster_binary_to_2d_array (Output_dir+"%s_Marsh.bil" % (site), site)
    print "Loading reference marsh"
    Reference, post_Reference, envidata_Reference =  ENVI_raster_binary_to_2d_array (Input_dir+"%s_ref_DEM_clip.bil" % (site), site)
    print "Evaluating the performance of the detection"
    Confusion_matrix, Performance, Metrix = Confusion (Platform_work, Reference, Nodata_value)
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_Platform, Output_dir+"%s_Confusion_DEM.bil" % (site), post_Platform, Confusion_matrix)
    cPickle.dump(Performance,open(Output_dir+"%s_Performance.pkl" % (site), "wb"))
    cPickle.dump(Metrix,open(Output_dir+"%s_Metrix.pkl" % (site), "wb"))


    
    
# Comment these 2 lines if you don't want to know how long the script run for.    
Stop = timeit.default_timer()
print 'Runtime = ', Stop - Start , 's'



