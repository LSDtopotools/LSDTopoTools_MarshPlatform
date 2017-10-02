"""
Example_Marsh_ID.py

This is the script where you enter your input. It is also where the calculations happen. Please read the README and the instructions in this script before you run it.

"""

#------------------------------------------------------------------
#0. Set up display environment in putty if you are working on a terminal with no graphical interface.
import matplotlib
matplotlib.use('Agg')

#------------------------------------------------------------------
#1. Load useful Python packages
import numpy as np
import cPickle
import timeit


#------------------------------------------------------------------
# Import the costum-made marsh-finding functions
from Example_functions import MARSH_ID
from Example_functions import Confusion
from Example_functions import ENVI_raster_binary_to_2d_array
from Example_functions import ENVI_raster_binary_from_2d_array


#------------------------------------------------------------------
#2. Set up the important variables

#Select site names. Simply add a site in the list to analyse multiple sites simultaneously.
Sites = ["FEL"]

# Set the value for empty DEM cells
Nodata_value = -9999

# Set the value for optimised detection parameters. Values are currently set to the published optimum.
opt1 = -2.0
opt2 = 0.85
opt3 = 8.0


#------------------------------------------------------------------
#3. Run the Marsh identification script

# Use this if you want to time the execution
Start = timeit.default_timer()

for site in Sites:
    print "Loading input data"
    print " Loading DEM"
    DEM, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array ("Input/%s_DEM_clip.bil" % (site), site)
    print " Loading Slopes"
    Slope, post_Slope, envidata_Slope =  ENVI_raster_binary_to_2d_array ("Input/%s_slope_clip.bil" % (site), site)
    print " Loading Curvature"
    Curvature, post_Curvature, envidata_Curvature =  ENVI_raster_binary_to_2d_array ("Input/%s_curvature_clip.bil" % (site), site)

    
    print "Identifying the platform and scarps"
    DEM_work = np.copy(DEM)
    
    # This is where the magic happens
    Search_space, Scarps, Platform = MARSH_ID(DEM, Slope, Curvature, Nodata_value, opt1, opt2, opt3)
    
    Platform_work = np.copy(Platform)
    Scarps[Scarps == 0] = Nodata_value


    print "Saving marsh features"
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "Output/%s_Search_space.bil" % (site), post_DEM, Search_space)
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "Output/%s_Scarps.bil" % (site), post_DEM, Scarps)
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "Output/%s_Marsh.bil" % (site), post_DEM, Platform)



    
    # Disable this section if you have no reference marsh.
    print " Loading Detected Marsh"
    Platform_work, post_Platform, envidata_Platform =  ENVI_raster_binary_to_2d_array ("Output/%s_Marsh.bil" % (site), site)

    print "Measuring performances"
    Reference, post_Reference, envidata_Reference =  ENVI_raster_binary_to_2d_array ("Input/%s_ref_DEM_clip.bil" % (site), site)
    
    Confusion_matrix, Performance, Metrix = Confusion (Platform_work, Reference, Nodata_value)
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_Platform, "Output/%s_Confusion_DEM.bil" % (site), post_Platform, Confusion_matrix)

    cPickle.dump(Performance,open("Output/%s_Performance.pkl" % (site), "wb"))
    cPickle.dump(Metrix,open("Output/%s_Metrix.pkl" % (site), "wb"))

    
    
Stop = timeit.default_timer()

print 'Runtime = ', Stop - Start , 's'



