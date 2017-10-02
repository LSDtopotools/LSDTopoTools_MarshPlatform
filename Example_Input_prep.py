"""
Example_Input_prep.py 

This is a Python script to prepare DEM files for the automated topographic detection of saltmarshes. This file is particularly important and performs the following actions:
1. Convert your DEM into the ENVI format used by LSDTopoTools
2. Apply a Wiener filter (optional)
3. Calculate basic topographic analysis (hillshade, slope and curvature)
4. Prepare a reference raster.
5. Clip the inputs to a specified domain
"""

#------------------------------------------------------------------
#0. Set up display environment in putty if you are working on a terminal with no graphical interface.
import matplotlib
matplotlib.use('Agg')


#------------------------------------------------------------------
#1. Load useful Python packages
import os
import sys


#------------------------------------------------------------------
#2. Set up the important variables

#Select site names. Simply add a site in the list to analyse multiple sites simultaneously.
Sites = ["FEL"]

# Set up value for empty DEM cells
Nodata_value = -9999

#------------------------------------------------------------------
#3. Run the preparation script.
# This part of the script runs in or three stages.
# For each stage, you will need to copy the file "Example_Input_pre.py" to the indicated location and comment out the other stages, unless the analysis pqckages are in the same directory.
# Each stage is therefore run individually.


for site in Sites:
    # Stage 1: If the original DEM of your site is not in ENVI .bil format or does not have a CRS, please run this stage.
    #          Run this stage from the directory where you stored your DEM.
    #          For me, this is : /home/s1563094/Datastore/Software/LSDTopoTools/LSDTopoTools_MarshExtraction/Example/Input
    
    """print "Converting DEM files to ENVI format"
    sourcefile = "%s_trans.asc" % (site) # This is the filename for your original DEM
    destinationfile = "%s_DEM.bil" % (site) # This is the filename for your processed DEM
    os.system("gdal_translate -of ENVI -a_srs EPSG:27700 " + sourcefile + " " +  destinationfile)"""
    
    # NOTE: The FEL_trans.asc file is too heavy to upload on GitHub. This stage has therefore already been taken care of in the example.
    # NOTE: The CRS is here set to be the British National GRid to correspond to the example data.
    # NOTE: Useful webpage for further use of GDAL: https://pcjericks.github.io/py-gdalogr-cookbook/
    
    
    
    # STAGE 2: If you wish to apply a Wiener filter, please run this stage.
    #          Run this stage from the directory where you stored the Wiener filter.
    #          For me, this is : /home/s1563094/Datastore/Software/LSDTopoTools/LSDTopoTools_ChannelExtraction/driver_functions_ChannelExtraction
    
    """print "Applying wiener filter"
    sourcedir = "/home/s1563094/Datastore/Software/LSDTopoTools/LSDTopoTools_MarshExtraction/Example/Input/ " # This is the directory where you store your DEM
    sourcefile = "%sdm_DEM " % (site) # This is the filename for your original DEM
    sourceformat = "bil" # DO NOT change the output format!
    os.system("./Wiener_filter.out " + sourcedir + sourcefile +  sourceformat)"""
    
    # NOTE: Applying a Wiener filter does not always yield better results. It is recommended mostly if your DEM is very noisy (e.g. it is not based on ground returns)


    
    
    # Stage 3: This stage is not optional. In this stage, the slope, curvature and hillshade rasters are computed.
    #          Run this stage from the directory where you store the TopoTools analysis pack.
    #          For me, this is: /home/s1563094/Datastore/Software/LSDTopoTools/Git_projects/LSDTopoTools_AnalysisDriver/Analysis_driver
    
    """print "Calculating slopes, curvatures and hillshade"
    sourcedir = "/home/s1563094/Datastore/Software/LSDTopoTools/LSDTopoTools_MarshExtraction/Input/Topography/HIN_degraded/ "
    driverfile = "%s_slope_curv.LSDTT_driver " % (site)
    os.system("./LSDTT_analysis_from_paramfile.out " + sourcedir + driverfile)"""


    
#STOP # Comment this and the rest of point 3. when you have finished preparing your data.
  

    
    
#------------------------------------------------------------------
#4. Prepare a reference raster
# Once you have commented out the previous stages, copy this file back in its original location and run again.
# This part of the script will transform a .shp file (polygon) into a raster which you can use to compare the results of the unsupervised classification to your own digitisation.

for site in Sites:  
    print " Converting reference marsh outline to a raster"
    directory = "Input/" # The directory where you store your refence. By default, it will be in the Input directory
    srcfile = "%s_marsh_DEM.shp" % (site)
    dstfile = "%s_marsh_DEM.bil"  % (site)
    os.system("gdal_rasterize -a Raster_val -of ENVI -a_srs EPSG:27700 -tr 1 1 " + directory+srcfile + " " +  directory+dstfile)
    


    
    
#------------------------------------------------------------------
#5. Clip your topographic data
# This part of the script will clip all your files to a specified domain boundary.

for site in Sites: 
    cutfile = "%s_domain.shp" % (site) # Use the same cutfile ( = domain) for all rasters if you don't want any trouble. 
    directory = "Input/" # The directory where you store your input files. By default, it will be in the Input directory
    
    print " Clipping reference marsh outline raster"
    srcfile = "%s_marsh_DEM.bil" % (site)
    dstfile = "%s_ref_DEM_clip.bil"  % (site)
    os.system("gdalwarp -overwrite -of ENVI -cutline " + directory+cutfile + " -crop_to_cutline " + directory+srcfile + " " +  directory+dstfile)

    
    STOP
    
    print " Clipping DEM raster"
    srcfile = "%s_DEM.bil" % (site)
    dstfile = "%s_DEM_clip.bil" % (site)
    os.system("gdalwarp -of ENVI -t_srs EPSG:27700 -cutline " + cutfile + " -crop_to_cutline " + directory+srcfile + " " +  directory+dstfile)
    
    print " Clipping slope raster"
    srcfile = "%s_slope.bil" % (site)
    dstfile = "%s_slope_clip.bil" % (site)
    os.system("gdalwarp -of ENVI -t_srs EPSG:27700 -cutline " + cutfile + " -crop_to_cutline " + directory+srcfile + " " +  directory+dstfile)
    
    print " Clipping curvature raster"
    srcfile = "%s_curvature.bil" % (site)
    dstfile = "%s_curvature_clip.bil" % (site)
    os.system("gdalwarp -of ENVI -t_srs EPSG:27700 -cutline " + cutfile + " -crop_to_cutline " + directory+srcfile + " " +  directory+dstfile)
    
    print " Clipping hillshade raster"
    srcfile = "%s_hs.bil" % (site)
    dstfile = "%s_hs_clip.bil" % (site)
    os.system("gdalwarp -of ENVI -t_srs EPSG:27700 -cutline " + cutfile + " -crop_to_cutline " + directory+srcfile + " " +  directory+dstfile)

