# Load useful Python packages
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal, osr, gdalconst
from osgeo.gdalconst import *
import cPickle



#---------------------------------------------------------------
def ENVI_raster_binary_to_2d_array(file_name, gauge):
    """
    This function transforms a raster into a numpy array.
    
    Args:
        file_name (ENVI raster): the raster you want to work on.
        gauge (string): a name for your file
    
    Returns:
        image_array (2-D numpy array): the array corresponding to the raster you loaded
        pixelWidth (geotransform, inDs) (float): the size of the pixel corresponding to an element in the output array.
    
    Source: http://chris35wills.github.io/python-gdal-raster-io/
    """ 
    
    
    print 'Opening %s' % (gauge)

    driver = gdal.GetDriverByName('ENVI')

    driver.Register()

    inDs = gdal.Open(file_name, GA_ReadOnly)

    if inDs is None:
        print "Couldn't open this file: " + file_name
        print "Perhaps you need an ENVI .hdr file? "
        sys.exit("Try again!")
    else:
        print "%s opened successfully" %file_name

        #print '~~~~~~~~~~~~~~'
        #print 'Get image size'
        #print '~~~~~~~~~~~~~~'
        cols = inDs.RasterXSize
        rows = inDs.RasterYSize
        bands = inDs.RasterCount

        #print "columns: %i" %cols
        #print "rows: %i" %rows
        #print "bands: %i" %bands

        #print '~~~~~~~~~~~~~~'
        #print 'Get georeference information'
        #print '~~~~~~~~~~~~~~'
        geotransform = inDs.GetGeoTransform()
        originX = geotransform[0]
        originY = geotransform[3]
        pixelWidth = geotransform[1]
        pixelHeight = geotransform[5]

        #print "origin x: %i" %originX
        #print "origin y: %i" %originY
        #print "width: %2.2f" %pixelWidth
        #print "height: %2.2f" %pixelHeight

        # Set pixel offset.....
        print '~~~~~~~~~~~~~~'
        print 'Convert image to 2D array'
        print '~~~~~~~~~~~~~~'
        band = inDs.GetRasterBand(1)
        image_array = band.ReadAsArray(0, 0, cols, rows)
        image_array_name = file_name
        print type(image_array)
        print image_array.shape

        return image_array, pixelWidth, (geotransform, inDs)




#---------------------------------------------------------------------
def ENVI_raster_binary_from_2d_array(envidata, file_out, post, image_array):
    """
    This function transforms a numpy array into a raster.
    
    Args:
        envidata: the geospatial data needed to create your raster
        file_out (string): the name of the output file
        post: coordinates for the goegraphical transformation
        image_array (2-D numpy array): the input raster
    
    Returns:
        new_geotransform
        new_projection: the projection in which the raster
        file_out (ENVI raster): the raster you wanted
    
    Source: http://chris35wills.github.io/python-gdal-raster-io/
    """ 
    driver = gdal.GetDriverByName('ENVI')

    original_geotransform, inDs = envidata

    rows, cols = image_array.shape
    bands = 1

    # Creates a new raster data source
    outDs = driver.Create(file_out, cols, rows, bands, gdal.GDT_Float32)

    # Write metadata
    originX = original_geotransform[0]
    originY = original_geotransform[3]

    outDs.SetGeoTransform([originX, post, 0.0, originY, 0.0, -post])
    outDs.SetProjection(inDs.GetProjection())

    #Write raster datasets
    outBand = outDs.GetRasterBand(1)
    outBand.WriteArray(image_array)

    new_geotransform = outDs.GetGeoTransform()
    new_projection = outDs.GetProjection()

    print "Output binary saved: ", file_out

    return new_geotransform,new_projection,file_out



#-----------------------------------------------------------------------------------------------------------
def Distribution(Data2D, Nodata_value):
    """
    This simple function takes a 2-D array (Data2D) and makes a probability distribution of its values. It is set to ignore elements with a specific value (Nodata_value).
    
    Args:
        Data2D (2D numpy array): the 2D array you want a distribution for
        Nodata_value (float): The value for ignored elements
    
    Returns:
        bins [1D numpy array]: the value bins
        hist [1D numpy array]: the probability associated to the bins
    
    Author: GCHG
    """ 

    Data1D = Data2D.ravel()

    Max_distribution = max(Data1D)    
    if len(Data1D[Data1D>Nodata_value]) == 0:
        Min_distribution = -1
    else:
        Min_distribution = min(Data1D[Data1D>Nodata_value])
    
    bin_size = (Max_distribution - Min_distribution) / 100
    
    X_values = np.arange(Min_distribution, Max_distribution, bin_size)
    
    
    hist, bins = np.histogram (Data1D, X_values, density=True)
    hist=hist/sum(hist)
    bins=bins[:-1]
   
    
    return bins,hist



#-----------------------------------------------------------------------------------------------------------


def Outline (Raster, Outline_value, Nodata_value):
    """
    This simple function takes a 2-D array (Raster) and attributes a specific value (Outline value) to elements at the limit of a bloc of elements with identical values. Effectively, it draws an outline around a group of elements with the same value. It is set to ignore elements with a specific value (Nodata_value).
    
    Args:
        Raster (2D numpy array): the 2-D array
        Outline_value (float): The value associated to the outline. Be smart and select a different value from those already in your 2-D array.
        Nodata_value (float): The value for ignored elements
    
    Returns:
        Raster (2D numpy array): the 2-D array, with the outlines given their own value.
    
    Author: GCHG
    """ 
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




#-----------------------------------------------------------------------------------------------------
def define_search_space (DEM, Slope, Nodata_value, opt): 
    """
   This function defines a search space (Search_space) within a 2-D array, based on the combined values of 2 2-D arrays (DEM and Slope) of the same dimensions. It defines the threshold for the selection of the search space according to a threshold value (opt). It is set to ignore elements with a specific value (Nodata_value).
    Args:
        DEM (2D numpy array): a 2-D array (here a DEM) used as a first condition for the definition of the search space
        Slope (2D numpy array): a 2-D array (here a DEM) used as a second condition for the definition of the search space
        Nodata_value (float): The value for ignored elements
        opt (float): the value of the threshold for the selection of the search space
    
    Returns:
        Search_space (2D numpy array): The resulting search space array. Search_space has a value of 0 for non-selected elements and 1 for selected elements.
        Crossover (2D numpy array): The array resulting of the multiplication of relative slope and relative relief.
        bins (1D array): the value bins for the Crossover array
        hist (1D array): the value hist for the Crossover array
        Inflecion_point(float): the value of the threshold for the search space selection.
    
    Author: GCHG
    """ 
    
    
    print 'Choosing a holiday destination ...'
    Height = len(DEM); Width = len(DEM[0,:])
    Search_space = np.zeros((Height,Width), dtype=np.float)

    # We calculate the relative relief of the DEM to have values of elevation between 0 and 1
    Relief = DEM-np.amin(DEM[DEM > Nodata_value])
    Rel_relief = Relief/np.amax(Relief)
    Rel_relief[DEM == Nodata_value] = Nodata_value

    # We then do the same thing for slope
    Rel_slope = Slope/np.amax(Slope)
    Rel_slope[Slope == Nodata_value] = Nodata_value

    # We then multiply these new relative relief and slope arrays and biologically name them "Crossover"
    Crossover = Rel_relief * Rel_slope
    Crossover[DEM == Nodata_value] = Nodata_value

    # We make a curve of the frequency of values in this Crossover
    # That curve should look like a decreasing exponential function
    data = Crossover.ravel(); data = data[data>0]
    step = (max(data) - min(data)) / 100
    value = np.arange(min(data), max(data), step)
    hist, bins = np.histogram (data, value, density=True)
    hist=hist/sum(hist); bins=bins[:-1]

    # We now find the slope of that curve
    hist_der = np.zeros(len(hist), dtype = np.float)
    for j in range(1, len(hist), 1):
        hist_der[j] = (hist[j]-hist[j-1])/step

    # If the slope gets above the -1 threshold, now that we have hit the closest point to the origin.
    # We call it the inflexion point even though it's not really an inflexion point.
    for j in range(1, len(hist)-1, 1):
        if hist_der[j] < opt and hist_der[j+1] >= opt:
            Inflexion_point = bins[j]

    # Points within the search space should have a Crossover value above the inflexion point
    Search = np.where(Crossover > Inflexion_point)
    Search_space[Search] = 1

    # We get rid of the borders of the DEM because otherwise it will be difficult to work with the smaller slope array
    Search_space[0,:] = 0; Search_space[Height-1,:] = 0; Search_space[:,0] = 0; Search_space[:,Width-1] = 0

    # And update the search locations for the shaved edges
    Search = np.where(Search_space == 1)

    # If this happens, your landscape is weird
    if np.amax(Search_space) == 0:
        print
        print " ... Your search space is empty! Are you sure there's a marsh platform here?"
        print
        STOP

    return Search_space, Crossover, bins, hist, Inflexion_point


#-----------------------------------------------------------------------------------------------------
def kernel (array, kernel_size, x_centre, y_centre):
    """
    This function defines a square kernel within an array (array), centred on (x_centre, y_centre). The is of a width of kernel_size.
    Args:
        array (2D numpy array): a 2-D array.
        kernel_size (float): the width of the square defining the size of the kernel. kernel_size MUST be an ODD number to account for the central element.
        x_centre (int): The index of the element in the 1st dimension.
        y_centre (int): The index of the element in the 2nd dimension.

    Returns:
        kernel (2D numpy array): The kernel of selected elements.

    Author: GCHG
    """
 
    if (-1)**kernel_size < 0:
        X_to_0 = x_centre
        X_to_End = len(array)-x_centre
        Y_to_0 = y_centre
        Y_to_End = len(array[0,:])-y_centre

        Lim_left = x_centre - min(np.floor(kernel_size/2), X_to_0)
        Lim_right = x_centre + min(np.floor(kernel_size/2)+1, X_to_End)
        Lim_top = y_centre - min(np.floor(kernel_size/2), Y_to_0)
        Lim_bottom = y_centre + min(np.floor(kernel_size/2)+1, Y_to_End)

        kernel = array [int(Lim_left):int(Lim_right), int(Lim_top):int(Lim_bottom)]

    else:
        print
        print " ... WARNING: you need to choose an odd kernel size, buddy"
        print
        pass

    return kernel


#-----------------------------------------------------------------------------------------------------
def peak_flag (Slope, Search_space, Order):
    """
    This function is the first stage of a routing process used to identify lines of maximum slopes.
    This function identifies multiple local maxima in an array (Slope), within a predefined search space (Search_space). The identified maxima are given a value of Order.
    
    Args:
        Slope (2D numpy array): the input 2-D array, here issued from a slope raster.
        Search_space (2D numpy array): the search space array in which to look for local maxima.
        Order (int): the value given to the local maxima points.

    Returns:
        Peaks (2D numpy array): a 2-D array where the local maxima have a value of Order and other elements are null.
        Slope_copy (2D numpy array): a copy of the input array where the value of the selected local maxima has been set to 0.

    Author: GCHG
    """
    
    print 'Finding local slope maxima ...'
    Slope_copy = np.copy(Slope) # the copy of the initial data array
    Search = np.where(Search_space == 1) # the searched locations
    Peaks = np.zeros((len(Slope),len(Slope[0,:])),dtype = np.float)

    for i in range(len(Search[0])):
        x=Search[0][i]; y=Search[1][i] # coordinates of the kernel's centre
        Kernel_slope = kernel (Slope, 3, x, y)
        Kernel_search = kernel(Search_space, 3, x, y)

        # if the centre of the kernel is its maximum and is not an isolated point
        if Kernel_slope[1,1] == np.amax(Kernel_slope) and np.amax(Kernel_search[Kernel_search<=Kernel_search[1,1]] > 0):
            Peaks[x,y] = Order # The kernel centre becomes a local peak
            Slope_copy[x,y] = 0 # The slope of the modified data array drops to 0

    return Peaks, Slope_copy



#-----------------------------------------------------------------------------------------------------
def initiate_ridge (Slope, Search_space, Peaks, Order):
    """
    This function is the second stage of a routing process used to identify lines of maximum slopes.
    This function identifies multiple duplets of elements in an array (Slope), within a predefined search space (Search_space) and within the neighbourhood of the local maxima identified in a second input array (Peaks). The identified elements are given a value of Order. To make this function work, the input array Slope should be the output array Slope_copy of the function peak_flag.
    
    Args:
        Slope (2D numpy array): the input 2-D array, here issued from a slope raster where the local maximal values have been replaced by 0. 
        Search_space (2D numpy array): the search space array.
        Peaks (2D numpy array): A 2-D array containing elements with a value of 1. These elements have the same indices as the elements with a value of 0 in Slope.
        Order (int): the value given to the identified elements. it should be superior by 1 to the value of Order in the function peak_flag.

    Returns:
        Ridges (2D numpy array): a 2-D array where the identified elements have a value of Order. This array is modified from the Peaks array and therefore also contains elements of a value equal to the Order in the function peak_flag.
        Slope_copy (2D numpy array): a copy of the input array where the value of the selected elements has been set to 0.

    Author: GCHG
    """
    
    print ' ... Starting ridges ...'
    Slope_copy = np.copy(Slope) # the copy of the initial data array
    Search = np.where(Search_space == 1) # the searched locations
    Search_peaks = np.where(Peaks == Order-1) # the searched locations where the peaks are
    Ridges = np.copy(Peaks)

    # Define Kernels
    for i in range(len(Search_peaks[0])):
        x=Search_peaks[0][i]; y=Search_peaks[1][i] # coordinates of the kernel's centre
        Kernel_slope = kernel (Slope, 3, x, y)
        Kernel_slope_copy = kernel (Slope_copy, 3, x, y)
        Kernel_ridges = kernel (Ridges, 3, x, y)
        Kernel_search = kernel (Search_space, 3, x, y)

        # 1/ If there are no other peaks, we have two ridge starters
        if np.count_nonzero(Kernel_ridges) == 1:
            Ridge_starter1 = np.where (Kernel_slope_copy == np.amax (Kernel_slope_copy))
            X1=Ridge_starter1[0][0]; Y1=Ridge_starter1[1][0]

            # if it is within the initial search space
            if Search_space[x+X1-1, y+Y1-1] != 0:
                Ridges[x+X1-1, y+Y1-1] = Order
                Slope_copy[x+X1-1, y+Y1-1] = 0

                # Look for a second ridge starter
                Ridge_starter2 = np.where (Kernel_slope_copy == np.amax (Kernel_slope_copy))
                X2=Ridge_starter2[0][0]; Y2=Ridge_starter2[1][0]
                Distance = np.sqrt((X2-X1)**2+(Y2-Y1)**2)

                # if it is within the initial search space AND not next to the first ridge starter
                if Search_space[x+X2-1, y+Y2-1] != 0 and Distance > np.sqrt(2):
                    Ridges[x+X2-1, y+Y2-1] = Order
                    Slope_copy[x+X2-1, y+Y2-1] = 0

                # Otherwise, look for second ridge starter elsewhere in the kernel
                elif Search_space[x+X2-1, y+Y2-1] != 0 and Distance <= np.sqrt(2):
                    for j in np.arange(0,9,1):
                        Kernel_slope_copy[X2, Y2] = 0

                        Ridge_starter2 = np.where (Kernel_slope_copy == np.amax (Kernel_slope_copy))
                        X2=Ridge_starter2[0][0]; Y2=Ridge_starter2[1][0]
                        Distance = np.sqrt((X2-X1)**2+(Y2-Y1)**2)

                        if Search_space[x+X2-1, y+Y2-1] != 0 and Distance > np.sqrt(2):
                            Ridges[x+X2-1, y+Y2-1] = Order
                            Slope_copy[x+X2-1, y+Y2-1] = 0
                            break


        # 2/ If there are two peaks, we have one ridge starter
        elif np.count_nonzero(Kernel_ridges) == 2:
            Ridge_starter1 = np.where (Kernel_slope_copy == np.amax (Kernel_slope_copy))
            X1=Ridge_starter1[0][0]; Y1=Ridge_starter1[1][0]

            # if it is within the initial search space
            if Search_space[x+X1-1, y+Y1-1] != 0:
                Ridges[x+X1-1, y+Y1-1] = Order
                Slope_copy[x+X1-1, y+Y1-1] = 0

    return Ridges, Slope_copy






#-----------------------------------------------------------------------------------------------------
def Continue_ridge (Slope, Search_space, Peaks, Order):
    """
    This function is the third and final stage of a routing process used to identify lines of maximum slopes.
    IMPORTANT: this function is meant to be run several times! It requires the incrementation of the Order value with each iteration.
    This function identifies multiple elements in an array (Slope), within a predefined search space (Search_space) and within the neighbourhood of the local maxima identified in a second input array (Peaks).  The identified elements are given a value of Order. To make this function work, the input array Slope should be the output array Slope_copy of the function initiate_ridge.

    Args:
        Slope (2D numpy array): the input 2-D array, here issued from a slope raster where the elements selected in the initiate_ridge function have been replaced by 0. 
        Search_space (2D numpy array): the search space array.
        Peaks (2D numpy array): A 2-D array containing elements with a value of 1. These elements have the same indices as the elements with a value of 0 in Slope.
        Order (int): the value given to the identified elements. On the first iteration it should be superior by 1 to the value of Order in the function initiate_ridge. the value of Order then needs to be incremented with every iteration.

    Returns:
        Ridges (2D numpy array): a 2-D array where the identified elements have a value of Order. This array is modified from the Peaks array and therefore also contains elements of a value equal to the Order in the functions peak_flag and initiate_ridge.
        Slope_copy (2D numpy array): a copy of the input array where the value of the selected elements has been set to 0.

    Author: GCHG
    """

    print ' ... Prolongating ridges ...'
    Slope_copy = np.copy(Slope) # the copy of the initial slope array
    Search = np.where(Search_space == 1) # the searched locations
    Search_peaks = np.where(Peaks == Order-1) # the searched locations where the peaks are
    Ridges = np.copy(Peaks)

    # Define Kernels
    for i in range(len(Search_peaks[0])):
        x=Search_peaks[0][i]; y=Search_peaks[1][i] # coordinates of the kernel's centre

        Kernel_slope = kernel (Slope, 3, x, y)
        Kernel_slope_copy = kernel (Slope_copy, 3, x, y)
        Kernel_ridges = kernel (Ridges, 3, x, y)
        Kernel_search = kernel (Search_space, 3, x, y)

        # Count the number of nonzero points in the kernel of the ridge array
        Ridge_count = np.count_nonzero(Kernel_ridges)

        # If there are only the 2 previous ridge points, draw a third point that is far enough from the previous point
        if Ridge_count == 2:
            New_point = np.where (Kernel_slope_copy == np.amax (Kernel_slope_copy))
            X=New_point[0][0]; Y=New_point[1][0]
            Grandad_point = np.where (Kernel_ridges == Order-2)
            Xgd=Grandad_point[0][0]; Ygd=Grandad_point[1][0]
            Distance = np.sqrt((X-Xgd)**2+(Y-Ygd)**2)

            if Search_space[x+X-1, y+Y-1] != 0 and Distance > np.sqrt(2):
                Ridges[x+X-1, y+Y-1] = Order
                Slope_copy[x+X-1, y+Y-1] = 0

            elif Search_space[x+X-1, y+Y-1] != 0 and Distance <= np.sqrt(2):
                for j in np.arange(0,9,1):
                    Kernel_slope_copy[X, Y] = 0

                    New_point = np.where (Kernel_slope_copy == np.amax (Kernel_slope_copy))
                    X=New_point[0][0]; Y=New_point[1][0]
                    Distance = np.sqrt((X-Xgd)**2+(Y-Ygd)**2)

                    if Search_space[x+X-1, y+Y-1] != 0 and Distance > np.sqrt(2):
                        Ridges[x+X-1, y+Y-1] = Order
                        Slope_copy[x+X-1, y+Y-1] = 0
                        break

    return Ridges, Slope_copy




#-----------------------------------------------------------------------------------------------------
def Clean_ridges (Peaks, DEM, Nodata_value, opt):
    """
    This function eliminates some of the ridges (Peaks) identified by the trio of functions (peak_flag, initiate_ridge and continue_ridge). The elimination process depends on local relief, which uses a DEM (DEM) and a threshold value (opt). It is set to ignore elements with a value of  Nodata_value.

    Args:
        Peaks (2D numpy array): the input 2-D arraym which is the output of the ridge identification process.
        DEM (2D numpy array): the DEM array used as a base for the elimination of unnecessary ridges.
        Nodata_value (float): The value for ignored elements.
        opt (float): The value of the threshold to eliminate unnecessary ridges.

    Returns:
        Peaks (2D numpy array): a 2-D array much like the input Peaks array, but the unnecessary elemets have been reset to 0.
        
    Author: GCHG
    """

    print "Cleaning up ridges ..."
    DEM_copy = np.copy(DEM)
    DEM_copy[DEM_copy==Nodata_value] = 0
    Search_ridge = np.where (Peaks != 0)

    Cutoff = np.percentile(DEM_copy,75)
    Threshold = np.amax(DEM_copy[DEM_copy<Cutoff])
    DEM_copy[DEM_copy>Threshold]=Threshold

    for i in range(len(Search_ridge[0])):
        x=Search_ridge[0][i]; y=Search_ridge[1][i] # coordinates of the kernel's centre
        Kernel_DEM = kernel (DEM_copy, 9, x, y)
        Kernel_DEM[Kernel_DEM==Nodata_value]=0

        if np.amax(Kernel_DEM)/Threshold < opt:
            Peaks[x,y] = 0

    Search_ridge = np.where (Peaks != 0)
    for i in range(len(Search_ridge[0])):
        x=Search_ridge[0][i]; y=Search_ridge[1][i] # coordinates of the kernel's centre
        Kernel_ridges = kernel (Peaks, 9, x, y)
        # If there aren't at least 8 ridge points in the neighbourhood of 10 by 10
        if np.count_nonzero(Kernel_ridges) < 8:
            Peaks[x,y] = 0

    return Peaks




#-----------------------------------------------------------------------------------------------------
def Fill_marsh (DEM, Peaks, Nodata_value, opt):
    """
    This function builds a marsh platform array by using the Peaks array as a starting point. It uses the DEM array to establish conditions on the elements to select. the opt parameter sets a threshold value to eliminate superfluous elements. It is set to ignore elements with a value of Nodata_value.

    Args:
        DEM (2D numpy array): the DEM array.
        Peaks (2D numpy array): the 2-D array of ridge elements, which is the output of the ridge identification and cleaning process.
        Nodata_value (float): The value for ignored elements.
        opt (float): The value of the threshold to eliminate unnecessary elements.

    Returns:
        Marsh (2D numpy array): a 2-D array where the marsh platform elements are identified by strictly positive values. Other elements have a valuof 0 or Nodata_value.
        
    Author: GCHG
    """
    
    print "Initiate platform ..."
    DEM_copy = np.copy(DEM)
    Marsh = np.zeros((len(DEM), len(DEM[0,:])), dtype = np.float)

    Counter = 1
    Search_ridges = np.where (Peaks > 0)
    for i in range(len(Search_ridges[0])):
        x=Search_ridges[0][i]; y=Search_ridges[1][i]
        Kernel_ridges = kernel (Peaks, 3, x, y)
        Kernel_DEM = kernel (DEM, 3, x, y)

        Marsh_point = np.where (np.logical_and (Kernel_DEM >= Kernel_DEM[1,1], Kernel_ridges == 0))
        for j in range(len(Marsh_point[0])):
            X=Marsh_point[0][j]; Y=Marsh_point[1][j]
            Marsh[x+X-1, y+Y-1] = Counter

    Search_marsh_start = np.where (Marsh == 1)
    for i in range(len(Search_marsh_start[0])):
        x=Search_marsh_start[0][i]; y=Search_marsh_start[1][i]
        Kernel_marsh = kernel (Marsh, 3, x, y)
        Kernel_ridges = kernel (Peaks, 3, x, y)
        if np.count_nonzero(Kernel_marsh) <=2:
            Marsh[x,y] = 0

    print ' ... Build platform ...'
    while Counter < 100:
        Counter = Counter+1
        Search_marsh = np.where (Marsh == Counter-1)
        for i in range(len(Search_marsh[0])):
            x = Search_marsh[0][i]; y = Search_marsh[1][i]
            Kernel_DEM = kernel (DEM, 3, x, y)
            Kernel_DEM_copy = kernel (DEM_copy, 3, x, y)
            Kernel_ridges = kernel (Peaks, 3, x, y)
            Kernel_marsh = kernel (Marsh, 3, x, y)
            Big_Kernel_DEM = kernel (DEM, 11, x, y)
            Big_Kernel_DEM_copy = kernel (DEM_copy, 11, x, y)
            

            Conditions = np.zeros((len(Kernel_DEM), len(Kernel_DEM[0,:])), dtype = np.float)
            # 1: free space
            Condition_1 = np.where (np.logical_and(Kernel_ridges == 0, Kernel_marsh == 0)); Conditions[Condition_1] = 1
            # 2: not topped
            Condition_2 = np.where (np.logical_and(Kernel_DEM_copy > np.amax(Big_Kernel_DEM_copy)-0.2, Conditions == 1)); Conditions[Condition_2] = 2
            
            
            #This is a distance thing to make sure you don't cross the ridges agin
            Here_be_ridges = np.where (Kernel_ridges != 0)
            Here_be_parents = np.where (Kernel_marsh == Counter-1)

            for j in range(len(Condition_2[0])):
                X=Condition_2[0][j]; Y=Condition_2[1][j]
                Distance_to_ridges = []
                Distance_to_parents = []

                for k in range(len(Here_be_ridges[0])):
                    Xr=Here_be_ridges[0][k]; Yr=Here_be_ridges[1][k]
                    Distance = np.sqrt((X-Xr)**2+(Y-Yr)**2)
                    Distance_to_ridges.append(Distance)

                for k in range(len(Here_be_parents[0])):
                    Xp=Here_be_parents[0][k]; Yp=Here_be_parents[1][k]
                    Distance = np.sqrt((X-Xp)**2+(Y-Yp)**2)
                    Distance_to_parents.append(Distance)

                if len(Distance_to_ridges)>0:
                    if min(Distance_to_ridges) > min(Distance_to_parents):
                        Marsh[x+X-1, y+Y-1] = Counter
                else:
                    Marsh[x+X-1, y+Y-1] = Counter
                    DEM_copy[x+X-1, y+Y-1] = 0

    
    print ' ... defining the elimination of low platforms ...'
    Platform = np.copy(Marsh)
    Platform[Platform > 0] = DEM [Platform > 0]
    Platform_bins, Platform_hist = Distribution(Platform,0)

    #1. Find the highest and biggest local maximum of frequency distribution
    # Initialize Index
    Index = len(Platform_hist)-1
    # Initiate Cutoff_Z value
    Cutoff_Z = 0
    
    for j in range(1,len(Platform_hist)-1):
        if Platform_hist[j]>0.9*max(Platform_hist) and Platform_hist[j]>Platform_hist[j-1] and Platform_hist[j]>Platform_hist[j+1]:
            Index  = j

    #2. Now run a loop from there toward lower elevations.
    Counter = 0
    for j in range(Index,0,-1):
        # See if you cross the mean value of frequency. Count for how many indices you are under.
        if Platform_hist[j] < np.mean(Platform_hist):
            Counter = Counter + 1
        # Reset the counter value if you go above average again
        else:
            Counter = 0 
    
        #If you stay long enough under (10 is arbitrary for now), initiate cutoff and stop the search
        if Counter > opt:
            Cutoff = j
            Cutoff_Z = Platform_bins[Cutoff]
            break
        
    # If you stay under for more than 5, set a Cutoff_Z value but keep searching    
    if Counter > opt/2:
        Cutoff = j
        Cutoff_Z = Platform_bins[Cutoff]

    Marsh[Platform<Cutoff_Z] = 0

 
    print " ... Fill high areas left blank ..."
    Search_marsh_condition = np.zeros((len(DEM), len(DEM[0,:])), dtype = np.float)
    Search_marsh = np.where (DEM >= Platform_bins[Index])
    Search_marsh_condition [Search_marsh] = 1
    Search_marsh_2 = np.where (np.logical_and(Marsh == 0, Search_marsh_condition == 1))
    Marsh[Search_marsh_2] = 3

    print ' ... Fill the interior of pools ...'
    for Iteration in np.arange(0,10,1):
        Counter = 100
        while Counter > 2:
            Counter = Counter-1
            Search_marsh = np.where (Marsh == Counter+1)
            Non_filled = 0
            for i in range(len(Search_marsh[0])):
                x = Search_marsh[0][i]; y = Search_marsh[1][i]
                Kernel_DEM = kernel (DEM, 3, x, y)
                Kernel_ridges = kernel (Peaks, 3, x, y)
                Kernel_marsh = kernel (Marsh, 3, x, y)

                if Non_filled <len(Search_marsh[0]):
                    if np.count_nonzero(Kernel_marsh) > 6:
                        Condition = np.where (np.logical_and(Kernel_marsh == 0, Kernel_ridges == 0))
                        for j in range(len(Condition[0])):
                            X=Condition[0][j]; Y=Condition[1][j]
                            Marsh[x+X-1, y+Y-1] = Counter
                    else:
                        Non_filled = Non_filled + 1
                        
    # Reapply the cutoff because the straight line thing is ugly
    Platform = np.copy(Marsh)
    Platform[Platform > 0] = DEM [Platform > 0]
    Marsh[Platform<Cutoff_Z] = 0

 

    # We fill in the wee holes
    Search_marsh = np.where (np.logical_and(Marsh == 0, Peaks == 0))
    for i in range(len(Search_marsh[0])):
        x = Search_marsh[0][i]; y = Search_marsh[1][i]
        Kernel_marsh = kernel (Marsh, 3, x, y)
        if np.count_nonzero(Kernel_marsh) == 8:
            Marsh[x,y] = 105

            
            
    print ' ... Adding the ridges'
    # We get rid of scarps that do not have a marsh next to them
    Search_false_scarp = np.where (Peaks > 0)
    for i in range(len(Search_false_scarp[0])):
        x = Search_false_scarp[0][i]; y = Search_false_scarp[1][i]
        Kernel_marsh = kernel (Marsh, 3, x, y)
        if np.count_nonzero (Kernel_marsh) == 0:
            Peaks[x, y] = 0

    # We get rid of the sticky-outy bits
    Search_ridge = np.where (Peaks > 0)
    for i in range(len(Search_ridge[0])):
        x=Search_ridge[0][i]; y=Search_ridge[1][i]
        Kernel_ridges = kernel (Peaks, 9, x, y)
        if np.count_nonzero(Kernel_ridges) < 8:
            Peaks[x,y] = 0
   
    # We put the scarps in the platform
    Search_side = np.where (Peaks > 0)
    Marsh[Search_side] = 110

    print " ... eliminate patches of empty elements ..."
    Search_marsh_condition = np.zeros((len(DEM), len(DEM[0,:])), dtype = np.float)
    Search_marsh = np.where (DEM >= Platform_bins[Index])
    Search_marsh_condition [Search_marsh] = 1
    Search_marsh_2 = np.where (np.logical_and(Marsh == 0, Search_marsh_condition == 1))
    Marsh[Search_marsh_2] = 3
    
    print ' ... Fill the interior of pools ...'
    for Iteration in np.arange(0,10,1):
        Counter = 110
        while Counter > 2:
            Counter = Counter-1
            Search_marsh = np.where (Marsh == Counter+1)
            Non_filled = 0
            for i in range(len(Search_marsh[0])):
                x = Search_marsh[0][i]; y = Search_marsh[1][i]
                Kernel_DEM = kernel (DEM, 3, x, y)
                Kernel_ridges = kernel (Peaks, 3, x, y)
                Kernel_marsh = kernel (Marsh, 3, x, y)

                if Non_filled <len(Search_marsh[0]):
                    if np.count_nonzero(Kernel_marsh) > 6:
                        Condition = np.where (np.logical_and(Kernel_marsh == 0, Kernel_ridges == 0))
                        for j in range(len(Condition[0])):
                            X=Condition[0][j]; Y=Condition[1][j]
                            Marsh[x+X-1, y+Y-1] = Counter
                    else:
                        Non_filled = Non_filled + 1
                        
    print ' ... defining the elimination of low platforms ...'
    Platform = np.copy(Marsh)
    Platform[Platform > 0] = DEM [Platform > 0]
    Marsh[Platform<Cutoff_Z] = 0

    Marsh[DEM == Nodata_value] = Nodata_value

    return Marsh

 



#---------------------------------------------------------------
def MARSH_ID (DEM, Slope, Nodata_value, opt1, opt2, opt3):
    """
    This is the master function for marsh identification. It defines in which order the functions define_search_space, peak_flag, initiate_ridge, Continue_ridge, Clean_ridges, Fill_marsh are executed. It is set to repeat the iteration of the Continue_ridge function 50 times.

    Args:
        DEM (2D numpy array): the input DEM array.
        Slope (2D numpy array): the input Slope array.
        Nodata_value (float): The value for ignored elements.
        opt1 (float): The value of the threshold used in the define_search_space function.
        opt2 (float): The value of the threshold used in the Clean_ridges function.
        opt3 (float): The value of the threshold used in the Fill_marsh function.

    Returns:
        Search_space (2D numpy array): The output search space of the define_search_space function.
        Ridge (2D numpy array): The output ridges of the peak_flag, initiate_ridge, Continue_ridge, Clean_ridges functions.
        Marsh (2D numpy array): The output marsh platform of the Fill_marsh function.
        
    Author: GCHG
    """
 
    DEM_work = np.copy(DEM); Slope_work = np.copy(Slope);

    Platform = np.copy(DEM_work)
    Ridge = np.copy(DEM_work)
    Marsh = np.copy(DEM_work)

    Platform[Platform != Nodata_value] = 0
    Summit = np.where (Platform==np.amax(Platform))
    Platform[Summit] = 1


    Search_space, Crossover, bins, hist, Inflexion_point = define_search_space (DEM_work, Slope_work, Nodata_value,opt1)

    Order = 1
    Ridge, Slope_temp = peak_flag (Slope_work, Search_space, Order)

    Order = Order+1
    Ridge, Slope_temp = initiate_ridge (Slope_temp, Search_space, Ridge, Order)

    while Order < 50:
        Order = Order+1
        Ridge, Slope_temp = Continue_ridge (Slope_temp, Search_space, Ridge, Order)

    Ridge = Clean_ridges (Ridge, DEM_work, Nodata_value, opt2)

    Marsh = Fill_marsh (DEM_work, Ridge, Nodata_value, opt3)


    print "My hovercraft is full of eels!"
    print


    return Search_space, Ridge, Marsh




#-----------------------------------------------------------------------------------------------------
def Confusion (Subject, Reference, Nodata_value):
    """
    This function compares a Subject 2-D array to a Reference 2-D array and returns an array of differences, which we call a confusion array or confusion map if it look like a map. It then calculates a number of metrics relative to the adequation between the subject and the reference. It is set to ignore elements with a value of Nodata_value.
    
    To learn more about confusion matrices and their associated metrics, please visit the Wikipedia page: https://en.wikipedia.org/wiki/Confusion_matrix
    
    Args:
        Subject (2D numpy array): the input array. This is the one you want to test
        Reference (2D numpy array): the reference array. This one is supposed to contain correct information
        Nodata_value (float): The value for ignored elements.

    Returns:
        Confusion_matrix (2D numpy array): an array containing the values 1 (True Positive), 2 (True Negative), -1 (False Positive) and -2 (False Negative).
        Performance (1D numpy array): the number of (respectively) True Positives, True Negatives, False Positives and False Negatives in Confusion_matrix.
        Metrix (1D numpy array): The values of (respectively) Accuracy, Reliability, Sensitivity, F1 derived from the Performance array.
        
    Author: GCHG
    """
    
    Height = len(Subject[:,0]); Width = len(Subject[0,:])
    Height_R = len(Reference[:,0]); Width_R = len(Reference[0,:])

    print Height, Width
    print Height_R, Width_R
   
    H = min (Height, Height_R)
    W = min (Width, Width_R)

    Confusion_matrix = Nodata_value*np.ones((Height, Width), dtype = np.float)

    Subject_marsh = np.where (np.logical_and(Subject != 0, Subject != Nodata_value))
    Reference_marsh = np.where (np.logical_and(Reference != 0, Reference != Nodata_value))

    Subject[Subject_marsh] = 1.
    Reference[Reference_marsh] = 1.

    for i in range (H):
        for j in range (W):
            if Subject[i,j] == 1 and Reference[i,j] == 1: # TRUE POSITIVE
                Confusion_matrix[i,j] = 1
            elif Subject[i,j] == 0 and Reference[i,j] == 0: # TRUE NEGATIVE
                Confusion_matrix[i,j] = 2
            elif Subject[i,j] == 1 and Reference[i,j] == 0: # FALSE POSITIVE
                Confusion_matrix[i,j] = -1
            elif Subject[i,j] == 0 and Reference[i,j] == 1: # FALSE NEGATIVE
                Confusion_matrix[i,j] = -2

    True_positive = np.sum(Confusion_matrix[Confusion_matrix == 1])
    True_negative = np.sum(Confusion_matrix[Confusion_matrix == 2])/2
    False_positive = -np.sum(Confusion_matrix[Confusion_matrix == -1])
    False_negative = -np.sum(Confusion_matrix[Confusion_matrix == -2])/2

    Reliability = True_positive / (True_positive+False_positive)
    Sensitivity = True_positive / (True_positive+False_negative)
    Accuracy = (True_positive+True_negative) / (True_positive+True_negative+False_positive+False_negative)
    F1 = 2*True_positive/(2*True_positive+False_positive+False_negative)

    Performance = np.array([True_positive,True_negative,False_positive,False_negative])
    Metrix = np.array([Accuracy, Reliability, Sensitivity, F1])


    return Confusion_matrix, Performance, Metrix

