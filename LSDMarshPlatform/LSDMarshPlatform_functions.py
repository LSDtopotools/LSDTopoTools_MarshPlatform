


#----------------------------------------------------------------
#1. Load useful Python packages
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal, osr, gdalconst
from osgeo.gdalconst import *
import cPickle



#---------------------------------------------------------------
# This function generates data distributions out of a 2D raster
def Distribution(Data2D, Nodata_value):
    

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





#-----------------------------------------------------------------------------------------------------
# This functions calculates local slope using the maximum slope method
# http://www.onlinegeographer.com/slope/Dunn_hickey_slope.pdf, equation 7

# It takes as input:
# 1/ The DEM array

# It returns:
# 1/ A Slope array


def maximum_slope (DEM):
    Height = len(DEM_work); Width = len(DEM_work[0,:])
    Slope_max = np.zeros((Height,Width), dtype=np.float)

    for i in range(1,len(DEM_work)-1): # rows
        for j in range(1,len(DEM_work[0,:])-1):  # cols
            # Here are the cells we consider
            Cells = np.zeros(9, dtype=np.float)
            Slopes = np.zeros(9, dtype=np.float)
            Cells[0]=DEM_work[i,j]; Cells[1]=DEM_work[i-1,j]; Cells[2]=DEM_work[i-1,j+1]; Cells[3]=DEM_work[i,j+1]; Cells[4]=DEM_work[i+1,j+1]; Cells[5]=DEM_work[i+1,j]; Cells[6]=DEM_work[i+1,j-1]; Cells[7]=DEM_work[i,j-1]; Cells[8]=DEM_work[i-1,j-1]

            # Calculate Slopes
            for a in range(1,len(Cells)):
                Slopes[a]= (Cells[a]-Cells[0])/1

            # Largest Slope is the slope value for i,j
            Slope_max[i,j]=max(Slopes)


    return Slope_max





#-----------------------------------------------------------------------------------------------------
# This functions defines a search space within an input raster

# It takes as input:
# 1/ The data array
# 2/ The slope array

# It returns:
# 1/ A search space within the data array

def define_search_space (DEM, Slope, Nodata_value,opt):
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
        
    
    
    
  
    #fig=plt.figure(3, facecolor='White',figsize=[4.7,4])
    #ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1, axisbg='white')
    #ax1.plot(bins, hist)
    #plt.savefig('Output/Paper/0_Fig3.png')
    #STOP

    return Search_space, Crossover, bins, hist, Inflexion_point


#-----------------------------------------------------------------------------------------------------
# This functions makes a kernel in an input raster

# It takes as input:
# 1/ The data array
# 2/ The desired kernel size (the total length or height of it!)
# 3/ The coordinates of the kernel centre

# It returns:
# 1/ A kernel, which is a lower ranking officer than Hollywood would have you believe

def kernel (array, kernel_size, x_centre, y_centre):
    if (-1)**kernel_size < 0:
        X_to_0 = x_centre
        X_to_End = len(array)-x_centre
        Y_to_0 = y_centre
        Y_to_End = len(array[0,:])-y_centre

        Lim_left = x_centre - min(np.floor(kernel_size/2), X_to_0)
        Lim_right = x_centre + min(np.floor(kernel_size/2)+1, X_to_End)
        Lim_top = y_centre - min(np.floor(kernel_size/2), Y_to_0)
        Lim_bottom = y_centre + min(np.floor(kernel_size/2)+1, Y_to_End)

        kernel = array [Lim_left:Lim_right, Lim_top:Lim_bottom]

    else:
        print
        print " ... WARNING: you need to choose an odd kernel size, buddy"
        print
        pass

    return kernel


#-----------------------------------------------------------------------------------------------------
# This functions detects and flags local peaks of a data array

# It takes as input:
# 1/ The data array
# 2/ The desired search space (0 = don't search; 1 = search)

# It returns:
# 1/ An array where local peaks have a value of 1 (all other values are 0)
# 2/ A copy of the data array where the locations of the peaks are 0.

def peak_flag (Slope, Search_space, Order):
    print 'Preparing the Kilimanjaro expedition ...'
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
# This functions detects the 2 or more points from which ridges are initiated. We call them ridge starters

# It takes as input:
# 1/ The data array
# 2/ The desired search space (0 = don't search; 1 = search)
# 3/ The array containing the location of local peaks

# It returns:
# 1/ An array where local peaks are 1 and ridge starters are 2 (all other values are 0)
# 2/ A copy of the data array where the locations of the peaks and ridge starters are 0.


def initiate_ridge (Slope, Search_space, Peaks, Order):
    print ' ... Rolling off Mt. Kilimanjaro ...'
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
# This functions continues the ridges from the ridge starters. It incorporates tidal range to filter out noise

# It takes as input:
# 1/ The data array
# 2/ The desired search space (0 = don't search; 1 = search)
# 3/ The array containing the location of local peaks
# 4/ The Tidal ranges metric
# 5/ The ridge Order

# It returns:
# 1/ An array where local peaks are 1 and ridge starters are 2 (all other values are 0)
# 2/ A copy of the data array where the locations of the peaks and ridge starters are 0.


def Continue_ridge (DEM, Slope, Search_space, Peaks, Order):
    #print ' ... Rolling down ...'

    DEM_copy = np.copy(DEM) # the copy of the initial DEM array
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
# This functions cleans up the ridges

# It takes as input:
# 1/ The ridge array
# 2/ The elevation array
# 3/ The slope array
# $/ the tidal range

# It returns:
# 1/ A ridge array cleansed of unnecessary ridges


def Clean_ridges (Peaks, DEM, Slope, Nodata_value,opt):
    print "Cleaning up: I want to break free ..."
    DEM_copy = np.copy(DEM)
    DEM_copy[DEM_copy==Nodata_value] = 0
    Search_ridge = np.where (Peaks != 0)

    print " ... What a relief ..."

    Cutoff = np.percentile(DEM_copy,75)
    Threshold = np.amax(DEM_copy[DEM_copy<Cutoff])
    DEM_copy[DEM_copy>Threshold]=Threshold

    for i in range(len(Search_ridge[0])):
        x=Search_ridge[0][i]; y=Search_ridge[1][i] # coordinates of the kernel's centre
        Kernel_DEM = kernel (DEM_copy, 9, x, y)
        Kernel_DEM[Kernel_DEM==Nodata_value]=0
        Kernel_relief = Kernel_DEM - np.amin(Kernel_DEM)
        Kernel_slope = kernel(Slope, 9, x, y)

        if np.amax(Kernel_DEM)/Threshold < opt:
            Peaks[x,y] = 0



    print " ... Shave the stubble ..."
    Search_ridge = np.where (Peaks != 0)
    for i in range(len(Search_ridge[0])):
        x=Search_ridge[0][i]; y=Search_ridge[1][i] # coordinates of the kernel's centre
        Kernel_ridges = kernel (Peaks, 9, x, y)
        # If there aren't at least 8 ridge points in the neighbourhood of 10 by 10
        if np.count_nonzero(Kernel_ridges) < 8:
            Peaks[x,y] = 0


    return Peaks






#-----------------------------------------------------------------------------------------------------
# This functions fills areas above the steep bits up the raster

# It takes as input:
# 1/ The ridge array
# 2/ The DEM
# 3/ Tidal properties

# It returns:
# 1/ An array with the marsh bits

def Fill_high_ground (DEM, Peaks, Nodata_value,opt):
    print "Paint me a platform ..."
    DEM_copy = np.copy(DEM)
    Marsh = np.zeros((len(DEM), len(DEM[0,:])), dtype = np.float)

    print " ... Start close to your sketch lines ..."
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

    print " ... Erase when you've overstepped the line ..."
    Search_marsh_start = np.where (Marsh == 1)
    for i in range(len(Search_marsh_start[0])):
        x=Search_marsh_start[0][i]; y=Search_marsh_start[1][i]
        Kernel_marsh = kernel (Marsh, 3, x, y)
        Kernel_ridges = kernel (Peaks, 3, x, y)
        if np.count_nonzero(Kernel_marsh) <=2:
            Marsh[x,y] = 0


    while Counter < 100:
        Counter = Counter+1
        #print ' ... Filling ... ', Counter
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

            """These conditions work! They generate a distribution where you can claerly see if there is some colouring on the tidal flat because you will see a peak at the lowest elevations. All you need to do now is identify that peak and cut it off!"""
            

    
    #This is where you define the cutoff spot!
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

   
    
    print " ... Fill high gaps ..."
    Search_marsh_condition = np.zeros((len(DEM), len(DEM[0,:])), dtype = np.float)
    Search_marsh = np.where (DEM >= Platform_bins[Index])
    Search_marsh_condition [Search_marsh] = 1
    Search_marsh_2 = np.where (np.logical_and(Marsh == 0, Search_marsh_condition == 1))
    Marsh[Search_marsh_2] = 3



        
    for Iteration in np.arange(0,10,1):
        Counter = 100
        while Counter > 2:
            Counter = Counter-1
            #print " ... Reverse filling ... ", Counter
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
    print " ... Filling ISOs ... "
    Search_marsh = np.where (np.logical_and(Marsh == 0, Peaks == 0))
    for i in range(len(Search_marsh[0])):
        x = Search_marsh[0][i]; y = Search_marsh[1][i]
        Kernel_marsh = kernel (Marsh, 3, x, y)
        if np.count_nonzero(Kernel_marsh) == 8:
            Marsh[x,y] = 105



    # We get rid of scarps that do not have a marsh next to them
    print " ... Eliminating false scarps ..."
    Search_false_scarp = np.where (Peaks > 0)
    for i in range(len(Search_false_scarp[0])):
        x = Search_false_scarp[0][i]; y = Search_false_scarp[1][i]
        Kernel_marsh = kernel (Marsh, 3, x, y)
        if np.count_nonzero (Kernel_marsh) == 0:
            Peaks[x, y] = 0


    # We get rid of the sticky-outy bits
    print " ... Shave the stubble ..."
    Search_ridge = np.where (Peaks > 0)
    for i in range(len(Search_ridge[0])):
        x=Search_ridge[0][i]; y=Search_ridge[1][i]
        Kernel_ridges = kernel (Peaks, 9, x, y)
        if np.count_nonzero(Kernel_ridges) < 8:
            Peaks[x,y] = 0


    
    
    # We put the scarps in the platform
    print " ... Filling ridges ..."
    Search_side = np.where (Peaks > 0)
    Marsh[Search_side] = 110
    
    
    
    
    
    
    # Some of our platforms are patchy. Try filling them now that we have added the scarps
    
    
    print " ... Fill high gaps ..."
    Search_marsh_condition = np.zeros((len(DEM), len(DEM[0,:])), dtype = np.float)
    Search_marsh = np.where (DEM >= Platform_bins[Index])
    Search_marsh_condition [Search_marsh] = 1
    Search_marsh_2 = np.where (np.logical_and(Marsh == 0, Search_marsh_condition == 1))
    Marsh[Search_marsh_2] = 3
    
    for Iteration in np.arange(0,10,1):
        Counter = 110
        while Counter > 2:
            Counter = Counter-1
            #print " ... Reverse filling ... ", Counter
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
    

    
    
    
    Marsh[DEM == Nodata_value] = Nodata_value

    return Marsh

 



#---------------------------------------------------------------
# This function is the MarshFinder: it finds marsh platforms, scarps and pioneer zones
# Here's how the MarshFinder works:
# STEP 0: Identify a sesarch space where elevation is above the Neap Low Tide and slope is in the higher 50%
# STEP 1: Identify scarps by:
#                          identifying local slope maxima (peaks)
#                          starting ridge lines along the highest slope values stemming from these peaks
#                          continuing the ridges along the shortest path of high slopes
# STEP 2: Clean the ridges if:
#                           they are too short
#                           the local relief is too small
#                           they are lower than the most frequent elevation of ridges minus a constant depending on tidal range
# STEP 3: Fill ground above ridges by:
#                                   filling ground directly in contact with ridges and above the ridge
#
# STEP 4:
# STEP 5:
# STEP 6:
# STEP 7:
# STEP 8:
# It takes as input:
# 1/ the DEM
# 2/ the Slopes
# 3/ the Curvature
# 4/ the Channels
# 5/ the Tidalstatistix
# It returns:
# 1/ the marsh platforms
# 2/ the marsh scarps
# 3/ the marsh channels
# 4/ the pioneer zones



#def MARSH_ID (DEM, Slope, Curvature, Metric2, Nodata_value):
#Added this line for optimisation purposes
def MARSH_ID (DEM, Slope, Curvature, Nodata_value, opt1, opt2, opt3):
    DEM_work = np.copy(DEM); Slope_work = np.copy(Slope); Curvature_work = np.copy(Curvature) #; Channels_work = np.copy(Channels)

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
        Ridge, Slope_temp = Continue_ridge (DEM, Slope_temp, Search_space, Ridge, Order)

    Ridge = Clean_ridges (Ridge, DEM_work, Slope_work, Nodata_value, opt2)

    Marsh = Fill_high_ground (DEM_work, Ridge, Nodata_value, opt3)


    print "My hovercraft is full of eels!"
    print


    return Search_space, Ridge, Marsh




#-----------------------------------------------------------------------------------------------------
# This functions compares the marsh identified automatically to a reference marsh, usually digitised from a picture
# It then builds a confusion matrix to see if the MarshFinder has performed as expected.

# It takes as input:
# 1/ The subject marsh array
# 2/ The reference marsh array


# It returns:
# 1/ An array of the confusion matrix and its associated metrix (in a raster and a figure)

def Confusion (Subject, Reference, Nodata_value):
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






#---------------------------------------------------------------
def ENVI_raster_binary_to_2d_array(file_name, gauge):
    print 'Opening %s' % (gauge)
    #Source : http://chris35wills.github.io/python-gdal-raster-io/
    '''
    Converts a binary file of ENVI type to a numpy array.
    Lack of an ENVI .hdr file will cause this to crash.
    '''

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
    #Source : http://chris35wills.github.io/python-gdal-raster-io/
    #util.check_output_dir(file_out)

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


#-------------------------------------------------------------------
#7.
# Some subaxes
#http://stackoverflow.com/questions/17458580/embedding-small-plots-inside-subplots-in-matplotlib
def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax
