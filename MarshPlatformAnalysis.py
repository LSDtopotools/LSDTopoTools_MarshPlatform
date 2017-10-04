"""
MarshPlatformAnalysis.py

This function drives marsh platform analysis

Authors: Guillaume CH Goodwin, Simon M. Mudd and Fiona J. Clubb, University of Edinburgh


"""

# First import the mecessary modules
import os
import sys
import LSDMarshPlatform as MP


#=============================================================================
# This is just a welcome screen that is displayed if no arguments are provided.
#=============================================================================
def print_welcome():

    print("\n\n=======================================================================")
    print("Hello! I'm going to run a marsh platform analysis.")
    print("You will need to tell me which directory to look in.")
    print("Use the -dir flag to define the working directory.")
    print("If you don't do this I will assume the data is in the same directory as this script.")
    print("For help type:")
    print("   python MarshPlatformAnalysis.py -h\n")
    print("=======================================================================\n\n ")

#=============================================================================
# This is the main function that runs the whole thing
#=============================================================================
def main(argv):

    # print("On some windows systems you need to set an environment variable GDAL_DATA")
    # print("If the code crashes here it means the environment variable is not set")
    # print("Let me check gdal enviroment for you. Currently is is:")
    # print(os.environ['GDAL_DATA'])
    #os.environ['GDAL_DATA'] = os.popen('gdal-config --datadir').read().rstrip()
    #print("Now I am going to get the updated version:")
    #print(os.environ['GDAL_DATA'])

    # If there are no arguments, send to the welcome screen
    if not len(sys.argv) > 1:
        full_paramfile = print_welcome()
        sys.exit()

    # Get the arguments
    import argparse
    parser = argparse.ArgumentParser()
    # The location of the data files
    parser.add_argument("-dir", "--base_directory", type=str, help="The base directory with the DEMs for the marsh analysis. If this isn't defined I'll assume it's the same as the current directory.")
    parser.add_argument("-sites", "--sites",type=str,default = "", help = "This is a comma delimited string that gets the list of sites you want for analysis and plotting. This is a prefix that preceeds all the other DEM extensions. Default = no sites")

    # What sort of analyses you want
    parser.add_argument("-MID", "--MarshID", type=bool, default=False, help="If this is true, this will run the marsh ID algorithm")

    # What sort of plots you want
    parser.add_argument("-MIDP", "--MarshID_plots", type=bool, default=False, help="If this is true I'll plot all the platform plots.")

    args = parser.parse_args()

    sites = []
    if not args.sites:
        print("WARNING! You haven't supplied your site names. Please specify this with the flag '-sites'")
        sys.exit()
    else:
        print("The sites you want to analyse are: ")
        sites = [str(item) for item in args.sites.split(',')]
        print(sites)

    # get the base directory
    if args.base_directory:

        this_dir = args.base_directory
        print("You gave me the base directory:")
        print(this_dir)
    else:
        this_dir = os.getcwd()
        print("You didn't give me a directory. I am using the current working directory:")
        print(this_dir)

    # Run the analysis if you want it
    if args.MarshID:
        MP.MarshID(Input_dir = this_dir, Output_dir = this_dir,Sites=sites)

    # make the plots depending on your choices
    if args.MarshID_plots:
        MP.Plot_platform_on_hillshade(Input_dir = this_dir, Output_dir = this_dir,Sites=sites)
        MP.Plot_marsh_outline_on_hillshade(Input_dir = this_dir, Output_dir = this_dir,Sites=sites)
        MP.Plot_Elevation_PDF(Input_dir = this_dir, Output_dir = this_dir,Sites=sites)


#=============================================================================
if __name__ == "__main__":
    main(sys.argv[1:])
