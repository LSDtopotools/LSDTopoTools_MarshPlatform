"""
Example_Motherscript.py

This function calls both platform extraction and plotting routines for the Marsh Platform finder. 

Authors: Guillaume CH Goodwin and Simon Mudd, University of Edinburgh


This Python script runs all the scripts necessary to:
(0. Prepare your data. Please read the README file before you enable this option.)
1. Extract marsh platforms and outlines from an input DEM
2. Save these features as .bil files
3. Plot the figures of our paper, including a comparison with a reference DEM
"""

# First import the mecessary modules
import os
import sys
import LSDMarshPlatform as MP

# Then run the scripts
print("First I will identify the marsh platform for you")
input_dir = "T:\\Git_projects\\LSDTopoTools_MarshPlatform\\Example_data\\"
output_dir = "T:\\Git_projects\\LSDTopoTools_MarshPlatform\\Example_data\\"
MP.Example_Marsh_ID(Input_dir = input_dir, Output_dir = output_dir)

print("Now I will generate some plots for you")
MP.Example_Plots