"""
Motherscript.py

This Python script runs all the scripts necessary to:
(0. Prepare your data. Please read the README file before you enable this option.)
1. Extract marsh platforms and outlines from an input DEM
2. Save these features as .bil files
3. Plot the figures of our paper, including a comparison with a reference DEM
"""

# First import the mecessary modules
import os
import sys

# Then run the scripts
#os.system("python Example_Input_prep.py")
os.system("python Example_Marsh_ID.py")
os.system("python Example_Plots.py")