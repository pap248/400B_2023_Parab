#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: paarth

This program will take the total masses of each galaxy given their type,
Halo, Disk, or Bulge, in units of 10^12 Msun.
A table is produced showing the masses of each of those types for each galaxy,
and also will show the total mass of the galaxy as well as the baryon fraction.
Also gives all the values for the entire local group, MW, M31 , and M33.
"""

from Readfile import *
from Readfile import Read
import numpy as np
import astropy.units as u
from astropy import units as u
from tabulate import tabulate
# the above is importing all the packages needed for the program

# the function below will take in the filename and the particle type (1,2, or 3)
# 1 is Halo, 2 is Disk , 3 is Bulge
def ComponentMass(filename, parttype):
    
    file = Read(filename)      # using the Read function from Readfile.py to take in the data
    data = file[2]             # grabbing data array
    index = np.where(data['type'] == parttype)      # indexing for particle data
    totmass = 0                            # this sets the total mass to 0 and 
                                           # and will be adding to this value in the loop
    for partnum in index:                  # iterating through each mass and adding to totmass
        allmass = data['m'][partnum-1]
        for i in allmass:
            totmass += i
    totmass = np.round(totmass*(1e-2),3)     # rounding each total mass to 3 decimal places
    return totmass                           # now returning the total mass based on filename and type





MWHalo = ComponentMass('MW_000.txt', 1)
MWDisk = ComponentMass('MW_000.txt', 2)
MWBulge = ComponentMass('MW_000.txt', 3)
MWtotal = MWHalo + MWDisk + MWBulge
# the above 4 lines is the different masses for Milky Way
# Halo, Disk, Bulge, and total mass of galaxy


M31Halo = ComponentMass('M31_000.txt', 1)
M31Disk = ComponentMass('M31_000.txt', 2)
M31Bulge = ComponentMass('M31_000.txt', 3)
M31total = M31Halo + M31Disk + M31Bulge
# the above 4 lines is the different masses for Andromeda Galaxy
# Halo, Disk, Bulge, and total mass of galaxy


M33Halo = ComponentMass('M33_000.txt', 1)
M33Disk = ComponentMass('M33_000.txt', 2)
M33Bulge = ComponentMass('M33_000.txt', 3)
M33total = M33Halo + M33Disk + M33Bulge
# the above 4 lines is the different masses for M33
# Halo, Disk, Bulge, and total mass of galaxy


LocalHalo = MWHalo + M31Halo + M33Halo
LocalDisk = MWDisk + M31Disk + M33Disk
LocalBulge = MWBulge + M31Bulge + M33Bulge
Localtotal = MWtotal + M31total + M33total
# the above 4 lines is the different masses for the Local Group
# Halo, Disk, Bulge, and total masses of Local Group

MWstellar = MWDisk+MWBulge
M31stellar = M31Disk+M31Bulge
M33stellar = M33Disk+M33Bulge
Localstellar = LocalDisk+LocalBulge
# the above 4 lines is the stellar values for the galaxies and the Local Group
# will be used later


fbarMW = ((MWstellar) / (MWtotal))
fbarM31 = ((M31stellar) / (M31total))
fbarM33 = ((M33stellar) / (M33total))
fbarLocal = ((Localstellar) / (Localtotal))
# the above 4 lines is Baryon fractions for the galaxies and the Local Group
# using Disk + Bulge masses and dividing by total mass




# the below is making the table with the values and column names
# using tabulate (method found online)

table = [['Milky Way',MWHalo,MWDisk,MWBulge,MWtotal,fbarMW], ['M31',M31Halo, M31Disk, M31Bulge,\
                                                              M31total, fbarM31],
          ['M33',M33Halo,M33Disk,M33Bulge,M33total, fbarM33], ['Local Group', LocalHalo, LocalDisk\
                                                      ,LocalBulge, Localtotal, fbarLocal ]]


print(tabulate(table, headers = ['Galaxy name\nUnits', 'Halo Mass\n10^12Msun', 'Disk Mass\n10^12Msun'\
                                 , 'Bulge Mass\n10^12Msun',\
                                 'Total mass\n10^12Msun' , 'fbar'] ))

# the above is printing the table with the appropriate headers
# sorry if the units row looks weird I couldn't make the 12 in 10^12 go to exponent or 
# make the sun for Msun go to subscript

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
   