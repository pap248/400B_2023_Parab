#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: paarth

This program will use the data from MW_000.txt to give the needed particle
values such as the magnitude of the distance in both kpc and light years. 
The magnitude of the velocity in km/s and the mass of the particle in 
terms of M_sun. 
"""

import numpy as np
import astropy.units as u
from Readfile import *
from astropy import units as u
from Readfile import Read
# importing needed packages

'''
filename = input('give filename:')
parttype = input('give particle type:')
partnum = int(input('give particle number:'))
data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
'''

# the above is test code, I originally thought to give inputs but wasn't sure


def Particleinfo(filename,parttype,partnum):      # making function 
    # parameters for function are the filename, particle type, and particle number
    file = Read(filename)      # using the Read function from Readfile.py to take in the data
    data = file[2]     # grabbing the data array from Readfile.py
    
    index = np.where(data['type'] == parttype )       # this is the indexing for the particle data
    # this is to sepearte the 3 different types from each other, 1 = Dark matter
    # 2 = Disk Stars and 3 = Bulge stars
    
    # for the below code, I use partnum-1 to make the sure that the particles are given the right number
    # so the first particle is actually particle number 1 instead of being particle 0
    # this also helps for organizing with the other types and their particles
    
    x = data['x'][index][partnum-1]       # this is grabbing the x value based on the particle number
    vx = data['vx'][index][partnum-1]     # this is grabbing velocity in x direction of the particle
    x = x*u.kpc                           # making units kpc
    vx = vx*(u.km/u.s)                    # making units km/s
        
    y = data['y'][index][partnum-1]       # this is grabbing the y value based on the particle number
    vy = data['vy'][index][partnum-1]     # this is grabbing velocity in y direction of the particle
    y = y*u.kpc                           # making units kpc
    vy = vy*(u.km/u.s)                    # making units km/s
        
    z = data['z'][index][partnum-1]       # this is grabbing the z value based on the particle number
    vz = data['vz'][index][partnum-1]     # this is grabbing velocity in z direction of the particle
    z = z*u.kpc                           # making units kpc
    vz = vz*(u.km/u.s)                    # making units km/s
    
    magdis = (np.sqrt((x**2)+(y**2)+(z**2)))       # finding magnitude for distance
    magdis = np.around(magdis,3)                   # rounding distance to 3 digits
    magvel = (np.sqrt((vx**2)+(vy**2)+(vz**2)))    # finding magnitude for velocity
    magvel = np.around(magvel,3)                   # rounding velocity to 3 digits
    
    mass = data['m'][index][partnum-1]             # taking the mass of the particle
    mass = mass / (10**10)                        # this is putting the mass back to terms of M_sun
    
    return magdis, magvel, mass                  # return all the values now for the particle
    
    
func = Particleinfo("MW_000.txt", 2, 100)          # now setting variable for function and taking the
                                                 # proper parameters, and calling the function

magdislyr = func[0].to(u.lyr)                     # now changing distance from kpc to light years

print("Distance is ",func[0], ",velocity is ",func[1], ",mass is ",func[2], \
      ",distance in lyr", magdislyr) # printing the distance of the particle in both kpc and
    # light years, the velocity in km/s and the mass of the particle