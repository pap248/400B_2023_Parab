#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: paarth

The fate of Sun analogs in Andromeda Galaxy
This program will find the fate of the sun analogs that are in M31, 8kpc
from galactic center. 
"""

import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from Readfile import Read

# modified CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease rmax instead of a factor of 2
from CenterOfMass_Solution import CenterOfMass
from io import StringIO


''' First, I plan on reading in the first snap shot of M31 and then looking at
how many sun analogs are there in total, seeing how many of them are 
about 8 kpc from the center of the galaxy'''

''' will look at every 50 - 100 snapshots. looking through the file, 
finding the solar particles that are distanced 7 - 9 kpc away from center of mass of M31.'''

''' my first step is to read in the first snapshot ,time 0, of the M31 galaxy 
and make a loop to see how many of the solar particles are similar to our sun's distance
from center of mass of galaxy'''

def __init__(self, filename, ptype):    
    ''' do type 2'''
    ''' Class to calculate the 6-D phase-space position of a galaxy's center of mass using
    a specified particle type. 
        
        PARAMETERS
        ----------
        filename : `str`
            snapshot file
        ptype : `int; 1, 2, or 3`
            particle type to use for COM calculations
    '''
 
    # read data in the given file using Read
    self.time, self.total, self.data = Read(filename)                                                                                             

    #create an array to store indexes of particles of desired Ptype                                
    self.index = np.where(self.data['type'] == ptype)

    
    # mass not needed but keeping here in case
    self.m = self.data['m'][self.index]
    
    #positions for the suns
    self.x = self.data['x'][self.index]
    self.y = self.data['y'][self.index]
    self.z = self.data['z'][self.index]
    
    # velocities not needed but keeping here in case
    self.vx = self.data['vx'][self.index]
    self.vy = self.data['vy'][self.index]
    self.vz = self.data['vz'][self.index]


''' now that the positions have been read in, my next step is to find the magnitude, 
the radial distance '''

def SunAnalogs(self, delta):
    '''Method to compute the position of the stars from the center of galaxy 

    PARAMETERS
    ----------
    delta : `float, optional`
        error tolerance in kpc. Default is 0.1 kpc
    
    RETURNS
    ----------
    position : `np.ndarray of astropy.Quantity'
        3-D position in kpc
    '''                                                                     
                                                                                                                       
    x_pos, y_pos, z_pos = self.COMdefine(self.x, self.y, self.z, self.m)
    
    # compute the magnitude of the position vector
    positions = np.sqrt(x_pos**2 + y_pos**2 + z_pos**2)
    
    ''' now that the radial positions have been found , I will make a loop to
    run through the array of positions and add the ones that are about 8 kpc 
    to a list '''

    
    sun_analogs = []
    for i in positions:
        if i > 7 and i < 9:
            np.append(sun_analogs, i)

    
    return sun_analogs


SunAnalogs(0.1)


''' in a previous homework we found where the merger takes place approximately, so 
I will make sure to use a snapshot from when the merger finishes also'''

''' for the plot, I want to show the radial postions as a function of time for M31.
x axis would be time, varying every 100 snapshots or so, then y is the radial postions
'''

''' for another plot I want to show a face on view of the sun analogs at the beginning
vs the end of merger, specifically for those within 10 kpc
not sure if that is too complex'''


''' Mainly I want to look at the stars that are at snapshot 0 then the final snapshot
and compare how many of them are still within 7 - 9 kpc from center'''











