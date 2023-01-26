#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: paarth

This program reads the file MW_000.txt and returns
the snapcshot time, the total number of particles,
and the data of the whole file


"""

import numpy as np
import astropy.units as u
# importing needed packages above

def Read(filename): # making function with filename parameter
    file = open("MW_000.txt", "r")       # open the file using read mode
    line1 = file.readline()              # read the first line
    label, value = line1.split()         # splits values 
    time = float(value)*u.Myr            # time is the value saved as float
    
    line2 = file.readline()              # reads next line
    label, value = line2.split()         # splits
    totpart = float(value)*u.Myr         # total number of particles
    
    file.close()                    # closing file
    
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    # above is all the data for the file, skipping the headers, 
    # and making the right headers
    
    
    #print(data['x'][2])       # this is a test
    #print(time)               # so is this
    
    return time, totpart, data      # now returning the time, total particles
# and the data

#print(Read("MW_000.txt"))  # also a test

