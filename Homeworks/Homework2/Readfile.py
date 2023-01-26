#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: paarth

"""
import numpy as np
import astropy.units as u


def Read(filename):
    file = open("MW_000.txt", "r")
    line1 = file.readline()
    label, value = line1.split()
    time = float(value)*u.Myr
    
    line2 = file.readline()
    label, value = line2.split()
    totpart = float(value)*u.Myr
    
    file.close()
    
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    
    #print(data['x'][2])
    #print(time)
    
    return time, totpart, data

#print(Read("MW_000.txt"))

