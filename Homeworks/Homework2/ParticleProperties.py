#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: paarth


"""

import numpy as np
import astropy.units as u
from Readfile import *
from astropy import units as u
from Readfile import Read

'''
filename = input('give filename:')
parttype = input('give particle type:')
partnum = int(input('give particle number:'))
'''
#data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)


def Particleinfo(filename,parttype,partnum):
    
    file = Read(filename)
    data = file[2]     # grabbing the data array from Readfile
    
    index = np.where(data['type'] == parttype )
    
    
    x = data['x'][index][partnum-1]
    vx = data['vx'][index][partnum-1]
    x = x*u.kpc
    vx = vx*(u.km/u.s)
        
    y = data['y'][index][partnum-1]
    vy = data['vy'][index][partnum-1]
    y = y*u.kpc
    vy = vy*(u.km/u.s)
        
    z = data['z'][index][partnum-1]
    vz = data['vz'][index][partnum-1]
    z = z*u.kpc
    vz = vz*(u.km/u.s)
    
    magdis = (np.sqrt((x**2)+(y**2)+(z**2)))
    magdis = np.around(magdis,3)
    magvel = (np.sqrt((vx**2)+(vy**2)+(vz**2)))
    magvel = np.around(magvel,3)
    
    mass = data['m'][index][partnum-1]
    mass = mass / (10**10)
    
    return magdis, magvel, mass
    
    
func = Particleinfo("MW_000.txt", 2, 100)

magdislyr = func[0].to(u.lyr)

print("Distance is ",func[0], ",velocity is ",func[1], ",mass is ",func[2], \
      ",distance in lyr", magdislyr)