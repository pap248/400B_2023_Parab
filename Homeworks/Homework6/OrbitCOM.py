#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

@author: Paarth

This program will store the separation and relative velocities
of the center of mass of the simulated Milky Way, Andromeda, and M33
over the entire simulation and plot the orbits.

"""
'''

Worked with Rey

'''


# import modules
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



def OrbitCOM(galaxy, start, end, n):
    """function that loops over all the desired snapshots to compute the COM 
    pos and vel as a function of time.
    
    inputs:
          galaxy = name of galaxy
          start = number of first snapshot to be read in
          end = number of the last snapshot to be read in
          n = integer for intervals over returning COM
    outputs: 
        files for Orbit of galaxy
    """
    
    # filenames for outputs
    fileout = 'orbit_{}.txt'.format(galaxy)
    
    # making if statement for the value of volDec 
    # M33 will be 4 for volDec due to it being tidally stripped severly
    if galaxy == 'M33':
        volDec = 4
    else:
        volDec = 2
    
    #  set tolerance and volDec for calculating COM_P in CenterOfMass
    delta = 0.1
    
    # testing
    #volDec = 2
    
    
    # generate the snapshot id sequence 
    snap_ids = np.arange(start,end+1, n)
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros([len(snap_ids),7])
    
    # a for loop 
    # loop over files
    
    # testing
    #print(orbit)
    for i in range(len(snap_ids)):
        
        # add a string of the filenumber to the value '000'
        ilbl = '000' + str(snap_ids[i])
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        filename='%s_'%(galaxy) + ilbl + '.txt'
        
        # Initialize an instance of CenterOfMass class, using disk particles
        COM = CenterOfMass(filename, 2)
        
        # testing
        #print(COM.time.to_value())
        #print(orbit[i])
        #print((COM.time).value)
        
        # changing time to store to Gyr
        orbit[i,0] = (COM.time).value / 1000
        
        # Store the COM pos and vel
        COM_P = (COM.COM_P(delta,volDec)).value
        orbit[i,1],orbit[i,2], orbit[i,3] = COM_P
        
        # to make easier to use
        x, y ,z = COM.COM_P(delta,volDec)
        COM_V = (COM.COM_V(x,y,z)).value
        
        orbit[i,4], orbit[i,5], orbit[i,6]= COM_V
        
        # testing
        #print(orbit[i])
        #print(orbit)
        
        # store the time, pos, vel in ith element of the orbit array,  
        # without units (.value) 

        
        # print snap_id to see the progress
        print(snap_ids[i])
    
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
# making of files 
#(OrbitCOM('MW', 0, 800, 5))
#(OrbitCOM('M31', 0, 800, 5))
#(OrbitCOM('M33', 0, 800, 5))

# reading in the files just made
# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt

fileMW = open('orbit_MW.txt','r')
MWlines = fileMW.readlines()
print(len(MWlines))
MWa = np.genfromtxt('orbit_MW.txt', delimiter="", skip_header = 1)


fileM31 = open('orbit_M31.txt','r')
M31lines = fileM31.readlines()
print(len(M31lines))
M31a = np.genfromtxt('orbit_M31.txt', delimiter="", skip_header = 1)


fileM33 = open('orbit_M33.txt','r')
M33lines = fileM33.readlines()
print(len(M33lines))
M33a = np.genfromtxt('orbit_M33.txt', delimiter="", skip_header = 1)

# making function to compute difference between vectors
# for both position and velocity 
# inputs are the vectors and output is the magnitude

def vectorsdiff(vec1, vec2):
    
    diff = vec1 - vec2
    mag = np.sqrt(diff[0]**2 + diff[1]**2 + diff[2]**2)
    return mag

# looping through the values to get the vectors and adding to array from
# text file
# first is MW and M31 then M33 and M31
MW_M31_rel = np.zeros([len(MWa), 3])
for i in range(len(MW_M31_rel)):
    MW_M31_rel[i, 0] = MWa[i][0]
    MW_M31_rel[i, 1] = vectorsdiff(np.array(MWa.tolist()[i][1:4]), \
                                   np.array(M31a.tolist()[i][1:4]))
    MW_M31_rel[i, 2] = vectorsdiff(np.array(MWa.tolist()[i][4:7]), \
                                   np.array(M31a.tolist()[i][4:7]))
    #print(MW_M31_rel)

M31_M33_rel = np.zeros([len(M31a), 3])
for i in range(len(MW_M31_rel)):
    M31_M33_rel[i, 0] = MWa[i][0]
    M31_M33_rel[i, 1] = vectorsdiff(np.array(M31a.tolist()[i][1:4]), \
                                   np.array(M33a.tolist()[i][1:4]))
    M31_M33_rel[i, 2] = vectorsdiff(np.array(M31a.tolist()[i][4:7]), \
                                   np.array(M33a.tolist()[i][4:7]))

# now making the plots for the separations and velocities
plt.figure()
plt.plot(MW_M31_rel[:,0], MW_M31_rel[:,1])
plt.xlabel('Time (Gyr)')
plt.ylabel('Distance (kpc)')
plt.title('Milky Way and Andromeda Separation')


plt.figure()
plt.plot(M31_M33_rel[:,0], M31_M33_rel[:,1])
plt.xlabel('Time (Gyr)')
plt.ylabel('Distance (kpc)')
plt.title('Andromeda and M33 Separation')

# now velocities
plt.figure()
plt.plot(MW_M31_rel[:,0], MW_M31_rel[:,2])
plt.xlabel('Time (Gyr)')
plt.ylabel('Velocity (km / s)')
plt.title('Milky Way and Andromeda Velo Separation')

plt.figure()
plt.plot(M31_M33_rel[:,0], M31_M33_rel[:,2])
plt.xlabel('Time (Gyr)')
plt.ylabel('Velocity (km / s)')
plt.title('Andromeda and M33 Velo Separation')







# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative 
# velocity for two galaxies over the entire orbit  

#def vectors(pos,vel):
    


# Determine the magnitude of the relative position and velocities 

# of MW and M31

# of M33 and M31




# Plot the Orbit of the galaxies 
#################################




# Plot the orbital velocities of the galaxies 
#################################

