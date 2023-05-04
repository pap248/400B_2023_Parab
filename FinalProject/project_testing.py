#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: paarth

This program will take the number of solar particles from each snap, and 
see how many of them are within 7.4 to 8.7 kpc from the center of the M31.
Changing for each snapshot to get the amount for each snapshot, not following
the same solar particles, but getting new ones each time that are within that distance.

"""
# importing modules

import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


from Readfile import Read
#from CenterOfMass_Solution import CenterOfMass
from io import StringIO


'''taking some code from class'''
class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot

    def __init__(self, filename, ptype):    
        ''' doing type 2 for disk'''
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
        
        # velocities 
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]

    ''' now that the positions have been read in, my next step is to find the magnitude, 
    the radial distance '''

    # fix to account for COM

    def COMdefine(self,a,b,c,m):
        ''' Method to compute the COM of a generic vector quantity by direct weighted averaging.
        
        PARAMETERS
        ----------
        a : `float or np.ndarray of floats`
            first vector component
        b : `float or np.ndarray of floats`
            second vector component
        c : `float or np.ndarray of floats`
            third vector component
        m : `float or np.ndarray of floats`
            particle masses
        
        RETURNS
        -------
        a_com : `float`
            first component on the COM vector
        b_com : `float`
            second component on the COM vector
        c_com : `float`
            third component on the COM vector
        '''
        
        
        # xcomponent Center of mass
        a_com = np.sum(a*m)/np.sum(m)
        # ycomponent Center of mass
        b_com = np.sum(b*m)/np.sum(m)
        # zcomponent Center of mass
        c_com = np.sum(c*m)/np.sum(m)
        
        # return the 3 components separately
        return a_com, b_com, c_com


    def COM_P(self, delta, volDec):
        '''Method to compute the position of the center of mass of the galaxy 
        using the shrinking-sphere method.

        PARAMETERS
        ----------
        delta : `float, optional`
            error tolerance in kpc. Default is 0.1 kpc
        
        RETURNS
        ----------
        p_COM : `np.ndarray of astropy.Quantity'
            3-D position of the center of mass in kpc
        '''                                                                     

        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        x_COM, y_COM, z_COM = self.COMdefine(self.x, self.y, self.z, self.m)

        # compute the magnitude of the COM position vector.
        r_COM = np.sqrt(x_COM**2 + y_COM**2 + z_COM**2)


        # iterative process to determine the center of mass
        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position
        
        x_new = self.x - x_COM
        y_new = self.y - y_COM
        z_new = self.z - z_COM
        r_new = np.sqrt(x_new**2.0 + y_new**2.0 +z_new**2.0)

        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at half that radius (reduced radius)                                                           
        r_max = max(r_new)/volDec
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        change = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (change > delta):
            # selecting all particles within the reduced radius
            index2 = np.where(r_new < r_max)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]
          

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius
            
            x_COM2, y_COM2, z_COM2 = self.COMdefine(x2,y2,z2,m2)
            
            # compute the new 3D COM position
            
            r_COM2 = np.sqrt(x_COM2**2 + y_COM2**2 + z_COM2**2) 

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            change = np.abs(r_COM - r_COM2)
                                                                                        

            # Before loop continues, reset : r_max, particle separations and COM                                        

            # reduce the volume by a factor of 2 again                                                                 
            r_max /= volDec
                                                                                            

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
           
            x_new = self.x - x_COM2
            y_new = self.y - y_COM2
            z_new = self.z - z_COM2
            r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2)
          
            # create an array  to store the COM position                                                                                                                                                       
            p_COM = np.array([x_COM, y_COM, z_COM])

        # set the correct units using astropy and round all values
        # and then return the COM positon vector
        return np.around(p_COM, 2)*u.kpc
        
    def COM_V(self, x_COM, y_COM, z_COM):
          ''' Method to compute the center of mass velocity based on the center of mass
          position.

          PARAMETERS
          ----------
          x_COM : 'astropy quantity'
              The x component of the center of mass in kpc
          y_COM : 'astropy quantity'
              The y component of the center of mass in kpc
          z_COM : 'astropy quantity'
              The z component of the center of mass in kpc
              
          RETURNS
          -------
          v_COM : `np.ndarray of astropy.Quantity'
              3-D velocity of the center of mass in km/s
          '''
          
          # the max distance from the center that we will use to determine 
          #the center of mass velocity                   
          rv_max = 15.0*u.kpc

          # determine the position of all particles 
          # relative to the center of mass position
          
          xV = self.x[:]*u.kpc - x_COM
          yV = self.y[:]*u.kpc - y_COM
          zV = self.z[:]*u.kpc - z_COM
          rV = np.sqrt(xV**2 + yV**2 + zV**2)
          
          # determine the index for those particles within the max radius
    
          indexV = np.where(rV < rv_max)
          
          # determine the velocity and mass of those particles within the mass radius
    
          vx_new = self.vx[indexV]
          vy_new = self.vy[indexV]
          vz_new = self.vz[indexV]
          m_new =  self.m[indexV]
         
          
          # compute the center of mass velocity using those particles
          
          vx_COM, vy_COM, vz_COM =   self.COMdefine(vx_new,vy_new,vz_new, m_new)
          
          # create an array to store the COM velocity
          
          v_COM = np.array([vx_COM,vy_COM,vz_COM])

          # return the COM vector
          # set the correct units using astropy
          # round all values
          return np.round(v_COM, 2)*u.km/u.s
      

''' now setting up the functions to make the arrays for the histogram'''

import scipy.optimize as so


'''apologies for the hard code, it was difficult to do the loop for the specific
values in each snap for each array'''

'''beginning my code for the graph'''


# Create a COM of object for M31 Disk Using Code from Assignment 4
COMD = CenterOfMass("M31_000.txt",2)
COMD1 = CenterOfMass("M31_100.txt",2)
COMD2 = CenterOfMass("M31_200.txt",2)
COMD3 = CenterOfMass("M31_300.txt",2)
COMD35 = CenterOfMass("M31_350.txt",2)
COMD4 = CenterOfMass("M31_400.txt",2)
COMD45 = CenterOfMass("M31_450.txt",2)
COMD5 = CenterOfMass("M31_500.txt",2)
COMD55 = CenterOfMass("M31_550.txt",2)
COMD6 = CenterOfMass("M31_600.txt",2)
COMD7 = CenterOfMass("M31_700.txt",2)
COMD8 = CenterOfMass("M31_800.txt",2)

# Compute COM of M31 using disk particles
COMP = COMD.COM_P(0.1,2)
COMV = COMD.COM_V(COMP[0],COMP[1],COMP[2])
COMP1 = COMD1.COM_P(0.1,2)
COMV1 = COMD1.COM_V(COMP1[0],COMP1[1],COMP1[2])
COMP2 = COMD2.COM_P(0.1,2)
COMV2 = COMD2.COM_V(COMP2[0],COMP2[1],COMP2[2])
COMP3 = COMD3.COM_P(0.1,2)
COMV3 = COMD3.COM_V(COMP3[0],COMP3[1],COMP3[2])
COMP35 = COMD35.COM_P(0.1,2)
COMV35 = COMD35.COM_V(COMP35[0],COMP35[1],COMP35[2])
COMP4 = COMD4.COM_P(0.1,2)
COMV4 = COMD4.COM_V(COMP4[0],COMP4[1],COMP4[2])
COMP45 = COMD45.COM_P(0.1,2)
COMV45 = COMD45.COM_V(COMP45[0],COMP45[1],COMP45[2])
COMP5 = COMD5.COM_P(0.1,2)
COMV5 = COMD5.COM_V(COMP5[0],COMP5[1],COMP5[2])
COMP55 = COMD55.COM_P(0.1,2)
COMV55 = COMD55.COM_V(COMP55[0],COMP55[1],COMP55[2])
COMP6 = COMD6.COM_P(0.1,2)
COMV6 = COMD6.COM_V(COMP6[0],COMP6[1],COMP6[2])
COMP7 = COMD7.COM_P(0.1,2)
COMV7 = COMD7.COM_V(COMP7[0],COMP7[1],COMP7[2])
COMP8 = COMD8.COM_P(0.1,2)
COMV8 = COMD8.COM_V(COMP8[0],COMP8[1],COMP8[2])

# Determine positions of disk particles relative to COM 
xD = COMD.x - COMP[0].value 
yD = COMD.y - COMP[1].value 
zD = COMD.z - COMP[2].value 
xD1 = COMD1.x - COMP1[0].value 
yD1 = COMD1.y - COMP1[1].value 
zD1 = COMD1.z - COMP1[2].value 
xD2 = COMD2.x - COMP2[0].value 
yD2 = COMD2.y - COMP2[1].value 
zD2 = COMD2.z - COMP2[2].value 
xD3 = COMD3.x - COMP3[0].value 
yD3 = COMD3.y - COMP3[1].value 
zD3 = COMD3.z - COMP3[2].value 

xD35 = COMD35.x - COMP35[0].value 
yD35 = COMD35.y - COMP35[1].value 
zD35 = COMD35.z - COMP35[2].value 

xD4 = COMD4.x - COMP4[0].value 
yD4 = COMD4.y - COMP4[1].value 
zD4 = COMD4.z - COMP4[2].value 

xD45 = COMD45.x - COMP45[0].value 
yD45 = COMD45.y - COMP45[1].value 
zD45 = COMD45.z - COMP45[2].value 

xD5 = COMD5.x - COMP5[0].value 
yD5 = COMD5.y - COMP5[1].value 
zD5 = COMD5.z - COMP5[2].value 

xD55 = COMD55.x - COMP55[0].value 
yD55 = COMD55.y - COMP55[1].value 
zD55 = COMD55.z - COMP55[2].value 

xD6 = COMD6.x - COMP6[0].value 
yD6 = COMD6.y - COMP6[1].value 
zD6 = COMD6.z - COMP6[2].value 
xD7 = COMD7.x - COMP7[0].value 
yD7 = COMD7.y - COMP7[1].value 
zD7 = COMD7.z - COMP7[2].value 
xD8 = COMD8.x - COMP8[0].value 
yD8 = COMD8.y - COMP8[1].value 
zD8 = COMD8.z - COMP8[2].value 
# total magnitude
rtot = np.sqrt(xD**2 + yD**2 + zD**2)
rtot1 = np.sqrt(xD1**2 + yD1**2 + zD1**2)
rtot2 = np.sqrt(xD2**2 + yD2**2 + zD2**2)
rtot3 = np.sqrt(xD3**2 + yD3**2 + zD3**2)
rtot4 = np.sqrt(xD4**2 + yD4**2 + zD4**2)
rtot5 = np.sqrt(xD5**2 + yD5**2 + zD5**2)

rtot35 = np.sqrt(xD35**2 + yD35**2 + zD35**2)
rtot45 = np.sqrt(xD45**2 + yD45**2 + zD45**2)
rtot55 = np.sqrt(xD55**2 + yD55**2 + zD55**2)

rtot6 = np.sqrt(xD6**2 + yD6**2 + zD6**2)
rtot7 = np.sqrt(xD7**2 + yD7**2 + zD7**2)
rtot8 = np.sqrt(xD8**2 + yD8**2 + zD8**2)
# getting the sun analogs location indexes radially, and using these indexes to plot
# for x , y , and z
rtotindex0 = np.where((rtot > 7.4) & (rtot < 8.7))
rtotindex0 = rtotindex0[0]
rtotindex1 = np.where((rtot1 > 7.4) & (rtot1 < 8.7))
rtotindex1 = rtotindex1[0]
rtotindex2 = np.where((rtot2 > 7.4) & (rtot2 < 8.7))
rtotindex2 = rtotindex2[0]
rtotindex3 = np.where((rtot3 > 7.4) & (rtot3 < 8.7))
rtotindex3 = rtotindex3[0]
rtotindex4 = np.where((rtot4 > 7.4) & (rtot4 < 8.7))
rtotindex4 = rtotindex4[0]
rtotindex5 = np.where((rtot5 > 7.4) & (rtot5 < 8.7))
rtotindex5 = rtotindex5[0]

rtotindex35 = np.where((rtot35 > 7.4) & (rtot35 < 8.7))
rtotindex35 = rtotindex35[0]
rtotindex45 = np.where((rtot45 > 7.4) & (rtot45 < 8.7))
rtotindex45 = rtotindex45[0]
rtotindex55 = np.where((rtot55 > 7.4) & (rtot55 < 8.7))
rtotindex55 = rtotindex55[0]

rtotindex6 = np.where((rtot6 > 7.4) & (rtot6 < 8.7))
rtotindex6 = rtotindex6[0]
rtotindex7 = np.where((rtot7 > 7.4) & (rtot7 < 8.7))
rtotindex7 = rtotindex7[0]
rtotindex8 = np.where((rtot8 > 7.4) & (rtot8 < 8.7))
rtotindex8 = rtotindex8[0]
# getting the sun analogs for x, y and z

snap0 = len(rtotindex0)
snap1 = len(rtotindex1)
snap2 = len(rtotindex2)
snap3 = len(rtotindex3)
snap4 = len(rtotindex4)
snap5 = len(rtotindex5)

snap35 = len(rtotindex35)
snap45 = len(rtotindex45)
snap55 = len(rtotindex55)

snap6 = len(rtotindex6)
snap7 = len(rtotindex7)
snap8 = len(rtotindex8)


snaps_array = np.array([snap0 , snap1,snap2,snap3,snap4,snap5,snap35,snap45,\
                        snap55,snap6,snap7,snap8])
#x = np.sort(snaps_array)
x = [0.0,1.43,2.86,4.23,5,5.71,6.43,7.14,7.86,8.57,10,11.43]
'''making the plots'''


plt.plot(x, snaps_array)
plt.ylabel('Number of solar particles within 7.4 - 8.7 kpc of M31')
plt.xlabel('Time in future (Gyr)')
plt.xticks(fontsize=14)
plt.yticks(fontsize=12)
plt.savefig('graphofsolarparticles.png')
plt.show()

