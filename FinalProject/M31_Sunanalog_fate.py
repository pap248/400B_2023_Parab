#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: paarth

This program takes the snapshot of the M31 galaxy and makes an image
of the Sun Analogs within the galaxy. Solar particles that are 7.4 - 8.7
kpc away from the center of the galaxy. It does this over 800 snapshots,
focusing on the merger of the Milky Way Galaxy, and Andromeda Galaxy. 

The fate of Sun analogs in Andromeda Galaxy


How do the postitions of the Sun analogs change ? 
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


''' First, I plan on reading in the first snap shot of M31 and then looking at
how many sun analogs are there in total, seeing how many of them are 
about 8 kpc from the center of the galaxy'''

''' will look at every 50 - 100 snapshots. looking through the file, 
finding the solar particles that are distanced 7.4 - 8.7 kpc away from center of M31.'''

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
      

''' now setting up the functions to make the contour plots for the snapshots'''

import scipy.optimize as so

def find_confidence_interval(x, pdf, confidence_level):
    return pdf[pdf > x].sum() - confidence_level

def density_contour(xdata, ydata, nbins_x, nbins_y, ax=None, **contour_kwargs):
    """ Create a density contour plot.
    Parameters
    ----------
    xdata : numpy.ndarray
    ydata : numpy.ndarray
    nbins_x : int
        Number of bins along x dimension
    nbins_y : int
        Number of bins along y dimension
    ax : matplotlib.Axes (optional)
        If supplied, plot the contour to this axis. Otherwise, open a new figure
    contour_kwargs : dict
        kwargs to be passed to pyplot.contour()
        i.e. unknown number of keywords 
    
    """

    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), normed=True, \
                                       density=True)
    
    # instead of normed=True, use density=True
    
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

    pdf = (H*(x_bin_sizes*y_bin_sizes))  
    
    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T  # transpose of distribution fxn
    fmt = {}
    
    
    if ax == None:
        contour = plt.contour(X, Y, Z, origin="lower", **contour_kwargs)

    else:
        contour = ax.contour(X, Y, Z, origin="lower", **contour_kwargs)
        
    return contour

''' now taking in the arrays that have been made of the position to find the 
radial values and find the sun analogs radially, then use those indexes for
the following snapshots to follow those particles'''

'''beginning my process for the images'''


# Create a COM of object for M31 Disk Using Code from Assignment 4
COMD = CenterOfMass("M31_000.txt",2)
#creating the base file to get the indexes

# Compute COM of M31 using disk particles
COMP = COMD.COM_P(0.1,1)
COMV = COMD.COM_V(COMP[0],COMP[1],COMP[2])
# Determine positions of disk particles relative to COM 
xD = COMD.x - COMP[0].value 
yD = COMD.y - COMP[1].value 
zD = COMD.z - COMP[2].value 
# total magnitude
rtot = np.sqrt(xD**2 + yD**2 + zD**2)
# getting the sun analogs location indexes radially, and using these indexes to plot
# for x , y , and z
rtotindex0 = np.where((rtot > 7.4) & (rtot < 8.7))
rtotindex0 = rtotindex0[0]
# getting the sun analogs for x, y and z
xD0 = xD[rtotindex0]
yD0 = yD[rtotindex0]
zD0 = zD[rtotindex0]

# getting the velocity values
vxD = COMD.vx - COMV[0].value 
vyD = COMD.vy - COMV[1].value 
vzD = COMD.vz - COMV[2].value 
# getting the velocity index locations
vxD0 = vxD[rtotindex0]
vyD0 = vyD[rtotindex0]
vzD0 = vzD[rtotindex0]

'''snapshot changing for the plots, then saving each one'''
COMD1 = CenterOfMass("M31_800.txt",2)
# Compute COM of M31 using disk particles
COMP1 = COMD1.COM_P(0.1,1)
COMV1 = COMD1.COM_V(COMP1[0],COMP1[1],COMP1[2])
# Determine positions of disk particles relative to COM 400 snapshot
xD1 = COMD1.x - COMP1[0].value 
yD1 = COMD1.y - COMP1[1].value 
zD1 = COMD1.z - COMP1[2].value 

xDsnap = xD1[rtotindex0]
yDsnap = yD1[rtotindex0]
zDsnap = zD1[rtotindex0]

# Determine velocities of disk particles relative to COM motion
vxD1 = COMD1.vx - COMV1[0].value 
vyD1 = COMD1.vy - COMV1[1].value 
vzD1 = COMD1.vz - COMV1[2].value 

# total velocity not needed
# vtot = np.sqrt(vxD**2 + vyD**2 + vzD**2)

# getting the velocities for those sun analogs
vxDsnap = vxD1[rtotindex0]
vyDsnap = vyD1[rtotindex0]
vzDsnap = vzD1[rtotindex0]

'''making the plots'''

# Vectors for r and v 
r = np.array([xDsnap,yDsnap,zDsnap]).T # transposed for the Rotate Function later
v = np.array([vxDsnap,vyDsnap,vzDsnap]).T
# M31 Disk Density 
fig, ax= plt.subplots(figsize=(10, 10))

# plot the particle density for M31 using a 2D historgram
# can modify bin number to make the plot smoother
plt.hist2d(xDsnap, yDsnap, bins=400, norm=LogNorm(), cmap='magma')
plt.colorbar()

# Add axis labels
plt.xlabel('x (kpc)', fontsize=22)
plt.ylabel('y (kpc)', fontsize=22)

#set axis limits
plt.ylim(-40,40)
plt.xlim(-40,40)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size


# Save to a file
plt.savefig('M31Disk800.png')

def RotateFrame(posI,velI):
    """a function that will rotate the position and velocity vectors
    so that the disk angular momentum is aligned with z axis. 
    
    PARAMETERS
    ----------
        posI : `array of floats`
             3D array of positions (x,y,z)
        velI : `array of floats`
             3D array of velocities (vx,vy,vz)
             
    RETURNS
    -------
        pos: `array of floats`
            rotated 3D array of positions (x,y,z) such that disk is in the XY plane
        vel: `array of floats`
            rotated 3D array of velocities (vx,vy,vz) such that disk angular momentum vector
            is in the +z direction 
    """
    
    # compute the angular momentum
    L = np.sum(np.cross(posI,velI), axis=0)
    # normalize the vector
    L_norm = L/np.sqrt(np.sum(L**2))


    # Set up rotation matrix to map L_norm to z unit vector (disk in xy-plane)
    
    # z unit vector
    z_norm = np.array([0, 0, 1])
    
    # cross product between L and z
    vv = np.cross(L_norm, z_norm)
    s = np.sqrt(np.sum(vv**2))
    
    # dot product between L and z 
    c = np.dot(L_norm, z_norm)
    
    # rotation matrix
    I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    v_x = np.array([[0, -vv[2], vv[1]], [vv[2], 0, -vv[0]], [-vv[1], vv[0], 0]])
    R = I + v_x + np.dot(v_x, v_x)*(1 - c)/s**2

    # Rotate coordinate system
    pos = np.dot(R, posI.T).T
    vel = np.dot(R, velI.T).T
    
    return pos, vel

# compute the rotated velocity vectors
rn, vn = RotateFrame(r,v)
# Rotated M31 Disk - EDGE ON (XZ)

# M31 Disk Density 
fig, ax= plt.subplots(figsize=(10, 10))

# plot the particle density for M31 
# can modify bin number (bin =100 smoothest)
plt.hist2d(rn[:,0], rn[:,2], bins=400, norm=LogNorm(), cmap='magma')
plt.colorbar()

# Add axis labels
plt.xlabel('x (kpc)', fontsize=22)
plt.ylabel('y (kpc)', fontsize=22)

#set axis limits
plt.ylim(-10,10)
plt.xlim(-40,40)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# make the contour plot
# x pos, y pos, contour res, contour res, axis, colors for contours.


plt.savefig('EdgeOn_Density800.png')

# Rotated M31 Disk - FACE ON

# M31 Disk Density 
fig, ax= plt.subplots(figsize=(10, 10))

# plot the particle density for M31 
# can modify bin number (bin =100 smoothest)
plt.hist2d(rn[:,0], rn[:,1], bins=400, norm=LogNorm(), cmap='magma')
plt.colorbar()

# Add axis labels
plt.xlabel('x (kpc)', fontsize=22)
plt.ylabel('y (kpc)', fontsize=22)

#set axis limits
plt.ylim(-40,40)
plt.xlim(-40,40)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# make the contour plot
# x pos, y pos, contour res, contour res, axis, colors for contours.


plt.savefig('FaceOn_Density800.png')








