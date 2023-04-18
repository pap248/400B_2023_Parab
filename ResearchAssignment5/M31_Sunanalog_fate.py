#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: paarth

The fate of Sun analogs in Andromeda Galaxy
This program will find the fate of the sun analogs that are in M31, 8kpc
from galactic center. 

How do the postitions of the Sun analogs change ? 
"""
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

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
    # write your own code to compute the generic COM 
    #using Eq. 1 in the homework instructions
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
    # write your own code below
    r_COM = np.sqrt(x_COM**2 + y_COM**2 + z_COM**2)


    # iterative process to determine the center of mass                                                            

    # change reference frame to COM frame                                                                          
    # compute the difference between particle coordinates                                                          
    # and the first guess at COM position
    # write your own code below
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
        # select all particles within the reduced radius (starting from original x,y,z, m)
        # write your own code below (hints, use np.where)
        index2 = np.where(r_new < r_max)
        x2 = self.x[index2]
        y2 = self.y[index2]
        z2 = self.z[index2]
        m2 = self.m[index2]
      

        # Refined COM position:                                                                                    
        # compute the center of mass position using                                                                
        # the particles in the reduced radius
        # write your own code below
        x_COM2, y_COM2, z_COM2 = self.COMdefine(x2,y2,z2,m2)
        
        # compute the new 3D COM position
        # write your own code below
        r_COM2 = np.sqrt(x_COM2**2 + y_COM2**2 + z_COM2**2) 

        # determine the difference between the previous center of mass position                                    
        # and the new one.                                                                                         
        change = np.abs(r_COM - r_COM2)
        # uncomment the following line if you want to check this                                                                                               
        # print ("CHANGE = ", CHANGE)                                                                                     

        # Before loop continues, reset : r_max, particle separations and COM                                        

        # reduce the volume by a factor of 2 again                                                                 
        r_max /= volDec
        # check this.                                                                                              
        #print ("maxR", r_max)                                                                                      

        # Change the frame of reference to the newly computed COM.                                                 
        # subtract the new COM
        # write your own code below
        x_new = self.x - x_COM2
        y_new = self.y - y_COM2
        z_new = self.z - z_COM2
        r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2)
      

        # set the center of mass positions to the refined values                                                   
        x_COM = x_COM2
        y_COM = y_COM2
        z_COM = z_COM2
        r_COM = r_COM2
        
        xSun = np.where(x_COM > 7 and x_COM < 9)
        ySun = np.where(y_COM > 7 and y_COM < 9)
        zSun = np.where(z_COM > 7 and z_COM < 9)

        # create an array  to store the COM position                                                                                                                                                       
        p_COM = np.array([x_COM, y_COM, z_COM])
        sun_analogs_pos = np.array([xSun, ySun, zSun])

    # set the correct units using astropy and round all values
    # and then return the COM positon vector
    return np.around(sun_analogs_pos, 2)*u.kpc
    
def COM_V(self, xSun, ySun, zSun):
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
      
      xV = self.x[:]*u.kpc - xSun
      yV = self.y[:]*u.kpc - ySun
      zV = self.z[:]*u.kpc - zSun
      rV = np.sqrt(xV**2 + yV**2 + zV**2)
      
      # determine the index for those particles within the max radius
      # write your own code below
      indexV = np.where(rV < rv_max)
      
      # determine the velocity and mass of those particles within the mas radius
      # write your own code below
      vx_new = self.vx[indexV]
      vy_new = self.vy[indexV]
      vz_new = self.vz[indexV]
      m_new =  self.m[indexV]
     
      
      # compute the center of mass velocity using those particles
      # write your own code below
      vx_COM, vy_COM, vz_COM =   self.COMdefine(vx_new,vy_new,vz_new, m_new)
      
      # create an array to store the COM velocity
      # write your own code below
      v_COM = np.array([vx_COM,vy_COM,vz_COM])
      #sun_analogs_vel = np.where(v_COM > 7 and v_COM < 9)

      # return the COM vector
      # set the correct units using astropy
      # round all values
      return np.round(v_COM, 2)*u.km/u.s
  
  

def SunAnalogs(self, delta):
    '''Method to compute the position of the stars from the center of galaxy 

    PARAMETERS
    ----------
    delta : `float, optional`
        error tolerance in kpc. Default is 0.1 kpc
    
    RETURNS
    ----------
    sun_analogs_pos : `np.ndarray of astropy.Quantity'
        3-D position in kpc of sun analogs
    '''                                                                     
                                                                                                                       
    x_pos, y_pos, z_pos = self.COMdefine(self.x, self.y, self.z, self.m)
    
    # compute the magnitude of the position vector
    positions = np.sqrt(x_pos**2 + y_pos**2 + z_pos**2)
    
    ''' now that the radial positions have been found , I will make a loop to
    run through the array of positions and add the ones that are about 8 kpc 
    to a list '''

    # locations of the sun analogs
    sun_analogs_pos = np.where(positions > 7 and positions < 9)
    

    
    return sun_analogs_pos


import scipy.optimize as so

def find_confidence_interval(x, pdf, confidence_level):
    return pdf[pdf > x].sum() - confidence_level

# Info about **kwargs, *args 
#https://book.pythontips.com/en/latest/args_and_kwargs.html

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
        
    Example Usage
    -------------
     density_contour(x pos, y pos, contour res, contour res, axis, colors for contours)
     e.g.:
     density_contour(xD, yD, 80, 80, ax=ax, 
             colors=['red','orange', 'yellow', 'orange', 'yellow'])

    """

    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), normed=True)
    # NOTE : if you are using the latest version of python, in the above: 
    # instead of normed=True, use density=True
    
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

    pdf = (H*(x_bin_sizes*y_bin_sizes))  
    
    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T  # transpose of distribution fxn
    fmt = {}
    # Contour Levels Definitions
    #one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
    #brentq is root finding method
    #two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
    #three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
    
    # You might need to add a few levels
    # I added a few between 1 and 2 sigma to better highlight the spiral arm
    #one_sigma1 = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.80))
    #one_sigma2 = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.90))

    # Array of Contour levels. Adjust according to the above
    #levels = [one_sigma, one_sigma1, one_sigma2, two_sigma, three_sigma][::-1]
    
    # contour level labels  Adjust accoding to the above.
    #strs = ['0.68', '0.8','0.9','0.95', '0.99'][::-1]

    
    ###### 
    
    if ax == None:
        contour = plt.contour(X, Y, Z, origin="lower", **contour_kwargs)

    else:
        contour = ax.contour(X, Y, Z, origin="lower", **contour_kwargs)
        
    return contour


# Create a COM of object for M31 Disk Using Code from Assignment 4
COMD = CenterOfMass("M31_000.txt",2)
# Compute COM of M31 using disk particles
COMP = COMD.COM_P(0.1,1)
COMV = COMD.COM_V(COMP[0],COMP[1],COMP[2])
# Determine positions of disk particles relative to COM 
xD = COMD.x - COMP[0].value 
yD = COMD.y - COMP[1].value 
zD = COMD.z - COMP[2].value 

# total magnitude
rtot = np.sqrt(xD**2 + yD**2 + zD**2)

# Determine velocities of disk particles relatiev to COM motion
vxD = COMD.vx - COMV[0].value 
vyD = COMD.vy - COMV[1].value 
vzD = COMD.vz - COMV[2].value 

# total velocity 
vtot = np.sqrt(vxD**2 + vyD**2 + vzD**2)

# Vectors for r and v 
r = np.array([xD,yD,zD]).T # transposed for the Rotate Function later
v = np.array([vxD,vyD,vzD]).T
# M31 Disk Density 
fig, ax= plt.subplots(figsize=(10, 10))

# plot the particle density for M31 using a 2D historgram
# can modify bin number to make the plot smoother
plt.hist2d(xD, yD, bins=150, norm=LogNorm(), cmap='magma')
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
#density_contour(xD, yD, 80, 80, ax=ax)


# Save to a file
plt.savefig('M31Disk.png')

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
plt.hist2d(rn[:,0], rn[:,2], bins=150, norm=LogNorm(), cmap='magma')
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
#density_contour(rn[:,0], rn[:,2], 80, 80, ax=ax)

plt.savefig('EdgeOn_Density.png')

# Rotated M31 Disk - FACE ON

# M31 Disk Density 
fig, ax= plt.subplots(figsize=(10, 10))

# plot the particle density for M31 
# can modify bin number (bin =100 smoothest)
plt.hist2d(rn[:,0], rn[:,1], bins=150, norm=LogNorm(), cmap='magma')
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
#density_contour(rn[:,0], rn[:,1], 80, 80, ax=ax)

plt.savefig('FaceOn_Density.png')




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

'''plot 2 is the sun analogs as a function of time or snapshots, amount of
sun analogs at each snapshot '''









