
# # Homework 7 Template
# 
# Rixin Li & G . Besla
# 

"""
@author: Paarth

This program will predict the future trajectory of M33 in the frame of M31.

"""


# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
# import Latex module so we can display the results with symbols
from IPython.display import Latex

# import CenterOfMass to determine the COM pos/vel of M33
from CenterOfMass_Solution import CenterOfMass

# import the GalaxyMass to determine the mass of M31 for each component
from GalaxyMass import ComponentMass
from Readfile import Read
# # M33AnalyticOrbit




class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self, filename): 
        """ This class will create functions that determine 
        the acceleration M33 feels from M31 and integrate its current 
        position and velocity forwards in time.
        The input is the filename for the file to store the integrated orbit"""

        ### gravitational constant (4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        
        ### store the output file name
        self.filename = filename
        
        ### get the current pos/vel of M33 
        # create an instance of the  CenterOfMass class for M33 
        M33_COM = CenterOfMass("M33_000.txt", 2)

        # store the position VECTOR of the M33 COM
        M33_COM_p = M33_COM.COM_P(0.1).value

        # store the velocity VECTOR of the M33 COM 
        M33_COM_v = M33_COM.COM_V(M33_COM_p[0],M33_COM_p[1],M33_COM_p[2]).value
        
        
        ### get the current pos/vel of M31 
        # create an instance of the  CenterOfMass class for M31 
        M31_COM = CenterOfMass("M31_000.txt", 2)

        # store the position VECTOR of the M31 COM 
        M31_COM_p = M31_COM.COM_P(0.1).value

        # store the velocity VECTOR of the M31 COM 
        M31_COM_v = M31_COM.COM_V(M31_COM_p[0],M31_COM_p[1],M31_COM_p[2])
        
        ### store the DIFFERENCE between the vectors posM33 - posM31
        # create two VECTORs self.r0 and self.v0 and have them be the
        # relative position and velocity VECTORS of M33
        self.r0 = M33_COM_p - M31_COM_p
        self.v0 = M33_COM_v - M31_COM_v
        
        ### get the mass of each component in M31 disk
        self.rdisk = 5

        self.Mdisk = 0.120 # from assignment 3
        
        ### bulge
        self.rbulge = 1

        self.Mbulge = 0.019 # assignment 3
        
        # Halo
        self.rhalo = 62 # scale length from assignment 5 for M31

        self.Mhalo = 1.921 # assignment 3
    
    def HernquistAccel(self, M , r_a , r):
        # r is vector
        # inputs are mass , scale length and position vector
        """ This function finds the gravitational acceleration induced
        by a Hernquist profile
        It will return the acceleration vector from a Hernquist potential."""
        
        ### store the magnitude of the position vector
        rmag = np.sqrt(r.dot(r)) # easiest way to get magnitude
        
        ### store the Acceleration
        Hern = - (self.G * M / (rmag * (r_a + rmag) ** 2) )* r
        # the last r is a vector
        
        return Hern
    
    
    
    def MiyamotoNagaiAccel(self, M, r_d, r):
        # inputs are mass , scale length and position vector
        """ This approximation mimics the exponential disk profile at distances
        far from the disk. It is called a Miyamoto-Nagai 1975 profile.
        Will return the acceleration vector from a Miyamoto-Nagai profile"""

        
        ### Acceleration
        # this is for the potential for this profile
        R = np.sqrt(r[0]**2 + r[1]**2) # for x^2 + y^2
        
        self.zdisk = self.rdisk/5
        B = self.rdisk + np.sqrt(r[2]**2 + self.zdisk**2)
        
        MNA = - self.G * M * r / ((R**2 + B**2) ** 1.5) * \
            np.array([1, 1, B/np.sqrt(r[2]**2 + self.zdisk**2)])
        # MNA is the acceleration
       
        return MNA
        # the np.array allows for a different value for the z component of the acceleration
     
    
    def M31Accel(self, r):
        """ This function sums all acceleration vectors from each galaxy component"""
        
        # below is the accelerations for the hernquist profile for each ptype
        halo_a = HernquistAccel(self.Mhalo, self.rhalo, r)
        disk_a = HernquistAccel(self.Mdisk, self.rdisk, r)
        bulge_a = HernquistAccel(self.Mbulge, self.rbulge, r)
        
        sum_a = np.sum(np.array([halo_a, disk_a, bulge_a]))
        
        # return the SUM of the output of the acceleration functions - a vector
        return sum_a
    
    
    
    def LeapFrog(self, dt , r , v): 
        """  This fucntion treats M33 as a point mass and  will adopt a variant 
        of the “Leap Frog” integration scheme given that a is a pure function of r
        The inputs are the time interval dt
        starting position vector r for the M33 COM position relative to the M31
        starting velocity vector v for the M33 relative to M31"""
        
        # predict the position at the next half timestep
        rhalf = r + v * dt/2 # from eq 7 in assignment 7
        
        # predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        vnew = v + (self.M31Accel(rhalf) * dt) # from eq 8
        
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don't know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        rnew = rnew = rhalf + (vnew * dt/2) # from eq 10
        
        # now return the new position and velcoity vectors
        return rnew, vnew
    
    
    
    def OrbitIntegration(self, t0, dt, tmax):
         """ This function loops over the LeapFrog integrator to solve the equations
         of motion and compute the future orbit of M33 for 10 Gyr into the future
         Inputs are starting time t0, time interval dt and final time tmax"""

         # initialize the time to the input starting time
         t = t0
        
        # initialize an empty array of size :  rows int(tmax/dt)+2  , columns 7
         orbit = np.zeros(int(tmax/dt)+2, 7)
        
        # initialize the first row of the orbit
         orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
        # this above is equivalent to 
        # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        
        
        # initialize a counter for the orbit.  
         i = 1 # since we already set the 0th values, we start the counter at 1
        
        # start the integration (advancing in time steps and computing LeapFrog at each step)
         while (t < tmax):
            # as long as t has not exceeded the maximal time 
         
            # advance the time by one timestep, dt
            t += dt
            # store the new time in the first column of the ith row
            orbit[i , 0] = t
            
            # advance the position and velocity using the LeapFrog scheme
            
            p, v = self.LeapFrog(dt, orbit[i-1, 1:4], orbit[i-1, 4:])
            
         
    
            # store the new position vector into the columns with indexes 1,2,3 
            # of the ith row of orbit
            
            # update counter i
            i += 1 
        
        
        # write the data to a file
            np.savetxt(self.filename, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        
        # there is no return function

