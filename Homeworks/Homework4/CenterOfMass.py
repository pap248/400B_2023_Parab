
# Homework 4
# Center of Mass Position and Velocity
# Solutions: G.Besla, R. Li, H. Foote


"""
@author: paarth

This program will show how each galaxy in our Local group moves accounting
for the forces that they show on each other. We will be computing the 
center of mass position and velocity vectors of each galaxy at any point in time.
Using Cartesian coordinates, and the data files of MW, M31, and M33. 

"""

# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl

from Readfile import Read




class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot

    def __init__(self, filename, ptype):
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

        # store the mass, positions, velocities of only the particles of the given type
        # shown below
        self.m = self.data['m'][self.index]
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]
        self.r = np.sqrt(((self.x)**2) + ((self.y)**2) + ((self.z)**2))


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
        # computing the generic COM below
        # using Eq. 1 in the homework instructions
        # xcomponent Center of mass
        x_sum = np.sum(np.multiply(self.x,self.m))
        a_com = (x_sum) / np.sum(self.m)
        # ycomponent Center of mass
        y_sum = np.sum(np.multiply(self.y,self.m))
        b_com = (y_sum) / np.sum(self.m)
        # zcomponent Center of mass
        z_sum = np.sum(np.multiply(self.z,self.m))
        c_com = (z_sum) / np.sum(self.m)
        
        # returning the 3 components separately
        return a_com, b_com, c_com
    
    
    def COM_P(self, delta):
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

        # Trying a first guess at the COM position by calling COMdefine                                                   
        x_COM, y_COM, z_COM = self.COMdefine(self.x, self.y, self.z, self.m)
        
        # computing the magnitude of the COM position vector
        r_COM = np.sqrt(((x_COM)**2) + ((y_COM)**2) + ((z_COM)**2))


        # iterative process to determine the center of mass                                                            

        # changing reference frame to COM frame                                                                          
        # computing the difference between particle coordinates                                                          
        # and the first guess at COM position
        x_new = (self.x) - x_COM
        y_new = (self.y) - y_COM
        z_new = (self.z) - z_COM
        r_new = np.sqrt(((x_new)**2) + ((y_new)**2) + ((z_new)**2))

        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at half that radius (reduced radius)                                                           
        r_max = max(r_new)/2.0
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        change = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (change > delta):
            # select all particles within the reduced radius 
            # starting from original x, y, z, m
            index2 = np.where(r_new < r_max)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]
            r2 = np.sqrt(((x2)**2) + ((y2)**2) + ((z2)**2))

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius
            x_COM2, y_COM2, z_COM2 = self.COMdefine(x2, y2, z2, m2)
            
            # compute the new 3D COM position
            r_COM2 = np.sqrt(((x_COM2)**2) + ((y_COM2)**2) + ((z_COM2)**2))

            # determine the difference between the previous center of mass position                                    
            # and the new one                                                                                       
            change = np.abs(r_COM - r_COM2)
            
            # following line is to check                                                                                                
            #print ("CHANGE = ", change)                                                                                     

            # Before loop continues, reset : r_max, particle separations and COM                                        

            # reduce the volume by a factor of 2 again                                                                 
            r_max /= 2.0
            # checking                                                                                              
            #print ("maxR", r_max)                                                                                      

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            x_new = x2 - x_COM2
            y_new = y2 - y_COM2
            z_new = z2 - z_COM2
            r_new = np.sqrt(((x_new**2) + ((y_new)**2) + ((z_new)**2)))

            # set the center of mass positions to the refined values                                                   
            x_COM = x_COM2
            y_COM = y_COM2
            z_COM = z_COM2
            r_COM = r_COM2

            # creating an array (np.array) to store the COM position                                                                                                                                                       
            p_COM = np.array([x_COM, y_COM, z_COM])
            p_COM = np.round(p_COM * u.kpc, 2 )

        # setting the correct units using astropy and round all values
        # and then return the COM positon vector
        return p_COM
        
        
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
        # the center of mass velocity                   
        rv_max = 15.0*u.kpc

        # determine the position of all particles relative to the 
        # center of mass position (x_COM, y_COM, z_COM)
        xV = (self.x * u.kpc) - x_COM
        yV = (self.y * u.kpc) - y_COM
        zV = (self.z * u.kpc) - z_COM
        rV = np.sqrt(((xV**2) + ((yV)**2) + ((zV)**2)))
        
        # determine the index for those particles within the max radius
        indexV = np.where(rV < rv_max)
        
        # determine the velocity and mass of those particles within the mas radius
        vx_new = self.vx[indexV]
        vy_new = self.vy[indexV]
        vz_new = self.vz[indexV]
        m_new =  self.m[indexV]
        
        # compute the center of mass velocity using those particles
        vx_COM, vy_COM, vz_COM = self.COMdefine(vx_new, vy_new, vz_new, m_new)
        
        # create an array to store the COM velocity
        v_COM = np.array([vx_COM, vy_COM, vz_COM])
        v_COM = np.round(v_COM * u.km / u.s, 2 )

        # return the COM vector
        # set the correct units using astropy
        # round all values
        return v_COM
    

# ANSWERING QUESTIONS
#######################
if __name__ == '__main__' : 

    # Create a Center of mass object for the MW, M31 and M33
    # below is an example of using the class for MW
    # doing type 2 for disk
    MW_COM = CenterOfMass("MW_000.txt", 2)
    
    # setting for other galaxies
    M31_COM = CenterOfMass("M31_000.txt", 2)
    M33_COM = CenterOfMass("M33_000.txt", 2)
    
    # below gives an example of calling the class's functions
    # MW:   store the position and velocity COM
    MW_COM_p = MW_COM.COM_P(0.1)
    # print(MW_COM_p)
    
    # below is for M31
    M31_COM_p = M31_COM.COM_P(0.1)
    
    # below is for M33
    M33_COM_p = M33_COM.COM_P(0.1)

    # now velocity vectors
    MW_COM_v = MW_COM.COM_V(MW_COM_p[0], MW_COM_p[1], MW_COM_p[2])

    # below is for M31
    M31_COM_v = M31_COM.COM_V(M31_COM_p[0], M31_COM_p[1], M31_COM_p[2])
    
    # below is for M33
    M33_COM_v = M33_COM.COM_V(M33_COM_p[0], M33_COM_p[1], M33_COM_p[2])
    
    
    
    
    '''
    I got my position functions to work but I couldn't figure 
    out why my velocity function was not working so the below 
    answers to the questions are only using the position vectors.
    '''
    
    
    # question 1
    
    print('\nQuestion 1')

    print('COM position vector for MW is', MW_COM_p)
    print('COM position vector for M31 is', M31_COM_p)
    print('COM position vector for M33 is', M33_COM_p)
    
    # question 2
    
    print('\nQuestion 2')
    
    MW_M31_sep_x = MW_COM_p[0] - M31_COM_p[0]
    MW_M31_sep_y = MW_COM_p[1] - M31_COM_p[1]
    MW_M31_sep_z = MW_COM_p[2] - M31_COM_p[2]
    sep_p = np.sqrt(((MW_M31_sep_x**2) + ((MW_M31_sep_y)**2) + ((MW_M31_sep_z)**2)))
    print('The magnitude of current separation between MW and M31 is'\
          , np.round(sep_p,3))
    
    # MW_M31_sep_vx = MW_COM_v[0] - M31_COM_v[0]
    # MW_M31_sep_vy = MW_COM_v[1] - M31_COM_v[1]
    # MW_M31_sep_vz = MW_COM_v[2] - M31_COM_v[2]
    # sep_v = np.sqrt(((MW_M31_sep_vx**2) + ((MW_M31_sep_vy)**2) + ((MW_M31_sep_vz)**2)))
    # print(sep_v)
    
    # question 3
    
    print('\nQuestion 3')
    
    M33_M31_sep_x = M33_COM_p[0] - M31_COM_p[0]
    M33_M31_sep_y = M33_COM_p[1] - M31_COM_p[1]
    M33_M31_sep_z = M33_COM_p[2] - M31_COM_p[2]
    sep_p_M33_M31 = np.sqrt(((M33_M31_sep_x**2) + ((M33_M31_sep_y)**2) + ((M33_M31_sep_z)**2)))
    print('The magnitude of current separation between M33 and M31 is'\
          , np.round(sep_p_M33_M31,3))
    
    # question 4
    
    print('\nQuestion 4')
    
    print('''
          The process to determine the center of mass is very important for 
          many reasons. Such as when the merger will initially take place, how 
          long it could last, and as well as how close M31 is to MW at all times, 
          etc. We can run many simulations and be prepared on what to expect.
          ''')
    





















