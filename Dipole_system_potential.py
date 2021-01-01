#!/usr/bin/python3

### ADD PATH TO ProbeParticleModel ###
ppafm = '/home/matyas/Scripts/Python/pypath/ProbeParticleModel'

import numpy as np
import matplotlib.pyplot as plt
import sys
import re
import argparse

sys.path.append(ppafm)
import pyProbeParticle.GridUtils as GU


    
def ReadGeo(geo_input):
    """
    ReadGeo(geo_input)
    Reads geometry file and returns Dipole positions and momentums and Lattice vectors
    
    Parameters
    ----------
    geo_input : str
       path to the geometry file 
    
    Returns
    -------
    Dip_pos : array 
       positions of dipoles
    Dip_moment : array
       moments of dipoles
    Lattice : 2D or 3D numpy.ndarray 
       Array containing the lattice vectors of the system
    """
    
    with open(geo_input) as f:
        #read Latice vectors
        Latice = np.loadtxt((x.replace(',',' ') for x in f),comments='#', max_rows=1)
        
        #read dipoles 
        Dipoles = np.loadtxt((x.replace(',',' ') for x in f), comments='#',  skiprows=1)
    
    
    Latice = Latice.reshape(2,3)

    try:
        Dip_pos = Dipoles[:,:3]
        Dip_moment = Dipoles[:,3:]
    except:
        Dip_pos = Dipoles[:3].reshape(1,3)
        Dip_moment = Dipoles[3:].reshape(1,3)

    return Dip_pos, Dip_moment, Latice
    

def CreateGrid(Lattice, point_grid, eval_height):
    """
    CreateGrid(Lattice, point_grid, eval_height)
    Creates a grid over sample space with base defined by lattice vectors and height as a eval height
   
    Parameters
    ----------
    Lattice: array with shape (2, 3) 
       Lattice vectors of the dipole system
    point_grid: list or 1d array size 3
       number of points to sample in each dimension [x,y,z]
    eval_height: list or 1d array size 2
       lower and upper z boundary of the sample saple 
    
    Returns
    -------
    PointGrid : 3d array 
       An Array representing each point point of the sampling space coordinates of each point are PointGrid[z,y,x]
    """
    
   
    agrid = np.linspace(0,Lattice[0], point_grid[0])
    bgrid = np.linspace(0,Lattice[1], point_grid[1])
    cgrid = np.linspace([0,0,eval_height[0]], [0,0,eval_height[1]], point_grid[2])
    
    PointGrid = np.zeros((point_grid[2], point_grid[1], point_grid[0], 3))
    
    for z, c in enumerate(cgrid):
        for y, b in enumerate(bgrid):
            for x, a in enumerate(agrid):
                
                PointGrid[z, y, x] =  a + b + c

    return PointGrid


def FindMaxLatice(Lattice, cutoff):
    """
    FindMaxLatice(Lattice, cutoff)
    Finds the furtherst cell which could possibe be within the cutoff
    
    Parameters
    ----------
    Lattice : array
       Lattice vectors of the dipole system
    cutoff : float 
       cutoff distance 

    Returns
    -------
    maxLat : list
       maximum possible number of cells in each direction 
    """
    
    return [int(cutoff/np.linalg.norm(x)) + 2 for x in Lattice]


def getELongRange(Dip_moments, Lattice, cutoff, long_range=True):
    """
    getELongRange(Dipoles, Lattice, cutoff, long_range=True):
    Aproximates the constant Electric field resulting from an infinite 2D system of Dipoles

    E =  (4/3)*pi*M_cell/(r*A_cell), r - cutoff, M_cell - Net dipole moment of the cell, A_cell - Area of the cell 
    
    Parameters
    ----------
    Dip_moments: array
       array containing the dipole momentum of all dipoles 
    Lattice : array
       Lattice vectors of the dipole system
    cutoff : float 
       cutoff distance 
    long_range : bool, optional 
       if False then function will only return array of zeros
    
    Returns
    -------
    E_long_range : array of shape (3,)
       Aproximation of the long range electric field from the lattice 
    """
    
    if long_range:
        
        cell_area = np.linalg.norm(np.cross(Lattice[0], Lattice[1]))
        Net_moment = Dip_moments.sum(axis=0)
        
        E_long_range = 4*np.pi*Net_moment/(3*cutoff*cell_area)
        
    else:
        E_long_range = np.zeros([3])
    
    return E_long_range


def CreateAllDip(Dip_pos, Dip_moments, Lattice, maxLat):
    """
    CreateAllDip(Dip_pos, Dip_mom, Lattice, maxLat)
    From periodic representation of the system creates a single array for all Dipole positions and the same for momenta
    
    Parameters
    ----------
    Dip_pos : array 
       positions of dipoles
    Dip_moment : arrayle
       moments of dipoles
    Lattice : array
       Lattice vectors of the dipole system
    maxLat : list
       maximum possible number of cells in each direction 

    Returns
    -------
    All_dip_pos : array
       positions of all dipoles in all unit cell 
    All_dip_mom : array
      moments of all dipoles in all unit cell 
    """

    #create posible unit cells in both directions
    alat = np.arange(-maxLat[0], maxLat[0]+1)
    blat = np.arange(-maxLat[1], maxLat[1]+1)

    alat = np.tile(alat,(3,1)).T
    blat = np.tile(blat,(3,1)).T

    alat = alat*Lattice[0]
    blat = blat*Lattice[1]

    #create all dipoles
    All_dip_pos = np.zeros((alat.shape[0], blat.shape[0], *Dip_pos.shape))

    for x, a in enumerate(alat):
        for y, b in enumerate(blat):
            for d, p in enumerate(Dip_pos):
                All_dip_pos[x, y, d] = a + b + p

    All_dip_pos = All_dip_pos.reshape(alat.shape[0]*blat.shape[0]*Dip_pos.shape[0], 3 )
    All_dip_mom = np.tile(Dip_moments, (alat.shape[0]*blat.shape[0], 1))
    
    return All_dip_pos, All_dip_mom


def CalculateGrid(Dip_pos, Dip_moments, PointGrid, Lattice, cutoff, long_range=False, K=14.3996):
    """
    CalculateGrid(Dip_pos, Dip_moments, PointGrid, Lattice, cutoff, long_range, K)
    Calculates the potential at each point of the Point from the system of Dipoles 
    To remove periodicity set cutoff = 0 
    
    Parameters
    ----------
    Dip_pos : array 
       positions of dipoles
    Dip_moment : arrayle
       moments of dipoles
    PointGrid : 3d array 
       points in space to be evaluated, coordinates of each point are PointGrid[z, y, x]
    Lattice : array
       Lattice vectors of the dipole system
    cutoff : float 
       cutoff distance if zero only one cell will be used 
    long_range : bool, optional 
       if True the calculation includes the aproxamiation of long range electric field
    K : float, optional 
       Columbic constant default is 14.3996
    
    Returns
    -------
    Potential : 3d array
       point [z, y, x] gives potential coresponding to coordinates from PointGrid[z, y, x]
    """
    
    #mamximum number of unit cell in one direction from origin being evaluated [+-x,=-y]
    maxLat = FindMaxLatice(Lattice, cutoff)
    cutoff_cube = abs(cutoff)**3
    
    if cutoff == 0:
        #create just one unit cell 
        maxLat = [0,0]
        Lattice = np.zeros([2,3])
        E_long_range = np.zeros([3])
        cutoff_cube = np.inf
    else:
        E_long_range = getELongRange(Dip_moments, Lattice, cutoff, long_range)

    #create array of all acceptable dipole positions and corespondings moments   
    Dip_pos, Dip_moments = CreateAllDip(Dip_pos, Dip_moments, Lattice, maxLat)

    Potential = np.zeros(PointGrid[:,:,:,0].shape)

    #loop over points
    for z, Z in enumerate(PointGrid):
        for y, Y in enumerate(Z):
            for x, pos in enumerate(Y):

                #add long range correction
                Potential[z, y, x] += np.dot(pos, E_long_range)

                #position from the dipoles frame of reference 
                rel_pos = pos - Dip_pos                   
                rel_pos_cube = np.sqrt( (rel_pos**2).sum(axis=1))**3
        
                #filter out distaces over cutoff it's sad how proud I am of this 😢️😭️   
                filt = (abs(rel_pos_cube) < cutoff_cube).astype(int) 
                
                #ignore to long distances and divide the by the cube of distance
                filt = filt/rel_pos_cube
                rel_pos *= np.tile(filt,(3,1)).T

                #now its just a sum of contributions of all dipoles in all directions
                Potential[z, y, x] += (rel_pos*Dip_moments).sum()

    Potential*=K
    
    # I am hungry
    return Potential


def DrawDipoles(Potential, Dip_pos, Dip_mom, Lattice, plotDip, axis='on', delta=0.8):
    #Draws a 2D repesentantion of the Dipoles
    #just make your own I don't know how to make this general
    
    size = max(Lattice[0][0],Lattice[1][1])
    plt.xlim(0, 1.1*size)
    plt.ylim(0, 1.1*size)

    if plotDip:
        Dip_mom *= delta
        Dip_pos -=  Dip_mom/2
        arrows = np.hstack((Dip_pos, Dip_mom))
        
    for i, Z in enumerate(Potential):

        plt.axis(axis)
        plt.imshow(Z)
        plt.draw()        

        if plotDip:
                
            for arr in arrows:
                mom = arr[3:]
                pos = arr[:3]
                
                plt.arrow(pos[0],pos[1],mom[0],mom[1], color='black', head_width=.5)

        plt.savefig('%s_pot.png' % i)
        plt.close()
        

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Draw potential from a dipole system')

    parser.add_argument('geo_input', help='file containing possitions and moments of all dipoles')
    parser.add_argument('-p', '--point_grid', dest='point_grid', type=int, nargs=3,
                        help='number of points along x y z lines in point grid')
    parser.add_argument('-c', '--cut_off', dest='cutoff', type=float, default=0,
                        help='cutoff of the dipoles, if zero or not selected system will be treated as non-periodic')
    parser.add_argument('-o', '--outfile', dest='outfile', default='Dipoles.xsf',
                        help='where to store result the final xsf geometry default is Dipoles.xsf')
    parser.add_argument('-e', '--evaluation_height', dest='eval_height', type=float, nargs=2,
                        help='the maximum and minimum height of the point grid')
    parser.add_argument('-l', '--long_range', dest='long_range', action='store_true', help='add an infinite system correction')
    parser.add_argument('-k', '--columbic_constat', dest = 'K', type=float, default = 14.3996, help='default is 14.3996')
        
    

    args = parser.parse_args()
    
    Dip_pos, Dip_moment, Lattice = ReadGeo(args.geo_input)
    PointGrid = CreateGrid(Lattice, args.point_grid, args.eval_height)
                    
    Potential = CalculateGrid(Dip_pos, Dip_moment, PointGrid, Lattice, args.cutoff, args.long_range, args.K)
    
    
    
    #create a lattice of the sample space  
    c = np.array([0,0, args.eval_height[1] - args.eval_height[0]])
    Lat  = np.vstack((c, Lattice))

    #create an xsf file as a vector representation
    GU.saveXSF(args.geo_input.split('.')[0], Potential, Lat)
    
            
        

