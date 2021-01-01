# Dipole-system-potential
Python script that creates a map of potential on a 3D grid caused by a 2 dimensional periodic system of dipoles and saves it in a .xsf file.

## Libraries
Python 3.6.9
* numpy 1.18.5
* argparse 1.4.0
* sys 
* re 
* ProbeParticleModel https://github.com/ProkopHapala/ProbeParticleModel
  * for creating the .xsf file if you wish 
  * add the path to the library in the script as ppafm right at the top
  * I found a similar function in ase but did not get it to work 😪️

## Usage 
### Terminal 
<pre>
Dip_potential_newVersion2.py [-h] [-p POINT_GRID POINT_GRID POINT_GRID]  
                                  [-c CUTOFF] [-o OUTFILE]  
                                  [-e EVAL_HEIGHT EVAL_HEIGHT] [-l]  
                                  [-k--columbic_constat K]  
                                  geo_input 

positional arguments:
  geo_input             file containing positions and moments of all dipoles

optional arguments:
  -h, --help            show this help message and exit
  -p POINT_GRID POINT_GRID POINT_GRID, --point_grid POINT_GRID POINT_GRID POINT_GRID
                        number of points along x y z lines in point grid
  -c CUTOFF, --cut_off CUTOFF
                        cutoff of the dipoles, if zero or not selected system
                        will be treated as non-periodic
  -o OUTFILE, --outfile OUTFILE
                        where to store the final xsf geometry default
                        is Dipoles.xsf
  -e EVAL_HEIGHT EVAL_HEIGHT, --evaluation_height EVAL_HEIGHT EVAL_HEIGHT
                        the maximum and minimum height of the point grid
  -l, --long_range      add an infinite system correction
  -k, --columbic_constat K
                        default is 14.3996
</pre>

#### Returns
outfile - .xsf file with the potential grid


### Python 
```python
import numpy as np
import Dipole_potential_map as Dip

Lattice = np.array([[3,0,0],
                    [0,5,0]])

eval_height = [10,15] # minimum and maximum height above the system                                                     
npoints = [12,20,5] # number sample points in x y z 

dip_pos = np.array([[0,1,0],
                   [2,3,0]]) # position of dipoles                                     
                   
dip_mom = np.array([[1,1,0],
                   [-1,-1,0]]) # their momentum 
                   
cutoff = 100

grid = Dip.CreateGrid(Lattice, npoints, eval_height) # points to be evaluated                                            
Potential = Dip.CalculateGrid(dip_pos,dip_mom, grid, Lattice, cutoff) # V(x,y,z) = Potential[z, y, x]  
```

## Example 

dip.in contains following

<pre>
10 0 0, 0 10 0 # Lattice vectors x y z, x y z 

#dipole positions and momentum
#x y z, x y z
5 5 0, 1 0 0
</pre>

running the following command on terminal 

```bash 
Dipole_potential_map.py dip.in -p 30 30 5 -c 200 -e 5 10
```
will create file Dipoles.xsf which I visualized with Vesta
