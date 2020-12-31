# Dipole-system-potential
Creates a map of potential on a 3D grid caused by a 2 dimensional periodic system of dipoles and saves it in a .xsf file.
usage: Dip_potential_newVersion2.py [-h] [-p POINT_GRID POINT_GRID POINT_GRID]
                                    [-c CUTOFF] [-o OUTFILE]
                                    [-e EVAL_HEIGHT EVAL_HEIGHT] [-l]
                                    [-k--columbic_constat K]
                                    geo_input

Draw potential of a dipole system

positional arguments:
  geo_input             file containing possitions and moments of all dipoles

optional arguments:
  -h, --help            show this help message and exit
  -p POINT_GRID POINT_GRID POINT_GRID, --point_grid POINT_GRID POINT_GRID POINT_GRID
                        number of points along x y z lines in point grid
  -c CUTOFF, --cut_off CUTOFF
                        cutoff of the dipoles, if zero or not selected system
                        will be treated as non-periodic
  -o OUTFILE, --outfile OUTFILE
                        where to store result the final xsf geometry default
                        is Dipoles.xsf
  -e EVAL_HEIGHT EVAL_HEIGHT, --evaluation_height EVAL_HEIGHT EVAL_HEIGHT
                        the maximum and minimum height of the point grid
  -l, --long_range      add an infinite system correction
  -k, --columbic_constat K
                        default is 14.3996

Libraries:
ProbeParticleModel https://github.com/ProkopHapala/ProbeParticleModel
                   for creating the .xsf file
                   add the path to the library in script as ppafm right at the top
                   I found a similiar function in ase but did not get it to work 😪️

Returns:
outfile - .xsf file with the potential grid
