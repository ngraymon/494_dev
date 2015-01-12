from __future__ import print_function


from MMTK import *
from PIGSHarmonicOscillatorFF           import PIGSHarmonicOscillatorForceField

from MMTK.ForceFields.ForceFieldTest    import gradientTest, forceConstantTest
from MMTK_PIGSNormalModeIntegrator      import PIGSLangevinNormalModeIntegrator
from MMTK.Environment                   import PathIntegrals
from MMTK.Trajectory                    import Trajectory, TrajectoryOutput
from MMTK.Dynamics                      import VelocityVerletIntegrator


# neil
import numpy as np
import sys


outdir  = "/scratch/$USER/"
outfile = "defaultfilename"

# handle arguments
import argparse
parser = argparse.ArgumentParser(description="Take a filename and directory in which to write the file")
parser.add_argument("-f", "--output_file_name", action="store", help="Provide a name for the file that output is written to")
parser.add_argument("-d", "--output_directory_name", action="store", help="Provide a name for the directory that the output file is written to")
parser.add_argument("-b", "--beta", type=float, action="store", help="Provide a beta value to replace the default of ((P-1)*tau)")
parser.add_argument("-t", "--tau", type=float, action="store", help="Provide a tau value to replace the default of 0.0625/Units.K")
parser_args = parser.parse_args() # prepare the argument parser

# get the input paramaters
if parser_args.output_file_name:
    outfile = parser_args.output_file_name
if parser_args.output_directory_name:
    outdir = parser_args.output_directory_name
if (parser_args.beta and parser_args.tau):
    print("You cannot provide both beta: " + str(parser_args.beta) + " and tau: " + str(parser_args.tau))
    sys.exit()



########################### CHOOSE SYSTEM PARAMETERS #########################
omega               = ((1.0 * Units.K * Units.k_B)/Units.hbar) # in THz
w_2                 = np.power(omega, 2)  # we use this often enough


P    = 129
tau  = float(0.0625/Units.K)
beta = 8.0

if parser_args.tau:
    P    = int(parser_args.tau)
    tau  = (beta / (P-1))
elif parser_args.beta:
    P    =  int(parser_args.beta)
    beta = ((P-1) * tau)




temperature         = 1.0/(Units.K*beta)             # in k-1
centroid_friction   = 0.03 / Units.ps                # in ps-1
delta_t             = 0.1 * Units.ps                 # in ps
skip                = int(10*Units.ps/delta_t)
num_steps           = int(10.0*Units.ns/delta_t)



print(  "tau                :" + str(tau)
    + "\nnumber of beads    :" + str(P)
    + "\nbeta               :" + str(beta)
    + "\ntemperature        :" + str(temperature)
    + "\ncentroid friction  :" + str(centroid_friction)
    + "\ndelta t            :" + str(delta_t)
    + "\nskip               :" + str(skip)
    + "\nomega              :" + str(omega)
    + "\nk_B                :" + str(Units.k_B)
    )

############################## SETUP UNIVERSE ################################

# create the universe, the first atom, and add the Langevin thermometer?
universe = InfiniteUniverse()
universe.atom1 = Atom('e', nbeads=P, position=Vector(0., 0., 0.), name='C1')
universe.addObject(Environment.PathIntegrals(temperature))

# create the forcefield and impose it on the universe
ff = PIGSHarmonicOscillatorForceField(universe.atomList()[0], Vector(0., 0., 0.), universe.atomList()[0].mass() * w_2)
universe.setForceField(ff)

# initialize the velocities, make the integrator the trajectory
universe.initializeVelocitiesToTemperature(temperature)
integrator = PIGSLangevinNormalModeIntegrator(universe, delta_t = delta_t, centroid_friction = centroid_friction)
trajectoryNVT = Trajectory(universe, (outdir + outfile + ".nc"), "w", "file_descriptor")

# integrate
integrator(steps = num_steps, actions = [TrajectoryOutput(trajectory = trajectoryNVT, data = ('configuration','time','velocities'), first = 0, last = None, skip = skip)] )

# collect information of energies
plotting_file = open(outdir+outfile+".plt","w")
potential_array = np.empty(len(trajectoryNVT))
semi_classical_array = np.empty((len(trajectoryNVT), 3))


# the energy estimator
def estimator_E(u_verse, i):
    def V(atom_index, bead_index): # use the norm squared, divided by three to give us the average position in 3dimensions
        return (0.5 * (u_verse.atomList()[atom_index].mass() * w_2) * (pow(u_verse.atomList()[atom_index].beadPositions()[bead_index].length(), 2) / 3.0 ))

    potential_array[i] = (V(atom_index = 0, bead_index = 0) + V(atom_index = 0, bead_index = -1)) * 0.5
    # estimator_E(u_verse, index)


# the semi classical estimator
def estimator_semi_classical(u_verse, i, dimension_index):
    # for now we only care about one particle
    atom = 0

    # define the middle bead
    M = (P-1) / 2

    # fetches one dimension of the position of a specific bead
    left    = lambda atom_index: u_verse.atomList()[atom_index].beadPositions()[M-1][dimension_index]
    middle  = lambda atom_index: u_verse.atomList()[atom_index].beadPositions()[M][dimension_index]
    right   = lambda atom_index: u_verse.atomList()[atom_index].beadPositions()[M+1][dimension_index]

    # this function 'fetches' one dimension of the momentum of an atom 
    m       = lambda atom_index: u_verse.atomList()[atom_index].mass()
    p       = lambda atom_index: (u_verse.velocities()[M][dimension_index] * m(atom_index)) / (P-1)

    # conversion factor to make units work out
    factor = 1.0 / Units.k_B
    #factor = 0

    # meat and potatoes

    semi_classical_array[i][dimension_index] = (
            (np.sqrt(2.0) / (P-1))
            * left(atom) 
            * right(atom) 
            * np.exp( factor 
                      * (( (beta * (P-1) * np.power(p(atom), 2)) / (2.0 * m(atom)))  
                        + (( m(atom) * tau * w_2 * 0.25) * ( np.power(left(atom), 2) + (2 * np.power(middle(atom), 2)) + np.power(right(atom), 2) )) 
                        )
                    ) 
            * np.cos( (p(atom) * (right(atom) - left(atom))) / Units.hbar ))
    '''

    semi_classical_array[i][dimension_index] = (
            (np.sqrt(tau / (np.pi * m(atom))))
            * left(atom) 
            * right(atom) 
            * np.exp( factor 
                      * (( (beta * (P-1) * np.power(p(atom), 2)) / (2.0 * m(atom)))  +
                        (( m(atom) * tau * w_2 * 0.25) * ( np.power(left(atom), 2) + (2 * np.power(middle(atom), 2)) + np.power(right(atom), 2) )) 
                        )
                    ) 
            * np.cos( (p(atom) * (right(atom) - left(atom))) / Units.hbar ))
    '''
    # function end


# operate over the whole trajectory
for index in xrange(len(trajectoryNVT)):
    universe.setFromTrajectory(trajectoryNVT, index)
    #estimator_E(universe, index)
    estimator_semi_classical(universe, index, 0) # x
    estimator_semi_classical(universe, index, 1) # y
    estimator_semi_classical(universe, index, 2) # z

print(str(semi_classical_array.shape))

# calculate E, std_dev, std_err
Energy  = np.mean(semi_classical_array, axis = 0)
std_dev = np.std(semi_classical_array, axis = 0, ddof = 1)
std_err = np.multiply(np.divide(std_dev, np.sqrt(len(trajectoryNVT))), 1.96)

print(str(Energy.shape))
print(str(std_dev.shape))
print(str(std_err.shape))

# write them to a file
if parser_args.tau:
    print(str(tau), file=plotting_file)
elif parser_args.beta:
    print(str(beta), file=plotting_file)
print(str(Energy[0]), file=plotting_file)
print(str(Energy[1]), file=plotting_file)
print(str(Energy[2]), file=plotting_file)
print(str(std_err[0]), file=plotting_file)
print(str(std_err[1]), file=plotting_file)
print(str(std_err[2]), file=plotting_file)

# close files
trajectoryNVT.close()
plotting_file.close()

