from MMTK import *
from PIGSHarmonicOscillatorFF           import PIGSHarmonicOscillatorForceField

from MMTK.ForceFields.ForceFieldTest    import gradientTest, forceConstantTest
from MMTK_PIGSNormalModeIntegrator      import PIGSLangevinNormalModeIntegrator
from MMTK.Environment                   import PathIntegrals
from MMTK.Trajectory                    import Trajectory, TrajectoryOutput


# neil
import estimator
import numpy as np
import sys


# handle arguments
import argparse
parser = argparse.ArgumentParser(prog='calculate_trajectories.py', description='calulate a PIGS trajectory', 
    usage = '%(prog)s [-h] [-d DIRECTORYNAME] [-f FILENAME] [-e] [-b BETAVALUE | -t TAUVALUE]' 
        + '\n%(prog)s [--help] [--outputdirectory DIRECTORYNAME] [--outputfile FILENAME] [--estimate] [--beta BETAVALUE | --tau TAUVALUE]',
    epilog='Extra information can be displayed here, such as a riddle: "What walks on four legs in the morning, two legs in the afternoon, and three legs in the evening" ')
parser.add_argument("-d", "--outputdirectory",  action="store",      help="Provide a name for the directory that the output file is written to")
parser.add_argument("-f", "--outputfile",       action="store",      help="Provide a name for the file that output is written to")
parser.add_argument("-e", "--estimate",         action="store_true", help="Set this flag to run the estimator after the trajectory has been calculated")
parser.add_argument("-b", "--beta", type=float, action="store",      help="Provide a beta value to replace the default of ((P-1)*tau)")
parser.add_argument("-t", "--tau", type=float,  action="store",      help="Provide a tau value to replace the default of 0.0625/Units.K")
parser_args = parser.parse_args() # prepare the argument parser

# get the input paramaters
outdir  = (parser_args.outputdirectory) if (parser_args.outputdirectory) else ("")
outfile = (parser_args.outputfile)      if (parser_args.outputfile)      else ("defaultfilename")


if (parser_args.beta and parser_args.tau):
    print("You cannot provide both beta: " + str(parser_args.beta) + " and tau: " + str(parser_args.tau))
    sys.exit()

########################### CHOOSE SYSTEM PARAMETERS #########################
# P  = 512 # tau  = float(0.015625/Units.K) # beta = 8.0

P    = 129                                      # P is the # of beads
tau  = float(0.0625/Units.K)  # inverse kelvin  # tau is the length of the springs
beta = 8.0                    # inverse kelvin  # beta is the length of the whole path

if(parser_args.tau):
    P = int(parser_args.tau)
    tau = (beta / (P-1))
elif(parser_args.beta):
    P = int(parser_args.beta)
    beta = ((P-1) * tau)

omega               = ((1.0 * Units.K * Units.k_B)/Units.hbar)  # in THz
w_2                 = np.power(omega, 2)
temperature         = 1.0  / beta                               # in K
centroid_friction   = 0.03 / Units.ps                           # in ps-1
delta_t             = 1.0 * Units.ps                            # in ps        # this is the distance in time between trajectory calculations
skip                = int(10.0 * Units.ps/delta_t)              # has no units # this is the distance in time between sampling points, this is chosen so that our data is decorrelated
# we need to give integrator skip step size in units of delta_t, so giving the integrator 5 means the skip time is the delta t * 5
num_steps           = int(1000.0 * (Units.ns/delta_t))    # number of steps is always x * (10-9 / 10-12)  which is equivalent to (x * 1000)

print(  "tau                :" + str(tau)
    + "\nnumber of beads    :" + str(P)
    + "\nbeta               :" + str(beta)
    + "\ntemperature        :" + str(temperature)
    + "\ncentroid friction  :" + str(centroid_friction)
    + "\ndelta t            :" + str(delta_t)
    + "\nskip size          :" + str(skip)
    + "\nNumber of steps    :" + str(num_steps)
    + "\nomega              :" + str(omega)
    )

############################## SETUP UNIVERSE ################################

# create the universe, the first atom, and add the Langevin thermometer?
universe = InfiniteUniverse()
universe.atom1 = Atom('e', nbeads=P, position=Vector(0., 0., 0.), name='C1')
universe.addObject(Environment.PathIntegrals(temperature))

# create the forcefield and impose it on the universe
ff = PIGSHarmonicOscillatorForceField(atom=universe.atomList()[0], center=Vector(0., 0., 0.), force_constant=(universe.atomList()[0].mass() * w_2), UsingPIGS = True)
universe.setForceField(ff)

# initialize the velocities, make the integrator the trajectory
universe.initializeVelocitiesToTemperature(temperature)
integrator = PIGSLangevinNormalModeIntegrator(universe, delta_t = delta_t, centroid_friction = centroid_friction)
trajectoryNVT = Trajectory(universe, (outdir + outfile + ".nc"), "w", "file_descriptor")

############################## CREATE TRAJECTORY ################################

# integrate
integrator(steps = num_steps, actions = [TrajectoryOutput(trajectory = trajectoryNVT, data = ('configuration','time','velocities'), first = 0, last = None, skip = skip)] )

# if we passed in the estimate flag then run the estimator
if parser_args.estimate:
    estimator.run(P, beta, tau, outdir, outfile)






