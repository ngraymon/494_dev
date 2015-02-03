# batch job script for estimator.py
 
# MMTK to read trajectory
from MMTK import *
#from MMTK                               import Units
from PIGSHarmonicOscillatorFF           import PIGSHarmonicOscillatorForceField
from MMTK.ForceFields.ForceFieldTest    import gradientTest, forceConstantTest
from MMTK.Dynamics                      import VelocityVerletIntegrator

from MMTK.Trajectory                    import Trajectory, TrajectoryOutput
from MMTK.Environment                   import PathIntegrals
# used for running unix commands, specifically to submit jobs to the server
import subprocess
import numpy as np

# handle arguments
import argparse
parser = argparse.ArgumentParser(prog='test.py', description='run classical trajectories', 
    usage = '%(prog)s [-h] [-d DIRECTORYNAME] [-f FILENAME] [-e]' 
        + '\n%(prog)s [--help] [--inputdirectory DIRECTORYNAME] [--inputfile FILENAME] [--estimate] ',
    epilog='Extra information can be displayed here, such as a riddle: "What walks on four legs in the morning, two legs in the afternoon, and three legs in the evening" ')
parser.add_argument("-d", "--inputdirectory",   action="store",      help="Provide a name for the directory that the output file is written to")
parser.add_argument("-f", "--inputfile",        action="store",      help="Provide a name for the file that output is written to")
parser.add_argument("-e", "--estimate",         action="store_true", help="Set this flag to run the estimator after the trajectory has been calculated")
parser_args = parser.parse_args() # prepare the argument parser

input_directory = parser_args.inputdirectory if parser_args.inputdirectory else "/workspace/"
input_filename  = parser_args.inputfile      if parser_args.inputfile      else "defaultfilename"

output_directory = "/scratch/ngraymon/"
output_filename  = "classical" + str(input_filename)


# try to open tracjectory file, and then prepare to launch classical's
try:
    trajectory_file = Trajectory(None, input_directory + input_filename + ".nc")
    universe = trajectory_file.universe
    universe.setFromTrajectory(trajectory_file, 0)
    P = universe.atomList()[0].numberOfBeads()      # retrive the number of beads
    M = (P-1) / 2                                   # redefine the middle bead
    #plotting_file        = open(outfile +".plt","w") if (outdir == "/workspace/") else open(outdir + outfile +".plt","w")
except IOError as e:
    sys.exit("ERROR could not open file: \'" + str(input_filename) + "\'")
else:
    for index in xrange(len(trajectory_file)):
        if (index == 10000):
            break
        universe.setFromTrajectory(trajectory_file, index)

        leftmiddlebead_position  = universe.atomList()[0].beadPositions()[M-1]   # M-1 bead's  x,y,z values, object is a Vector
        rightmiddlebead_position = universe.atomList()[0].beadPositions()[M+1]   # M+1 bead's  x,y,z values, object is a Vector
        initial_position = universe.atomList()[0].beadPositions()[M]             # M bead's  x,y,z values,     object is a Vector

        ############################## CREATE TRAJECTORY
        omega               = ((1.0 * Units.K * Units.k_B)/Units.hbar)
        w_2                 = np.power(omega, 2)
        delta_t             = 0.1  * Units.ps
        num_steps           = int(1000)

        
        new_universe = InfiniteUniverse()
        new_universe.atom1 = Atom('e', position=initial_position, name='C1')
        ff = PIGSHarmonicOscillatorForceField(atom=new_universe.atomList()[0], center=Vector(0., 0., 0.), force_constant=(new_universe.atomList()[0].mass() * w_2), UsingPIGS = False)
        new_universe.setForceField(ff)
        # we artifically set the momentum to zero, p=v*m, so the velocity is also zero
        new_universe.setVelocities(velocities = Vector(0., 0., 0.))
        integrator = VelocityVerletIntegrator(new_universe, delta_t = delta_t)
        trajectoryNVE = Trajectory(new_universe, (output_directory + output_filename + "_" + str(index) + ".nc"), "w", "file_descriptor")

        # integrate
        integrator(steps = num_steps, actions = [TrajectoryOutput(trajectory = trajectoryNVE, data = ('configuration','time','velocities'), first = 0, last = None)] )



        '''
        classic_integrator.integrate(leftmiddlebead_position, rightmiddlebead_position, initial_position, initial_momentum, estimate = parser_args.estimate)


        unix_cmd = 'qsub -N "tauho{0}" -v "FILENAME=ho{1}","TAU={2}" classical_integrator.sh'.format(quote())
        p = subprocess.Popen(unix_cmd, shell=True, stderr=subprocess.PIPE)
        while True:
            out = p.stderr.read(1)
            if out == '' and p.poll() != None:
                break
            if out != '':
                sys.stdout.write(out)
                sys.stdout.flush()
        '''


    # for loop


