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

############################# IMPORT THE TRAJECTORY ##########################
def run(P = 129.0, beta = 8.0, tau = (0.0625 / Units.K), outdir = None, outfile = None):
    if (outdir == None) or (outfile == None) :
        sys.exit("outdir or outfile not provided, can't run() estimator.py")

    beta = ((P-1) * tau) if (P != 129) else float(8.0)

    outdir = "./workspace/"
    outfile = "betaho"
    for P in range(21, 201, 2):
        beta = ((P-1) * tau) if (P != 129) else float(8.0)
        # gather our information and prepare files to write to
        try:
            print(outdir + outfile + str(P) + ".nc")
            trajectory_file = Trajectory(None, outdir + outfile + str(P) + ".nc")
            universe = trajectory_file.universe
            plotting_file        = open("./workspace/"+ outfile + str(P) + ".plt","w")
        except IOError as e:
            print("ERROR could not open file: \'" + str(outfile) + "\'", file = sys.stderr)
            continue
        else:
            #semi_classical_array    = np.empty(len(trajectory_file))  # store the results in this array
            #middle_point            = (len(trajectory_file) / 2)      # calculate the 'middle point' of our data sample, 
            potential_array = np.empty(len(trajectory_file))
            # this value is ACTUALLY half the size of the data sample, it is refered to as the 'middle point' because
            # it is used as an index when claculating values from either the first or second half of semi_classical_array


            # operate over the whole trajectory
            for index in xrange(len(trajectory_file)):
                universe.setFromTrajectory(trajectory_file, index)
                # for the first half of our data we need are calculating the numerator              
                # for the second half of our data we need are calculating the denominator
                #denominator = (False) if (index < middle_point) else (True)                
                #semi_classical_array[index] = estimator_semi_classical(universe, index, denominator, P, beta, tau)  # this is where the estimation happens
                potential_array[index] = estimator_semi_classical(universe, index, False, P, beta, tau)
                #potential_array[index] = estimator_E(universe, index)
                            
            '''
            # find the mean of the numerator and denominator
            numerator_mean     = np.mean(semi_classical_array[0:middle_point:1], axis = 0)
            denominator_mean   = np.mean(semi_classical_array[middle_point::1],  axis = 0)
            print(str(numerator_mean) + " " + str(denominator_mean))

            # the actual result of our estimator
            estimator_result    = np.divide(numerator_mean, denominator_mean)

            # calculate the standard error of the numerator, denominator and the combined standard error
            numerator_std_err   = np.std(semi_classical_array[0:middle_point:1], axis = 0, ddof = 1) / np.sqrt(middle_point)
            denominator_std_err = np.std(semi_classical_array[middle_point::1], axis = 0, ddof = 1) / np.sqrt(middle_point)
            combined_std_err    = (numerator_std_err / denominator_mean) + ((numerator_mean * denominator_std_err) / np.power(denominator_mean, 2))
            '''
            estimator_result = np.mean(potential_array)
            numerator_mean=0 
            enominator_mean=0
            numerator_std_err=0
            denominator_std_err=0
            combined_std_err = np.std(potential_array, axis = 0, ddof = 1) / np.sqrt(len(trajectory_file))

            print(str(beta), file=plotting_file)
            print(str(estimator_result),    file=plotting_file)
            print(str(numerator_mean),      file=plotting_file)
            print(str(enominator_mean),     file=plotting_file)
            print(str(numerator_std_err),   file=plotting_file)
            print(str(denominator_std_err), file=plotting_file)
            print(str(combined_std_err),    file=plotting_file)
            
            '''
            # print for debugging purposes
            print(estimator_result)
            print(numerator_std_err)
            print(denominator_std_err)
            print(combined_std_err)

            # write them to a file
            print(str(beta), file=plotting_file)
            print(str(estimator_result),    file=plotting_file)
            print(str(numerator_mean),      file=plotting_file)
            print(str(enominator_mean),     file=plotting_file)
            print(str(numerator_std_err),   file=plotting_file)
            print(str(denominator_std_err), file=plotting_file)
            print(str(combined_std_err),    file=plotting_file)

            # carry out the estimation
            estimator(outfile, trajectory_file, plotting_file)
            '''
                

        # done


# the energy estimator
def estimator_E(u_verse, i):
    omega = ((1.0 * Units.K * Units.k_B)/Units.hbar) # in THz
    w_2   = np.power(omega, 2)  # we use this often enough

    def V(atom_index, bead_index): # use the norm squared, divided by three to give us the average position in 3dimensions
        return (0.5 * (u_verse.atomList()[atom_index].mass() * w_2) * (pow(u_verse.atomList()[atom_index].beadPositions()[bead_index].length(), 2) / 3.0 ))

    return (V(atom_index = 0, bead_index = 0) + V(atom_index = 0, bead_index = -1)) * 0.5
   # end

######################### SEMI CLASSICAL ESTIMATOR ##########################
def estimator_semi_classical(u_verse, i, denominator, P, beta, tau):

    # we use these two parameters in our equation
    omega = ((1.0 * Units.K * Units.k_B)/Units.hbar) # in THz
    w_2   = np.power(omega, 2)  # we use this often enough

    # for now we only care about one particle
    atom = 0

    # define the middle bead
    M = (P-1) / 2

    # fetches one dimension of the position of a specific bead
    left    = lambda dimension_index: u_verse.atomList()[0].beadPositions()[M-1][dimension_index]
    middle  = lambda dimension_index: u_verse.atomList()[0].beadPositions()[M][dimension_index]
    right   = lambda dimension_index: u_verse.atomList()[0].beadPositions()[M+1][dimension_index]

    # this function 'fetches' one dimension of the momentum of an atom 
    m       = u_verse.atomList()[0].mass()
    p       = lambda dimension_index: u_verse.velocities()[M][dimension_index] * m / (P-1)

    # the q's
    q       = lambda dimension_index, denominator: 1 if (denominator is True) else (left(dimension_index) * right(dimension_index))

    # meat and potatoes, equation version three
    fn = lambda d, denominator: (
        (np.sqrt(2.0) / (P-1))
        * q(d, denominator)
        * np.exp(  (m * Units.k_B * np.power((right(d) - left(d)), 2) / (4.0 * Units.hbar * Units.hbar * tau))
                +  ((beta * (P-1) * np.power(p(d),2)) / (2.0 * Units.k_B * m))
                ) 
        * np.cos( (p(d) * (left(d) - right(d))) / Units.hbar ))

    return (fn(0, denominator) + fn(1, denominator) + fn(2, denominator)) / 3

######################### RUNNING ESTIMATOR DIRECTLY ##########################
if(__name__ == "__main__"):
    import argparse
    parser = argparse.ArgumentParser(prog='estimator.py', description='Preforms black magic rituals developed by Dimitry', 
    usage = '%(prog)s [-h] [-d DIRECTORYNAME] [-f FILENAME] [-p BEADSVALUE' 
        + '\n%(prog)s [--help] [--outputdirectory DIRECTORYNAME] [--outputfile FILENAME] [--beads BEADSVALUE]',
    epilog='Extra information can be displayed here, such as a riddle: "What walks on four legs in the morning, two legs in the afternoon, and three legs in the evening" ')
    parser.add_argument("-d", "--outputdirectory",  action="store", help="Provide a name for the directory that the output file is written to")
    parser.add_argument("-f", "--outputfile",       action="store", help="Provide a name for the file that output is written to")
    parser.add_argument("-p", "--beads", type=int,  action="store", help="Provide the number of beads to replace the default of 129")
    parser_args = parser.parse_args() # prepare the argument parser

    # get the input paramaters
    if (parser_args.outputfile == False) or (parser_args.beads == False):
        print("You must provide a trajectory filename and a number of beads")
        sys.exit()

    directoryname = (parser_args.outputdirectory) if (parser_args.outputdirectory) else ("./workspace/")

    # run the job
    run(P = int(parser_args.outputfile[-2] + parser_args.outputfile[-1]) , outdir = directoryname, outfile = parser_args.outputfile)




