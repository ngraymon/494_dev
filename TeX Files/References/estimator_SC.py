from __future__ import print_function

from PIGSHarmonicOscillatorFF           import PIGSHarmonicOscillatorForceField
from MMTK.ForceFields.ForceFieldTest    import gradientTest, forceConstantTest
from MMTK_PIGSNormalModeIntegrator      import PIGSLangevinNormalModeIntegrator
from MMTK.Environment                   import PathIntegrals
from MMTK.Trajectory                    import Trajectory, TrajectoryOutput
from MMTK.Dynamics                      import VelocityVerletIntegrator
from MMTK                               import Units

# neil
import numpy as np
import warnings
import sys

############################ CLASSICAL ESTIMATOR ############################
def estimator_classical(u_verse, t, P, beta, tau, m, QA, QM, QB):   
    atom = 0
    p       = lambda dimension_index:    u_verse.velocities()[0][dimension_index] * m
    q       = lambda dimension_index:    u_verse.atomList()[0].position()[dimension_index]
    qA      = lambda dimension_index:    QA[dimension_index]
    qM      = lambda dimension_index:    QM[dimension_index]
    qB      = lambda dimension_index:    QB[dimension_index]
    # hamiltonian_constant
    S_pq    = lambda dimension_index:    0.5 * (p(dimension_index) * q(dimension_index)) #- (pM(dimension_index) * qM(dimension_index)) )

    fn = lambda d: (qA(d) * qB(d)
        * np.exp(1j * p(d) * (qA(d) - q(d)) / Units.hbar)
        * np.exp(  (-m * Units.k_B / (2 * Units.hbar * Units.hbar * tau))
                *  ( np.power( qA(d) - q(d), 2) - np.power( qA(d) - qM(d), 2) )
                )
        * np.exp(1j * S_pq(d) / Units.hbar) 
       ) # fn

    return np.mean([fn(0), fn(1), fn(2)], dtype = complex) #, dtype=np.float64

'''
with warnings.catch_warnings(record=True) as w:
    # Cause all warnings to always be triggered.
    warnings.simplefilter("always")   
    if len(w) > 0:
        print(w)
        print("p ~ " + str(p(0)) + " " + str(p(1)) + " " + str(p(2)))
        print("q ~ " + str(q(0)) + " " + str(q(1)) + " " + str(q(2)))
        print("1j * S / hbar ~ " + str(1j * S_pq(0) / Units.hbar))
        print("1j * S / hbar ~ " + str(1j * S_pq(1) / Units.hbar))
        print("1j * S / hbar ~ " + str(1j * S_pq(2) / Units.hbar))
    try:  
        return np.mean([fn(0), fn(1), fn(2)]) #, dtype=np.float64
    except(RuntimeError, OverflowError) as exception_instance:
        print(w)
        print(type(exception_instance))
        print(exception_instance.args)
        print(exception_instance)

    return np.mean([fn(0), fn(1), fn(2)]) #, dtype=np.float64
'''

''' third version of classical function
# meat and potatoes, equation version three
   fn = lambda d: (
        (np.sqrt(2.0) / (P-1))
        * qA(d) * qB(d)
        * np.exp( 1j * pM(d) * (qA(d) - qM(d)) / Units.hbar)
        * np.exp(  (-m * Units.k_B/ (4 * Units.hbar * Units.hbar * tau))
                *  ( np.power( qA(d) - qM(d) , 2) - np.power( qA(d) - qM(d), 2) )
                )
        * np.exp(1j * S_pq(d) / Units.hbar) 
       ) # fn

    return (fn(0) + fn(1) + fn(2)) / 3
'''
''' second version of classical function
# meat and potatoes, equation version three
   fn = lambda d: (
        (np.sqrt(2.0) / (P-1))
        * qA(d) * qB(d)
        * np.exp( 1j * pM(d) * (qA(d) - qM(d)) / Units.hbar)
        * np.exp(  (-m / (4 * Units.hbar * Units.hbar * tau))
                *  ( np.power( qA(d) - qM(d) , 2) - np.power( qA(d) - qM(d), 2) )
                )
        * np.exp(1j * S_pq(d) / Units.hbar) 
       ) # fn

    return (fn(0) + fn(1) + fn(2)) / 3
'''


''' First version of classical function
# meat and potatoes, equation version three
fn = lambda d: (
    (np.sqrt(2.0) / (P-1))
    * qA(d) * qB(d)
    * np.exp(  (-2 * tau * E0)
            +  (1j * S_pq(d) / Units.hbar) 
            +  (    (m / (2 * Units.hbar * Units.hbar * tau))  
                *   (np.power(qA(d) - q(d), 2) - np.power(qB(d) - qM(d), 2)) 
                )
            )
    * np.exp(  1j * ((p(d) * (qA(d) - q(d)))  -  (pM(d) * (qB(d) - qM(d)))) / Units.hbar )
   ) # fn
'''
############################ CLASSICAL ESTIMATOR ############################

#[3,5,9,17,33,65,129,257,513,1025]
P = 1025 #17
M = P/2
tau = (0.0625 / Units.K)
beta = ((P-1) * tau)
number_of_classical_trajectories = 10000
number_of_classical_steps        = 1000 + 1

initial_traj     = Trajectory(None, "./workspace/betaho" + str(P) + ".nc")
initial_universe = initial_traj.universe

real_array      = np.zeros((number_of_classical_trajectories, number_of_classical_steps))
imaginary_array = np.zeros((number_of_classical_trajectories, number_of_classical_steps))
sval            = np.zeros(number_of_classical_steps, dtype = complex)

left    = np.zeros((number_of_classical_trajectories, 3))
middle  = np.zeros((number_of_classical_trajectories, 3))
right   = np.zeros((number_of_classical_trajectories, 3))


for step in range(number_of_classical_trajectories):
    initial_universe.setFromTrajectory(initial_traj, step)
    left[step]      = initial_universe.atomList()[0].beadPositions()[M-1] 
    middle[step]    = initial_universe.atomList()[0].beadPositions()[M] 
    right[step]     = initial_universe.atomList()[0].beadPositions()[M+1]
del initial_traj
del initial_universe
print("Collected initial values, P = " + str(P))

#for temp, step in enumerate([10, 980, 2505, 5043, 8922, 9500]):
for step in range(number_of_classical_trajectories):
    if (step in range(500, 10001, 500)):
        print(str(step))

    try:
        trajectory_file = Trajectory(None, "./workspace/classicalbetaho" + str(P) + "_" + str(step) + ".nc")
    except:
        print("./workspace/classicalbetaho" + str(P) + "_" + str(step) + ".nc")
        continue
    else:
        universe = trajectory_file.universe
        
        #print("Current classical trajectory " + str(step))
        for index in xrange(number_of_classical_steps):
            universe.setFromTrajectory(trajectory_file, index)
            
            complex_number = estimator_classical(universe, 
                                                    index, P, beta, tau,
                                                    universe.atomList()[0].mass(),
                                                    left[step], 
                                                    middle[step], 
                                                    right[step])  # this is where the estimation happens

            #real_array[step, index]         = complex_number.imag
            #imaginary_array[step, index]    = complex_number.imag
            
            real_array[step, index]         = complex_number.real
            imaginary_array[step, index]    = complex_number.imag
            #sval[index]         = left[step][0] - middle[step][0]
            #sval[index]         = np.exp(1j * universe.velocities()[0][0] * universe.atomList()[0].mass() * (left[step][0] - universe.atomList()[0].position()[0]) / Units.hbar)
            #sval[index]         = np.exp(  (-universe.atomList()[0].mass() / (2 * Units.hbar * Units.hbar * tau))
            #                    *  ( np.power( left[step][0] - universe.atomList()[0].position()[0] , 2) - np.power( left[step][0] - middle[step][0], 2) )
            #                    )
            #sval[index]         = 1j * 0.5 * universe.velocities()[0][0] * universe.atomList()[0].mass() * universe.atomList()[0].position()[0] / Units.hbar




    #print("Saving data to file")
    #np.savetxt("./" + str(P) + "_" + str(step) + ".sval", sval)

    del universe
    del trajectory_file
    # end of for loop



print("Calculating averaged C(t) real values")
averaged_real  = np.mean(real_array, axis = 0)

print("Calculating averaged C(t) imaginary values")
averaged_imag  = np.mean(imaginary_array, axis = 0)

print("Calculating standard error of the real values")
real_std_err        = np.std(real_array, axis = 0, ddof = 1) / np.sqrt(number_of_classical_steps)

print("Calculating standard error of the imaginary values")
imag_std_err        = np.std(imaginary_array, axis = 0, ddof = 1) / np.sqrt(number_of_classical_steps)


print("Saving data to file")
#for temp, step in enumerate([10, 980, 2505, 5043, 8922, 9500]):
#np.savetxt("./workspace/betaho" + str(P) + str(step) + ".real_val", real_array)
#np.savetxt("./workspace/betaho" + str(P) + str(step) + ".imag_val", real_array)
np.savetxt("./workspace/betaho" + str(P) + ".real_avg", averaged_real)
np.savetxt("./workspace/betaho" + str(P) + ".imag_avg", averaged_imag)
np.savetxt("./workspace/betaho" + str(P) + ".real_std_err", real_std_err)
np.savetxt("./workspace/betaho" + str(P) + ".imag_std_err", imag_std_err)



