#
# pigs_plotter.py
#
"""Plot the Energy vs (Beta or Tau)"""
from __future__ import print_function

# SciPy packages
import MMTK.Units as u
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
matplotlib.use('Agg')
from pylab import *
import os




# get the input paramaters
try:
    input_file = open("plotdata.batchjoboutput")
except OSError as err:
    print("Could not find input file: plotdata.batchjoboutput\n Catastrophic Error! Please find the nearest table and hide underneath it\n")
    raise OSError(err)
else:
    with input_file:
        start       = float(input_file.readline().strip())
        stop        = float(input_file.readline().strip())
        step        = float(input_file.readline().strip())
        parameter   = input_file.readline().strip()
        input_file.close()


# Prepare the plot
figure_object = figure(figsize=(18, 10), facecolor='w', edgecolor='k') # figsize = (width, heigh)



beta_tau            = []
avg_estimated_val   = []
avg_err             = []

# get the data and draw the plot
#print("Beads|   \\beta   |   est_val   | numerator mean | srd_err numerator | denominator mean | std_err denominator |  combined std_err \n" +
      #"-----------------------------------------------------------")

color_array = ['b', 'g', 'r', 'c', 'm', 'y', '#FF9666']
errorbar_color = ['c', 'y', '#00FFCC', '#FF9933', '#996633', '#0066CC', '#FF9666']

def ft(c_real, c_imag, delta_t):
    size = 10000
    weta_min = 0.0
    weta_max = 1.0
    t_j = np.multiply(np.arange(len(c_real)), delta_t)

    #w_trans  = np.transpose(np.array([(x * (weta_max - weta_min) / size) for x in range(size)]))
    #test     = np.sum((np.cos(np.multiply(w_trans, t_j)) * r) - (np.sin(np.multiply(w_trans, t_j)) * i), axis=0)

    w_value  = np.array([(x * (weta_max - weta_min) / size) for x in range(size)])
    fourier = lambda w, r, i: np.sum((np.cos(w * t_j) * r) - (np.sin(w * t_j) * i))

    spectrum = np.multiply(np.array([fourier(w, c_real, c_imag) for w in w_value]), 2*delta_t)

    #print(test.shape)
    #print(spectrum.shape)
    #print(test)
    #print(spectrum)
    #sample = range(0, 1001, 2)

    return w_value, spectrum#, w_half, spectrum_half


for index, num in enumerate([1025]): #513
#for index, num in enumerate([3,5,9,17,33,65,129,257,513,1025]):
#for index, num in enumerate(np.arange(start, stop+step, step, dtype = int)):
    #print ("index " + str(index) + " num {0:.2f}" ).format(num)
    #template_string = 'workspace/betaho{0:.2f}.plt'
    #filename = template_string.format(num)
    if(num == 0):
        continue

    filename = "./workspace/betaho" + str(num) + ".real_avg"   
    #filename = "workspace/" + "betaho" + str(num) + ".plt"
    try:
        if (os.stat(filename).st_size > 0):
            #bt, val, mean_top, mean_bot, err_top, err_bot, err_total = np.loadtxt(filename, unpack=True)
            tau = (0.0625 / u.K)
            beta = ((num-1) * tau)

            #for temp, step in enumerate([10, 980, 2505, 5043, 8922, 9500]):

            #sub = figure_object.add_subplot(3,3, temp+1, title ="$\\langle C(t)\\rangle$ vs " + "$t$ \n I am" + str(step), xlabel = "$time$", ylabel = "$\\langle C(t)\\rangle$")
            #sub = figure_object.add_subplot(111, title ="$\\langle C(t)\\rangle$ vs " + "$t$ \n I am" + str(step), xlabel = "$time$", ylabel = "$\\langle C(t)\\rangle$")
            #sub.axhline(y = ( u.hbar * u.hbar  / ( 2.0 * 0.000548579909 * u.K * u.k_B)) , xmin=0, xmax=1, linewidth=2, color = 'k')
            #sub.axhline(y = 0.0, xmin=0, xmax=1, linewidth=2, color = 'k')



            real_Corr = np.loadtxt("./workspace/betaho" + str(num) + ".real_avg")
            imag_Corr = np.loadtxt("./workspace/betaho" + str(num) + ".imag_avg")
            #long_real_Corr = np.loadtxt("./workspace/betaho" + str(1025) + ".real_avg")
            #long_imag_Corr = np.loadtxt("./workspace/betaho" + str(1025) + ".imag_avg")
            real_std_err              = np.loadtxt("./workspace/betaho" + str(num) + ".real_std_err") / np.sqrt(1001)
            imag_std_err              = np.loadtxt("./workspace/betaho" + str(num) + ".imag_std_err") / np.sqrt(1001)
            #index = temp
            #real_Corr[real_Corr >  1600] = 0
            #real_Corr[real_Corr < 0] = 0
            #std_err[std_err >   1600] = 0
            #std_err[std_err <  0] = 0

            omega    = ((1.0 * u.K * u.k_B)/u.hbar)  # in THz
            constant = u.hbar / (2 * 0.000548579909 * omega)
            delta_t  = 0.1  * u.ps



            real_compare        = [cos(omega * time * delta_t ) * constant for time in range(0, 1001)]
            imaginary_compare   = [-sin(omega * time * delta_t ) * constant for time in range(0, 1001)]
            sample = range(0, 1001, 4)

            #sub.errorbar(sample, real_Corr[sample], label = "$\\langle C(t).real\\rangle$" + " \n$P = " + str(num) + "$\n$\\beta = " + str(beta) + "$\n$\\tau = " + str(tau) + "$", 
            #    yerr=real_std_err[sample], fmt='-o', color='k', ecolor=errorbar_color[index], elinewidth=1, linewidth=4, linestyle="-",clip_on=False)

            #sub.errorbar(sample, imag_Corr[sample], label = "$\\langle C(t).imag\\rangle$" + " \n$P = " + str(num) + "$\n$\\beta = " + str(beta) + "$\n$\\tau = " + str(tau) + "$", 
            #    yerr=imag_std_err[sample], fmt='o', color=color_array[index+1], ecolor=errorbar_color[index+1], elinewidth=1, linewidth=4, linestyle="-", clip_on=False)

            #sub.plot(range(0, 1001), real_compare, label = "$\\langle C(t).real\\rangle$" + " \n$P = " + str(num) + "$\n$\\beta = " + str(beta) + "$\n$\\tau = " + str(tau) + "$", color='k')
            #sub.plot(range(0, 1001), imaginary_compare, label = "$\\langle C(t).imag\\rangle$" + " \n$P = " + str(num) + "$\n$\\beta = " + str(beta) + "$\n$\\tau = " + str(tau) + "$", color=color_array[index+1])


            #############  fourier stuff  ###############
            fx, fy = ft(real_Corr, imag_Corr, delta_t)
            scaling_factor = 2
            half = range(0, 1000/scaling_factor, scaling_factor)
            fx_half, fy_half = ft(real_Corr[half], imag_Corr[half], 2*delta_t)

            fx = fx * u.hbar / u.k_B
            fx_half = fx_half * u.hbar / u.k_B

            fy = fy / fy[0]    
            fy_half = fy_half / fy_half[0]      
            
            sub = figure_object.add_subplot(111, title ="Harmonic Oscillator spectrum", xlabel = "$\\omega$  $(THz)$", ylabel = "$I(\\omega)$")
            #sub2 = figure_object.add_subplot(212, title ="Harmonic Oscillator spectrum", xlabel = "$\\omega (THz)$", ylabel = "$I(\\omega)$")
            sub.axvline(x = 1, ymin=0, linewidth=2, color = 'k') # 0.131 THz
            plabel = "Spectrum" + " \n$P$ = " + str(num) + "\n$\\beta$ = $" + str(beta) + " \\mathrm{K^{-1}}$\n$\\tau$  = $" + str(tau) + " \\mathrm{K^{-1}}$"
            sub.plot(fx, fy, label = plabel + "\n$\Delta t$ = $" + str(0.1) + " ps$", color='g')
            sub.plot(fx_half, fy_half, label = plabel + "\n$\Delta t$ = $" + str(0.2) + " ps$", color='r')

            xticks, xticklabels = matplotlib.pyplot.xticks()
            yticks, yticklabels = matplotlib.pyplot.yticks()

            from scipy.signal import argrelextrema # find the maximums and label them
            maximum_point = argrelextrema(fy, np.greater)[0][1]
            y_shift = yticks[-1] / 13
            x_shift = xticks[-1] / 60
            sub.annotate("x = {0:f}".format(fx[maximum_point]), xy=(fx[maximum_point], fy[maximum_point]), xycoords='data', xytext=(fx[maximum_point] + x_shift, fy[maximum_point] + y_shift), arrowprops=dict(arrowstyle="->"))
            
            maximum_point = argrelextrema(fy_half, np.greater)[0][0]
            y_shift = yticks[-1] / 20
            x_shift = xticks[-1] / 60
            sub.annotate("x_half = {0:f}".format(fx_half[maximum_point]), xy=(fx_half[maximum_point], fy_half[maximum_point]), xycoords='data', xytext=(fx_half[maximum_point] + x_shift, fy_half[maximum_point] + y_shift), arrowprops=dict(arrowstyle="->"))
 
            #############  fourier stuff  ###############


            sub.legend()

            print("Plotted trajectory with P value of: " + str(num))
            #beta_tau.append(bt)
            #avg_estimated_val.append(val)
            #avg_err.append(err_total)
            #output_string = "  {0:} |  {1:.3f} | {2: .4e} | {3: .4e} | {4: .4e} | {5: .4e} | {6: .4e} | {7: .4e} "
            #print(output_string.format(num, bt, val, mean_top, mean_bot, err_top, err_bot, err_total))
        else:
            print("Could not open file: " + filename)
            
    except OSError as err:
        print("Could not find file: " + filename)
        #raise OSError(err)
    




#sub.plot(beta_tau, avg_estimated_val, label = str(parameter) + " averaged over all dimensions", color = 'r')
#sub.plot(beta_tau, avg_estimated_val, label = str(parameter) + " averaged over all dimensions", color = 'r')
#sub.plot(beta_tau, avg_estimated_val, label = str(parameter) + " averaged over all dimensions", color = 'r')
#sub.errorbar(beta_tau, avg_estimated_val, label = "$\\beta$" + " averaged over all dimensions", yerr=avg_err, xerr=None, fmt='-o', color='g', ecolor='r', elinewidth=1, linewidth=1, linestyle="-", marker="o", markersize=2, clip_on=False)
#sub.errorbar(beta_tau, y_dimension, label = str(parameter) + " y dimension", yerr=y_err, xerr=None, fmt='-o', ecolor='r', elinewidth=1, linewidth=1, linestyle="-", marker=" ")
#sub.errorbar(beta_tau, z_dimension, label = str(parameter) + " z dimension", yerr=z_err, xerr=None, fmt='-o', ecolor='b', elinewidth=1, linewidth=1, linestyle="-", marker=" ")
#sub.legend(prop=matplotlib.font_manager.FontProperties(size=24))


xticks, xticklabels = matplotlib.pyplot.xticks()
yticks, yticklabels = matplotlib.pyplot.yticks()
xmin = (3*xticks[0] - xticks[1])/2.
xmax = (3*xticks[-1] - xticks[-2])/2.
ymin = (3*yticks[0] - yticks[1])/2.
ymax = (3*yticks[-1] - yticks[-2])/2.
matplotlib.pyplot.xlim(xmin, xmax)
matplotlib.pyplot.xticks(xticks)
matplotlib.pyplot.ylim(ymin, ymax)
matplotlib.pyplot.yticks(yticks)
matplotlib.pyplot.tight_layout(pad=0.0, h_pad=None, w_pad=None, rect=None)


# save plot
pp = PdfPages("p" + str(parameter) + ".pdf")
pp.savefig(figure_object, )
pp.close()
#figure_object.savefig("p" + str(parameter) + ".png", dpi = 300, bbox_inches = 'tight')
#close(figure_object)






