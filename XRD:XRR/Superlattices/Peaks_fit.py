"""
author: Annabella Drewanowski
date: 30.6.22

fit peaks in XRD signals of superlattices
1. data path + file & change the directory path to write the data into csv
2. plot to see the peaks by eye and then write them in params list
3. by running again the peaks should be fitted

"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import csv

data_path = '/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/Analysis/XRD:XRR/SCO/'
file = '22167_SL(CO)L/22167_SL(CO)L_002_XRD.csv'
values_header = ['two_theta', 'yobs', 'ycal', 'bkg', 'diff']
df = pd.read_csv(data_path+file, names = values_header, header = 0)


def gaussian(x, A1, sigma1, mu1): # A2, sigma2, mu2
    return np.abs(A1)*np.exp(-(x-np.abs(mu1))**2/(2*np.abs(sigma1)**2)) + 2 #+ np.abs(A2)*np.exp(-(x-np.abs(mu2))**2/(2*np.abs(sigma2)**2))

def gaussian_fit(x, y, mean, sigma):
    A_1 = 3E2
    popt, pcov = curve_fit(gaussian, x, y, p0=[A_1, sigma, mean])
    return gaussian(x, *popt), popt[2], pcov[2, 2]


"""Data"""
# plot original diffraction spectrum to make estimates for the gaussian fit
x = df.two_theta
y = df.yobs
plt.plot(x, y)        # uncomment once the estimates are made
plt.yscale('log')     # uncomment """
plt.show()            # uncomment """

"""Plot"""
# plug in the parameters for gaussian fit of the peaks: two_theta_max, two_theta_min, mean, sigma
params = [[44.45, 43.40, 44.15, 0.1], [45.66, 44.69, 45.3, 0.1], [46.9,46.1, 46.5, 0.1]]
# empty list for writing csv file with peak position and std deviation
mu = []
mu_dev = []

plt.figure(figsize=(10,8))
plt.yscale('log')
plt.title('Positions of peaks {}'.format(file[:5]))
plt.xlabel(r'2$theta$(°)')
plt.ylabel('I (a.u.)')
plt.tight_layout()

# iterate through nested parameter list and fit each peak using the parameters
for param in params:
    x_input_fit = df.two_theta[(df.two_theta >= param[1]) & (df.two_theta <= param[0])].to_numpy()
    y_input_fit = df.yobs[(df.two_theta >= param[1]) & (df.two_theta <= param[0])].to_numpy()
    y_gauss, mean, var = gaussian_fit(x_input_fit, y_input_fit, param[2], param[3])
    mu.append(mean)
    mu_dev.append(np.sqrt(var))
    plt.plot(x_input_fit, y_gauss, label=r'position {:.4f}$\pm${:.4f}'.format(mean, np.sqrt(var)), lw=1.5)

plt.plot(x, y, lw = 0.5)
plt.legend()
plt.show()


"""write data into csv file"""
header = ['peak positions 2theta (deg)', ' std. deviation']
with open('SCO:LCO/{}/SL_data_{}_{}.csv'.format(file[:5], file[:5], file[-11:-8]), 'w') as f:
    writer = csv.writer(f)
    counter = 0
    writer.writerow(header)
    for i, item in enumerate(mu):
        writer.writerow(['{:.4f}'.format(mu[i]), '{:.4f}'.format(mu_dev[i])])




