import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def linear(x,m,q):
    # linear function
    return m*x+q

# modify range according to your data to input into the fit
x_range = np.arange(0,800,1)
y_range = np.arange(0,40,1)

# input coordinates
x_i = np.array([375, 500, 750])        # x_i = np.array([100, 120, 190])
y_i = np.array([15, 17, 35])        # y_i = np.array([23.04, 31.98, 39.05])
#y_err = np.array([1.7, 3.26, 8.3])        # y_err = np.array([1.7, 3.26, 8.3])
x_units = 'sec'
y_units = 'uc' # r'$\AA$'

# fit
a_guess = np.abs(x_i[0]-x_i[1])/np.abs(y_i[0]-y_i[1])
popt, pcov = curve_fit(linear, x_i, y_i)# sigma = y_err, absolute_sigma=True
slope = popt[0] # uc/s
dev_slope = np.sqrt(pcov[0][0])
delta_slope = dev_slope/slope**2
print(1/slope)

plt.figure(figsize=(10,8))
plt.tight_layout()
y_fit = linear(x_range, *popt)
plt.ylabel(r'$\Lambda$ ({})'.format(y_units))
plt.xlabel('sputtering time ({})'.format(x_units))
plt.plot(x_range, y_fit)
#plt.errorbar(x_i, y_i, yerr = y_err, fmt = 'r.', label = 'data')
plt.plot(x_i, y_i, 'ro')
plt.title(r'SCO/LCO//LAO 0.015 sccm, growth rate = {} $\pm$ {} {}/{}, y offset = {} {}'.format(round((1/slope),3), round(delta_slope), x_units, y_units, round(popt[1]), y_units))
#plt.savefig('growth_rate_LCOL_0.015.png')
plt.show()