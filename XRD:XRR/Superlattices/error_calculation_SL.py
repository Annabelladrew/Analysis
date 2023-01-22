import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import csv

"""import csv"""
data_path = '/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/Analysis/XRD:XRR/Superlattices/SCO:LCO/22169/'
file = 'SL_data_22169_002.csv'
values_header = ['x_i', 'peak_pos', 'delta_2theta', 'theta', 'delta_pos', 'y_i', 'y_err'] # values
df = pd.read_csv(data_path+file, names = values_header, header = 0)
print(df.head())
def linear_func (x, intercept, slope):
    y = slope*x + intercept
    return y

popt, pcov = curve_fit(linear_func, df.x_i, df.y_i, sigma = df.y_err, absolute_sigma=True)
intercept = popt[0]
slope = popt[1]
dev_inter = np.sqrt(pcov[0][0])
dev_slope = np.sqrt(pcov[1][1])
y_fit = linear_func(df.x_i, *popt)

plt.figure(figsize=(10,8))
plt.title(r'{}: y = ({}$\pm${})*x'.format(file[8:-4], slope, dev_slope))
plt.errorbar(df.x_i, df.y_i, yerr = df.y_err, fmt = 'r.', label = 'data')
plt.plot(df.x_i, y_fit, label = 'fit')
plt.tight_layout()
plt.legend()
plt.show()

"""add column with std. deviation of slope"""
df_result = pd.DataFrame({'slope': [slope],'delta_slope': [dev_slope]})
print(df_result)
df_result.to_csv(data_path+file[:-4]+'_result.csv')
