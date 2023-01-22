"""Created Friday 29. April 2022
author: Annabella Drewanowski"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import pandas as pd


data_path = '/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/Analysis/XRD:XRR/LCO/22163_LCOT/'
file = '22163_002_XRD.csv'
values_header = ['two_theta', 'yobs', 'ycal', 'bkg', 'diff']
df = pd.read_csv(data_path+file, names = values_header, header = 0)

data_path2 = '/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/Analysis/XRD:XRR/LCO/22157b_LCOT/'
file2 = '22157b_LCOT_XRD.csv'
values_header = ['two_theta', 'yobs', 'ycal', 'bkg', 'diff']
df2 = pd.read_csv(data_path2+file2, names = values_header, header = 0)



"""uncomment the following if you wanna plot multiple xrds"""
# files = [f for f in glob.glob('CSV/*.csv')]
# files.sort()
# dataframes = [pd.read_csv(file) for file in files ]
#
# colors = ['red', 'blue','green','purple','orange']
# substrates = ['STO', 'LSAT', 'NGO', 'LAO', 'NAO'] #substrates = ['LAO', 'LAO annealed']
# fig, ax = plt.subplots(figsize=(10,8))
# offset = 1

# range (in case there is signal for XRR -  then there will be an annoying line)
two_theta_min = 42
two_theta_max = 52


# for i,data in enumerate(dataframes):
#     data = data.loc[data.iloc[:,1] != 0] # select values that are not zero
#     data = data.loc[(data.iloc[:,0] >= two_theta_min) & (data.iloc[:,0] <= two_theta_max)]
#     plt.plot(data.iloc[:,0], data.iloc[:,1]*offset, label = substrates[i], color = colors[i], lw = 1.5) #data.iloc[:,0], data.iloc[:,1]*offset  label = files[i].split('_')[2]
#     offset = offset*1000


# set general font size
#plt.rcParams.update({'font.size': 22})
#plt.rcParams['font.size'] = '26'

# Set tick font size
# for label in (ax.get_xticklabels() + ax.get_yticklabels()):
#     label.set_fontsize(26)
plt.figure(figsize=(10, 8))
plt.yscale('log')
plt.plot(df.two_theta, df.yobs, label = '22163 LCOT')
plt.plot(df2.two_theta, df2.yobs, label = '22157b LCOT')
#plt.xlim((42, 52))
#plt.plot(x, y, lw = 0.5)
plt.xlabel(r'2$\Theta (°)$')
plt.ylabel('Intensity (arb. units)')
#plt.title('22102 002 Peak (annealed)')
plt.legend()
plt.tight_layout()
plt.savefig('LCO_163_157b.png', dpi = 100)
plt.show()