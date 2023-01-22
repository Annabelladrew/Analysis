import matplotlib.pyplot as plt
import pandas as pd
import os

os.chdir("/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/Analysis/XRD:XRR")
values = ['Substrate', 'c_axis', 'strain_value', 'c_axis_strained', 'c_axis_strained_poisson1', 'strain_value_2', 'c_axis_strained_poisson2', 'c_axis_measured', 'nothing']
df = pd.read_excel('Expected_c_axis.xlsx', skiprows = 7, names = values) #
print(df.head())

y = df.c_axis_measured
y_expected = df.c_axis_strained_poisson1
x = -df.strain_value

print(y)
c_bulk = 3.885



plt.figure(figsize=(10,8))
plt.xlabel('strain (%)')
plt.ylabel(r'c-axis ($\AA$)')
plt.axhline(y = c_bulk, linestyle = '-', color = 'lightblue')
plt.axvline(x = 0, linestyle = '-', color = 'lightblue')
plt.scatter(x, y, label = 'measured')
plt.scatter(x, y_expected, color = 'red', label = 'expected')
plt.legend()
plt.savefig('LCO_0,11sccm_c_axis.png')
plt.show()



# plt.xlabel('strain in %')
# plt.ylabel('film thickness (nm)')
# #plt.savefig('STO_thickness_vs_strain.pdf', dpi = 100)

