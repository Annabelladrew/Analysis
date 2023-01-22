import matplotlib.pyplot as plt
import pandas as pd
import os

data_path = '/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/RIXS'
data_file = 'c_axis_strain.csv'
#os.chdir("D:/01Coding/Python/data_sets/myowndata")
print(data_path+data_file)
df = pd.read_csv(data_path+data_file)

print(df)







# plt.figure(figsize=(10,8))
# plt.plot(strain, c_axis_thickness, 'bo')
# plt.xlabel('strain in %')
# plt.ylabel('film thickness (nm)')
# #plt.savefig('STO_thickness_vs_strain.pdf', dpi = 100)
# plt.show()

