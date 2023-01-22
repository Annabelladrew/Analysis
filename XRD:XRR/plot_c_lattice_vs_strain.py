import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv('')



plt.figure(figsize=(10,8))
plt.plot(strain, c_axis_thickness, 'bo')
plt.xlabel('strain in %')
plt.ylabel('film thickness (nm)')
#plt.savefig('STO_thickness_vs_strain.pdf', dpi = 100)
plt.show()

