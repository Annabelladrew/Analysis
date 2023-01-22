import pandas as pd


data_path = '/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/RIXS'
data_file = 'c_axis_strain.csv'

print(data_path+data_file)
df = pd.read_csv(data_path+data_file)

print(df)