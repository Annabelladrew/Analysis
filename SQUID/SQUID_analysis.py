import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def linear(x, a, b=0):
    return a * x + b

def polynom(x, A, B, C):
    return A*x + B*x**2+ C*x**3

class Measurement:
    def __init__(self, name, reference, label, areas = [1, 1]):
        data_path = '/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/Analysis/SQUID/'
        data_path_ref = '/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/Analysis/SQUID/'
        self.filename = name
        self.filename_ref = reference
        self.area_sample = areas[0]
        self.area_substrate = areas[1]
        self.label = label
        self.label_sample = 'sample (film + substrate) '
        self.label_substrate = 'substrate '
        self.label_temp = label

        # Sample measurement data: path file & names
        df = pd.read_csv(data_path+self.filename)
        self.df = df.drop(df.loc[:, 'Transport Action': 'Map 16'], axis=1)
        names = ['comment', 'Time', 'Temperature', 'Magnetic_field', 'Moment', 'Std_err']
        self.df.columns = names
        # corresponding reference measurement data: path file & names
        df_ref = pd.read_csv(data_path_ref+self.filename_ref)
        self.df_ref = df_ref.drop(df_ref.loc[:, 'Transport Action': 'Map 16'], axis=1)
        self.df_ref.columns = names
        # initializing objects corresponding to columns
        self.df = self.df.loc[(df.iloc[:,4] != '--') & (df_ref.iloc[:,4] != '--')].reset_index()
        self.df_ref = self.df_ref.loc[(df.iloc[:,4] != '--') & (df_ref.iloc[:,4] != '--')].reset_index()
        self.time = self.df.Time
        self.time_ref = self.df_ref.Time
        self.applied_field = self.df.Magnetic_field
        self.moment = self.df['Moment'].astype('float64')/self.area_sample
        self.applied_field_ref = self.df_ref.Magnetic_field
        self.moment_ref = self.df_ref['Moment'].astype('float64')/self.area_substrate
        self.linear_fit_ref = self.linear_fit()
        self.ref_linear_subtraction = self.subtract(self.moment_ref, self.linear_fit_ref)
        self.moment_linear_subtraction = self.subtract(self.moment, self.linear_fit_ref)
        self.moment_linear_subtraction_moment_ref = self.subtract(self.moment_linear_subtraction, self.ref_linear_subtraction)



    # specifically subtracts reference
    # def subtract_reference(self):
    #     return self.moment-self.moment_ref

    # linear fit: fitting slope to reference measurement
    def linear_fit(self):
        x_range = self.applied_field_ref[(self.applied_field <= -30000) & (self.applied_field>= -70000)].to_numpy()
        y_range = self.moment_ref[(self.applied_field <= -30000) & (self.applied_field>= -70000)].to_numpy()
        a_guess = (y_range[0] - y_range[-1]) / np.abs(x_range[-1] - x_range[0])
        p0 = np.array([a_guess, 0.1])
        popt, pcov = curve_fit(linear, x_range, y_range, p0)
        x_range_pos = self.applied_field[(self.applied_field >= 30000) & (self.applied_field <= 75000)].to_numpy()
        y_range_pos = self.moment_ref[(self.applied_field >= 30000) & (self.applied_field <= 75000)].to_numpy()
        #a_guess_pos = (y_range_pos[-1] - y_range_pos[0]) / np.abs(x_range_pos[-1] - x_range_pos[0])
        #print('guess', a_guess, a_guess_pos)
        p0_pos = np.array([a_guess, 0.1])
        popt_pos, pcov_pos = curve_fit(linear, x_range_pos, y_range_pos, p0)
        #print('a_pos', popt_pos[0], 'a', popt[0])
        # evaluate linear function with parameters from fit
        a_average = (popt[0] + popt_pos[0])/2
        #print('average', a_average)

        linear_ref = linear(self.applied_field, a_average)

        return linear_ref

    # def linear_fit_2(self):
    #     y_range = self.moment_linear_subtraction[(self.applied_field <= -20000) & (self.applied_field >= -60000)].to_numpy()
    #     x_range = self.applied_field[(self.applied_field <=- 20000) & (self.applied_field >= -60000)].to_numpy()
    #     a_guess = (y_range[-1]-y_range[0])/np.abs(x_range[-1]-x_range[0])
    #     p0 = np.array([a_guess, 0.1])
    #     popt, pcov = curve_fit(linear, x_range, y_range, p0)
    #     b = popt[-1]
    #     return linear(self.applied_field, *popt) - b

    # def polynom_fit(self):
    #     popt, pcov = curve_fit(polynom, self.applied_field_ref, self.ref_linear_subtraction)
    #     return polynom(self.applied_field, *popt)

    # subtract any function from another
    def subtract(self, f1, f2):
        # self.corrected = True
        return f1 - f2

class MT:
    def __init__(self, name, name_ref):
        path = '/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/Analysis/SQUID/'
        self.name = name
        path_ref = '/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/Analysis/SQUID/'
        self.name_ref = name_ref
        df1 = pd.read_csv(path+self.name)
        self.MT = df1.drop(df1.loc[:, 'Transport Action': 'Map 16'], axis=1)
        self.MT.columns = ['comment', 'Time', 'Temperature', 'Magnetic_field', 'Moment', 'Std_err']
        df2 = pd.read_csv(path_ref + self.name_ref)
        self.MT_ref = df2.drop(df2.loc[:, 'Transport Action': 'Map 16'], axis=1)
        self.MT_ref.columns = ['comment', 'Time', 'Temperature', 'Magnetic_field', 'Moment', 'Std_err']

        self.temperature = self.MT.Temperature
        self.moment = self.MT.Moment
        self.temperature_ref = self.MT_ref.Temperature
        self.moment_ref = self.MT_ref.Moment



def main():
    area_substrate = 25e-6  # area of large substrate/film [m^2]
    area_sample = 5.876e-6  # mm**2
    areas = [area_sample, area_substrate]

    # M1 = Measurement('22163_LCOT/'+'LCOT_MH_5K_FC_p1T.csv', 'Reference/STO/A6/'+'STO_MH_5K_FC_1T.csv', label = '5K')
    M2 = Measurement('22163_LCOT/' + 'LCOT_MH_10K_FC+1T.csv', 'Reference/STO/A6/'+'STO_MH_10K_FC_1T.csv', '10K')
    # M3 = Measurement('22163_LCOT/' + 'LCOT_MH_20K_FC+1T.csv', 'Reference/STO/A6/'+'STO_MH_20K_FC_1T.csv')
    # M4 = Measurement('22163_LCOT/' + 'LCOT_MH_50K_FC+1T.csv', 'Reference/STO/A6/'+'STO_MH_50K_FC_1T.csv')
    # M5 = Measurement('22163_LCOT/' + 'LCOT_MH_100K_FC+1T.csv', 'Reference/STO/A6/'+'STO_MH_100K_FC_1T.csv')
    # M6 = Measurement('22163_LCOT/' + 'LCOT_MH_200K_FC+1T.csv', 'Reference/STO/A6/'+'STO_MH_200K_FC_1T.csv')
    # M7 = Measurement('22163_LCOT/' + 'LCOT_MH_300K_FC+1T.csv', 'Reference/STO/A6/'+'STO_MH_300K_FC_1T.csv')
    # M8 = Measurement('22163_LCOT/' + 'LCOT_MH_350K_FC+1T.csv', 'Reference/STO/A6/'+'STO_MH_350K_FC_1T.csv')
    # M = MT('22163_LCOT/' + 'LCOT_MT_FC+1T.csv', 'Reference/STO/A6/' + 'STO_MT_FC_1T.csv')
    #
    #D1 = Measurement('22157b_LCOT/' + 'LCOT_MH_5K_FC+1T_2.csv', 'STO_MH_5K_FC_1T.csv', areas)
    # D2 = Measurement('22157b_LCOT/' + 'LCOT_MH_350K_FC+1T_2.csv', 'STO_MH_350K_FC_1T.csv', cut_off)
    # D3 = Measurement('22157b_LCOT/' + 'LCOT_MH_350K_FC+1T_2.csv', 'STO_MH_350K_FC_1T.csv', cut_off)
    # N1 = Measurement('22169_LCOT/' + '22169_LCOT_MH_5K_FC_0d5T.csv', 'Reference/STO/B1/'+'STOsubs_MH_5K_FC_0d5T.csv', '5K')
    N2 = Measurement('22169_LCOT/' + '22169_LCOT_MH_10K_FC_0d5T.csv', 'Reference/STO/B1/'+'STOsubs_MH_10K_FC_0d5T.csv', '10K')
    # N3 = Measurement('22169_LCOT/' + '22169_LCOT_MH_50K_FC_0d5T.csv', 'Reference/STO/B1/'+'STOsubs_MH_50K_FC_0d5T.csv', '50K')
    # N4 = Measurement('22169_LCOT/' + '22169_LCOT_MH_100K_FC_0d5T.csv', 'Reference/STO/B1/'+'STOsubs_MH_100K_FC_0d5T.csv', '100K')
    N5 = Measurement('22169_LCOT/' + '22169_LCOT_MH_200K_FC_0d5T.csv', 'Reference/STO/B1/'+'STOsubs_MH_200K_FC_0d5T.csv', '200K')
    # N6 = Measurement('22169_LCOT/' + '22169_LCOT_MH_250K_FC_0d5T.csv', '22169_LCOT/' + '22169_LCOT_MH_250K_FC_0d5T.csv', '250K') #no substrate reference
    # N7 = Measurement('22169_LCOT/' + '22169_LCOT_MH_275K_FC_0d5T.csv', '22169_LCOT/' + '22169_LCOT_MH_275K_FC_0d5T.csv', '275K') #no substrate reference
    # N8 = Measurement('22169_LCOT/' + '22169_LCOT_MH_280K_FC_0d5T.csv', '22169_LCOT/' + '22169_LCOT_MH_280K_FC_0d5T.csv', '280K') #no substrate reference
    # N9 = Measurement('22169_LCOT/' + '22169_LCOT_MH_300K_FC_0d5T.csv', '22169_LCOT/' + '22169_LCOT_MH_300K_FC_0d5T.csv', '300K') #no substrate reference
    # N10 = Measurement('22169_LCOT/' + '22169_LCOT_MH_350K_FC_0d5T.csv', 'Reference/STO/B1/' + 'STOsubs_MH_350K_FC_0d5T.csv', '350K')
    # N = MT('22169_LCOT/'+'22169_LCOT_MT_FC_0d5T.csv', 'Reference/STO/B1/'+'STOsubs_MT_FC_0d5T.csv')

    O1 = Measurement('22172_LCOL/' + '22172_LCOL_MH_10K_FC_p0d5T_2.csv', 'Reference/LAO/C3/'+'LAOsubs_22172_MH_10K_FC_p0d5T.csv', '10K')
    O2 = Measurement('22172_LCOL/' + '22172_LCOL_MH_50K_FC_p0d5T.csv', 'Reference/LAO/C3/' + 'LAOsubs_22172_MH_10K_FC_p0d5T.csv', '50K') #no substrate reference
    O3 = Measurement('22172_LCOL/' + '22172_LCOL_MH_200K_FC_p0d5T.csv', 'Reference/LAO/C3/' + 'LAOsubs_22172_MH_200K_FC_p0d5T.csv', '200K')  # no substrate reference
    O4 = Measurement('22172_LCOL/' + '22172_LCOL_MH_300K_FC_p0d5T.csv', 'Reference/LAO/C3/' + 'LAOsubs_22172_MH_300K_FC_p0d5T.csv', '300K')  # no substrate reference
    O = MT('22172_LCOL/' + '22172_LCOL_MT_FC_p0.5T.csv', 'Reference/LAO/C3/' + 'LAOsubstrate_22172_MT_F_p0.5T.csv')
    # conversion: m_B = 9.274099e-24 Am^2
    # 1 emu = 1e-3 Am**2 ->  Am**2 = 1e3 m_B / 9.274099e-24
    # numbers of unit cells in volume of film
    thickness_STO_thin = 35
    thickness_STO = 70 # thickness in uc of film
    thickness_LAO = 50
    height = 0.5e-3 # height of substrate: [m]
    STO_lattice_param = 3.905e-10 # [m] in plane
    LAO_lattice_param = 3.789e-10 # [m] in plane
    # number_uc = area_substrate*thickness/(substrate_lattice_param**2)
    number_uc_sample_STO = area_substrate*thickness_STO/(STO_lattice_param**2)
    number_uc_sample_STO_thin = area_substrate * thickness_STO_thin / (STO_lattice_param ** 2)
    number_uc_sample_LAO = area_substrate*thickness_LAO / (LAO_lattice_param ** 2)
    # number_uc_substrate = area_substrate*height/substrate_lattice_param**3
    # oested in tesla factor
    field_factor = 1/10000
    #print("{:e}".format(number_uc))
    # calculating magnetic moment (emu) times conversion factor gives Bohr magneton dividing by number_uc gives Bohr Magneton/uc
    conversion_factor = 1e21/(9.274099)

    print("{:e}".format(number_uc_sample_LAO))

    """plot"""
    """MT"""

    """22169"""
    # plt.figure(figsize=(10, 8))
    # plt.xlabel('applied Field (T)')
    # plt.ylabel('Magnetic moment (emu)')
    # plt.title('22172_LCOL (FC+0.5T)')
    # plt.plot(O.temperature, O.moment)
    # plt.tight_layout()
    # plt.plot(N.temperature, (N.moment-N.moment_ref)*factor/number_uc_sample, label = 'raw film on substrate - substrate')
    # plt.scatter(N.temperature_ref, N.moment_ref+4.48E-5, label = 'raw substrate', color = 'orange')
    # plt.scatter(N.temperature, N.moment, label = 'raw film + substrate')
    # plt.legend()

    #plt.yscale('log')

    """MH"""
    plt.figure(figsize=(10, 8))
    plt.xlabel('applied Field (T)')
    plt.ylabel(r'Moment $\mu_b$/uc') # 'Moment (emu)'
    #plt.title('22172_LCOL M(H) subtracted by substrate at {} (FC+0.5T)'.format(O3.label_temp))
    plt.title('Comparison')
    plt.plot(M2.applied_field * field_factor, (M2.moment - M2.moment_ref) * conversion_factor / number_uc_sample_STO_thin,'-o', label='STO 1', markersize = 1)
    plt.plot(O1.applied_field*field_factor, (O1.moment-O1.moment_ref)*conversion_factor/number_uc_sample_LAO, '-o', label = 'LAO', markersize = 4)
    plt.plot(N2.applied_field * field_factor, (N2.moment - N2.moment_ref) * conversion_factor / number_uc_sample_STO, '-o', label='STO', markersize = 1)
    # plt.plot(O2.applied_field*field_factor, O2.moment, label = O2.label_temp)
    # plt.plot(O3.applied_field*field_factor, O3.moment-O3.moment_ref, '-o', label = O3.label_temp)
    # plt.plot(O4.applied_field * field_factor, O4.moment - O4.moment_ref, '-o', label=O4.label_temp)
    #plt.plot(O4.applied_field*field_factor, O4.moment, label = O4.label_temp)

    # plt.xlim(-0.75,0.75)
    # plt.ylim(-0.0003, 0.0003)
    # plt.plot(O1.applied_field_ref*field_factor, O1.moment_ref, label = O1.label_substrate)
    # plt.plot(O1.applied_field_ref*field_factor, O1.linear_fit_ref)
   # plt.legend()



    """22169"""
    # #plt.plot(N6.applied_field * field_factor, N6.moment, label = N6.label_substrate)
    # plt.plot(N1.applied_field * field_factor, N1.moment-N1.moment_ref, label = N1.label_sample)
    # plt.plot(N2.applied_field * field_factor, N2.moment-N2.moment_ref, label=N2.label)
    # plt.plot(N3.applied_field * field_factor, N3.moment_ref, label=N3.label_substrate)
    # plt.plot(N4.applied_field * field_factor, N4.moment_ref, label=N4.label_substrate)
    # plt.plot(N5.applied_field * field_factor, N5.moment_ref, label=N5.label_substrate)
    # plt.plot(N6.applied_field * field_factor, N6.moment, label=N6.label_sample)
    # plt.plot(N7.applied_field * field_factor, N7.moment, label=N7.label_sample)
    # plt.plot(N8.applied_field * field_factor, N8.moment, label=N8.label_sample)
    # plt.plot(N9.applied_field * field_factor, N9.moment, label=N9.label_sample)
    # plt.plot(N10.applied_field * field_factor, N10.moment_ref, label=N10.label_substrate)
    # plt.vlines(0, -2E-5, 2E-5, color = 'grey', alpha = 0.5)
    # plt.xlim(-0.5, 0.5)
    # plt.ylim(-4E-5, 4E-5)
    # aplt.legend()


    # plt.figure(figsize=(10, 8))
    # plt.xlabel('applied Field (T)')
    # plt.ylabel(r'Moment $\mu_b$/uc')
    # plt.title(r'STO B1 raw signal & linear correction at 350K converted to $\mu_B$/uc assuming thickness of film')
    # plt.plot(N5.applied_field * field_factor, (N5.moment-N5.moment_ref)*conversion_factor/number_uc_sample_STO, label= '22168 200K')  # *factor/number_uc
    # plt.plot(N6.applied_field * field_factor, (N6.ref_linear_subtraction*factor/number_uc_sample), label=N6.label_substrate)
    # plt.legend()
    # plt.xlim(-0.5,0.5)
    # plt.ylim(-0.2, 0.2)
    """22163"""
    # plt.plot(M1.applied_field * field_factor, M1.moment, label = 'raw film + substrate')
    # plt.plot(M1.applied_field_ref * field_factor, M1.moment_ref, label = 'raw substrate')
    # plt.plot(M1.applied_field * field_factor, M1.moment-M1.moment_ref, label = 'film')
    # plt.plot(M1.applied_field_ref * field_factor, M1.ref_linear_subtraction, label = 'substrate linear sub')

    # plt.hlines(0, -3, 3, color = 'grey', alpha = 0.3)
    # plt.vlines(0,-0.5,0.5, color = 'grey', alpha = 0.3)
    # plt.xlim(-0.2,0.2)
    # plt.ylim(-0.1,0.1)
    plt.tight_layout()
    plt.legend()
   # plt.savefig('LCOT_linear_sub_10K.png')
    #plt.savefig('22163_LCOT_MH_5K_FC1T_all.png')
    plt.show()


main()