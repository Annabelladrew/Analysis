import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def linear(x, a, b):
    return a * x + b

class Measurement:
    def __init__(self, name):
        data_path = '/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/Analysis/XTreme/processed data/XMCD grazing/'
        self.filename = name
        data_names = ['Energy', 'minus1', 'plus1', 'plus2', 'minus2', 'minus3', 'plus3', 'I_cp', 'I_cm', 'cp_Hm', 'cm_Hm',
                      'XMCD_Hp', 'XMCD_Hm', 'XAS', 'XMCD', 'XAS_normalized_by_jump', 'XMCD_normalized_by_XAS_jump']
        self.df = pd.read_csv(data_path+self.filename, names = data_names, header = 0, skiprows = 1)
        self.energy = self.df.Energy
        self.E_max = 572.5 # keV defines range for linear fit pre peak range
        self.E_cutoff = 600 # keV
        self.E_min = 589 # keV defines range after peaks for Fermi step h
        self.epsilon = 0 # defines position of Fermi step relative to peaks
        self.delta = 0.1 # defines width of Fermi step
        self.I_p = self.df.I_cp
        self.I_m = self.df.I_cm
        self.xmcd = False
        self.subtracted_I_p = self.subtract_linear_background(self.I_p)
        self.subtracted_I_m = self.subtract_linear_background(self.I_m)
        self.subtracted_Fermi_I_p = self.Subtract_Fermi(self.subtracted_I_p)
        self.subtracted_Fermi_I_m = self.Subtract_Fermi(self.subtracted_I_m)
        self.xmcd = self.XMCD(self.subtracted_Fermi_I_m, self.subtracted_Fermi_I_p)


    def plot_linear_background(self):
        plt.figure(figsize=(10, 8))
        plt.plot(self.energy, self.subtracted_I_p, color='r')
        plt.plot(self.energy, self.I_p, lw=3, label='cp')
        plt.legend()
        plt.hlines(0, 565, 605, color='g')
        plt.title('linear fit pre edge range')
        plt.tight_layout()
        plt.show()
        plt.figure(figsize=(10, 8))
        plt.plot(self.energy, self.subtracted_I_m, color='r')
        plt.plot(self.energy, self.I_m, lw=3, label='cm', color='orange')
        plt.legend()
        plt.hlines(0, 565, 605, color='g')
        plt.title('linear fit pre edge range')
        plt.tight_layout()
        plt.show()

    def plot_xmcd(self):
        plt.figure(figsize=(10, 8))
        plt.plot(self.energy, self.xmcd, lw=2)
        plt.hlines(0, 565, 605)
        plt.show()
        plt.tight_layout()

    def XMCD(self, I_m, I_p):
        return I_m-I_p

    def subtract_linear_background(self, I):
        linear_background = self.linear_background_fit(I)
        subtracted_linear = I - linear_background
        return subtracted_linear

    def linear_background_fit(self, I):
        x_range = self.energy[self.energy <= self.E_max].to_numpy()
        y_range = I[self.energy <= self.E_max].to_numpy()
        # fitting parameter estimates : it might be enough to estimate slope of one polarization (they are very similar)
        a_guess = (y_range[-1]-y_range[0])/np.abs(x_range[-1]-x_range[0])
        p0 = np.array([a_guess, 0.1])
        popt, pcov = curve_fit(linear, x_range, y_range, p0)
        # evaluate linear function with parameters from fit
        linear_background_function = linear(self.energy, *popt)
        return linear_background_function

    def Subtract_Fermi(self, I):
        return I - self.Fermi_step(I)

    def plot_Fermi(self):
        plt.figure(figsize=(10, 8))
        plt.plot(self.energy, self.subtracted_Fermi_I_p, color='r')
        plt.plot(self.energy, self.subtracted_I_p, lw=3, label='cp')
        plt.legend()
        plt.hlines(0, 565, 605, color='g')
        #plt.title('Fermi step subtracted')
        plt.tight_layout()
        plt.show()
        plt.figure(figsize=(10, 8))
        plt.plot(self.energy, self.subtracted_Fermi_I_m, color='r')
        plt.plot(self.energy, self.subtracted_I_m, lw=3, label='cm', color='orange')
        plt.legend()
        plt.hlines(0, 565, 605, color='g')
        #plt.title('Fermi step subtracted')
        plt.tight_layout()
        plt.show()

    def Fermi_step(self, I):
        def mu_step(E_l3, E_l2, height, a):
            b = 1 - a
            return height * (1 - a * (1 / (1 + np.exp((self.energy - E_l3 - self.epsilon) / self.delta))) - b * (1 / (1 + np.exp((self.energy - E_l2 - self.epsilon) / self.delta))))
        # select energy range for step size h and find mean as an estimate
        energy_post_range = I[(self.energy >= self.E_min) & (self.energy <= self.E_cutoff)]
        h = energy_post_range.mean()
        print('h = ', h)
        # find peak positions
        l3_index = I[(575 <= self.energy) & (self.energy <= 580)].idxmax()
        peak_l3 = self.energy[l3_index]
        print('peak L3: ', peak_l3)
        l2_index = I[(583 <= self.energy) & (self.energy <= 588)].idxmax()
        peak_l2 = self.energy[l2_index]
        print('peak L2: ', peak_l2)
        # plug in parameters into step function and return function
        return mu_step(peak_l3, peak_l2, h, 2/3)




def main():
    # write cut-off
    m1 = Measurement('157b_1642-1659_2T_75K.csv')
    m2 = Measurement('157b_1710-1727_6T_75K.csv')
    # m3 = Measurement('157b_1743-1800_minus_6T_75K.csv') # nicely transform all .dat into csv.!
    plt.figure(figsize=(10, 8))
    plt.title('grazing incidence, 75K')
    plt.plot(m1.energy, m1.xmcd, lw = 2, color = 'blue', label = '2T')
    plt.plot(m2.energy, m2.xmcd, lw = 2, color = 'orange', label = '6T')
    plt.legend()
    plt.tight_layout()
    plt.show()

main()
