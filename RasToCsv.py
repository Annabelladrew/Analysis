# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 14:21:58 2020

@author: simon
"""

import tkinter as tk
import tkinter.filedialog
import easygui



class MyFrame(tk.Frame):
    def __init__(self):
        self.transformData()
    
    def transformData(self):

        filename = easygui.fileopenbox()
        # Tk().withdraw()

        b = filename[:-3]

        raw = open(b+'ras', 'r',errors='ignore')
        data_raw = [x for x in raw]
        raw.close()
        data_only = []
        append_it = False
        for count, i in enumerate(data_raw):
            if append_it:
                data_only.append(i)
            if count == 372:
                a = 10
                continue
            if i == '*RAS_INT_START\n':
                append_it = True
            elif i == '*RAS_INT_END\n':
                append_it = False
                data_only = data_only[:-1]
            else:
                continue
        data_new = []
        for i in data_only:
            sub_data = i.split(' ')
            new_sub = []
            for j in sub_data:
                new_sub.append(float(j))
            data_new.append(new_sub)
    

        new_txt = open(b+'csv', 'w')
        new_txt.write('#twotheta, yobs, ycal, bkg, diff\n')
        for i in data_new:
            # line = str(i[0]) + ',' + str(i[1])
            line = str(i[0]) + ',' + str(i[1]) + ',' + str(i[2]) + ',' + str(i[2]) + ',' + str(i[2]) + '\n'
            print(line)
            new_txt.write(line)
        new_txt.close()
        
MyFrame()


        


