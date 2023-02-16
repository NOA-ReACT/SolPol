# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 13:59:33 2021

@author: Lilly Daskalopoulou for the PhD thesis titled "The impact of triboelectrification on desert dust flow dynamics"
         under GPL3 license
"""

import os
import numpy as np
import numpy.ma as ma
import math
import scipy.io
import pandas as pd

# Retrieve the file date from pathname and convert to datetime
def filedate(key_file):
    splt_filename = key_file.split('\\')
    date_str = splt_filename[1].split('_')
    temp1 = pd.to_datetime(date_str[1],format='%d%m%Y')
    temp1 = temp1.date()
    
    return temp1

def darkread(datafile):  
    darktime = datafile[datafile['Data_Type']==7]
    darktime['DateTime'] = pd.to_datetime(darktime['Values'])
    darktime = darktime.reset_index()
    darktime['DateTime'] = darktime['DateTime'].apply(set_date2)
    
    # DC current
    I_dark = datafile[datafile['Data_Type']==9]
    I_dark['Values'] = I_dark['Values'].astype(float)
    I_dark = I_dark.reset_index()
    
    dark_data = pd.concat([I_dark['Values'],darktime['DateTime']],axis=1).set_index(darktime['DateTime'])
    
    # Data reduction
    dark_data_rdc = dark_data.iloc[4::5,:]
    dark_data_rdc = dark_data_rdc.drop(columns='DateTime')
    
    I_dark_max = max(dark_data_rdc['Values'])
    print(I_dark_max)
    
    I_dark_min = min(dark_data_rdc['Values'])
    print(I_dark_min)
    
    I_dark_mean = np.mean(dark_data_rdc['Values'])
    print(I_dark_mean)
    
    return(I_dark_mean,dark_data_rdc)

# Read SolPol dark text file from folder
import glob

solpol_dark_raw = {}
for dark_filename in glob.glob('C:/Users/Lilly Daskalopoulou/Desktop/PhD/Papers/Dust_Orientation/Data_Antikythera/txt_dark/*.txt'):
    solpol_dark_raw[dark_filename[:-4]] = pd.read_csv(dark_filename, sep=' ', header=12)
    print(dark_filename)

## =========================== SolPol Dark Data retrieval ================================= ##
## =================================================================================== ##
solpol_dark_raw1 = pd.DataFrame()
I_dark_mean_list = []
dark_list = []

for dkey,dvalue in solpol_dark_raw.items():
    solpol_dark_raw1 = dvalue
    solpol_dark_dates = dkey
    darkdate = filedate(dkey)
    
    print(darkdate)

    import datetime as dt
    
    def set_date2(dt):
        return dt.replace(day=darkdate.day, month=darkdate.month, year=darkdate.year)
    
    solpol_dark_raw1 = solpol_dark_raw1.drop(columns=['2w'])
    solpol_dark_raw1.columns = ['Data_Type','Values']
    solpol_dark_raw1['Data_Type'] = solpol_dark_raw1['Data_Type'].astype(int)
    
    I_dark_av,dark_data = darkread(solpol_dark_raw1)
    
    I_dark_mean_list.append(I_dark_av) 
    
    dark_list.append(dark_data)
    
I_dark_mean = np.mean(I_dark_mean_list) #total mean for all days

dark_df = pd.DataFrame()
for darkdc in dark_list:
    dark_df = dark_df.append(darkdc)  


