### ==================================== ###
###    SolPol data retrieval algorithm   ###
### ==================================== ###

"""
Created on Thu Jun 18 18:32:02 2020

@author: Lilly Daskalopoulou for the PhD thesis titled "The impact of triboelectrification on desert dust flow dynamics"
         under GPL3 license
"""

import os
import numpy as np
import numpy.ma as ma
import math
import scipy.io
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cmaps
import pytz
import seaborn as sns

import warnings
import pvlib

from pvlib import atmosphere
from pvlib.tools import datetime_to_djd, djd_to_datetime


# Retrieve the file date from pathname and convert to datetime
def filedate(key_file):
    splt_filename = key_file.split('\\')
    date_str = splt_filename[1].split('_')
    temp1 = pd.to_datetime(date_str[1],format='%d%m%Y')
    temp1 = temp1.date()
    
    iris_str = date_str[0] # input iris size from filename
    return temp1, iris_str


# Phase complex averaging
def phase_ave(s):
    a = []
    for i in range(0,np.max(A1_Ch1_phs.index),20):
        
        
        if (i+4) >= np.max(A1_Ch1_phs.index):
            break
        elif (i+9) >= np.max(A1_Ch1_phs.index) and np.max(A1_Ch1_phs.index) > (i+4):
            a.append(-(- s[i+4]))
        elif (i+14) >= np.max(A1_Ch1_phs.index) and np.max(A1_Ch1_phs.index) > (i+9):
            a.append(-(- s[i+4]+ s[i+9])/2)
        elif (i+19) >= np.max(A1_Ch1_phs.index) and np.max(A1_Ch1_phs.index) > (i+14):
            a.append(-(- s[i+4]+ s[i+9] - s[i+14])/3)
        else:
            a.append(-(- s[i+4] + s[i+9] - s[i+14] + s[i+19])/4)
    
    print(a)
    return a


# Plot Polarization
def pol_plot(name,ii,C,D,xdate):  
    plt.figure(ii)
    
    myFmt = mdates.DateFormatter('%d-%m %H:%M')
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(myFmt)
    ax.xaxis_date()
    
    plt.plot(C,'bo',label='Lin_Pol')
    plt.plot(D,'ro',label='Circ_Pol')
    plt.xlabel('Time (UTC)') #for Athens, Cyprus & several Antik.
    plt.ylabel('Polarization in ppm')
    plt.grid(True)
    plt.ylim((-200,200))
#    plt.xlim(xdate.replace(hour=4,minute=30,second=0),xdate.replace(hour=17,minute=30,second=0))
    plt.xlim(xdate.replace(hour=6,minute=30,second=0),xdate.replace(hour=19,minute=30,second=0)) #for Mindelo data
    plt.xticks(rotation=30)
    plt.tight_layout()
    plt.legend()
    
    # save plots
#    plt.savefig(name + '.png')
#    plt.savefig(name + '.eps',format ='eps')


# Plot EVPAs
def evpa_plot(name1,jj,E,F):
    plt.figure(jj)
    fig, ax = plt.subplots()
    
    plt.plot(E,F,'ko',label='EVPA')
    plt.plot(np.zeros(180),'r')
    plt.xlim((0,180))
    # plt.plot(np.zeros(math.ceil(max(EVPA_perday_df['zenith']))+int(min(EVPA_perday_df['zenith']))),'r')
    plt.xlabel('Solar Zenith Angle (degs)') 
    plt.ylabel('Electric Vector Polarization Angle (degs)')
    plt.grid(True)
    plt.title(name1)
    plt.xticks(rotation=30)
    plt.tight_layout()
    
    # save plots
#    plt.savefig(name1 + '.png')
#    plt.savefig(name1 + '.eps',format ='eps')
        

# Chauvenet criterion for outlier detection for 120 mins
def chauvenet(A,B):
    chauvenet_120m = []
    # mean value per 1 hour
    A_hrmean = A.resample('120T', label='right', closed='right').mean()
    A_hrmean = A_hrmean.reset_index(drop=True)
    
    # standard deviation
    A_std = A.resample('120T', label='right', closed='right').std()
    A_std = A_std.reset_index(drop=True)
    
    # counts resampled items
    A_cnt = A.resample('120T', label='right', closed='right').count()
    A_cnt = A_cnt.reset_index(drop=True)
    
    # Merge with 1-hr mean dataframe
    A_df120m_merge = pd.concat([A_hrmean,A_cnt,A_std],axis=1)
    A_df120m_merge.columns=['Mean','Count','Std']
    
    # Duplicate mean/std values per counts interval in the initial dataset
    mean_std_dup = A_df120m_merge.loc[np.repeat(A_df120m_merge.index.values,A_df120m_merge['Count'])]
    mean_std_dup = mean_std_dup.reset_index(drop=True)
    mean_std_dup = mean_std_dup.set_index([A.index])
    
    # Calculate the Chauvenet parameter
    #Zscore = abs(xi-μ)/σ , Pi = erfc(Zscore/sqrt(2)), Chauv_par = Pi*n < 1/2
    B['Zscore'] = abs(A-mean_std_dup['Mean'])/mean_std_dup['Std']
    B['ChauvenetPar1'] = scipy.special.erfc(B['Zscore']/np.sqrt(2))
    
    # Find entries which satisfy the first Chauvenet filtering, Dmax equals to 1/2 for 50% chance of survival for each measurement
    A_Ch = B['ChauvenetPar1']*mean_std_dup['Count'] < 0.5
    
    # Count the number of false/true arguments, true -> discarding the respective value
    A_Ch.value_counts()
    
    # Create dataframe with the smoothed dataset
    A_Ch = A_Ch.to_frame('Chauv_120min')
    chauvenet_120m_pre = pd.concat([A_Ch,A],axis=1)
    
    # After Chauvenet Dataframe
    chauvenet_120m = chauvenet_120m_pre.loc[chauvenet_120m_pre['Chauv_120min'] == False]
    
    print(chauvenet_120m)
    return chauvenet_120m


## ================ Specify data folders/path ===================== ##
## ================================================================ ##
# Specify data folder e.g. Data_Antikythera or Data_Athens, Data_Cyprus, Data_mindelo
from pathlib import Path
Data_folder = Path('C:/Users/Lilly Daskalopoulou/Desktop/PhD/Papers/Dust_Orientation/Data_Mindelo/txt')
str(Data_folder)


# Read SolPol raw texts file from folder
import glob

solpol_raw = {}
for filename in glob.glob('C:/Users/Lilly Daskalopoulou/Desktop/PhD/Papers/Dust_Orientation/Data_Mindelo/temp/*.txt'):
    solpol_raw[filename[:-4]] = pd.read_csv(filename, sep=' ', header=12)
    print(filename)

# Insert SolPol position
# lat, long, alt alter for each location
# @Antik: (35.86099, 23.30982, 193), @Athens:(37.966295, 23.710624, 60), @Cy:(35.14063, 33.38135, 181)
# @Mindelo: (35.86099, 23.30982, 20)
    
location = "Mindelo"

if location == 'Antikythera' :
    lat = 35.86099
    long = 23.30982
    alt = 193
elif location == 'Athens':
    lat = 37.966295
    long = 23.710624
    alt = 60
elif location == 'Cyprus':
    lat = 35.14063
    long = 33.38135
    alt = 181
elif location == 'Mindelo':
    lat = 16.87775
    long = -24.994889
    alt = 20

## ====================== SolPol Data retrieval ====================== ##
## =================================================================== ##
cnt = 0   
temp_lin = []
last_lst = []
b = []
c = []
k = []
l = []
f = []

solpol_list = [] # for the raw data plotting

## ===================== Pre-processing ============================== ##
## =================================================================== ##
jj =0 
for key,value in solpol_raw.items():
    solpol_raw1 = value
    solpol_dates = key
    date, aperture_str = filedate(key)
    print(date)
    print(aperture_str)
   
    # input aperture size from filename
    if aperture_str == 'poliris45':
        aperture = 4.5  # in mm diameter
    elif aperture_str == 'poliris55':
        aperture = 5.5
    elif aperture_str == 'poliris7':
        aperture = 7.0
    else: aperture = 5.5

    # Function to replace current date to filename date
    import datetime as dt
    
    def set_date(dt):
        return dt.replace(day=date.day, month=date.month, year=date.year)
    
    solpol_raw1 = solpol_raw1.drop(columns=['2w'])
    solpol_raw1.columns = ['Data_Type','Values']
    solpol_raw1['Data_Type'] = solpol_raw1['Data_Type'].astype(int)
    
    # Polarimeter position (degs)
    pol_pos = solpol_raw1[solpol_raw1['Data_Type']==1]
    pol_pos['Values'] = pol_pos['Values'].astype(int)
    pol_pos = pol_pos.reset_index()
    
    # Rotator/Polarizer Position (degs)
    rot_pos = solpol_raw1[solpol_raw1['Data_Type']==2]
    rot_pos['Values'] = rot_pos['Values'].astype(float)
    rot_pos = rot_pos.reset_index()
    
    # PEM setting (nm)
    PEM_set = solpol_raw1[solpol_raw1['Data_Type']==3]
    PEM_set['Values'] = PEM_set['Values'].astype(int)
    PEM_set = PEM_set.reset_index()
    
    # Retardation (waves)
    retard = solpol_raw1[solpol_raw1['Data_Type']==4]
    retard['Values'] = retard['Values'].astype(float)
    retard = retard.reset_index()
    
    # Wavelength Filter - Wavelength Bandwidth
    wave_flt = solpol_raw1[solpol_raw1['Data_Type']==5]
    wave_flt['Values'] = wave_flt['Values'].astype(int)
    wave_flt = wave_flt.reset_index()
    
    # ND Filter
    ND_flt = solpol_raw1[solpol_raw1['Data_Type']==6]
    ND_flt['Values'] = ND_flt['Values'].astype(float)
    ND_flt = ND_flt.reset_index()
    
    # Datetime
    time = solpol_raw1[solpol_raw1['Data_Type']==7]
    time['DateTime'] = pd.to_datetime(time['Values'])
    time = time.reset_index()
    time['DateTime'] = time['DateTime'].apply(set_date)
    time_ave = time.iloc[4::5,:]
    time_ave = time_ave.drop(['index','Values','Data_Type'],axis=1)
    
    # Bias Voltage on Diode
    bias_volt = solpol_raw1[solpol_raw1['Data_Type']==8]
    bias_volt['Values'] = bias_volt['Values'].astype(float)
    bias_volt = bias_volt.reset_index()
    
    # DC current
    I = pd.DataFrame(solpol_raw1[solpol_raw1['Data_Type']==9])
    I['Values'] = I['Values'].astype(float)
    
    ########## Subtract the mean dark DC value per day of measurement #########
    import dark
    
    dark_dr_gr = dark.dark_df.groupby(pd.Grouper(freq='D'))
    
    for date1, df in dark_dr_gr:
        dark_daily_df = df

        if date == date1.date():
            
            I_dark_daily_mean = np.mean(dark_daily_df['Values'])
            print(I_dark_daily_mean)
            
            I['Values'] = I['Values'] - I_dark_daily_mean
            break
    
    
    # date before the 2020-08-27 don't have darks, for those subtract the mean dark value from all the dark measurements
    if date < pd.to_datetime('2020-08-27').date():
            I['Values'] = I['Values'] - float(dark.I_dark_mean)
    ###########################################################################
    
    # Geometric correction
    geom = (aperture*10**(-3)/2)**2/(4*(53*10**(-2))**2) # for the normalization with solid angle, where geom = Ω/4pi
#    geom1 = math.pi/4*((4.5*10**(-3)/2)**2)
#    I = I*geom
    I = I.drop(columns='Data_Type')
    I = I.reset_index()
    
    # Labjack other output
    labjack = solpol_raw1[solpol_raw1['Data_Type']==10]
    labjack['Values'] =labjack['Values'].astype(float)
    labjack = labjack.reset_index()
    
    # AC lock-in amplifier output in 1ω = 47 KHz
    lockin1 = solpol_raw1[solpol_raw1['Data_Type']==11]
    lockin1[['Amplitude1','Phase1']] = lockin1.Values.str.split(pat=',',expand=True) 
    lockin1['Amplitude1'] = pd.to_numeric(lockin1['Amplitude1'])
    lockin1['Phase1'] = pd.to_numeric(lockin1['Phase1'])
    A1_init = 0
    lockin1['Amplitude1'] = lockin1['Amplitude1'] - A1_init
    P1_init = 0
    lockin1['Phase1'] = lockin1['Phase1'] - P1_init
    lockin1 = lockin1.reset_index()
    
    # AC lock-in amplifier output in 2ω = 94 KHz
    lockin2 = solpol_raw1[solpol_raw1['Data_Type']==12]
    lockin2[['Amplitude2','Phase2']] = lockin2.Values.str.split(pat=',',expand=True) 
    lockin2['Amplitude2'] = pd.to_numeric(lockin2['Amplitude2'])
    lockin2['Phase2'] = pd.to_numeric(lockin2['Phase2'])
    A2_init = 0
    lockin2['Amplitude2'] = lockin2['Amplitude2'] - A2_init
    P2_init = 0
    lockin2['Phase2'] = lockin2['Phase2'] - P2_init
    lockin2 = lockin2.reset_index()
    
    # Make new dataframe
    headers = ['Polarimeter position (degs)','Rotator/Polarizer Position (degs)','PEM Setting (nm)','Retardation (waves)','Wavelength Filter','ND Filter','DateTime (UTC)','Bias Voltage','mean DC','Labjack','Amplitude 1w','Phase 1w (degs)','Amplitude 2w','Phase 2w (degs)']
    temp_data = [pol_pos['Values'],rot_pos['Values'],PEM_set['Values'],retard['Values'],wave_flt['Values'],ND_flt['Values'],time['DateTime'],bias_volt['Values'],I['Values'],labjack['Values'],lockin1['Amplitude1'],lockin1['Phase1'],lockin2['Amplitude2'],lockin2['Phase2']]
    solpol_data = pd.concat(temp_data,axis=1,keys=headers)
    #### change to time['Values'] if problem with datetime ####
    
    # Set flag with aperture size
    if aperture == 4.5:
        solpol_data['Flag'] = np.ones((solpol_data.index.size), dtype=int)
    elif aperture == 5.5:
        solpol_data['Flag'] = 2*np.ones((solpol_data.index.size), dtype=int)
    elif aperture == 7:
        solpol_data['Flag'] = 3*np.ones((solpol_data.index.size), dtype=int)
    else: 
        solpol_data['Flag'] = 2*np.ones((solpol_data.index.size), dtype=int)
    
    # Plot the raw data before the initial reduction
    print(solpol_data)
    
    solpol_list.append(solpol_data)
    
    # Data reduction, select every fifth measurement for the same polarizer position
    solpol_data_rdc = solpol_data.iloc[4::5,:]
    
    ######### Plot raw DC signal per measurement in dif. colours ##############
    temp_I = pd.concat([solpol_data_rdc['DateTime (UTC)'].to_frame(),solpol_data_rdc['mean DC'].to_frame()],axis=1)
    temp_I = temp_I.set_index([temp_I['DateTime (UTC)']]).drop(columns=['DateTime (UTC)'])

    ### !! convert to UTC here, for older Antik. data !!
    if date < pd.to_datetime('2020-01-01').date():
        temp_I = temp_I.tz_localize('Europe/Athens')
        temp_I = temp_I.tz_convert(None)
    
    date_str = date.strftime("%Y-%m-%d")
    
    plt.plot(temp_I['mean DC'],'*',markersize=18)
    plt.xlabel('DateTime (UTC)')
    plt.ylabel('I (Volts)')
    plt.xlim([pd.to_datetime(date_str + ' 04:30:00'), pd.to_datetime(date_str + ' 17:30:00')])
    plt.ylim(0,6) # change limits if dark measurements
    plt.title('Raw DC voltage - Regular measurement')
    plt.xticks(rotation=30)
    plt.tight_layout()
    jj += 1 
    ###########################################################################

    # !! Save to csv
#    solpol_data_rdc.to_csv(str(jj) +'dark_volt_17102020.csv', sep=',')


    ############ Circular and Linear Polarization calculation #############

### 1w corrected output, no phase modulated (Circular)
#    ch1_eff = 0.7342 # pre-defined channel 1 modulation efficiency
#    A1_Ch1 = (solpol_data_rdc['Amplitude 1w']/solpol_data_rdc['mean DC'])/ch1_eff

    ### !! in order to later divide with mean I, use the following expressions
    J1_A = 0.519153 # Bessel function J1(A) for A = 2.4048
    ch1_eff = np.sqrt(2)/J1_A
    A1_Ch1 = solpol_data_rdc['Amplitude 1w']*ch1_eff
    ###
    
    # Averaging for all polarizer positions - for each PEM rotator position
    temp2 = A1_Ch1.index.values
    idx = np.array([temp2, temp2, temp2, temp2]).T.flatten()[:len(temp2)]
    
    A1_Ch1_ave = A1_Ch1.groupby(idx).mean()
    A1_Ch1_ave = A1_Ch1_ave.reset_index()
    A1_Ch1_ave = A1_Ch1_ave.drop(columns=['index'])
    A1_Ch1_ave.columns = ['Average_Values']
    
    time_ave_gr = time_ave.iloc[3::4,:]
    time_ave_gr = time_ave_gr.reset_index()
    time_ave_gr = time_ave_gr.drop(['index'],axis=1)
    
    A1_Ch1_ave = pd.concat([A1_Ch1_ave['Average_Values'],time_ave_gr['DateTime']],axis=1)
    print(A1_Ch1_ave)
    
    # Phase, multiplication by cos
    A1_Ch1_phs = A1_Ch1*np.cos(np.deg2rad(solpol_data_rdc['Phase 1w (degs)']))

    circ_pol = pd.DataFrame(phase_ave(A1_Ch1_phs))
    circ_pol = circ_pol * 0.99027 # sin2a = 0.99027

#    # 2w corrected output, no phase modulated (Linear)
#    ch2_eff = 0.6106 # pre-defined channel 2 modulation efficiency
#    A2_Ch1 = (solpol_data_rdc['Amplitude 2w']/solpol_data_rdc['mean DC'])/ch2_eff 

    ### !! In order to later divide with mean I, use the following expressions
    J2_A = 0.431751 # Bessel function J2(A) for A = 2.4048
    ch2_eff = np.sqrt(2)/J2_A
    A2_Ch1 = solpol_data_rdc['Amplitude 2w']*ch2_eff
    ### 
    
    A2_Ch1_ave = A2_Ch1.groupby(idx).mean()
    A2_Ch1_ave = A2_Ch1_ave.reset_index()
    A2_Ch1_ave = A2_Ch1_ave.drop(columns=['index'])
    A2_Ch1_ave.columns = ['Average_Values']
    
    A2_Ch1_ave = pd.concat([A2_Ch1_ave['Average_Values'],time_ave_gr['DateTime']],axis=1)
    print(A2_Ch1_ave)
    
    # Phase, multiplication by cos
    A2_Ch1_phs = A2_Ch1*np.cos(np.deg2rad(solpol_data_rdc['Phase 2w (degs)']))
    lin_pol = pd.DataFrame(phase_ave(A2_Ch1_phs))
    lin_pol = lin_pol * 0.99027 # sin2a = 0.99027
    
    # Last averaging per 2 pairs
    circ_pol_ave = ((circ_pol + circ_pol.shift(-1)) / 2)[::2]
    lin_pol_ave = ((lin_pol + lin_pol.shift(-1)) / 2)[::2]
    
    time_ave_ave = time_ave_gr.iloc[1::2,:]
    time_ave_ave = time_ave_ave.reset_index()
    time_ave_ave = time_ave_ave.drop(['index'],axis=1)

    if date < pd.to_datetime('2020-01-01').date():
        # UTC time convertion for DOLP plotting when data not in UTC
        time_ave_gr = time_ave_gr.set_index(['DateTime'])
        time_ave_gr = time_ave_gr.tz_localize('Europe/Athens')
        time_ave_gr = time_ave_gr.tz_convert(None)
        time_ave_gr = time_ave_gr.reset_index()
        
        # UTC time convert for the earlier Antikythera data
        time_ave_ave = time_ave_ave.set_index(['DateTime'])
        time_ave_ave = time_ave_ave.tz_localize('Europe/Athens')
        time_ave_ave = time_ave_ave.tz_convert(None)
        time_ave_ave = time_ave_ave.reset_index()
    
    #### Test plots
    # cnt += 1
    # plt.figure(cnt)
    temp_lin = lin_pol_ave.iloc[0:len(time_ave_ave)]
    temp_lin.rename(columns={0:'Lin_Pol'},inplace=True)
    temp_lin = temp_lin.reset_index()
    temp_lin = temp_lin.drop(['index'],axis=1)
    # plt.plot(time_ave_gr,temp_lin,'bo')
    
    temp_circ = circ_pol_ave.iloc[0:len(time_ave_ave)]
    temp_circ.rename(columns={0:'Circ_Pol'},inplace=True)
    temp_circ = temp_circ.reset_index()
    temp_circ = temp_circ.drop(['index'],axis=1)
    # plt.plot(time_ave_gr,temp_circ,'ro')
    
    # concat with time
    lin_pol_df = pd.concat([temp_lin['Lin_Pol'],time_ave_ave['DateTime']],axis=1)
    lin_pol_df = lin_pol_df.set_index('DateTime')
    
    circ_pol_df = pd.concat([temp_circ['Circ_Pol'],time_ave_ave['DateTime']],axis=1)
    circ_pol_df = circ_pol_df.set_index('DateTime')
    
    # lists of all the averaged lin_pol/circ_pol dataframes
    b.append(lin_pol_df)
    c.append(circ_pol_df)
    
    ################## Calculation of the SZA and EVPA ########################
    # Electric Vector Polarization Angle (EVPA) for Q/I and U/I
    #!!! recheck EVPA
    A2_Ch1_phs_df = pd.DataFrame(A1_Ch1_phs)
    A2_Ch1_phs_df = A2_Ch1_phs_df.rename(columns={0:'Stokes_Values'})
    temp_stokes = pd.concat([solpol_data_rdc['Polarimeter position (degs)'],A2_Ch1_phs_df['Stokes_Values'],time_ave['DateTime'],solpol_data_rdc['Flag']],axis=1)
    
    # Q/I, for the 0 degs PEM position
    zero_degs = temp_stokes.loc[temp_stokes['Polarimeter position (degs)']==0]
    zero_degs = zero_degs.reset_index()
    zero_degs = zero_degs.drop(columns='index')
    
    # U/I, for the 45 degs PEM position
    ffive_degs = temp_stokes.loc[temp_stokes['Polarimeter position (degs)']==45]
    ffive_degs = ffive_degs.reset_index()
    ffive_degs = ffive_degs.drop(columns='index')
    
    # UTC conversion of ffive_degs (for the older Antik. data)
    if date < pd.to_datetime('2020-01-01').date():
        ffive_degs = ffive_degs.set_index(['DateTime'])
        ffive_degs = ffive_degs.tz_localize('Europe/Athens')
        ffive_degs = ffive_degs.tz_convert(None)
        ffive_degs = ffive_degs.reset_index()
        
    # keep the same number of measurements in 0 and 45 degrees
    if len(zero_degs) > len(ffive_degs):
        zero_degs = zero_degs.iloc[0:len(ffive_degs)]
    elif len(zero_degs) < len(ffive_degs):
        ffive_degs = ffive_degs.iloc[0:len(zero_degs)]
    
    # EVPA calculation
    idx1 = zero_degs.loc[zero_degs['Stokes_Values'] == 0]
    
    EVPA_df = pd.DataFrame(index=ffive_degs['DateTime'],columns=range(1))
    EVPA_df.rename(columns={0:'EVPA'},inplace=True)
    # EVPA_df.where(zero_degs.loc[zero_degs['Stokes_Values'] == 0],'nan')
    
    # for index in EVPA_df.index:
    for i in pd.Series(range(1,len(zero_degs))):
        if zero_degs['Stokes_Values'][i] == 0:
            EVPA_df['EVPA'][i]= math.nan
        else:
            EVPA_df['EVPA'][i] = ma.arctan(ffive_degs['Stokes_Values'][i]/zero_degs['Stokes_Values'][i])/2
    EVPA_df['EVPA'] = EVPA_df['EVPA'].astype(float)
    EVPA_df['EVPA (degs)'] = np.rad2deg(EVPA_df['EVPA'])
    
    ##### Solar Zenith Angle (SZA) calcuation #########
    SZA_df1 = pvlib.solarposition.get_solarposition(time_ave_ave['DateTime'], lat, long, alt)
    # plt.plot(SZA_df1['zenith'])
    SZA_df2 = pvlib.solarposition.get_solarposition(ffive_degs['DateTime'], lat, long, alt)
    #############
    
    EVPA_df = pd.concat([EVPA_df,SZA_df2['zenith']],axis=1)
        
    EVPA_df2 = EVPA_df.iloc[3::4,:]
    EVPA_df2 = pd.concat([EVPA_df2,lin_pol_df['Lin_Pol']],axis=1)
    
    # plt.plot(SZA_df2['zenith'])
    # plt.plot(SZA_df2['zenith'],EVPA_df['EVPA (degs)'],'mo')
    k.append(EVPA_df)
    l.append(SZA_df1)
    f.append(EVPA_df2)
    
    # List containing phase modulated signal against Q/I and U/I, plus lin. pol
    temp_phs1 = pd.DataFrame(A1_Ch1_phs).set_index(solpol_data_rdc['DateTime (UTC)'])
    temp_phs1 = temp_phs1.rename(columns={0:'Phase cor. signal'})
    
    temp_phs2 = pd.concat([solpol_data_rdc['Polarimeter position (degs)']])
    temp_phs2 = pd.DataFrame(temp_phs2).set_index(solpol_data_rdc['DateTime (UTC)'])
    temp_phs = pd.concat([temp_phs1,temp_phs2],axis=1)
    
    if len(lin_pol) > len(time_ave_gr):
        lin_pol = lin_pol.iloc[0:len(time_ave_gr)] # edw kati paizei me to linpol argotera des to !!!
    
    stokes_df = lin_pol.set_index(time_ave_gr['DateTime'])
    stokes_df = stokes_df.rename(columns={0:'Stokes'})
    
    ### U -> 0 degrees, -Q -> 45 degrees (swapped)
    U = stokes_df.iloc[0::2,:]
    U = U.rename(columns={'Stokes':'U/I'})
    
    Q = -stokes_df.iloc[1::2,:]
    Q = Q.rename(columns={'Stokes':'Q/I'})
    
    # (swapped Q <-> U) if measurement sequence stops before finishing the entire 45 degs sequence then cut the previous Q value in order to match
    if len(U) > len(Q):
        U = U.iloc[0:len(Q)]
    
    last_df = pd.concat([temp_phs,U],axis=1)
    last_df = pd.concat([last_df,Q],axis=1)
#    updated_last_df = last_df.drop([U.index[-1],last_df.index[-1]])
    temp_flag = pd.DataFrame(temp_stokes['Flag']).set_index(temp_stokes['DateTime'])
    last_df = pd.concat([last_df,temp_I,temp_flag],axis=1)
    last_df = last_df.loc[last_df.index <= Q.index[-1]]
    
    ## !!!! careful: if data files don't have the same naming it doesn't sort over ascending index value
    last_lst.append(last_df)

# print max, min mean I values
I_df = pd.DataFrame()
for meandc in last_lst:
    I_df = I_df.append(meandc)

print(max(I_df['mean DC']))
print(min(I_df['mean DC']))
print(np.mean(I_df['mean DC']))

## ============ Group data per day & dataframes with all the Lin/Circ pol, EVPA ============= ##
## ========================================================================================== ##

# initialize dataframes
b_df = pd.DataFrame()
c_df = pd.DataFrame()
k_df = pd.DataFrame()
l_df = pd.DataFrame()
f_df = pd.DataFrame()

for df1 in b:
    b_df = b_df.append(df1)

for df2 in c:
    c_df = c_df.append(df2)
    
for df3 in k:
    k_df = k_df.append(df3)

for df4 in l:
    l_df = l_df.append(df4)

for df5 in f:
    f_df = f_df.append(df5)

#solpol_all_df = pd.concat([b_df['Lin_Pol'],c_df['Circ_Pol'],l_df['zenith']],axis=1)
solpol_all_df = pd.concat([b_df['Lin_Pol'],c_df['Circ_Pol']],axis=1)
solpol_all_df["Lin_Pol"] = solpol_all_df["Lin_Pol"]*(10**6)
solpol_all_df["Circ_Pol"] = solpol_all_df["Circ_Pol"]*(10**6)
solpol_all_df = solpol_all_df.sort_index()  # in order to subtract with the correct values each time

# grouped data
solpol_grouped = solpol_all_df.groupby(pd.Grouper(freq='D'))


## ============ Plot lin/circ smoothed dataset through Chauvenet criterion ===================== ##
## =============================================================================================== ##

# initialize parameters
solpol_perday_df = pd.DataFrame()
solpol_perday_smoothed1 = pd.DataFrame()
solpol_perday_smoothed2 = pd.DataFrame()
ch_lin = []
ch_circ =[]

cnt = 1
for date1, df in solpol_grouped:

    print(len(df))

    if len(df) > 0:
        solpol_perday_df = solpol_grouped.get_group(date1)
        
        ## write into file
        d = pd.to_datetime(date1).date()
        fname = "polirisT2_{curr_date}_Antik".format(curr_date = d) #change sufix according to location
#        df.to_csv(fname, sep=',')

        cnt += 1
        
        # plot the initial datasets
        pol_plot(fname,cnt,solpol_perday_df['Lin_Pol'],solpol_perday_df['Circ_Pol'],date1)
        
        # Apply the chauvenet criterion per 1 hr
        solpol_perday_smoothed1 = chauvenet(solpol_perday_df['Lin_Pol'],solpol_perday_df)
        solpol_perday_smoothed2 = chauvenet(solpol_perday_df['Circ_Pol'],solpol_perday_df)
        
        ch_lin.append(solpol_perday_smoothed1)
        ch_circ.append(solpol_perday_smoothed2)
        
        ## write into file
        fname2 = "polirisT2_{curr_date}_Antik_smoothed".format(curr_date = d)
#        solpol_perday_smoothed1.to_csv(fname2, sep=',')
        
        # plot the smoothed datasets
        pol_plot(fname2,cnt,solpol_perday_smoothed1['Lin_Pol'],solpol_perday_smoothed2['Circ_Pol'],date1)
        
    else: 
        print('Not existing measurements during this date')
# plt.close(fig=None)


### ============= Plot initial EVPA per SZA =============== ##

EVPA_grouped = k_df.groupby(pd.Grouper(freq='D'))

cnt1 =1
for date2,val in EVPA_grouped:
  
    if len(val) > 0:
        EVPA_perday_df = EVPA_grouped.get_group(date2)  
    
        ## write into file
        d1 = pd.to_datetime(date2).date()
        fname1 = "EVPA_{curr_date}_Antik".format(curr_date = d1) #change sufix according to location: (Ath, Antik, Cy)
#        val.to_csv(fname1, sep=',')
    
        cnt1 += 1
        
        # plot EVPA
        evpa_plot(fname1,cnt1,EVPA_perday_df['zenith'],EVPA_perday_df['EVPA (degs)'])

    else: 
        print('Not existing measurements during this date')


## =============== Filtered EVPA for Lin_pol values over 50ppm ============== ##

# initialize dataframes
ch_lin_df = pd.DataFrame()
ch_circ_df = pd.DataFrame()

for df6 in ch_lin:
    ch_lin_df = ch_lin_df.append(df6)
   
# Find common indexes between the EVPA and chauvenet smoothed lin datasets
temp_ch = f_df.loc[ch_lin_df.index]
temp_ch['Lin_Pol'] = temp_ch['Lin_Pol']*10**6

# EVPA for Lin_pol above 50ppm
final_evpa_50ppm = temp_ch.loc[temp_ch['Lin_Pol'] >= 50]

# Filtered EVPA per day
EVPA_flt_grouped = final_evpa_50ppm.groupby(pd.Grouper(freq='D'))

cnt2 = 1
for date3,val1 in EVPA_flt_grouped:
  
    if len(val1) > 0:
        EVPA_flt_perday_df = EVPA_flt_grouped.get_group(date3)  
    
    ## write into file
        d3 = pd.to_datetime(date3).date()
        fname3 = "EVPA_flt_{curr_date}_Antik".format(curr_date = d3) #change sufix according to location: (Ath, Antik, Cy)
#        val1.to_csv(fname3, sep=',')
        
        cnt2 += 1
        
        # plot filtered EVPA
        evpa_plot(fname3,cnt2,EVPA_flt_perday_df['zenith'],EVPA_flt_perday_df['EVPA (degs)'])

    else: 
        print('Not existing measurements during this date')


#### ================= DOLP calculation, norm. Q, U., DOLP plotting ============== ##
## =============================================================================================== ##

last_all_df = pd.DataFrame()
for df6 in last_lst:
    last_all_df = last_all_df.append(df6)
    last_all_df=last_all_df.sort_index()

# Appply the Chauvenet criterion
# df initialization
temp_last = pd.DataFrame()
temp_dolp = pd.DataFrame()
temp_Q = pd.DataFrame()
temp_U = pd.DataFrame()
temp_smooth = pd.DataFrame()
temp_dolp_smoothed = pd.DataFrame()
norm_Q_smoothed = pd.DataFrame()
norm_U_smoothed = pd.DataFrame()
error = pd.DataFrame()

# Mean I for every set of 4 measurements, should be independent of polarizer position
flag_df = pd.concat([last_all_df['mean DC'],last_all_df['Flag']],axis=1)
mean_I = flag_df.groupby(np.arange(len(last_all_df['mean DC']))//4).mean()

# if you don't wish to include flag enable the following
#mean_I = last_all_df['mean DC'].groupby(np.arange(len(last_all_df['mean DC']))//4).mean()
#mean_I = mean_I.to_frame()
mean_I.rename(columns={'mean DC':'mean I'},inplace=True)

temp_mean = pd.DataFrame(np.repeat(mean_I.values,4,axis=0))
temp_mean.rename(columns={0:'mean I'},inplace=True)

# Relative difference of measured I(DC) from mean I 
rel_dif_I = 100*abs(last_all_df['mean DC'].reset_index(drop=True) - temp_mean['mean I']) / temp_mean['mean I']
rel_dif_I = pd.DataFrame(rel_dif_I).set_index(last_all_df.index)
rel_dif_I.rename(columns={0:'Rel. Dif (%)'},inplace=True)
print(max(rel_dif_I.values))

# Normalized Q & mean I df
temp_Q = last_all_df['Q/I']
temp_Q = temp_Q.to_frame()
temp_Q = temp_Q.dropna()

# Normalized U & mean I df
temp_U = last_all_df['U/I']
temp_U = temp_U.to_frame()
temp_U = temp_U.dropna()

# !!! Group df with normalized Q, U over mean I
temp_last = temp_Q.join(temp_U, how='outer')
#temp_last = pd.concat([temp_Q,temp_U],axis=1) # with the unsmoothed Q, U
mean_I = mean_I.set_index(temp_last.index)

# find duplicated indexes, which should not exist
indx_dupl = temp_last.index.duplicated()

if indx_dupl == 'True':
    temp_last = temp_last[~temp_last.index.duplicated(keep='first')] # delete duplicate indexes
    mean_I = mean_I[~mean_I.index.duplicated(keep='first')]
    
    temp_last = pd.concat([temp_last,mean_I],axis=1)
    
else:
    temp_last = pd.concat([temp_last,mean_I],axis=1)

# divide by mean I
temp_last['Q/I'] = temp_last['Q/I'] / (2*temp_last['mean I'])
norm_Q = temp_last['Q/I'].to_frame().dropna().reset_index(drop=True)
norm_Q_indx = temp_last['Q/I'].to_frame().dropna()

temp_last['U/I'] = temp_last['U/I'] / (2*temp_last['mean I'])
norm_U = temp_last['U/I'].to_frame().dropna().reset_index(drop=True)

# chauvenet smoothed Q/I & U/I
norm_Q_smoothed = chauvenet(temp_last['Q/I'].dropna(),norm_Q_smoothed).drop(['Chauv_120min'],axis=1)
norm_U_smoothed = chauvenet(temp_last['U/I'].dropna(),norm_U_smoothed).drop(['Chauv_120min'],axis=1)
#temp_smooth = norm_Q_smoothed.merge(norm_U_smoothed)
temp_smooth = pd.concat([norm_Q_smoothed,norm_U_smoothed],axis=1)

###### Degree of Linear Polarization (DOLP) calculation #########
temp_dolp['DOLP'] = np.sqrt(norm_Q['Q/I']**2 + norm_U['U/I']**2) 
temp_dolp = temp_dolp.set_index(norm_Q_indx.index)
#temp_dolp = temp_dolp.set_index(temp_Q.index) 
#temp_dolp = temp_dolp.set_index(solpol_all_df.index)

# DOLP stadard error calculation, per 12 hrs
error['Stdev'] = temp_dolp['DOLP'].resample('2H', label='right', closed='right').std()
error['Cnt'] = temp_dolp['DOLP'].resample('2H', label='right', closed='right').count()
error['StdError'] = error['Stdev'] / np.sqrt(error['Cnt'])

std_error = error.loc[np.repeat(error.index.values,error['Cnt'])]
std_error = std_error.reset_index(drop=True)
std_error = std_error.set_index(temp_dolp.index)

# concat to DOLP dataframe
temp_dolp['StdError'] = std_error['StdError']

# smoothed dataset through Chauvenet
temp_dolp_smoothed = chauvenet(temp_dolp['DOLP'],temp_dolp)
temp_dolp_smoothed = temp_dolp_smoothed.drop(['Chauv_120min'],axis=1)
temp_dolp_smoothed = temp_dolp_smoothed.join(temp_dolp['StdError'])

# for the unsmoothed dataset enable this
#last_all_df1 = pd.concat([temp_last,temp_dolp_smoothed],axis=1)
#last_all_df_gr = last_all_df1.groupby(pd.Grouper(freq='D'))

# for the smoothed dataset enable this
last_all_df1_smooth = temp_smooth.join(temp_dolp_smoothed)
last_all_df1_smooth = pd.concat([last_all_df1_smooth, mean_I], axis=1)

SZA_all = pvlib.solarposition.get_solarposition(last_all_df1_smooth.index, lat, long, alt)

last_all_df1_smooth = pd.concat([last_all_df1_smooth, SZA_all['zenith']], axis=1)
last_all_df1_smooth.rename(columns={'zenith':'SZA'},inplace=True)


# Group by day & Plot
last_all_df_gr = last_all_df1_smooth.groupby(pd.Grouper(freq='D'))

for datey,dfy in last_all_df_gr:
    
    try:
        last_day_df = last_all_df_gr.get_group(datey)
    except KeyError:
        print(f"Date {datey} not found, skipping")
        continue
    
#    temp1000 = np.zeros(len(last_day_df.index))
#    temp1000 = pd.DataFrame(temp1000).set_index(last_day_df.index)
    
    plt.figure(i)
    myFmt = mdates.DateFormatter('%d-%m-%y %H:%M')
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(myFmt)
    ax.xaxis_date()
    ax.grid()
    
#    plt.plot(last_day_df['Phase cor. signal'],'ro', label = 'Phase cor. signal')
    plt.plot(last_day_df['DOLP']*10**6,'bo', label = 'DOLP') #/geom*10**6 this without the scaling factor,  last_day_df['DOLP']*10**6 - 2*
#    plt.errorbar(last_day_df['DOLP'].index,last_day_df['DOLP']*10**6, yerr = [2*last_day_df['StdError']*10**6,2*last_day_df['StdError']*10**6], ecolor=(0.8, 0.8, 0.8),elinewidth=2, capsize=5,capthick=2)

#    plt.errorbar(last_day_df['DOLP'].index,last_day_df['DOLP']*10**6, yerr = [last_day_df['StdError']*10**6, 2*(last_day_df['StdError']*10**6)], ecolor=(0.8, 0.8, 0.8),capthick=2)
    plt.plot(last_day_df['Q/I']*10**6,'ko',label = 'Q/I') #/geom*10**6
    plt.plot(last_day_df['U/I']*10**6,'mo',label = 'U/I') #/geom*10**6
    plt.xlabel('Time (UTC)') #for Athens, Cyprus & several Antik.
#    plt.plot(temp1000,'r')
    plt.ylim((-0.0002*10**6,0.0002*10**6))
#    plt.ylim((-0.001,0.001))
#    plt.xlim(datey.replace(hour=4,minute=30,second=0),datey.replace(hour=17,minute=30,second=0))
    plt.xlim(datey.replace(hour=6,minute=30,second=0),datey.replace(hour=19,minute=30,second=0))  # for Mindelo data
    plt.xticks(rotation=30)
    plt.tight_layout()
    plt.legend()
    i+=1
    
    # save plots
    fname4 = "dolp_{curr_date}_{curr_location}".format(curr_date = datey.strftime("%Y-%m-%d"),curr_location = location) #change sufix according to location   
    
    plt.savefig(fname4 + '.png')
    plt.savefig(fname4 + '.eps',format ='eps')
    
#    dfy.to_csv(fname4 + '.csv', sep=',') # save to file
    

###============ Plot multiple DOLP "instances" in one plot ============######       
#dolp_minus_dark = pd.DataFrame()
#dolp_minus_dark = last_day_df['DOLP']*10**6
#dolp_minus_dark = pd.DataFrame(dolp_minus_dark)
#dolp_minus_dark.rename(columns={'DOLP':'DOLP_regular'},inplace=True)

#dolp_minus_dark = pd.concat([last_day_df['DOLP']*10**6,dolp_minus_dark],axis=1)
#dolp_minus_dark.rename(columns={'DOLP':'DOLP_Imean'},inplace=True)

#plt.plot(dolp_minus_dark['DOLP_regular'],'ko',label = 'Regular DOLP')
#plt.plot(dolp_minus_dark['DOLP_Imax'],'mo',label = 'DOLP - max. DC dark')
#plt.plot(dolp_minus_dark['DOLP_Imean'],'bo',label = 'DOLP - mean DC dark')
#plt.plot(dolp_minus_dark['DOLP_Imin'],'ro',label = 'DOLP - min. DC dark')
#plt.ylim(0,200)
#plt.xticks(rotation=30)
#plt.legend()
    
