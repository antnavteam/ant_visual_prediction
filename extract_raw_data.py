#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 13:23:26 2024

@author: Oceane DAUZERE-PERES
"""

import csv
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.api as sm

# calibration values (in mm per sensor reading)
calib_0 = 0.309570687371
calib_1 = 0.256095062487
# ball diameter (in mm)
ball_diameter = 48.7
ball_perimeter = ball_diameter*np.pi

def clean(input):
    """clean the raw data in a temporary file by removing the coma (",") between each 3 digits of the numbers """
    tmpFile = "tmp.csv"
    with open(input, "r") as file, open(tmpFile, "w",newline='') as outFile:
        reader = csv.reader(file, delimiter=";")
        writer = csv.writer(outFile, delimiter=";")
        header = next(reader)
        writer.writerow(header)
        for row in reader:
            colValues=row
            for col in range(0,8):
                colValues[col]="".join([n for n in row[col] if n != ","])
            writer.writerow(colValues)
    outFile.close()
    return(tmpFile)


def trackball_data_cruncher (idt, trackball_raw_data_4column, time_stamps, ball_diameter, calib_0, calib_1):
    """transform the raw data detected by the lasers sensors into coordinates and velocity signals"""
    end=len(trackball_raw_data_4column)
    trackball_rel_data = [[trackball_raw_data_4column[i][j]-trackball_raw_data_4column[i-1][j] for i in range(2,end)]for j in range(4)]
    time_stamps_rel_sec =  [(time_stamps[i] - time_stamps[i-1])/1000 for i in range (1,end-1)] # put in second (was in ms)

    X0=[a*calib_0 for a in trackball_rel_data[:][0]]
    X1=[a*calib_1 for a in trackball_rel_data[:][1]]
    Xs= [(X0[i] + X1[i]) / 2 for i in range(len(X0))] # Make average of X0 and X1 (they should get same signal)
    a= np.corrcoef(X0,X1)[0,1]
    if a<0.9 :
        print('Correlation X0-X1 =', a, idt)

# Ang_speed based on time stamps
    ball_rot_frame = [a/ball_perimeter for a in Xs] # in ball rotation per frame
    rad_frame = [a*2*np.pi for a in ball_rot_frame]
    radsec = [rad_frame[i]/time_stamps_rel_sec[i] for i in range(len(rad_frame))]
    Ang_speed_degsec = [(a/np.pi)*180 for a in radsec] # Put them in ball_rotation /sec

# Deal with the Ys translation to get speed
    Y0=[a*calib_0 for a in trackball_rel_data[:][2]] # (calib in mm)
    Y1=[a*calib_1 for a in trackball_rel_data[:][3]]

    cm_frame = [np.sqrt(Y0[i]**2 + Y1[i]**2)/10 for i in range(len(Y0))] # cm/frame
    Fwd_speed_cm_sec = [cm_frame[i]/time_stamps_rel_sec[i] for i in range(len(cm_frame))] #cm/sec

# Get the Path
    Theta_path = np.cumsum(rad_frame) # in radian (per frame)

    dx = [cm_frame[i]*np.cos(Theta_path[i]) for i in range(len(cm_frame))] #theta be in radian
    dy = [cm_frame[i]*np.sin(Theta_path[i]) for i in range(len(cm_frame))] #dr (derivative of movement = speed)

    X=np.cumsum(dx)-dx[0]
    Y=np.cumsum(dy)-dy[0]

    return (Ang_speed_degsec, Fwd_speed_cm_sec ,X, Y, Theta_path)

### EXPERIMENT GAIN ########################################################################################################################################
# This part extract the raw data to obtain dictionnaries corresponding to our parameters of interest for each trial as well as the beginning and  the end of each condition inside each trial. 
# idt correspond to the name, or "identity", of a singular trial and is used as the index inside all the dictionnaries to access the data of this specific trial.
# at this end of this part a csv table will containing the values used for the statictical analysis will be created. 
# YOU  HAVE TO RUN THIS CODE LAST IF YOU WANT TO VIZUALIZE THE DATA CORRESPONDING TO THIS EXPERIMENT (such as mean oscillation cycles, trajectories or angular velocity signals)

mainfolder = '' ### WRITE THE PATH TO THE FILE WITH RAW DATA CORRESPONDING TO THE MODULATION OF THE GAIN EXPERIMENT
list_idt=[]
ant={}
eye={}
order={}
cookie={}
condition={}
Ang_speed_degsec={}
Ang_speed_smooth={}
Fwd_speed_cm_sec={}
X={}
Y={}
Theta_path={}
time_stamps_raw={}
time_stamps_sec={}
frequence={}
Ang_acceleration={}
Fwd_acceleration={}
events_start={}
events_end={}
events={}


for f in glob.glob(mainfolder+'*Ant*'):
    os.chdir(f)

    # COORDS
    for file in glob.glob(f+'/*coords*'):
        filename=file[-64:-43]
        idt=filename[3:8]
        list_idt+=[idt]
        eye[idt]=filename[6:8]
        ant[idt]=filename[3:5]
        order[idt]=filename[12:18]
        condition[idt]=filename[9:11]
        cookie[idt]=filename[19:21]
        fclean=clean(file)
        clean_file=f+"/"+fclean
        with open(clean_file,"r") as csvfile :
            reader=csv.reader(csvfile, delimiter=';')
            time_stamps=[]
            i=0
            for row in reader:
                if i==0:
                    trackball_data_raw = [row[3:7]]
                else:
                    trackball_data_raw += [[float(a) for a in row[3:7]]]
                    time_stamps +=  [float(row[7])]
                i+=1
        csvfile.close()
        os.remove(clean_file)
        
        #speed and position
        Ang_speed_degsec[idt], Fwd_speed_cm_sec[idt], X[idt], Y[idt], Theta_path[idt] = trackball_data_cruncher (idt,trackball_data_raw, time_stamps, ball_diameter, calib_0, calib_1)
        #time stamps
        time_stamps_raw[idt] =  [a/1000 for a in time_stamps]# in sec
        print(idt)
        # resample a data point every 33 milliseconds
        df = pd.DataFrame(np.array([Fwd_speed_cm_sec[idt], Ang_speed_degsec[idt],X[idt],Y[idt],time_stamps_raw[idt][:len(time_stamps_raw[idt])-1]]).T, columns=["Fwd_speed","Ang_speed","X","Y","Time"],index=pd.to_datetime(time_stamps_raw[idt][:len(time_stamps_raw[idt])-1], unit='s'))
        df2 = df.resample('33L').mean().reset_index().interpolate(method='ffill')
        time_stamps_sec[idt] = df2.Time.tolist()
        X[idt] = df2.X.tolist()
        Y[idt] = df2.Y.tolist()
        #smoothing whith a window of 20 frames (1 time fwd speed, 2 times angular speed)
        Fwd_speed_cm_sec[idt] = df2.Fwd_speed.rolling(20,min_periods=1).mean().tolist()
        Ang_speed_degsec[idt] = df2.Ang_speed.tolist()
        df2['Ang_speed_smooth'] = df2.Ang_speed.rolling(20,min_periods=1).mean().tolist()
        Ang_speed_smooth[idt] = df2.Ang_speed_smooth.rolling(20,min_periods=1).mean().tolist()
        #detection error for ant 05-RE --> removing the erronous part of the signal
        if idt=="05_RE":
            Ang_speed_smooth[idt]=Ang_speed_smooth[idt][:11091]+Ang_speed_smooth[idt][11121:]
            time_stamps_sec[idt]=time_stamps_sec[idt][:11091]+time_stamps_sec[idt][11121:]
            Fwd_speed_cm_sec[idt]=Fwd_speed_cm_sec[idt][:11091]+Fwd_speed_cm_sec[idt][11121:]
        #angular acceleration
        Ang_acceleration[idt]= [(Ang_speed_smooth[idt][i]-Ang_speed_smooth[idt][i-1])/(time_stamps_sec[idt][i]-time_stamps_sec[idt][i-1]) for i in range (1,len(time_stamps_sec[idt])-1)]
        #forward acceleration
        Fwd_acceleration[idt]= [(Fwd_speed_cm_sec[idt][i]-Fwd_speed_cm_sec[idt][i-1])/(time_stamps_sec[idt][i]-time_stamps_sec[idt][i-1]) for i in range (1,len(time_stamps_sec[idt]))]

    #EVENT (obtain start and end of the different conditions of a trial)
    for filename in glob.glob('*events*'):
        with open(os.path.join(os.getcwd(), filename)) as event_file:
            d=0
            idt=filename[3:8]
            events_start[idt]=[]
            events_end[idt]=[]
            events[idt]=[]
            reader = csv.reader(event_file, delimiter=";")
            header = next(reader)
            for row in reader:
                if row[1] == "view_cleared":
                    events_start[idt]+=["".join([n for n in row[0] if n !=","])]
                if row[1] == "view_filled_RGBA(0.000, 0.000, 0.000, 1.000)":
                    if len(events_start[idt]) > len(events_end[idt]):
                        events_end[idt]+=["".join([n for n in row[0] if n !=","])]
                    else:
                        events_start[idt]+=["".join([n for n in row[0] if n !=","])]
                        events[idt]+=["B"]
                if "rot_gain" in row[1]:
                    events[idt]+=["".join([n for n in row[1][16:]])]


### Writing the new file used for statistical analysis 
os.chdir(r"") ### WRITE THE PATH TO THE FILE WHERE YOU WANT THE NEW FILE TO BE CREATED
Data_file = "expgain_stat.csv"
Newfile = open(Data_file, "w",newline='')
writer = csv.writer(Newfile, delimiter=";")
writer.writerow(["Ant","Eye","Trial","Gain","Cookie","Order","Repetition","Mean angular speed","blocked oscillation","fourrier","nb_osci","nb_block"])


starts={}
ends={}
j=0

for idt in list_idt:
    j+=1
    starts[idt]=[]
    ends[idt]=[]

    for e in range (len(events[idt])):
        start = round(time_stamps_raw[idt][int(events_start[idt][e])]/0.033)+20
        end = round(time_stamps_raw[idt][int(events_end[idt][e])]/0.033)
        starts[idt]+=[start]
        ends[idt]+=[end]
        max=[]
        min=[]
        # Time spent turning left or right
        time_left= len([a for a in Ang_speed_smooth[idt][start:end]  if a>3])*0.033 #margin to avoid freezing period
        time_right= len([a for a in Ang_speed_smooth[idt][start:end] if a<-3])*0.033 #margin to avoid freezing period
    
        #frequency of direction switches
        nb_changes=0
        for i in range(start+1,end):
           if (Ang_speed_smooth[idt][i]<=1 and Ang_speed_smooth[idt][i-1]>=1): #baseline at 1 instead of 0 to avoid counting freezing as switches
                nb_changes+=1
        frequency0= (nb_changes)/(time_left+time_right)

        #frequency of oscillations
        nb_pics=0
        a=0
        slope= Ang_speed_smooth[idt][start+1]-Ang_speed_smooth[idt][start]
        for i in range(start+1,end):
            if (Ang_speed_smooth[idt][i+1]-Ang_speed_smooth[idt][i])*slope<0:
                if np.abs(Ang_speed_smooth[idt][i]-a)>10: #if oscillation amplitude is enough
                    nb_pics+=1
                a=Ang_speed_smooth[idt][i]
            slope= Ang_speed_smooth[idt][i+1]-Ang_speed_smooth[idt][i]
        frequency2=nb_pics/(2*(time_left+time_right))
        nb_osci=round(nb_pics/2)
        
        #FOURRIER
        autocor = sm.tsa.acf([a for a in Ang_speed_smooth[idt][start:end] if a!=0], nlags=len(Ang_speed_smooth[idt][start:end]))
        fourierTransform = np.fft.rfft(autocor)
        frequency = np.linspace(0, 1/(2*0.033), len(np.square(abs(fourierTransform))))
        argfourrier=np.argmax(np.square(abs(fourierTransform[12:])))+12
        fourrier=frequency[argfourrier]
        
        p_blocked=(frequency2-frequency0)/frequency2

        writer.writerow([ant[idt],eye[idt],j,events[idt][e],cookie[idt],e+1,condition[idt],
                         np.mean([np.abs(s) for s in Ang_speed_smooth[idt][start:end] if s!=0]),
                         p_blocked,fourrier,nb_changes,nb_osci-nb_changes])
Newfile.close()


### EXPERIMENT VISUAL STRUCTURE ########################################################################################################################################
# This part extract the raw data to obtain dictionnaries corresponding to our parameters of interest for each trial as well as the beginning and  the end of each condition inside each trial. 
# idt correspond to the name, or "identity", of a singular trial and is used as the index inside all the dictionnaries to access the data of this specific trial.
# at this end of this part a csv table will containing the values used for the statictical analysis will be created. 
# YOU  HAVE TO RUN THIS CODE LAST IF YOU WANT TO VIZUALIZE THE DATA CORRESPONDING TO THIS EXPERIMENT (such as mean oscillation cycles, trajectories or angular velocity signals)

mainfolder = '' ### WRITE THE PATH TO THE FILE WITH RAW DATA CORRESPONDING TO THE VISUAL STRUCTURE EXPERIMENT

list_idt=[]
ant={}
cookie={}
condition={}
Ang_speed_degsec={}
Ang_speed_smooth={}
Fwd_speed_cm_sec={}
X={}
Y={}
Theta_path={}
time_stamps_raw={}
time_stamps_sec={}
Ang_acceleration={}
Fwd_acceleration={}
events_start={}
events_end={}
events={}

for f in glob.glob(mainfolder+'*Ant*'):
    os.chdir(f)

    # COORDS
    for file in glob.glob(f+'/*coords*'):
        filename=file[-36:]
        idt=(filename[3:5])
        list_idt+=[idt]
        ant[idt]=filename[3:5]
        cookie[idt]=filename[13:15]
        fclean=clean(file)
        clean_file=f+"/"+fclean
        with open(clean_file,"r") as csvfile :
            reader=csv.reader(csvfile, delimiter=';')
            time_stamps=[]
            i=0
            for row in reader:
                if i==0:
                    trackball_data_raw = [row[3:7]]
                else:
                    trackball_data_raw += [[float(a) for a in row[3:7]]]
                    time_stamps +=  [float(row[7])]
                i+=1
        csvfile.close()
        os.remove(clean_file)

        #speed and position
        Ang_speed_degsec[idt], Fwd_speed_cm_sec[idt], X[idt], Y[idt], Theta_path[idt] = trackball_data_cruncher (idt,trackball_data_raw, time_stamps, ball_diameter, calib_0, calib_1)
        #time stamps
        time_stamps_raw[idt] =  [a/1000 for a in time_stamps]# in sec
        print(idt)
        # resample a data point every 33 milliseconds
        df = pd.DataFrame(np.array([Fwd_speed_cm_sec[idt], Ang_speed_degsec[idt],X[idt],Y[idt],time_stamps_raw[idt][:len(time_stamps_raw[idt])-1]]).T, columns=["Fwd_speed","Ang_speed","X","Y","Time"],index=pd.to_datetime(time_stamps_raw[idt][:len(time_stamps_raw[idt])-1], unit='s'))
        df2 = df.resample('33L').mean().reset_index().interpolate(method='ffill')
        time_stamps_sec[idt] = df2.Time.tolist()
        X[idt] = df2.X.tolist()
        Y[idt] = df2.Y.tolist()
        #smoothing whith a window of 20 frames (1 time fwd speed, 2 times angular speed)
        Fwd_speed_cm_sec[idt] = df2.Fwd_speed.rolling(20,min_periods=1).mean().tolist()
        Ang_speed_degsec[idt] = df2.Ang_speed.tolist()
        df2['Ang_speed_smooth'] = df2.Ang_speed.rolling(20,min_periods=1).mean().tolist()
        Ang_speed_smooth[idt] = df2.Ang_speed_smooth.rolling(20,min_periods=1).mean().tolist()
        #angular acceleration
        Ang_acceleration[idt]= [(Ang_speed_smooth[idt][i]-Ang_speed_smooth[idt][i-1])/(time_stamps_sec[idt][i]-time_stamps_sec[idt][i-1]) for i in range (1,len(time_stamps_sec[idt])-1)]
        #forward acceleration
        Fwd_acceleration[idt]= [(Fwd_speed_cm_sec[idt][i]-Fwd_speed_cm_sec[idt][i-1])/(time_stamps_sec[idt][i]-time_stamps_sec[idt][i-1]) for i in range (1,len(time_stamps_sec[idt]))]

    #EVENT (obtain start and end of the different conditions of a trial)
    for filename in glob.glob('*events*'):
        with open(os.path.join(os.getcwd(), filename)) as event_file:
            idt=(filename[3:5])
            events_start[idt]=[]
            events_end[idt]=[]
            events[idt]=[e for e in filename[6:12]]
            reader = csv.reader(event_file, delimiter=";")
            header = next(reader)
            for row in reader:
                if row[1] == "view_cleared":
                    events_start[idt]+=["".join([n for n in row[0] if n !=","])]
                if row[1] == "view_filled_RGBA(0.000, 0.000, 0.000, 1.000)":
                    events_end[idt]+=["".join([n for n in row[0] if n !=","])]
                    

### Writing the new file used for statistical analysis 
os.chdir(r"") ### WRITE THE PATH TO THE FILE WHERE YOU WANT THE NEW FILE TO BE CREATED
Data_file = "expstructure_stat.csv"
Newfile = open(Data_file, "w",newline='')
writer = csv.writer(Newfile, delimiter=";")
writer.writerow(["Ant","Condition","Cookie","Order","Mean angular speed","blocked oscillation","fourrier","nb_osci","nb_block"])


starts={}
ends={}


for idt in list_idt:
    starts[idt]=[]
    ends[idt]=[]

    for e in range (len(events[idt])):
        start = round(time_stamps_raw[idt][int(events_start[idt][e])]/0.033)+20
        end = round(time_stamps_raw[idt][int(events_end[idt][e])]/0.033)
        starts[idt]+=[start]
        ends[idt]+=[end]
        # Time spent turning left or right
        time_left= len([a for a in Ang_speed_smooth[idt][start:end]  if a>3])*0.033 #margin to avoid freezing period
        time_right= len([a for a in Ang_speed_smooth[idt][start:end] if a<-3])*0.033 #margin to avoid freezing period
    
        #frequency of direction switches
        nb_changes=0
        for i in range(start+1,end):
           if (Ang_speed_smooth[idt][i]<=1 and Ang_speed_smooth[idt][i-1]>=1): #baseline at 1 instead of 0 to avoid counting freezing as switches
                nb_changes+=1
        frequency0= (nb_changes)/(time_left+time_right)

        #frequency of oscillations
        nb_pics=0
        a=0
        slope= Ang_speed_smooth[idt][start+1]-Ang_speed_smooth[idt][start]
        for i in range(start+1,end):
            if (Ang_speed_smooth[idt][i+1]-Ang_speed_smooth[idt][i])*slope<0:
                if np.abs(Ang_speed_smooth[idt][i]-a)>10: #if oscillation amplitude is enough
                    nb_pics+=1
                a=Ang_speed_smooth[idt][i]
            slope= Ang_speed_smooth[idt][i+1]-Ang_speed_smooth[idt][i]
        frequency2=nb_pics/(2*(time_left+time_right))
        nb_osci=round(nb_pics/2)
        
        #FOURRIER analysis to obtain frequency
        autocor = sm.tsa.acf([a for a in Ang_speed_smooth[idt][start:end] if a!=0], nlags=len(Ang_speed_smooth[idt][start:end]))
        fourierTransform = np.fft.rfft(autocor)
        frequency = np.linspace(0, 1/(2*0.033), len(np.square(abs(fourierTransform))))
        argfourrier=np.argmax(np.square(abs(fourierTransform[12:])))+12
        fourrier=frequency[argfourrier]
        
        p_blocked=(frequency2-frequency0)/frequency2

        writer.writerow([ant[idt],events[idt][e],cookie[idt],e+1,
                         np.mean([np.abs(s) for s in Ang_speed_smooth[idt][start:end] if s!=0]),
                         p_blocked,fourrier,nb_changes,nb_osci-nb_changes])
Newfile.close()

### EXPERIMENT WEIGHT BALL ########################################################################################################################################
# This part extract the raw data to obtain dictionnaries corresponding to our parameters of interest for each trial as well as the beginning and  the end of each condition inside each trial. 
# idt correspond to the name, or "identity", of a singular trial and is used as the index inside all the dictionnaries to access the data of this specific trial.
# at this end of this part a csv table will containing the values used for the statictical analysis will be created. 
# YOU  HAVE TO RUN THIS CODE LAST IF YOU WANT TO VIZUALIZE THE DATA CORRESPONDING TO THIS EXPERIMENT (such as mean oscillation cycles, trajectories or angular velocity signals)

mainfolder = '' ### WRITE THE PATH TO THE FILE WITH RAW DATA CORRESPONDING TO THE HEAVY BALL EXPERIMENT

list_idt=[]
ant={}
size1={"01":0.89,"02":1.06,"03":0.76,"04":0.71,"05":0.8,"06":0.95,"07":0.85,"08":0.91,"09":0.76,"10":0.78,"11":0.87,"12":0.685,"13":0.685,"14":0.74,"15":0.70,"16":0.87,"17":0.755,"18":0.77,"19":0.8,"20":0.69,"21":0.76,"22":0.72,"23":0.815,"24":0.81}
size2={"01":0.388,"02":0.477,"03":0.339,"04":0.3,"05":0.36,"06":0.415,"07":0.37,"08":0.4,"09":0.3,"10":0.35,"11":0.39,"12":0.31,"13":0.33,"14":0.32,"15":0.32,"16":0.4,"17":0.35,"18":0.355,"19":0.38,"20":0.29,"21":0.344,"22":0.277,"23":0.374,"24":0.32}
order={}
cookie={}
condition={}
Ang_speed_degsec={}
Ang_speed_smooth={}
Fwd_speed_cm_sec={}
X={}
Y={}
Theta_path={}
time_stamps_raw={}
time_stamps_sec={}
Ang_acceleration={}
Fwd_acceleration={}
events_start={}
events_end={}
events={}

for f in glob.glob(mainfolder+'*Ant*'):
    os.chdir(f)

    # COORDS
    for file in glob.glob(f+'/*coords*'):
        #with open(os.path.join(os.getcwd(), filename)) as file: # open in read only mode
        filename=file[-38:]
        idt=filename[3:8]
        list_idt+=[idt]
        ant[idt]=filename[3:5]
        cookie[idt]=filename[15:17]
        order[idt]=filename[12:14]
        condition[idt]=filename[6:8]
        fclean=clean(file)
        clean_file=f+"/"+fclean
        with open(clean_file,"r") as csvfile :
            reader=csv.reader(csvfile, delimiter=';')
            time_stamps=[]
            i=0
            for row in reader:
                if i==0:
                    trackball_data_raw = [row[3:7]]
                else:
                    trackball_data_raw += [[float(a) for a in row[3:7]]]
                    time_stamps +=  [float(row[7])]
                i+=1
        csvfile.close()
        os.remove(clean_file)

        #speed and position
        Ang_speed_degsec[idt], Fwd_speed_cm_sec[idt], X[idt], Y[idt], Theta_path[idt] = trackball_data_cruncher (idt,trackball_data_raw, time_stamps, ball_diameter, calib_0, calib_1)
        #time stamps
        time_stamps_raw[idt] =  [a/1000 for a in time_stamps]# in sec
        print(idt)
        # resample a data point every 33 milliseconds
        df = pd.DataFrame(np.array([Fwd_speed_cm_sec[idt], Ang_speed_degsec[idt],X[idt],Y[idt],time_stamps_raw[idt][:len(time_stamps_raw[idt])-1]]).T, columns=["Fwd_speed","Ang_speed","X","Y","Time"],index=pd.to_datetime(time_stamps_raw[idt][:len(time_stamps_raw[idt])-1], unit='s'))
        df2 = df.resample('33L').mean().reset_index().interpolate(method='ffill')
        time_stamps_sec[idt] = df2.Time.tolist()
        X[idt] = df2.X.tolist()
        Y[idt] = df2.Y.tolist()
        #smoothing whith a window of 20 frames (1 time fwd speed, 2 times angular speed)
        Fwd_speed_cm_sec[idt] = df2.Fwd_speed.rolling(20,min_periods=1).mean().tolist()
        Ang_speed_degsec[idt] = df2.Ang_speed.tolist()
        df2['Ang_speed_smooth'] = df2.Ang_speed.rolling(20,min_periods=1).mean().tolist()
        Ang_speed_smooth[idt] = df2.Ang_speed_smooth.rolling(20,min_periods=1).mean().tolist()
        #angular acceleration
        Ang_acceleration[idt]= [(Ang_speed_smooth[idt][i]-Ang_speed_smooth[idt][i-1])/(time_stamps_sec[idt][i]-time_stamps_sec[idt][i-1]) for i in range (1,len(time_stamps_sec[idt])-1)]
        #forward acceleration
        Fwd_acceleration[idt]= [(Fwd_speed_cm_sec[idt][i]-Fwd_speed_cm_sec[idt][i-1])/(time_stamps_sec[idt][i]-time_stamps_sec[idt][i-1]) for i in range (1,len(time_stamps_sec[idt]))]

    #EVENT (obtain start and end of the different conditions of a trial)
    for filename in glob.glob('*events*'):
        with open(os.path.join(os.getcwd(), filename)) as event_file:
            d=0
            idt=filename[3:8]
            events_start[idt]=[]
            events_end[idt]=[]
            events[idt]=[e for e in filename[9:11]]
            reader = csv.reader(event_file, delimiter=";")
            header = next(reader)
            for row in reader:
                if row[1] == "view_cleared":
                    events_start[idt]+=["".join([n for n in row[0] if n !=","])]
                if row[1] == "view_filled_RGBA(0.000, 0.000, 0.000, 1.000)":
                    events_end[idt]+=["".join([n for n in row[0] if n !=","])]
                    

### Writing the new file used for statistical analysis 
os.chdir(r"") ### WRITE THE PATH TO THE FILE WHERE YOU WANT THE NEW FILE TO BE CREATED
Data_file = "expballs_stat.csv"
Newfile = open(Data_file, "w",newline='')
writer = csv.writer(Newfile, delimiter=";")
writer.writerow(["Ant","Size1","Size2","Ball","Trial","Condition","Cookie","Order","Mean angular speed","blocked oscillation","fourrier","nb_osci","nb_block"])

starts={}
ends={}

starts={}
ends={}
j=0

for idt in list_idt:
    j+=1
    starts[idt]=[]
    ends[idt]=[]

    for e in range (len(events[idt])):
        start = round(time_stamps_raw[idt][int(events_start[idt][e])]/0.033)+20
        end = round(time_stamps_raw[idt][int(events_end[idt][e])]/0.033)
        starts[idt]+=[start]
        ends[idt]+=[end]
        
        # Time spent turning left or right
        time_left= len([a for a in Ang_speed_smooth[idt][start:end]  if a>3])*0.033 #margin to avoid freezing period
        time_right= len([a for a in Ang_speed_smooth[idt][start:end] if a<-3])*0.033 #margin to avoid freezing period
    
        #frequency of direction switches
        nb_changes=0
        for i in range(start+1,end):
           if (Ang_speed_smooth[idt][i]<=1 and Ang_speed_smooth[idt][i-1]>=1): #baseline at 1 instead of 0 to avoid counting freezing as switches
                nb_changes+=1
        frequency0= (nb_changes)/(time_left+time_right)

        #frequency of oscillations
        nb_pics=0
        a=0
        slope= Ang_speed_smooth[idt][start+1]-Ang_speed_smooth[idt][start]
        for i in range(start+1,end):
            if (Ang_speed_smooth[idt][i+1]-Ang_speed_smooth[idt][i])*slope<0:
                if np.abs(Ang_speed_smooth[idt][i]-a)>10: #if oscillation amplitude is enough
                    nb_pics+=1
                a=Ang_speed_smooth[idt][i]
            slope= Ang_speed_smooth[idt][i+1]-Ang_speed_smooth[idt][i]
        frequency2=nb_pics/(2*(time_left+time_right))
        nb_osci=round(nb_pics/2)
        
        #FOURRIER analysis to obtain frequency
        autocor = sm.tsa.acf([a for a in Ang_speed_smooth[idt][start:end] if a!=0], nlags=len(Ang_speed_smooth[idt][start:end]))
        fourierTransform = np.fft.rfft(autocor)
        frequency = np.linspace(0, 1/(2*0.033), len(np.square(abs(fourierTransform))))
        argfourrier=np.argmax(np.square(abs(fourierTransform[12:])))+12
        fourrier=frequency[argfourrier]
        
        p_blocked=(frequency2-frequency0)/frequency2

        writer.writerow([ant[idt],size1[ant[idt]],size2[ant[idt]],condition[idt],j,events[idt][e],cookie[idt],e+1,
                         np.mean([np.abs(s) for s in Ang_speed_smooth[idt][start:end] if s!=0]),
                         p_blocked,fourrier,nb_changes,nb_osci-nb_changes])
Newfile.close()

### EXPERIMENT NB BARS : SUPPLEMENTARY ########################################################################################################################################
# This part extract the raw data to obtain dictionnaries corresponding to our parameters of interest for each trial as well as the beginning and  the end of each condition inside each trial.
# idt correspond to the name, or "identity", of a singular trial and is used as the index inside all the dictionnaries to access the data of this specific trial.
# at this end of this part a csv table will containing the values used for the statictical analysis will be created.
# YOU  HAVE TO RUN THIS CODE LAST IF YOU WANT TO VIZUALIZE THE DATA CORRESPONDING TO THIS EXPERIMENT (such as mean oscillation cycles, trajectories or angular velocity signals)  

mainfolder = '' ### WRITE THE PATH TO THE FILE WITH RAW DATA CORRESPONDING TO THE VERTICAL BARS EXPERIMENT

list_idt=[]
ant={}
cookie={}
condition={}
Ang_speed_degsec={}
Ang_speed_smooth={}
Fwd_speed_cm_sec={}
X={}
Y={}
Theta_path={}
time_stamps_raw={}
time_stamps_sec={}
Ang_acceleration={}
Fwd_acceleration={}
events_start={}
events_end={}
events={}

for f in glob.glob(mainfolder+'*Ant*'):
    os.chdir(f)

    # COORDS
    for file in glob.glob(f+'/*coords*'):
        filename=file[-36:]
        idt=(filename[3:5])
        list_idt+=[idt]
        ant[idt]=filename[3:5]
        cookie[idt]=filename[13:15]
        fclean=clean(file)
        clean_file=f+"/"+fclean
        with open(clean_file,"r") as csvfile :
            reader=csv.reader(csvfile, delimiter=';')
            time_stamps=[]
            i=0
            for row in reader:
                if i==0:
                    trackball_data_raw = [row[3:7]]
                else:
                    trackball_data_raw += [[float(a) for a in row[3:7]]]
                    time_stamps +=  [float(row[7])]
                i+=1
        csvfile.close()
        os.remove(clean_file)

        #speed and position
        Ang_speed_degsec[idt], Fwd_speed_cm_sec[idt], X[idt], Y[idt], Theta_path[idt] = trackball_data_cruncher (idt,trackball_data_raw, time_stamps, ball_diameter, calib_0, calib_1)
        #time stamps
        time_stamps_raw[idt] =  [a/1000 for a in time_stamps]# in sec
        print(idt)
        # resample a data point every 33 milliseconds
        df = pd.DataFrame(np.array([Fwd_speed_cm_sec[idt], Ang_speed_degsec[idt],X[idt],Y[idt],time_stamps_raw[idt][:len(time_stamps_raw[idt])-1]]).T, columns=["Fwd_speed","Ang_speed","X","Y","Time"],index=pd.to_datetime(time_stamps_raw[idt][:len(time_stamps_raw[idt])-1], unit='s'))
        df2 = df.resample('33L').mean().reset_index().interpolate(method='ffill')
        time_stamps_sec[idt] = df2.Time.tolist()
        X[idt] = df2.X.tolist()
        Y[idt] = df2.Y.tolist()
        #smoothing whith a window of 20 frames (1 time fwd speed, 2 times angular speed)
        Fwd_speed_cm_sec[idt] = df2.Fwd_speed.rolling(20,min_periods=1).mean().tolist()
        Ang_speed_degsec[idt] = df2.Ang_speed.tolist()
        df2['Ang_speed_smooth'] = df2.Ang_speed.rolling(20,min_periods=1).mean().tolist()
        Ang_speed_smooth[idt] = df2.Ang_speed_smooth.rolling(20,min_periods=1).mean().tolist()
        #angular acceleration
        Ang_acceleration[idt]= [(Ang_speed_smooth[idt][i]-Ang_speed_smooth[idt][i-1])/(time_stamps_sec[idt][i]-time_stamps_sec[idt][i-1]) for i in range (1,len(time_stamps_sec[idt])-1)]
        #forward acceleration
        Fwd_acceleration[idt]= [(Fwd_speed_cm_sec[idt][i]-Fwd_speed_cm_sec[idt][i-1])/(time_stamps_sec[idt][i]-time_stamps_sec[idt][i-1]) for i in range (1,len(time_stamps_sec[idt]))]

    #EVENT (obtain start and end of the different conditions of a trial)
    for filename in glob.glob('*events*'):
        with open(os.path.join(os.getcwd(), filename)) as event_file:
            idt=(filename[3:5])
            events_start[idt]=[]
            events_end[idt]=[]
            events[idt]=[e for e in filename[6:12]]
            reader = csv.reader(event_file, delimiter=";")
            header = next(reader)
            for row in reader:
                if row[1] == "view_cleared":
                    events_start[idt]+=["".join([n for n in row[0] if n !=","])]
                if row[1] == "view_filled_RGBA(0.000, 0.000, 0.000, 1.000)":
                    events_end[idt]+=["".join([n for n in row[0] if n !=","])]
                    

### Writing the new file used for statistical analysis 
os.chdir(r"") ### WRITE THE PATH TO THE FILE WHERE YOU WANT THE NEW FILE TO BE CREATED
Data_file = "expbars_stat.csv"
Newfile = open(Data_file, "w",newline='')
writer = csv.writer(Newfile, delimiter=";")
writer.writerow(["Ant","Condition","Cookie","Order","Mean angular speed","blocked oscillation","fourrier","nb_osci","nb_block"])

    
starts={}
ends={}


for idt in list_idt:
    starts[idt]=[]
    ends[idt]=[]

    for e in range (len(events[idt])):
        start = round(time_stamps_raw[idt][int(events_start[idt][e])]/0.033)+20
        end = round(time_stamps_raw[idt][int(events_end[idt][e])]/0.033)
        starts[idt]+=[start]
        ends[idt]+=[end]
        
        # Time spent turning left or right
        time_left= len([a for a in Ang_speed_smooth[idt][start:end]  if a>3])*0.033 #margin to avoid freezing period
        time_right= len([a for a in Ang_speed_smooth[idt][start:end] if a<-3])*0.033 #margin to avoid freezing period
    
        #frequency of direction switches
        nb_changes=0
        for i in range(start+1,end):
           if (Ang_speed_smooth[idt][i]<=1 and Ang_speed_smooth[idt][i-1]>=1): #baseline at 1 instead of 0 to avoid counting freezing as switches
                nb_changes+=1
        frequency0= (nb_changes)/(time_left+time_right)

        #frequency of oscillations
        nb_pics=0
        a=0
        slope= Ang_speed_smooth[idt][start+1]-Ang_speed_smooth[idt][start]
        for i in range(start+1,end):
            if (Ang_speed_smooth[idt][i+1]-Ang_speed_smooth[idt][i])*slope<0:
                if np.abs(Ang_speed_smooth[idt][i]-a)>10: #if oscillation amplitude is enough
                    nb_pics+=1
                a=Ang_speed_smooth[idt][i]
            slope= Ang_speed_smooth[idt][i+1]-Ang_speed_smooth[idt][i]
        frequency2=nb_pics/(2*(time_left+time_right))
        nb_osci=round(nb_pics/2)
        
        #FOURRIER analysis to obtain frequency
        autocor = sm.tsa.acf([a for a in Ang_speed_smooth[idt][start:end] if a!=0], nlags=len(Ang_speed_smooth[idt][start:end]))
        fourierTransform = np.fft.rfft(autocor)
        frequency = np.linspace(0, 1/(2*0.033), len(np.square(abs(fourierTransform))))
        argfourrier=np.argmax(np.square(abs(fourierTransform[12:])))+12
        fourrier=frequency[argfourrier]
        
        p_blocked=(frequency2-frequency0)/frequency2

        writer.writerow([ant[idt],events[idt][e],cookie[idt],e+1,
                         np.mean([np.abs(s) for s in Ang_speed_smooth[idt][start:end] if s!=0]),
                         p_blocked,fourrier,nb_changes,nb_osci-nb_changes])
Newfile.close()


       
