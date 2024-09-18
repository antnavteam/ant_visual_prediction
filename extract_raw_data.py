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
        df2 = df.resample('33ms').mean().reset_index().ffill()
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
writer.writerow(["Ant","Eye","Trial","Gain","Order","Repetition","Mean angular speed","Mean angular velocity","blocked oscillation","fourrier","nb_osci","nb_block","Mean forward speed"])


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

        #extremums detection
        nb_pics=0
        a=0
        mins=[]
        maxs=[]
        slope= Ang_speed_smooth[idt][start+1]-Ang_speed_smooth[idt][start]
        for i in range(start+1,end):
            if (Ang_speed_smooth[idt][i+1]-Ang_speed_smooth[idt][i])*slope<0:
                if np.abs(Ang_speed_smooth[idt][i]-a)>3: #if oscillation amplitude is enough
                    if slope > 0 : #maximums
                        if len(mins)==len(maxs) :
                            maxs+=[Ang_speed_smooth[idt][i]]
                        else:
                            maxs[-1]=np.max([maxs[-1],Ang_speed_smooth[idt][i]])
                    elif slope < 0 : #minimums
                        if len(mins) < len(maxs) : 
                            mins+=[Ang_speed_smooth[idt][i]]
                        elif len(mins)!=0:
                            mins[-1]=np.min([mins[-1],Ang_speed_smooth[idt][i]])
                a=Ang_speed_smooth[idt][i]
            slope= Ang_speed_smooth[idt][i+1]-Ang_speed_smooth[idt][i]
        
        
        # Total number of oscillations (number of extremums detected divided by 2)
        nb_osci=round((len(maxs)+len(mins))/2)
        
        #Counting the numbers of blocked oscillations by counting the extremums that belongs to blocked oscillations (if a maximums is below 0 or a minimum above 0)
        nb_block=0
        for i in maxs:
            if i < 0 :
                nb_block+=1
        for i in mins:
            if i > 0:
                nb_block+=1
        #proportion of blocked oscillations
        p_blocked = nb_block/nb_osci
                
        #FOURRIER
        autocor = sm.tsa.acf([a for a in Ang_speed_smooth[idt][start:end] if a!=0], nlags=len(Ang_speed_smooth[idt][start:end]))
        fourierTransform = np.fft.rfft(autocor)
        frequency = np.linspace(0, 1/(2*0.033), len(np.square(abs(fourierTransform))))
        argfourrier=np.argmax(np.square(abs(fourierTransform[12:])))+12
        fourrier=frequency[argfourrier]

        writer.writerow([ant[idt],eye[idt],j,events[idt][e],e+1,condition[idt],
                         np.mean([np.abs(s) for s in Ang_speed_smooth[idt][start:end] if s!=0]),
                         np.mean([(s) for s in Ang_speed_smooth[idt][start:end] if s!=0]),
                         p_blocked,fourrier,nb_osci-nb_block,nb_block, 
                         np.mean([np.abs(s) for s in Fwd_speed_cm_sec[idt][start:end] if s!=0])])
Newfile.close()


### EXPERIMENT VISUAL STRUCTURE ########################################################################################################################################
# This part extract the raw data to obtain dictionnaries corresponding to our parameters of interest for each trial as well as the beginning and  the end of each condition inside each trial. 
# idt correspond to the name, or "identity", of a singular trial and is used as the index inside all the dictionnaries to access the data of this specific trial.
# at this end of this part a csv table will containing the values used for the statictical analysis will be created. 
# YOU  HAVE TO RUN THIS CODE LAST IF YOU WANT TO VIZUALIZE THE DATA CORRESPONDING TO THIS EXPERIMENT (such as mean oscillation cycles, trajectories or angular velocity signals)

mainfolder = '' ### WRITE THE PATH TO THE FILE WITH RAW DATA CORRESPONDING TO THE VISUAL STRUCTURE EXPERIMENT

list_idt=[]
ant={}
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
        df2 = df.resample('33ms').mean().reset_index().ffill()
        time_stamps_sec[idt] = df2.Time.tolist()
        X[idt] = df2.X.tolist()
        Y[idt] = df2.Y.tolist()
        #smoothing with a window of 20 frames (1 time fwd speed, 2 times angular speed)
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
writer.writerow(["Ant","Condition","Order","Mean angular speed","blocked oscillation","fourrier","nb_osci","nb_block"])


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

        #extremums detection
        nb_pics=0
        a=0
        mins=[]
        maxs=[]
        slope= Ang_speed_smooth[idt][start+1]-Ang_speed_smooth[idt][start]
        for i in range(start+1,end):
            if (Ang_speed_smooth[idt][i+1]-Ang_speed_smooth[idt][i])*slope<0:
                if np.abs(Ang_speed_smooth[idt][i]-a)>3: #if oscillation amplitude is enough
                    if slope > 0 : #maximums
                        if len(mins)==len(maxs) :
                            maxs+=[Ang_speed_smooth[idt][i]]
                        else:
                            maxs[-1]=np.max([maxs[-1],Ang_speed_smooth[idt][i]])
                    elif slope < 0 : #minimums
                        if len(mins) < len(maxs) : 
                            mins+=[Ang_speed_smooth[idt][i]]
                        elif len(mins)!=0:
                            mins[-1]=np.min([mins[-1],Ang_speed_smooth[idt][i]])
                a=Ang_speed_smooth[idt][i]
            slope= Ang_speed_smooth[idt][i+1]-Ang_speed_smooth[idt][i]
        
        
        # Total number of oscillations (number of extremums detected divided by 2)
        nb_osci=round((len(maxs)+len(mins))/2)
        
        #Counting the numbers of blocked oscillations by counting the extremums that belongs to blocked oscillations (if a maximums is below 0 or a minimum above 0)
        nb_block=0
        for i in maxs:
            if i < 0 :
                nb_block+=1
        for i in mins:
            if i > 0:
                nb_block+=1
        #proportion of blocked oscillations
        p_blocked = nb_block/nb_osci

        
        #FOURRIER analysis to obtain frequency
        autocor = sm.tsa.acf([a for a in Ang_speed_smooth[idt][start:end] if a!=0], nlags=len(Ang_speed_smooth[idt][start:end]))
        fourierTransform = np.fft.rfft(autocor)
        frequency = np.linspace(0, 1/(2*0.033), len(np.square(abs(fourierTransform))))
        argfourrier=np.argmax(np.square(abs(fourierTransform[12:])))+12
        fourrier=frequency[argfourrier]

        writer.writerow([ant[idt],events[idt][e],e+1,
                         np.mean([np.abs(s) for s in Ang_speed_smooth[idt][start:end] if s!=0]),
                         p_blocked,fourrier,nb_osci-nb_block,nb_block])
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
        df2 = df.resample('33ms').mean().reset_index().ffill()
        time_stamps_sec[idt] = df2.Time.tolist()
        X[idt] = df2.X.tolist()
        Y[idt] = df2.Y.tolist()
        #smoothing with a window of 20 frames (1 time fwd speed, 2 times angular speed)
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
writer.writerow(["Ant","Size1","Size2","Ball","Trial","Condition","Order","Mean angular speed","blocked oscillation","fourrier","nb_osci","nb_block"])

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

        #extremums detection
        nb_pics=0
        a=0
        mins=[]
        maxs=[]
        slope= Ang_speed_smooth[idt][start+1]-Ang_speed_smooth[idt][start]
        for i in range(start+1,end):
            if (Ang_speed_smooth[idt][i+1]-Ang_speed_smooth[idt][i])*slope<0:
                if np.abs(Ang_speed_smooth[idt][i]-a)>3: #if oscillation amplitude is enough
                    if slope > 0 : #maximums
                        if len(mins)==len(maxs) :
                            maxs+=[Ang_speed_smooth[idt][i]]
                        else:
                            maxs[-1]=np.max([maxs[-1],Ang_speed_smooth[idt][i]])
                    elif slope < 0 : #minimums
                        if len(mins) < len(maxs) : 
                            mins+=[Ang_speed_smooth[idt][i]]
                        elif len(mins)!=0:
                            mins[-1]=np.min([mins[-1],Ang_speed_smooth[idt][i]])
                a=Ang_speed_smooth[idt][i]
            slope= Ang_speed_smooth[idt][i+1]-Ang_speed_smooth[idt][i]
        
        
        # Total number of oscillations (number of extremums detected divided by 2)
        nb_osci=round((len(maxs)+len(mins))/2)
        
        #Counting the numbers of blocked oscillations by counting the extremums that belongs to blocked oscillations (if a maximums is below 0 or a minimum above 0)
        nb_block=0
        for i in maxs:
            if i < 0 :
                nb_block+=1
        for i in mins:
            if i > 0:
                nb_block+=1
        #proportion of blocked oscillations
        p_blocked = nb_block/nb_osci
        
        #FOURRIER analysis to obtain frequency
        autocor = sm.tsa.acf([a for a in Ang_speed_smooth[idt][start:end] if a!=0], nlags=len(Ang_speed_smooth[idt][start:end]))
        fourierTransform = np.fft.rfft(autocor)
        frequency = np.linspace(0, 1/(2*0.033), len(np.square(abs(fourierTransform))))
        argfourrier=np.argmax(np.square(abs(fourierTransform[12:])))+12
        fourrier=frequency[argfourrier]

        writer.writerow([ant[idt],size1[ant[idt]],size2[ant[idt]],condition[idt],j,events[idt][e],e+1,
                         np.mean([np.abs(s) for s in Ang_speed_smooth[idt][start:end] if s!=0]),
                         p_blocked,fourrier,nb_osci-nb_block,nb_block])
Newfile.close()

### EXPERIMENT NB BARS : SUPPLEMENTARY ########################################################################################################################################
# This part extract the raw data to obtain dictionnaries corresponding to our parameters of interest for each trial as well as the beginning and  the end of each condition inside each trial.
# idt correspond to the name, or "identity", of a singular trial and is used as the index inside all the dictionnaries to access the data of this specific trial.
# at this end of this part a csv table will containing the values used for the statictical analysis will be created.
# YOU  HAVE TO RUN THIS CODE LAST IF YOU WANT TO VIZUALIZE THE DATA CORRESPONDING TO THIS EXPERIMENT (such as mean oscillation cycles, trajectories or angular velocity signals)  
mainfolder = "" ### WRITE THE PATH TO THE FILE WITH RAW DATA CORRESPONDING TO THE VERTICAL BARS EXPERIMENT

list_idt=[]
ant={}
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
        df2 = df.resample('33ms').mean().reset_index().ffill()
        time_stamps_sec[idt] = df2.Time.tolist()
        X[idt] = df2.X.tolist()
        Y[idt] = df2.Y.tolist()
        #smoothing with a window of 20 frames (1 time fwd speed, 2 times angular speed)
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
writer.writerow(["Ant","Condition","Order","Mean angular speed","blocked oscillation","fourrier","nb_osci","nb_block"])

    
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

        #extremums detection
        nb_pics=0
        a=0
        mins=[]
        maxs=[]
        slope= Ang_speed_smooth[idt][start+1]-Ang_speed_smooth[idt][start]
        for i in range(start+1,end):
            if (Ang_speed_smooth[idt][i+1]-Ang_speed_smooth[idt][i])*slope<0:
                if np.abs(Ang_speed_smooth[idt][i]-a)>3: #if oscillation amplitude is enough
                    if slope > 0 : #maximums
                        if len(mins)==len(maxs) :
                            maxs+=[Ang_speed_smooth[idt][i]]
                        else:
                            maxs[-1]=np.max([maxs[-1],Ang_speed_smooth[idt][i]])
                    elif slope < 0 : #minimums
                        if len(mins) < len(maxs) : 
                            mins+=[Ang_speed_smooth[idt][i]]
                        elif len(mins)!=0:
                            mins[-1]=np.min([mins[-1],Ang_speed_smooth[idt][i]])
                a=Ang_speed_smooth[idt][i]
            slope= Ang_speed_smooth[idt][i+1]-Ang_speed_smooth[idt][i]
        
        
        # Total number of oscillations (number of extremums detected divided by 2)
        nb_osci=round((len(maxs)+len(mins))/2)
        
        #Counting the numbers of blocked oscillations by counting the extremums that belongs to blocked oscillations (if a maximums is below 0 or a minimum above 0)
        nb_block=0
        for i in maxs:
            if i < 0 :
                nb_block+=1
        for i in mins:
            if i > 0:
                nb_block+=1
        #proportion of blocked oscillations
        p_blocked = nb_block/nb_osci
        
        #FOURRIER analysis to obtain frequency
        autocor = sm.tsa.acf([a for a in Ang_speed_smooth[idt][start:end] if a!=0], nlags=len(Ang_speed_smooth[idt][start:end]))
        fourierTransform = np.fft.rfft(autocor)
        frequency = np.linspace(0, 1/(2*0.033), len(np.square(abs(fourierTransform))))
        argfourrier=np.argmax(np.square(abs(fourierTransform[12:])))+12
        fourrier=frequency[argfourrier]

        writer.writerow([ant[idt],events[idt][e],e+1,
                         np.mean([np.abs(s) for s in Ang_speed_smooth[idt][start:end] if s!=0]),
                         p_blocked,fourrier,nb_osci-nb_block,nb_block])
Newfile.close()

### EXPERIMENT FREE ANTS : SUPPLEMENTARY ########################################################################################################################################
# This part extract the raw data to obtain dictionnaries corresponding to our parameters of interest for each trial as well as the beginning and the end of each condition inside each trial. 
# idt correspond to the name, or "identity", of a singular trial and is used as the index inside all the dictionnaries to access the data of this specific trial.
# at this end of this part a csv table will containing the values used for the statictical analysis will be created. 
# YOU  HAVE TO RUN THIS CODE LAST IF YOU WANT TO VIZUALIZE THE DATA CORRESPONDING TO THIS EXPERIMENT (such as mean oscillation cycles, trajectories or angular velocity signals)

mainfolder = "" ### WRITE THE PATH TO THE FILE WITH RAW DATA CORRESPONDING TO THE MODULATION OF THE GAIN EXPERIMENT
list_idt=[]
ant={}
order={}
condition={}
Ang_speed_degsec={}
Ang_speed_smooth={}
Fwd_speed_cm_sec={}
angles={}
angles_cum={}
vectors={}
X={}
Y={}
X2={}
Y2={}
Theta_path={}
time_stamps_raw={}
time_stamps_sec={}
frequence={}
Ang_acceleration={}
Fwd_acceleration={}

def vector_direction(vector):
    """gives the orientation in dregree of a vector (a, b)"""
    a=vector[0]
    b=vector[1]
    if a >= 0 and b >= 0:
        if a==0:
            return(90)
        else:
            return (np.arctan(b/a)*(180/np.pi))
    if a <= 0 and b >= 0:
        if a==0:
            return(90)
        else:
            return (180 - np.abs(np.arctan(b/a)*(180/np.pi)))
    if a <= 0 and b <= 0:
        if a==0:
            return(270)
        else:
            return (180 + np.abs(np.arctan(b/a)*(180/np.pi)))
    if a >= 0 and b <= 0:
        if a==0:
            return(270)
        else:
            return (360 - np.abs(np.arctan(b/a)*(180/np.pi)))

for file in glob.glob(mainfolder+'*ant*'):

    filename=file[-24:]
    idt=filename[3:5]+'_'+filename[-5]
    list_idt+=[idt]
    ant[idt]=filename[3:5]
    order[idt]=filename[-7]
    condition[idt]=filename[-5]
    with open(file,"r") as csvfile :
        reader=csv.reader(csvfile, delimiter=';')
        time_stamps_raw[idt]=[]
        X[idt],Y[idt] = [],[]
        i=0
        for row in reader:
            if i > 0 :
                X[idt] += [float(row[0])]
                Y[idt] += [float(row[1])]
                time_stamps_raw[idt] +=  [float(row[3])]
            i+=1
    csvfile.close()
    
        
    #time stamps
    print(idt)
    # resample a data point every 100 milliseconds
    df = pd.DataFrame(np.array([X[idt][2:],Y[idt][2:],time_stamps_raw[idt][2:]]).T, columns=["X","Y","Time"],index=pd.to_datetime(time_stamps_raw[idt][2:], unit='s'))
    df2 = df.resample('100ms').mean().reset_index().ffill()
    #smoothing the trajectory (x, y) 2 times with a window of 20 frames 
    X2[idt] = df2.X.rolling(20,min_periods=1).mean().tolist()
    Y2[idt] = df2.Y.rolling(20,min_periods=1).mean().tolist()
    time_stamps_sec[idt] = df2.Time.tolist()
    
    #speed and position
    vectors[idt] = [(X2[idt][1]-X2[idt][0],Y2[idt][1]-Y2[idt][0])]
    angles[idt]=[vector_direction(vectors[idt][0])]
    angles_cum[idt]=[vector_direction(vectors[idt][0])]
    Ang_speed_degsec[idt] = []
    #for each point in the trajectory : calculating the orientation of the ant and making the difference with the previous orientation to obtain the angle of the ant's turn.
    #then computing a cumulative orientation (modulo 360 of the absolute orientation)
    for i in range (2,len(X2[idt])):
        vectors[idt] += [(X2[idt][i]-X2[idt][i-1],Y2[idt][i]-Y2[idt][i-1])]
        angles[idt]+= [vector_direction(vectors[idt][i-1])]
        angles_cum[idt] += [angles_cum[idt][i-2]+(angles[idt][i-1]-angles[idt][i-2]) if np.abs(angles[idt][i-1]-angles[idt][i-2]) < 180 else angles_cum[idt][i-2]-np.sign(angles[idt][i-1]-angles[idt][i-2])*(360-np.abs(angles[idt][i-1]-angles[idt][i-2]))]

    #smoothing the cumulative orientation with a window of 10 frames 
    df3 = pd.DataFrame(np.array([angles_cum[idt],time_stamps_sec[idt][1:]]).T, columns=["Angle","Time"],index=pd.to_datetime(time_stamps_sec[idt][1:], unit='s'))
    angles_cum[idt] = df3.Angle.rolling(10,min_periods=1).mean().tolist()
    
    # making the difference between the cumulative orientation angles at each time step in order to obtain the amplitude of the turn performed by the ant thus the corresponding angular velocity
    for i in range (1,len(angles[idt])):
        Ang_speed_degsec[idt] += [(angles_cum[idt][i]-angles_cum[idt][i-1])/0.05]
        
    df4 = pd.DataFrame(np.array([Ang_speed_degsec[idt],time_stamps_sec[idt][2:]]).T, columns=["Ang_speed","Time"],index=pd.to_datetime(time_stamps_sec[idt][2:], unit='s'))
    #smoothing the angular speed whith a window of 20 frames 2 times 
    Ang_speed_smooth[idt] = df4.Ang_speed.rolling(20,min_periods=1).mean().tolist()

### Writing the new file used for statistical analysis 
os.chdir(r"") ### WRITE THE PATH TO THE FILE WHERE YOU WANT THE NEW FILE TO BE CREATED
Data_file = "expfree_stat.csv"
Newfile = open(Data_file, "w",newline='')
writer = csv.writer(Newfile, delimiter=";")
writer.writerow(["Ant","Trial","Condition","Average angular velocity","blocked oscillation","nb_osci","nb_block"])


for idt in list_idt:
    
    #extremums detection
    maxs=[]
    mins=[]
    a=0
    slope= Ang_speed_smooth[idt][1]-Ang_speed_smooth[idt][0]
    for i in range(len(Ang_speed_smooth[idt])-1):
        if (Ang_speed_smooth[idt][i+1]-Ang_speed_smooth[idt][i])*slope<0:
            if np.abs(Ang_speed_smooth[idt][i]-a)>15: #if oscillation amplitude is enough
                if slope > 0 : #maximum
                    if len(maxs)==len(mins):
                        maxs+=[Ang_speed_smooth[idt][i]]
                    else :
                        maxs[-1]=np.max([maxs[-1],Ang_speed_smooth[idt][i]])
                if slope < 0 : #minimum
                    if len(mins)<len(maxs):
                        mins+=[Ang_speed_smooth[idt][i]]
                    elif len(mins)!=0:
                        mins[-1]=np.min([mins[-1],Ang_speed_smooth[idt][i]])     
            a=Ang_speed_smooth[idt][i]
        slope= Ang_speed_smooth[idt][i+1]-Ang_speed_smooth[idt][i]
    
    # Total number of oscillations (number of extremums detected divided by 2)
    nb_osci=round((len(maxs)+len(mins))/2)
    
    #Counting the numbers of blocked oscillations by counting the extremums that belongs to blocked oscillations (if a maximums is below 0 or a minimum above 0)
    nb_block=0
    for i in maxs:
        if i < 0 :
            nb_block+=1
    for i in mins:
        if i > 0:
            nb_block+=1
    #proportion of blocked oscillations
    p_blocked = nb_block/nb_osci
    
    writer.writerow([ant[idt],order[idt],condition[idt],
                         np.mean([np.abs(s) for s in Ang_speed_smooth[idt] if s!=0]),
                         p_blocked,nb_osci-nb_block,nb_block])
Newfile.close()
