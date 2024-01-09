#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 13:23:26 2024

@author: Oceane DAUZERE-PERES
"""
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

def trajectory(idt,event):
    """Display the trajectory correponding to the ant indicated by idt in the condition specified by event
    Yellow portions of the trajectory mean that the is turning left (positive angular velocity) and purple portions mean that the ant is trurning right (negative angular velocity)"""
    # finding the beginning and the end of the chosen event in the signal 
    j=0
    e=0
    for cond in events[idt]:
        if events[idt][j]==event:
            e=j
        j+=1
    start=starts[idt][e]
    end=ends[idt][e]
    
    # Plotting the associated trajectory
    X1=[x-X[idt][starts[idt][e]] for x in X[idt][start:end]]
    Y1=[y-Y[idt][starts[idt][e]] for y in Y[idt][start:end]]
    
    points = np.array([X1, Y1]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    fig, ax = plt.subplots()
    lc = LineCollection(segments, cmap='viridis')
    lc.set_array([np.sign(v) for  v in Ang_speed_smooth[idt][(start+20):(end+20)]])  #np.sign(v)
    lc.set_linewidth(1)
    line = ax.add_collection(lc)
    fig.colorbar(line,ax=ax)
    
    plt.scatter(X1,Y1,linewidth=0.6,alpha=0.7,c=[np.sign(v) for  v in Ang_speed_smooth[idt][start+20:end+20]],s=1,cmap='viridis')
    plt.show()
    


def angular_velocity_signal(idt,event):
    """Display the angular velocity signal correponding to the ant indicated by idt in the condition specified by event
    The blue curve correspond to the raw signal and the green curve to the smoothed signal"""
    # finding the beginning and the end of the chosen event in the signal 
    j=0
    e=0
    for cond in events[idt]:
        if events[idt][j]==event:
            e=j
        j+=1
    start=starts[idt][e]
    end=ends[idt][e]
    
    # Plotting the associated angular velocity signal
    plt.plot(np.linspace(0,len(time_stamps_sec[idt][start:end])*0.033,len(time_stamps_sec[idt][start:end])),(Ang_speed_degsec[idt][start:end]),color="lightblue")
    plt.plot(np.linspace(0,len(time_stamps_sec[idt][start:end-20])*0.033,len(time_stamps_sec[idt][start:end-20])), Ang_speed_smooth[idt][start+20:end],color="green")
    plt.plot(np.linspace(0,len(time_stamps_sec[idt][start:end])*0.033,1000),np.zeros(1000),color="black")
    plt.show()

##################################### RUN THIS CODE AFTER EXTRACTING THE DATA FROM THE DATASET OF THE EXPERIMENT YOU WANT TO LOOK AT #########################################

#examples of possible figures (here with the gain experiment dataset)
trajectory("19_RE","B")
angular_velocity_signal("19_RE","B")

# enter an "idt" as the first positional argument and an "event" as the second positional argument depending on the experiment you want to look at
# for experiment on the gain : gain -1 = "-1", gain 0 = "0", gain 1 = "1", gain 3 = "3", gain 5 = "5", in the dark = "B"
# for experiment on the visual structure : gain 2.5 = "2", gain 0 = "0", homogeneous black = "B", homogeneous white = "W", horizontal bars = "H", vertical bars = "V"
# for experiment on the weight of the ball : in the dark = "B", gain 0 = "0"
# for supplementary experiment on the number of bars : homogeneous black = "B", horizontal bars = "H", 1 bar verical white in the front = "W", 1 vertical bar black in the front = "1", 4 vertical bars = "4", 8 vertical bars = "8"
