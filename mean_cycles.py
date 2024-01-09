#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 13:23:26 2024

@author: Oceane DAUZERE-PERES
"""

mport pandas as pd
import numpy as np
import matplotlib.pyplot as plt

##############

def mean_individual_cycle(idt,event,list_cycle,list_conditions1,list_conditions2,pt_synchro,oscillations):
    """
    Compute the mean oscillation cycle for a given ant idt in the condition specified by event
    list_conditions1 and list_conditions2 correspond to the conditions we want mean oscillation cycle for at the population level. We are adding to this list the individual mean oscillation cycle computed.
    list_cycle is a list of the conditions of mean oscillations cycle that we are computing at the population level, it is only updated if the individual mean oscillation cycle computed correspond to new conditions
    pt_synchro correponds to a dictionnary of indexes correponding to the middle of the mean oscillation cyles for each cycle type in list_cycle. We are adding to the dictionnary the index of the middle of individual mean oscillation cycle computed.
    oscillations is a dictionnary containing lists corresponding to mean oscillation cycles for each cycle type in list_cycle. We are adding to the dictionnary the individual mean oscillation cycle computed.
    """"
    print(idt)
    e=event
    Ang_speed=Ang_speed_smooth[idt][starts[idt][e]:ends[idt][e]-20] 

    osci_middle=[]
    maxs=[]
    mins=[]
    
    slope = Ang_speed[1]-Ang_speed[0]
    a=Ang_speed[0]
#    for i in range (1,len(Ang_speed)) :
#        if (Ang_speed[i]-Ang_speed[i-1])*slope<0 :
#            if np.abs(Ang_speed[i]-a)>20: #enough amplitude to be considered an oscillation
#                if slope >0 and len(mins)==len(maxs):
#                    maxs+=[i-1]
#                    a=Ang_speed[i]
#                if slope <0  and len(mins)<len(maxs):
#                    mins+=[i-1]
#                    a=Ang_speed[i]
#                if len(maxs)==len(mins) and len(maxs)!=0: # new oscillation when a new local minimum and local maximum has been detected
#                    middle=np.mean(Ang_speed[maxs[-1]:mins[-1]]) 
#                    for j in range (maxs[-1],mins[-1]):
#                        if Ang_speed[j-1]>=middle and Ang_speed[j]<=middle:
#                            osci_middle+=[j]
#        slope=Ang_speed[i]-Ang_speed[i-1]
#
    
    ## Detection of the oscillations in the signal
    slope = Ang_speed[1]-Ang_speed[0]
    a=Ang_speed[0]
    baseline=np.mean(Ang_speed)
    osci_up=[]
    osci_down=[]
    for i in range (1,len(Ang_speed)) :
        if Ang_speed[i]>baseline and Ang_speed[i-1]<=baseline:
            osci_up+=[[i]]
            if len(osci_down)!=0:
                osci_down[len(osci_down)-1]+=[i-1]
        if Ang_speed[i]<baseline and Ang_speed[i-1]>=baseline:
            osci_down+=[[i]]
            osci_middle+=[i]
            if len(osci_up)!=0:
                osci_up[len(osci_up)-1]+=[i-1]

   # Averaging the different oscillations detected
    l=1
    center_index=0
    mean_osci= [np.mean([Ang_speed[osci_middle[i]] for i in range(len(osci_middle)-2)])]
    
    while l<80:
        mean_osci+=[np.mean([Ang_speed[osci_middle[i]+l] for i in range(len(osci_middle)-2)])] 
        mean_osci=[np.mean([Ang_speed[osci_middle[i]-l] for i in range(len(osci_middle)-2)])] + mean_osci 
        center_index+=1
        l+=1
        
    #updating the global variables at the population level with this new mean oscillation cycle
    if "1E" in list_conditions1 and (idt[3:]=="RE" or idt[3:]=="LE") :
        cycle=str(events[idt][e])+"_1E" 
    else :
        cycle=str(events[idt][e])+idt[-3:]
    if cycle not in list_cycle:
        list_cycle+=[cycle]
        pt_synchro[cycle]=[center_index]
        oscillations[cycle]=[mean_osci]
    else:
        pt_synchro[cycle]+=[center_index]
        oscillations[cycle]+=[mean_osci]
    
    #plotting the individual mean oscillation cycle in transparency
    line=(":" if "LE" in cycle or "LB" in cycle or "1E" in cycle else ("--" if "RE" in cycle else "-"))
    j= list_conditions2.index((cycle[0] if cycle[0]!="-" else "-1"))
    plt.plot(np.linspace(0,len(mean_osci)*0.033,len(mean_osci))-(center_index)*0.033,(mean_osci),color=colors[j],linestyle=line,alpha=0.1)
    plt.show()
    
    return (list_cycle,pt_synchro,oscillations)



def mean_population_cycle(list_idt,list_cycle,list_conditions1,list_conditions2):
    """
    Compute the mean oscillation cycles at the population level (using mean oscillation cycles at the individual level) corresponding to the combination of conditions in list_conditions1 and list_conditions2.
    list_cycle is a list of the conditions of mean oscillations cycle that we are computing at the population level (corresponding to combinationsof conditions in list_conditions1 and list_conditions2).
    """"
    list_cycle=[]
    pt_synchro={}
    oscillations={}
    #Compute the individual mean oscillation cycles corresponding to all the chosen combination of conditions
    for idt in list_idt:
        if idt[3:] in list_conditions1 or ("1E" in list_conditions1 and (idt[3:]=="RE" or idt[3:]=="LE")):
            for e in range (len(events[idt])):
                if events[idt][e] in list_conditions2:
                    list_cycle,pt_synchro,oscillations=mean_individual_cycle(idt,e,list_cycle,list_conditions1,list_conditions2,pt_synchro,oscillations)
                    
    #Averaging the individual mean oscillation cycles to obtain the mean oscillation cycle at the population level for each oscillation type in list_cycle (corresponding to all the chosen combination of conditions) 
    for cycle in list_cycle:
        l=1
        center_index=0
        osci_middle=pt_synchro[cycle]
        osci=oscillations[cycle]
        mean_osci= [np.mean([osci[j][osci_middle[j]] for j in range(len(osci))])]
        while l<80:
            mean_osci+=[np.mean([osci[j][osci_middle[j]+l] for j in range(len(osci))])]
            mean_osci= [np.mean([osci[j][osci_middle[j]-l] for j in range(len(osci))])] + mean_osci
            center_index+=1
            l+=1
       
        #plotting the mean oscillation cycle for each combination of conditions
        line=(":" if "LE" in cycle or "LB" in cycle or "1E" in cycle else ("--" if "RE" in cycle else "-"))
        print(cycle,list_cycle)
        j= list_conditions2.index((cycle[0] if cycle[0]!="-" else "-1"))
        
        plt.plot(np.linspace(0,len(mean_osci)*0.033,len(mean_osci))-(center_index)*0.033,(mean_osci),color=colors[j],linestyle=line,label=cycle)
        plt.legend(loc="upper left")
        plt.plot(np.linspace(-2.5,2.5,60),np.zeros(60),color="black")
        plt.show()

##################################### RUN THIS CODE AFTER EXTRACTING THE DATA FROM THE DATASET OF THE EXPERIMENT YOU WANT TO LOOK AT #########################################

colors=["purple","red","orange","green","b","black"]

## Figure corresponding to the ones in the article (RUN WITH THE DATASET OF THE GAIN EXPERIMENT)
mean_population_cycle(list_idt,list_cycle,["2E"],["-1","0","1","3","5","B"])
mean_population_cycle(list_idt,list_cycle,["2E","1E"],["B"])
mean_population_cycle(list_idt,list_cycle,["2E","1E"],["5"])
mean_population_cycle(list_idt,list_cycle,["2E","1E"],["0"])

# For list_conditions1 (3rth positional argument) write the list of condition between trials as specified in the name of the file :
# for experiment on the gain : two eyes="2E", right eye covered ="RE", left eye covered="LE"  /!\ if you want to combine RE and LE write "1E" (= 1 eye covered)
# for experiment on the weight of the ball : Heavy ball = "HB", Light ball ="LB"
# for other experiment write list_conditions1=[""]

# For list_conditions2 (4rth positional argument) Write the list of condition inside a trial as specified in the name of the file :
# for experiment on the gain : gain -1 = "-1", gain 0 = "0", gain 1 = "1", gain 3 = "3", gain 5 = "5", in the dark = "B"
# for experiment on the visual structure : gain 2.5 = "2", gain 0 = "0", homogeneous black = "B", homogeneous white = "W", horizontal bars = "H", vertical bars = "V"
# for experiment on the weight of the ball : in the dark = "B", gain 0 = "0"
# for supplementary experiment on the number of bars : homogeneous black = "B", horizontal bars = "H", 1 bar verical white in the front = "W", 1 vertical bar black in the front = "1", 4 vertical bars = "4", 8 vertical bars = "8"



