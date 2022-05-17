#!/usr/bin/env python
# coding: utf-8

# In[12]:


#####only know one order information, task location cen be repeat ##########################
from __future__ import division
import math
import random
import xlwt
import xlrd # import this moduel so that i can read data from a file.
import numpy
from datetime import datetime # To calculate solution time
from itertools import combinations
Vs=1.5  #maximum speed of shuttle
#As=1    #accelerarion and decelerarin of shuttle
VHc=2.5 #maximum speed of crane in horizontal direction
VVc=0.5 #maximum speed of crane in vertival direction
#AHc=0.5 #accelerarion and decelerarin of crane in horizontal direction
#AVc=0.5 #accelerarion and decelerarin of crane in verical direction
L=1.4 #length of a storage location (shuttle)
w=1.4 #width of a storage location
H=2 #height of a storage location, in z direction 
Ncx=16 #number of storage cells in x direction (shuttle)
Ncy= 80 #number of storage cells in y direction 
Ncz=3 #number of storage cells in z direction
V=3600 #storage capacity  
##############claculating time for shuttle
Piq=[] #time for shuttle to finish a task oin a lane
for i in range(Ncx):
    Piq.append(2*i*L/Vs)
####time for lift move from I/O to lane i
Pi=[]
for i in range(1,Ncy*Ncz+1):
    Pi.append(max((i//Ncy+1)*H/VVc,i%Ncy*w/VHc))    
#########time for life move from lane i to lane j
Pij=[[]for i in range(Ncy*Ncz)] 
for i in range(1,Ncy*Ncz+1):
    for j in range(1,Ncy*Ncz+1):
        Pij[i-1].append(max(abs(i%Ncy-j%Ncy)*w/VHc,abs(i//Ncy-j//Ncy)*H/VVc))
###generate tasks 
def simulation_randomshuttle(i):
    random.seed(i)
    Lanes=[i for i in range(1,Ncy*Ncz+1)]
    N0_1=random.sample(Lanes,NoS)
    N0_1=sorted(N0_1)
    N1_1=[item for item in Lanes if item not in N0_1]
    N1_1=sorted(N1_1) 
    ######generate tasks
    tasklist=[]
    while len(tasklist)<N_tasks:
        task=(random.choice(Lanes),random.randint(1,Ncx))
        #while task in tasklist: 
            #task=(random.choice(Lanes),random.randint(1,Ncx))
        tasklist.append(task)         
    movement_time=0
    TSC=0 #total single cycle time
    roundn=0
    movement_distance=0
    for task in tasklist:
        if roundn<warmup:
            if task[0] not in N0_1:
                shuttlelane=random.choice(N0_1)
                N0_1.remove(shuttlelane)
                N0_1.append(task[0])
                N1_1.append(shuttlelane)
                N1_1.remove(task[0])
        else:
            if task[0] in N0_1:
                TSC+=max(Pi[task[0]-1],Piq[task[1]-1])+Pi[task[0]-1]
            else:
                shuttlelane=random.choice(N0_1)
                N0_1.remove(shuttlelane)
                N0_1.append(task[0])
                N1_1.append(shuttlelane)
                N1_1.remove(task[0])
                TSC+=Pi[shuttlelane-1]+Pij[shuttlelane-1][task[0]-1]+Piq[task[1]-1]+Pi[task[0]-1]
                movement_time+=1
                #movement_distance+=abs(shuttlelane-task[0])
                movement_distance+=Pij[shuttlelane-1][task[0]-1]
        roundn+=1
    #print (TSC) 
    return TSC/(N_tasks-warmup)

import numpy as np
lis=[]
warmup=500
N_tasks=2000
Rs_l=0.0

NoL=Ncz*Ncy #number of lanes
for i in range(10):
    lis=[]
    Rs_l+=0.1
    
    NoS=int(math.ceil(Ncz*Ncy*Rs_l)) #number of shuttles 
    print (NoS)
    for i in range(50):  
        lis.append(simulation_randomshuttle(i))  

    T_horizontal=Ncy*w/VHc
    T_vertical=Ncz*H/VVc
    T_shuttle=2*Ncx*L/Vs
    T1=max(T_horizontal,T_vertical)
    T2=max(T1,T_shuttle)
    beta=min(T_horizontal/T1,T_vertical/T1)
    b=min(T_horizontal/T2,T_vertical/T2,T_shuttle/T2)
    for i in [T_horizontal/T2,T_vertical/T2,T_shuttle/T2]:
        if i!=1 and i!=b:
            a=i
    EC=NoS/(Ncy*Ncz)*(T1*(beta**2/6+1/2)+T2*(b**3/(12*a)+a**2/6+1/2))+(1-NoS
        /(Ncy*Ncz))*(T1*(beta**2/3+1)+T_shuttle/2+T1*(1/3+beta**2/6-beta**3/30))
    #print abs(np.mean(lis_1)-ES_R1)
    #print np.std(lis_1)
    
    print (abs(np.mean(lis)-EC),np.std(lis))


# In[ ]:




