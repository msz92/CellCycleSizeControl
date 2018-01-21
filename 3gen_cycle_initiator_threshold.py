#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 15:59:41 2018

@author: michael
"""

from EcoliCell import EcoliCellBCD_lw
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import time
import sys



start = time.time()


#########################################################PARAMETERS ############################################################ 


exp_lengths=np.loadtxt('exp_sorb_lengths')


exp_time=500
exp_timestep=1
cell_num=100#len(exp_lengths)
init_threshold_mean=1
init_threshold_cv=0.2

mean_start_length=2.5
cv_start_length=0.2
mean_start_width=0.5
cv_start_width=0.05

doubling_time=45

start_growth_rate=np.log(2)/(doubling_time/60)#np.log(2);abs(np.random.normal(0.59,0.01))
g_rate_var=0.105*start_growth_rate

c_phase_length=38
c_phase_cv=0.118
d_phase_length=37
d_phase_cv=0.073



############################################INITIALIZE STARTING CULTURE ################################################# 

###instantiate cell objects
cell_population = [EcoliCellBCD_lw(np.random.normal(mean_start_length,cv_start_length),np.random.normal(mean_start_width,cv_start_width),np.random.normal(start_growth_rate,g_rate_var),1) for count in range(cell_num)]
#cell_population = [EcoliCellBCD(birth_mean_volume,birth_volume_sigma,np.random.normal(start_growth_rate,start_growth_rate*0.15),1) for count in range(cell_num)]
#cell_population = [EcoliCellBCD_lw(exp_lengths[i],abs(np.random.normal(0.5,0.01)),np.random.normal(start_growth_rate,g_rate_var),1) for i in range(len(exp_lengths))]

#cell_population=np.ndarray.tolist(np.random.choice(cell_population,5000))


#
for i in range(len(cell_population)):
    cell_population[i].initialize_cycle_variation(0,20)



start_g_rates=[i.g_rate for i in cell_population]
start_volumes=[i.volume for i in cell_population]
start_lengths=[i.length for i in cell_population]
start_widths=[i.width for i in cell_population]


cell_number,total_mass,mean_volumes,std_volumes,generation_times,division_volume=[],[],[],[],[],[]
newborn_volume=[]
elongation=[]
growth_rates=[]
volumecell1=[]
lengthcell1=[]
widthcell1=[]
###output cell volumes for each cell

division_size=[]    ## size at division
division_g_rate=[] ##g_rate at division
division_counter=0


##################################### SIMULATE CELL GROWTH AND DIVISION  ##########################################

if (exp_time/exp_timestep).is_integer():
    time_iterator=int(exp_time/exp_timestep)
else:
    sys.exit('Error: time series is not int')

#iterate over 



for v in tqdm(range(time_iterator)):  
    #iterate over single cells present at t
    for i in range(len(cell_population)):
        #sizer_threshold=2.2#np.random.normal(1.6,0) ##division criterion if sizer model is chosen
        #timer_threshold=np.random.normal(5,5); ##division criterion if timer model is chosen
        #increment_threshold=abs(np.random.normal(2.2,0.2595))
        
        ##add volume by exponential growth in time-increment exp_timestep        
        #cell_age=cell_population[i].age
        #volume_increment=cell_population[i].volume-cell_population[i].newborn_volume
        #if (cell_population[i].volume)>sizer_threshold and hasattr(cell_population[i],'cycle_timer')==False:
        volume_per_ori=cell_population[i].volume/cell_population[i].num_oric
        time_c_phase=np.random.normal(c_phase_length,c_phase_cv*c_phase_length)
        time_d_phase=time_c_phase+np.random.normal(d_phase_length,d_phase_cv*d_phase_length)
        
        #cell_population[i].add_volume(exp_timestep)
        
        if hasattr(cell_population[i],'mother_cycle_timer') and cell_population[i].mother_cycle_timer>time_d_phase:
        ##if D-phase is finished
                division_size.append(cell_population[i].volume)
                division_g_rate.append(cell_population[i].g_rate)
                gen_time=cell_population[i].age
                generation_times.append(gen_time)
                newborn_volume.append(cell_population[i].newborn_volume) #newborn volume of cell to divide
                elongation.append(cell_population[i].volume/cell_population[i].newborn_volume)
                
                cell_population[i].divide() ### DIVIDE
                
                cell_population[i].g_rate=np.random.normal(start_growth_rate,g_rate_var)
            
                #cell is not in cell cycle anymore
                
                division_counter+=1
                
                if hasattr(cell_population[i],'daughter_cycle_timer') and hasattr(cell_population[i],'granddaughter_cycle_timer'):
                    num_oric_newborn=4
                    daughter_cycle_timer=cell_population[i].daughter_cycle_timer
                    granddaughter_cycle_timer=cell_population[i].granddaughter_cycle_timer
                    del cell_population[i].daughter_cycle_timer
                    del cell_population[i].granddaughter_cycle_timer
                    cell_population[i].initiate_mother_C_phase()
                    cell_population[i].mother_cycle_timer=daughter_cycle_timer+1
                    cell_population[i].initiate_daughter_C_phase()
                    cell_population[i].daughter_cycle_timer=granddaughter_cycle_timer+1
                    
                    newCell=EcoliCellBCD_lw(cell_population[i].newborn_length,cell_population[i].newborn_width,np.random.normal(start_growth_rate,g_rate_var),num_oric_newborn) ### a baby is born
                    newCell.initiate_mother_C_phase()
                    newCell.initiate_daughter_C_phase()
                    newCell.init_timer=cell_population[i].init_timer
                    newCell.mother_cycle_timer=daughter_cycle_timer+1
                    newCell.daughter_cycle_timer=granddaughter_cycle_timer+1
                    cell_population.append(newCell)
                
                elif hasattr(cell_population[i],'daughter_cycle_timer'):
                    num_oric_newborn=2
                    daughter_cycle_timer=cell_population[i].daughter_cycle_timer
                    del cell_population[i].daughter_cycle_timer
                    cell_population[i].initiate_mother_C_phase()
                    cell_population[i].mother_cycle_timer=daughter_cycle_timer+1
                    
                    newCell=EcoliCellBCD_lw(cell_population[i].newborn_length,cell_population[i].newborn_width,np.random.normal(start_growth_rate,g_rate_var),num_oric_newborn) ### a baby is born
                    newCell.initiate_mother_C_phase()
                    newCell.init_timer=cell_population[i].init_timer
                    newCell.mother_cycle_timer=daughter_cycle_timer+1
                    cell_population.append(newCell)
                   
                
                else:
                    num_oric_newborn=1
                    del cell_population[i].mother_cycle_timer
                    cell_population[i].num_oric=1
                    cell_population[i].num_terc=1
                    newCell=EcoliCellBCD_lw(cell_population[i].newborn_length,cell_population[i].newborn_width,np.random.normal(start_growth_rate,g_rate_var),num_oric_newborn) ### a baby is born
                    cell_population.append(newCell)
                
                
                #create and add daughter cell
                
            
                
                
                
                
        else: #if D-phase is not finished
            initiator_threshold=np.random.normal(init_threshold_mean,init_threshold_cv)#,0.095*1)
            #volume_increment_ori=(cell_population[i].volume/cell_population[i].num_oric)-(cell_population[i].newborn_volume/cell_population[i].num_oric)
            #if volume_increment_ori>initiator_threshold:
            if volume_per_ori>initiator_threshold: ##if volume has reached threshold
                
                
                
                if hasattr(cell_population[i],'mother_cycle_timer') and hasattr(cell_population[i],'daughter_cycle_timer'):
                    
                    cell_population[i].initiate_granddaughter_C_phase()
                    cell_population[i].init_timer=0
                
                
                    if cell_population[i].mother_cycle_timer>=time_c_phase: ##check if mother C-phase has ended
                        cell_population[i].initiate_D_phase()
                        cell_population[i].add_length(exp_timestep)
                        cell_population[i].init_timer+=1
                        cell_population[i].mother_cycle_timer+=exp_timestep
                        if hasattr(cell_population[i],'daughter_cycle_timer'):
                            cell_population[i].daughter_cycle_timer+=exp_timestep
                
                
                
                elif hasattr(cell_population[i],'mother_cycle_timer'): ##if cell is already in cycle
                    
                    cell_population[i].initiate_daughter_C_phase()
                    cell_population[i].init_timer=0
                   
                    
                    if cell_population[i].mother_cycle_timer>=time_c_phase: ##check if mother C-phase has ended
                        cell_population[i].initiate_D_phase()
                        cell_population[i].add_length(exp_timestep)
                        cell_population[i].init_timer+=1
                        
                        
                        cell_population[i].mother_cycle_timer+=exp_timestep
                        if hasattr(cell_population[i],'daughter_cycle_timer'):
                            cell_population[i].daughter_cycle_timer+=exp_timestep
                            
                        
                        
                        
                    else: ## if mother C-phase is ongoing
                        cell_population[i].add_length(exp_timestep)
                        cell_population[i].init_timer+=1
                        
                        if hasattr(cell_population[i],'daughter_cycle_timer'):
                            cell_population[i].daughter_cycle_timer+=exp_timestep
                            
                            cell_population[i].mother_cycle_timer+=exp_timestep
                        
                        else:
                            cell_population[i].mother_cycle_timer+=exp_timestep
                        
                        
                
                
                
                    
                
                
                
                
                else: ##if cell is not yet in cycle
                    cell_population[i].initiate_mother_C_phase()
                    cell_population[i].init_timer=0
                    cell_population[i].add_length(exp_timestep)
                    
                    cell_population[i].mother_cycle_timer+=exp_timestep
                    
            else:   ###if volume has not reached threshold
                
                
                if hasattr(cell_population[i],'mother_cycle_timer') and cell_population[i].mother_cycle_timer>=time_c_phase:
                    cell_population[i].initiate_D_phase()
                    cell_population[i].add_length(exp_timestep)
                    cell_population[i].init_timer+=1
                    
                    
                    if hasattr(cell_population[i],'daughter_cycle_timer') and hasattr(cell_population[i],'granddaughter_cycle_timer'):
                        cell_population[i].granddaughter_cycle_timer+=exp_timestep
                        cell_population[i].daughter_cycle_timer+=exp_timestep
                        cell_population[i].mother_cycle_timer+=exp_timestep
                    
                    
                    if hasattr(cell_population[i],'daughter_cycle_timer'):
                        cell_population[i].daughter_cycle_timer+=exp_timestep
                        cell_population[i].mother_cycle_timer+=exp_timestep
                    else:
                        cell_population[i].mother_cycle_timer+=exp_timestep
                
                
                else:
                    cell_population[i].add_length(exp_timestep)
                    cell_population[i].init_timer+=1
                    if hasattr(cell_population[i],'mother_cycle_timer'):
                        if hasattr(cell_population[i],'daughter_cycle_timer'):
                            cell_population[i].mother_cycle_timer+=exp_timestep
                            cell_population[i].daughter_cycle_timer+=exp_timestep
                        else:
                            cell_population[i].mother_cycle_timer+=exp_timestep
                    
            
            
            
            
            
                            
                    
                    
                    
                    
                    
                        
                    
#if volume_increment>increment_threshold and hasattr(cell_population[i],'cycle_timer')==False:
        #if cell_age>timer_threshold:
#        if volume_per_ori>initiator_threshold:
#        
#        
#            cell_population[i].initiate_mother_C_phase() #start mother's division cycle
#            
#            if hasattr(cell_population[i],'mother_cycle_timer'):   # check if cell is already in her cycle
#                cell_population[i].initiate_daughter_C_phase()       # initiate daughter cycle if true
#                
#        
#        
#        
#        elif cell_population[i].mother_cycle_timer>time_c_phase:
#            cell_population[i].initiate_D_phase()
        
    volumes=[i.volume for i in cell_population]
    g_rates=[i.g_rate for i in cell_population]
    cell_number.append(len(cell_population))
    total_mass.append(sum(volumes))
    mean_volumes.append(np.mean(volumes))
    std_volumes.append(np.std(volumes))
    growth_rates.append(np.mean(g_rates))
    volumecell1.append(cell_population[0].volume)
    lengthcell1.append(cell_population[0].length)
    widthcell1.append(cell_population[0].width)
        

num_oris=[i.num_oric for i in cell_population]
 
    
############################################### PLOT DATA #######################################################    

t=np.linspace(0,exp_time,time_iterator)        

fig, axs = plt.subplots(3, 2)
axs[0,0].plot(t,total_mass)
axs[0,0].set_ylabel(r'total volume [$\mu m^{3}$]')
axs[0,0].set_xlabel('time [min]')
axs[0,1].plot(t,cell_number)
axs[0,1].set_ylabel('cell number')
axs[0,1].set_xlabel('time [min]')
axs[1,0].plot(t,mean_volumes)
axs[1,0].set_ylim([0 ,2])
axs[1,0].set_ylabel(r'cell volume [$\mu m^{3}$]')
axs[1,0].set_xlabel('time [min]')
axs[1,1].plot(t,std_volumes)
axs[1,1].set_ylabel(r'cell volume $\sigma$ [$\mu m$]')
axs[1,1].set_ylim([0 ,0.8])
axs[1,1].set_xlabel('time [min]')
axs[2,0].hist(newborn_volume,100,density=0)
axs[2,0].set_title('newborn volumes')
axs[2,0].set_xlabel(r'cell volume [$\mu m^{3}$]')
axs[2,0].set_xlim([0 ,4])
axs[2,1].hist(generation_times,100,density=0)
axs[2,1].set_title('generation times')
axs[2,1].set_xlabel('time [min]')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig('ss_adaptation_45min.pdf',format='pdf')
plt.show()
#

sample_end=np.random.choice(volumes,cell_num)
plt.figure(2)
plt.hist(start_volumes,50,alpha=0.6,label='start')
plt.hist(sample_end,50,alpha=0.6,label='end')
plt.legend()
plt.show()



fig, ax = plt.subplots()
fit = np.polyfit(np.log(newborn_volume), np.log(elongation), deg=1)
ax.plot(np.log(newborn_volume), fit[0] * np.log(newborn_volume) + fit[1], color='red')
ax.scatter(np.log(newborn_volume), np.log(elongation))
fig.show()


plt.figure(5)
plt.plot(t,volumecell1)
#plt.hist(volumes,1000)
plt.ylabel(r'cell volume [$\mu m^{3}$]',fontsize=18)
plt.xlabel('time [min]',fontsize=18)
plt.savefig('single_cell_volume_ss.pdf',format='pdf')
plt.show()

    
#    fig, ax = plt.subplots()
#    ax.scatter(newborn_volume,division_size)
#    fit = np.polyfit(newborn_volume, division_size, deg=1)
#    ax.plot(newborn_volume, fit[0] * division_size + fit[1], color='red')
#    fig.show()
    
end = time.time()
    
total_time=end-start

print('doubling time: {}'.format(doubling_time))
print('Execution time: {} minutes and {} seconds'.format(int(np.floor(total_time/60)),int(np.round(total_time % 60,1))))
print('{} cell divisions took place'.format(division_counter))
print('mean start: {}, std start: {}'.format(np.around(np.mean(start_volumes),3),np.around(np.std(start_volumes),3)))
print('mean now: {}, std now: {}'.format(np.around(np.mean(volumes),3),np.around(np.std(volumes),3)))
print('fit slope: {}'.format(np.around(fit[0],3)))