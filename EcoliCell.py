#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 31 17:02:36 2017

@author: michael
"""


import numpy as np



#cell_volumes_t1=np.ones(cell_num)


class EcoliCell:
    
    def __init__(self,start_volume,sigma,g_rate):
        self.volume=abs(np.random.normal(start_volume,sigma))
        #self.volume=abs(np.random.lognormal(start_volume,sigma))
        #self.volume=start_volume
        self.newborn_volume=self.volume
        self.g_rate=g_rate
        self.age=0;
        
        return;

    def initialize_age_variation(self,age_min,age_max):
        self.age=np.random.uniform(age_min,age_max)
        #self.age=np.random.randint(age_min,age_max)
        return;

    #shows volume after delta_t assuming steady growth rate
    def add_volume(self,delta_t):
        self.volume=self.volume+(self.volume*((np.exp((self.g_rate)*(delta_t/60)))-1))
        self.age+=1
        return;
        
        
    def initiate_cell_cycle(self,cycle_offset):
        self.cycle_timer=cycle_offset


    
    def divide(self):#,sizer_threshold,timer_threshold):#,division_error):
        #self.division_error=division_error
        self.volume_predivision=self.volume
        self.volume=0.5*self.volume
        self.newborn_volume=self.volume
        #self.volume=np.random.normal(0.5,self.division_error)*self.volume
        self.age=0
        self.cycle_timer=0
        return;
        
class EcoliCellBCD:
    
    def __init__(self,start_volume,start_volume_sigma,g_rate,num_oric):
        self.volume=abs(np.random.normal(start_volume,start_volume_sigma))
        #self.volume=abs(np.random.lognormal(start_volume,sigma))
        #self.volume=start_volume
        self.newborn_volume=self.volume
        self.g_rate=g_rate
        self.age=0;
        self.num_oric=num_oric
        self.num_terc=1
        self.init_timer=0
        return;

    def initialize_cycle_variation(self,age_min,age_max):
        self.mother_cycle_timer=np.floor(np.random.uniform(age_min,age_max))
        #self.age=np.random.randint(age_min,age_max)
        return;

    #shows volume after delta_t assuming steady growth rate
    def add_volume(self,delta_t):
        self.volume=self.volume+(self.volume*((np.exp((self.g_rate)*(delta_t/60)))-1))
        self.age+=1
        return;
        
    def add_length(self,delta_t):
        self.length=self.volume+(self.length*((np.exp((self.g_rate)*(delta_t/60)))-1))
        self.volume=(2*np.pi*(self.length-self.width))+(4*np.pi*((self.width/2)**2))
        self.age+=1
        return;
        
        
    def initiate_mother_C_phase(self):
        self.mother_cycle_timer=0
        self.num_oric=2  ##number of replication origins
        self.num_terc=1  ##number of replication termini
            
    def initiate_daughter_C_phase(self):
        self.daughter_cycle_timer=0
        self.num_oric=4  ##number of replication origins
        self.num_terc=1  ##number of replication termini
    
    
    def initiate_granddaughter_C_phase(self):
        self.granddaughter_cycle_timer=0
        self.num_oric=8  ##number of replication origins
        self.num_terc=1
    
    def initiate_D_phase(self):##number of replication origins
        self.num_terc=2
    
    
    
    def divide(self):#,sizer_threshold,timer_threshold):#,division_error):
        #self.division_error=division_error
        self.volume_predivision=self.volume
        self.volume=0.5*self.volume
        self.newborn_volume=self.volume
        #self.volume=np.random.normal(0.5,self.division_error)*self.volume
        self.age=0
        self.cycle_timer=0
        return;              
        
class EcoliCellBCD_lw:
    
    def __init__(self,start_length,start_width,g_rate,num_oric):
        self.length=start_length
        self.width=start_width
        self.volume=(((self.width/2)**2)*np.pi*(self.length-self.width))+((4/3)*np.pi*((self.width/2)**3))
        self.surface=(2*np.pi*(self.length-self.width))+(4*np.pi*((self.width/2)**2))
        #self.volume=abs(np.random.lognormal(start_volume,sigma))
        #self.volume=start_volume
        self.newborn_volume=self.volume
        self.newborn_length=self.length
        self.newborn_width=self.width
        self.g_rate=g_rate
        self.age=0;
        self.num_oric=num_oric
        self.num_terc=1
        self.init_timer=0
        return;

    def initialize_cycle_variation(self,age_min,age_max):
        self.mother_cycle_timer=np.random.uniform(age_min,age_max)
        #self.age=np.random.randint(age_min,age_max)
        return;

    #shows volume after delta_t assuming steady growth rate
    def add_volume(self,delta_t):
        self.volume=self.volume+(self.volume*((np.exp((self.g_rate)*(delta_t/60)))-1))
        self.length=self.width+((self.volume-((4/3)*np.pi*((self.width/2)**3)))/(((self.width/2)**2)*np.pi))
        self.age+=1
        return;
        
    def add_length(self,delta_t):
        self.length=self.length+(self.length*((np.exp((self.g_rate)*(delta_t/60)))-1))
        self.width=np.random.normal(self.newborn_width,0.05*self.newborn_width)
        self.volume=(((self.width/2)**2)*np.pi*(self.length-self.width))+((4/3)*np.pi*((self.width/2)**3))
        self.age+=1
        return;
        
        
    def initiate_mother_C_phase(self):
        self.mother_cycle_timer=0
        self.num_oric=2  ##number of replication origins
        self.num_terc=1  ##number of replication termini
            
    def initiate_daughter_C_phase(self):
        self.daughter_cycle_timer=0
        self.num_oric=4  ##number of replication origins
        self.num_terc=1  ##number of replication termini
    
    def initiate_granddaughter_C_phase(self):
        self.granddaughter_cycle_timer=0
        self.num_oric=8  ##number of replication origins
        self.num_terc=1
    
    def initiate_D_phase(self):##number of replication origins
        self.num_terc=2
    
    
    
    def divide(self):#,sizer_threshold,timer_threshold):#,division_error):
        #self.division_error=division_error
        self.volume_predivision=self.volume
        self.volume=0.5*self.volume
        self.length=0.5*self.length
        self.newborn_volume=self.volume
        self.newborn_length=self.length
        self.newborn_width=self.width
        #self.volume=np.random.normal(0.5,self.division_error)*self.volume
        self.age=0
        self.cycle_timer=0
        return;       
    
     
        

   
        
    
        
    

        
        
        
        
                
     



