# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 17:21:28 2021

@author: wiese
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import math as ma
from pathlib import Path
from tqdm import tqdm
from stl import mesh
# import pandas as pd
        
def read_cloudPositions():
    """
    reads cloudPositions and stroes
    - head                         = Everything above the particles positions
    - tail                         = Everything after the particles positions
    - number of initial particles  = integer 
    """
    
    p = Path(__file__).parents[1]

    #read file
    with open(p.joinpath('2_TRACKING','constant','cloudPositions'),'r') as file:
        file = file.readlines()
        
        #get file head
        head = file[0:17]
        
        #get number of initial particles
        numParticles = 0
        for line in file:
            if line[0] == '(':
                if line[-2] == ')':
                    numParticles += 1
                    
        #get file tail
        tail = file[16+numParticles+1:]
    
    return numParticles, head, tail

def read_cloudProperties():
    """
    reads cloudProperties and stroes
    - head                         = Everything above the injectionmodels
    - model                        = InjectionSettings
    - tail                         = Everything after the injectionmodels
    """
    
    p = Path(__file__).parents[1]
    with open(p.joinpath('2_TRACKING','constant','cloudProperties'),'r') as file:
        file = file.readlines()
        # print(file)
        # print(file.index('    injectionModels\n'))
    
    #get head
    head = file[0:file.index('    injectionModels\n')+2]
    
    #get model
    model = []
    stop = 0
    for line in file[len(head):]:
        if stop == 3:
            break
        
        model.append(line)
        
        if '}' in line:
            stop += 1
    
    #get tail
    tail = file[len(head)+len(model):]
    
    return head, tail, model

def read_mapping(Map='AVE'):
    """
    reads mappings
    """
    p = Path(__file__).parents[1]
    
    if Map == 'SUM':
        with open(p.joinpath('1_MAPPING','Sum_Growth'),'r') as file:
            file = file.readlines()
    else:
        with open(p.joinpath('1_MAPPING','Ave_Growth'),'r') as file:
            file = file.readlines()
    
    mappings = []
    for i in range(int(file[20])):
        mappings.append(int(float(file[i+22].replace('\n',''))))

    return mappings

def read_probabilites(NumZones):
    p = Path(__file__).parents[1]
    v = np.genfromtxt(p.joinpath('2_TRACKING',f'long_term_steady_state{NumZones}'),delimiter=',')
    return v

def latest_time():
    return 125000

def read_postions():
    """
    reads positons given a C file
    """
    p = Path(__file__).parents[1]
    with open(p.joinpath('1_MAPPING',f'{latest_time()}','C'),'r') as file:
        file = file.readlines()
        
        c = 0
        for line in file:
            c += 1
            try:
                positions_in_mesh = int(line.replace('\n',''))
                c+=1
                break
            except:
                continue
    
    positions = []
    for i in range(positions_in_mesh):
        positions.append(file[c+i])    

    return positions

def calculate_cell_addition(cells, substrate, Volume, mu_max, Ks, uptake, dt=0.005):
    """
    Calculates when to add new particles to the particle cloud
    """
    results = {}
    prev_cells = cells
    t, c, key = 0, 0, 1
    
    while substrate*Volume > 0:
        
        #update time
        t += dt
        
        #update substrate
        for cell in range(prev_cells):
            substrate -= uptake*dt
        
        # update mu_t
        mu_t = mu_max*(substrate/(substrate+Ks))
        
        # determine neww number of cells
        new_cells = ma.floor(cells*np.exp(mu_t*t))
        
        if new_cells > prev_cells:
            key += 1
            # print(f'Add Cells at time step {round(t,4)}')
            results[key] = (round(t,4),new_cells-prev_cells)
        
        c += 1
        if c == 100_000:
            print('Still Running')
            c = 0
        prev_cells = new_cells
    
    #plot the results
    cells_now = cells
    plt.plot(0,cells,'g.')
    for entry in results:
        cells_now += results[entry][1]
        plt.plot(results[entry][0],cells_now,'r.')
    plt.title('Growth Prediction')
    plt.xlabel('Time [s]')
    plt.ylabel('Cells [-]')
    plt.savefig('Growth_plot.png')
    
    return results

def write_injection_file(PositionFile_Head, PositionFile_Tail, Properties_Head, 
                         Properties_Tail, Properties_Model, potPositions,
                         mapping,long_term_steady_state, particlesToAdd):
    """
    writes all cloudPostions files and the cloudProperties file
    
    Takes: 
        PostionFile_Head, PositionFile_Tail                --> list, list
        Properties_Head, Properties_Tail, Properties_Model --> list, list, list
        potPositions                                       --> list
        mapping                                            --> list
        long_term_steady_state                             --> np.array
        particlesToAdd                                     --> dict (tuple (float, int))
    Generates:
        n cloudPostions files
        cloudProperties file
    Returns:
        Nothing
    """
    
    selections =[]
    #get path 1 layer above
    p = Path(__file__).parents[0]
    
    #prepare file elements (no mod needed)
    Head_Pos = PositionFile_Head
    Tail_Pos = PositionFile_Tail
    
    #prepare file elements
    Head_Pro, Tail_Pro = Properties_Head,Properties_Tail
    Model = Properties_Model
    
    #generate choice options
    choice_options = []
    for i in range(1,len(long_term_steady_state)+1):
        choice_options.append(i)
    
    #correct long-term steady-state
    c = 0
    for probability in long_term_steady_state:
        if probability < 1e-8:
            long_term_steady_state[c] = 0
        c += 1
    
    #wirte for each entry a new model to the cloudPositionsfile
    print('adding particles\n')
    
    with open(p.joinpath('constant','cloudProperties'),'w') as cPro:
        for line in Head_Pro:
            cPro.write(line)
        for line in Model:
            cPro.write(line)
    
        pre_mod_num = '        model1\n'
        pre_time = '            SOI             0;\n'
        pre_cloud_pos = '            positionsFile   "cloudPositions";\n'
        
        for key in tqdm(particlesToAdd.keys()):
            
            sub_model = Model
            new_time = particlesToAdd[key][0]
            sub_model[0] = sub_model[0].replace(pre_mod_num,f'        model{key+1}\n')
            sub_model[6] = sub_model[6].replace(pre_time,f'            SOI             {new_time};\n')
            sub_model[7] = sub_model[7].replace(pre_cloud_pos,f'            positionsFile   "cloudPositions{key}";\n')
            
            with open(p.joinpath('constant',f'cloudPositions{key}'),'w') as cPos:
                
                for line in Head_Pos:
                    cPos.write(line)

                for particle in range(particlesToAdd[key][1]):
                    
                    field_selection = np.random.choice(choice_options, p = long_term_steady_state)
                    indices = [index for index, element in enumerate(mapping) if element == field_selection]
                    point_selection_index = np.random.randint(0,len(indices)-1)
                    point = potPositions[point_selection_index]
                    cPos.write(point)
                    selections.append(point)
                for line in Tail_Pos:
                    cPos.write(line)
            
            pre_mod_num = sub_model[0]
            pre_time = sub_model[6]
            pre_cloud_pos = sub_model[7]
            
            for line in sub_model:
                cPro.write(line)
            
        
        for line in Tail_Pro:
            cPro.write(line)     
    
    return selections

def check_injections(values):
    
    a = []
    for line in values: a.append(line.strip('\n').strip('(').strip(')').split())
    
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    
    # geometries=['BOTTOM.stl','WALL.stl']
    geometries = ['STIRRER.stl','HOLDER.stl']
   
    for geometry in geometries:
        ma_mesh = mesh.Mesh.from_file(geometry)
        ax.add_collection(mplot3d.art3d.Poly3DCollection(ma_mesh.vectors, color = 'gray'))
    
    for entry in a:
        ax.plot3D(float(entry[0]),float(entry[1]),float(entry[2]))
    
    plt.xlim((-0.1,0.1))
    plt.ylim((-0.1,0.1))
    ax.set_zlim((0.05,0.2))
    ax._axis3don = False
    return    

def t0_to_2t0(cells, substrate, Volume, mu_max, Ks, uptake, dt=0.005):
    """
    Calculates when to add new particles to the particle cloud
    """
    results = {}
    prev_cells = cells
    t, c, key = 0, 0, 0
    
    while substrate*Volume > 0:
        
        #update time
        t += dt
        
        #update substrate
        for cell in range(prev_cells):
            substrate -= uptake*dt
        
        # update mu_t
        mu_t = mu_max*(substrate/(substrate+Ks))
        
        # determine neww number of cells
        new_cells = ma.floor(cells*np.exp(mu_t*t))
        
        if new_cells > prev_cells:
            
            # print(f'Add Cells at time step {round(t,4)}')
            results[key] = (round(t,4),new_cells-prev_cells)
            key += 1
            
        c += 1
        if c == 100_000:
            print(f'Still Running t={round(t,3)} N={ma.floor(cells*np.exp(mu_t*t))}')
            c = 0
        prev_cells = new_cells
        
        if key > cells:
            break
    
    key = max(results.keys())
    
    target = results[key][0]/2
    # print(target)
    
    #starting for comparision
    close = results[int(key/2)][0]

    
    #compare
    c_key = int(cells/2) -1
    
    for i in range(int(cells/2)):
        if (target-results[c_key][0])**2 < (target- close)**2:
            close = results[c_key][0]
            f_key = c_key
        c_key -= 1    
    
    returndict = {}
    returndict[1] = (dt*.5, f_key)
    returndict[2] = (dt, cells-f_key)
    
    # print(returndict)
    
    return returndict

def check_probability(p):
    """
    Correct probabilites to sum up to 1.0 for later use 
    """
    if sum(p) != 1.0:
        print('p not equal to 1')
        c = 0
        for state in p:
            if state < 1E-10:
                p[c] = 0
            c+=1
        ps = sum(p)
        # print(ps)
        minval = np.min(p[np.nonzero(p)])
        
        c = 0
        for i in p:
            if i == minval:
                p[c] += 1-ps
        
    print(sum(p))
    return p
            
    


def three_day_sim(cells, substrate, Volume, mu_max, Ks, uptake, total_plot_steps = 3):
    """
    Calculates when to add new particles to the particle cloud
    """
    results = {}
    
    
    prev_cells = cells
    t, i, = 0, 0
    
    total_substrate = substrate*Volume
    
    total_time = 60*60*24*3
    
    step = total_time/total_plot_steps
    
    while i <= total_time:
        
        mu_t = mu_max*(total_substrate/(total_substrate+Ks))
        
        results[i] = cells*np.exp(mu_t*i)
        
        total_substrate -= cells*np.exp(mu_t*i)*uptake
        
        #increment index
        i += step

    print(results)
    
    return results

    
def main1(substrate,Volume,mu_max, Ks, uptake, dt):
    print(f'runing scipt:\n{__file__}\n')
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)

    print('get example cloudPostions\n')
    numParticles, PosHead, PosTail = read_cloudPositions()
      
    print('get cloudProperties and model basics\n')
    ProHead, ProTail, model = read_cloudProperties()
       
    print('get the mappings in a list\n')
    mappings = read_mapping()
        
    print('get long-term steady-state\n')
    long_term_steady_state = read_probabilites(3)
    long_term_steady_state = check_probability(long_term_steady_state)
    
    print('get positions\n')
    postions = read_postions()
        
        
    print('calculate important positions\n')
    cells_to_add = calculate_cell_addition(numParticles,substrate,Volume,mu_max, Ks, uptake, dt)

    
    print('write new file\n')
    selections = write_injection_file(PosHead, PosTail, ProHead, ProTail,
                                      model, postions, mappings, 
                                      long_term_steady_state, cells_to_add)
    
    print('checking')
    check_injections(selections)
    
    return selections


def main2(substrate,Volume,mu_max, Ks, uptake, dt):
    print(f'runing scipt:\n{__file__}\n')
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)

    print('get example cloudPostions\n')
    numParticles, PosHead, PosTail = read_cloudPositions()
      
    print('get cloudProperties and model basics\n')
    ProHead, ProTail, model = read_cloudProperties()
       
    print('get the mappings in a list\n')
    mappings = read_mapping()
        
    print('get long-term steady-state\n')
    long_term_steady_state = read_probabilites(3)
    long_term_steady_state = check_probability(long_term_steady_state)
    print(long_term_steady_state, sum(long_term_steady_state))
    print('get positions\n')
    postions = read_postions()
        
    print('calculate important positions\n')
    cells_to_add = t0_to_2t0(numParticles,substrate,Volume,mu_max, Ks, uptake, dt)
    
    print('write new file\n')
    selections = write_injection_file(PosHead, PosTail, ProHead, ProTail,
                                      model, postions, mappings, 
                                      long_term_steady_state, cells_to_add)
    
    return selections
    
    
if __name__ == '__main__':
    
    """
    Settings
    """
    
    substrate = 15.5          # g/L 
    Volume = 5.2              # read from V                         
    uptake = (0.65/3600)*0.45 # g/g
    dt = 0.00005
    mu_max = 0.65/3600        # 0.65 #h-1
    Ks = 0.009 
    
    """
    Runing the Code
    """
    # selections = main1(substrate,Volume,mu_max,Ks,uptake,dt)
    selections = main2(substrate,Volume,mu_max,Ks,uptake,dt)
    # selections = three_day_sim(10, substrate, Volume, mu_max, Ks, uptake)