# -*- coding: utf-8 -*-
"""
Generate the Particlelifelines and the state transion 
probability matricies for each particle as well as a 
long term steady state of the system.
"""
import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm

class Particle:
    """
    Position descibed by x, y, z
    Age described by timesteps which are currently updated later
    """
    def __init__(self,x,y,z,name=None):
        #location
        self.name = name
        self.x = x
        self.y = y
        self.z = z
        self.age = 0
        self.value_history = []
        self.value_index = []
        self.state_transistion = None
        self.times = []
        self.doubling = 0
        
    def __str__(self):
        if self.name != None:
            try:
                return f'{self.name}\nx: {self.x}\ny: {self.y}\nz: {self.z}\nage: {self.age}'
            except:
                return f'{self.name}\nx: {self.x}\ny: {self.y}\nz: {self.z}'
        else:
            try:
                return f'x: {self.x}\ny: {self.y}\nz: {self.z}\nage: {self.age}'
            except:
                return f'x: {self.x}\ny: {self.y}\nz: {self.z}'
    
    def update_age(self,dt):
        self.age = len(self.x)*dt
    
#Helper Functions      
def list2intList(a_list):
    """
    Function turns a list of strings to a list of integers
    """
    r_list = []
    for element in a_list:
        r_list.append(int(element))
    return r_list

def get_particle_positions(filename,ps = False):
    """
    Function generates list of Particles given a vtk datafile. Using the Keyword
    Argument ps it is posible to get the particle space.
    """
    
    #initialize lists and variables
    x,y,z, particle_length,particle_list = [],[],[],[],[]
    last, new, particle_number = 0,0,0
    
    #opening given vtk file
    print(f'\nOpening file: {filename:>20}')
    tracks = open(f'{filename}')
    
    #transform vtk to list of particles
    counter = 0
    loop = 0 
    for line in tracks:
        counter += 1
        if counter == 5: 
            loop = line.split(' ')[1]
            print(f'\nReading Datapoints: {loop:>9}')
        
        if counter > 5 and counter <= int(loop)+5:
            line = line.replace('\n','').split(' ')
            # print(line)
            x.append(float(line[0]))
            y.append(float(line[1]))
            z.append(float(line[2]))
            
        if counter == int(loop)+6:
            loop2 = line.split(' ')[1]
        
        if counter > int(loop)+6 and counter <= int(loop)+6+int(loop2):
            
            line = list2intList(line.replace('\n','').split(' '))
            particle_length.append(line[0])

    print(f'\nInitialize creation of {loop2} particles')
        
    for i in particle_length:
        new += i
        particle_list.append(Particle(x[last:new],y[last:new],z[last:new],f'P{particle_number}'))
        # print(x[last:new])
        last = new
        particle_number += 1
    
    def particle_space():
        """
        This defines the particle space given the initial vtk file as a list 
        of tuples of form [(min x, max x ), (min y, max y), (min z, max z)]
        """
        return [(min(x),max(x)),(min(y),max(y)),(min(z),max(z))]
    
    if ps == True:
        return particle_list, particle_space()
    else:
        return particle_list

def get_numbered_directories():
    all_dirs = next(os.walk('.'))[1]
    num_dirs = []
    for di in all_dirs:
        try:
            if '.' not in di:
                num_dirs.append(int(di))
            else:
                num_dirs.append(float(di))
        except:
            continue
    num_dirs.sort()
    return num_dirs

def get_latesttime():
    return get_numbered_directories()[-1]


def num_zones(dangerzones_file):
    """
    Determine the unique values in a list of dangerzones
    """
    values = []
    
    with open(dangerzones_file,'r') as file:
        
        # readfile
        a = file.read()
        b = a.split('\n')
        file.close()
        # loop until
        loop = int(b[20])+22
        
        for i in range(22,loop):
            if b[i] not in values:
                values.append(b[i])
    
    return len(values)

def read_zones(path):
    """
    Read information about cells form a file
    """
    values = []
    with open(path,'r') as file:
        
        # readfile
        a = file.read()
        b = a.split('\n')
        
        # loop until
        loop = int(b[20])+22
        
        
        for i in range(22,loop):
            values.append(int(float(b[i])))
                 
    return values

def read_locations(path):
    """
    Read coordinates of cell centres form a file
    """
    os.system('postProcess -func writeCellCentres -latestTime')  
    locations = []
    
    with open(path,'r') as file:
            
        # readfile
        a = file.read()
        b = a.split('\n')
        # print(b[19])
        # loop until
        loop = int(b[19])+21
        # print(loop)
        for i in range(21,loop):
            c = b[i].replace('(','').replace(')','').split(' ')
            # print(c)
            loc = [float(c[0]), float(c[1]), float(c[2])]
            locations.append(loc)
                 
    return locations

def read_particledata(particledata):
        with open('particle.data', 'rb') as filehandle:
        #read the data from binary data stream
            particles,particle_space = pickle.load(filehandle)
            filehandle.close()
            return particles, particle_space

def particle_lifelines(particle_list,cell_list,zone_list,num_zones):
    """
    Generates the particles history in reagards to the values it had to endure.
    
    Parameters:
    particle_list --> list of Particle
    cell_list     --> list of lists
    zone_list     --> list of integers
    num_zones     --> integer
    
    Output:
    end_particles --> list of Particle
                      These Particles have now stroed information about there
                      value_history which shows where they where exposed to
                      the different zones.
    """
    #for each particle get the closest point from the list of cell centeres.
    end_particles = []
    
    #for progressbar
    c, pb = 0,0
    print('Mapping particles to cell-centers for zone determination')
    #print(cell_list)
    for particle in tqdm(particle_list):
        c+= 1
        for i in range(len(particle.x)):
            points_ = []
            for point in cell_list:
                points_.append(((point[0]-particle.x[i])**2+(point[1]-particle.y[i])**2+(point[2]-particle.z[i])**2)**(1/2))
            #print(points_)
            particle.value_index.append(points_.index(min(points_)))
        end_particles.append(particle)
        
        if c == len(particle_list)//50:
                pb += 1
                sys.stdout.write('\r')
                sys.stdout.write("[%-50s] %d%%" % ('='*pb, 2*pb))
                sys.stdout.flush()
                c = 0
            
            
    #match indecies to zone values
    for particle in particle_list:
        for i in particle.value_index:
            particle.value_history.append(zone_list[i])
            
        # #might be used later
        # zone_dist = {}
        # for key in range(1,num_zones+1):
        #     zone_dist[key] = 0
        #     for element in particle.value_history:
        #         if element == key:
        #             zone_dist[key] += 1
        #     zone_dist[key] = zone_dist[key]/len(particle.x)
                    
        # particle.zone_dist = zone_dist
        
    return end_particles

def particle_state_transions(list_of_particles,num_states):
    """
    Given a list_of_particles and a number of available states (num_states) a state
    transion matrix can be caluclated and is added to each particle.
    """
    
    for particle in list_of_particles:
        #get the list of states for a particle
        list_of_states = particle.value_history
        
        #iniate a nxn zero matrix where n is equal to num_states
        out = np.zeros((num_states,num_states)) 
        
        #increment values in the nxn matrix 
        for i in range(len(list_of_states)-1):
            curr_ = list_of_states[i]
            next_ = list_of_states[i+1]
            out[curr_-1][next_-1] += 1
        
        # get percentages for each value in the system
        for n in range(len(out)):
            line_sum = sum(out[n])
            for m in range(len(out)):
                
                if out[n][m] != 0:
                    out[n][m] = out[n][m]/line_sum
        
        #add the state_transition_probability matrix to the particle
        particle.state_transions = out
    return list_of_particles

def system_state_transion(list_of_particles, Number_of_states):
    """
    Given a list_of_particles the systems average state-transion-probability-matrix (P) is calculated
    """
    P = np.zeros((len(list_of_particles[0].state_transions),len(list_of_particles[0].state_transions))) 
    for particle in list_of_particles:
        P+=particle.state_transions
    
    print(P)    
    for i in range(Number_of_states):
        if sum(P[i]) != 0:
            P[i] = P[i]/sum(P[i])
    
    return P
 
def long_term_steady_state(P):
    """
    Assuming Av=b where v is the long term steady-state-transion probability
    A is equal to the transposed state transion probability matrix with an added line to describe
    the n values where n is representing the number of states
    b is equal to a column vector of shape n+1 where all values are 0 but the last one which is equal to 1
    """
    
    A = np.append(np.transpose(P)-np.identity(len(P)),[np.ones(len(P))],axis=0)
    b = np.transpose(np.append(np.zeros(len(P)),1))
    
    v = np.linalg.solve(np.transpose(A).dot(A),np.transpose(A).dot(b))
    
    for state in range(len(P)):
        print(f'State {state+1}: {v[state]}')
    
    return v

def mc_process_plot(STPM,NumIter=50):
    """
    Given a STPM (state transition probability matrix), a number of cells (NumCells), 
    and a number of iterations (NumIter) the long therm steady state is visualized
    """
    num_states = len(STPM)
    lables = [f'State {i}' for i in range(num_states)]
    
    #state 0 is 100% of pop in the beginnnig
    state=np.array([np.append(1,np.zeros(num_states-1))])
    stateHist=state
    for x in range(NumIter):
        state=np.dot(state,STPM)
        stateHist=np.append(stateHist,state,axis=0)
        dfDistrHist0 = pd.DataFrame(stateHist,columns = lables)
    dfDistrHist0.plot(ylim=(0,1),
                      xlim=(-0.1,NumIter),
                      ylabel = 'Distribution [-]',
                      xlabel = 'State Changes [-]',
                      title = 'longterm steady state\nstarting from State 0')
    plt.show()
    
    #Even distribution of all states
    state=np.array([np.ones(num_states)/num_states])
    stateHist=state
    for x in range(NumIter):
        state=np.dot(state,STPM)
        stateHist=np.append(stateHist,state,axis=0)
        dfDistrHist1 = pd.DataFrame(stateHist,columns = lables)
    dfDistrHist1.plot(ylim=(0,1),
                      xlim=(-0.1,NumIter),
                      ylabel = 'Distribution [-]',
                      xlabel = 'State Changes [-]',
                      title = 'longterm steady state\nstarting from even distribution')
    plt.show()
    
    #state n is 100% of pop in the beginnnig
    state=np.array([np.append(np.zeros(num_states-1),1)])
    stateHist=state
    for x in range(NumIter):
        state=np.dot(state,STPM)
        stateHist=np.append(stateHist,state,axis=0)
        dfDistrHist2 = pd.DataFrame(stateHist,columns = lables)
    dfDistrHist2.plot(ylim=(0,1),
                      xlim=(-0.1,NumIter),
                      ylabel = 'Distribution [-]',
                      xlabel = 'State Changes [-]',
                      title = f'longterm steady state\nstarting from State {num_states-1}')
    plt.show()
    return dfDistrHist1

def main(Num_Zones):
    print(f'Running Script with {Num_Zones} states')
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    
    try: #to read the VTK/particleTracks.vtk file otherwise generate & read it
        particles,particle_space = get_particle_positions('VTK/particleTracks.vtk',ps = True)
    except:
        os.system('particleTracks')
        try:
            particles,particle_space = get_particle_positions('VTK/particleTracks.vtk',ps = True)
        except:
            particles,particle_space = get_particle_positions('VTK\\particleTracks.vtk',ps = True)
        
    #with open('particle.data', 'wb') as filehandle:
    #    # store the data as binary data stream
    #    pickle.dump([particles,particle_space], filehandle)
    
    print('\nReading files for zones, coordinates and particles')
    try:
        os.chdir(os.path.dirname(os.getcwd())+'\\1_MAPPING')
    except:
        os.chdir(os.path.dirname(os.getcwd())+'/1_MAPPING')
    
    if Num_Zones == 3:
        zones = read_zones('Ave_Growth')
        #number_of_zones = num_zones('Ave_Growth')
    else:
        zones = read_zones('Sum_Growth')
        #number_of_zones = num_zones('Sum_Growth')
 
    #print(zones)
    
    #get back to current working direcotry
    os.chdir(dname)
    #print(get_latesttime())    
    try:
        coordinates = read_locations(f'{get_latesttime()}\\C')
    except:
        coordinates = read_locations(f'{get_latesttime()}/C')
        
    # particles, particle_space = read_particledata('particle.data')
    #print('coordinates',coordinates) 
    print('\nDetermine lifelines of particles\n')
    end_particles = particle_lifelines(particles,coordinates,zones,Num_Zones)
    end_particles = particle_state_transions(end_particles,Num_Zones)
    
    P = system_state_transion(end_particles,Num_Zones)
    print(f'\n\nThe State transion matrix for this system looks as follows:\n\n{P}')
    
    print('\n\nWhere the long-term steady-state-solution vector is:\n')
    steady_state = long_term_steady_state(P)
    #dist = mc_process_plot(P)
    
    print('\n\n')
    return P
         
if __name__ == '__main__':
    main(3)    
    main(9)
