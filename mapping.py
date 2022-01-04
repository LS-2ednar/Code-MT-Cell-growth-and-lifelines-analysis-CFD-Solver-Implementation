"""
Calculating local energy dissipationrate
"""
import os
import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt

class Organism:
    """
    Organism class to handle the different field calculations
    """
    def __init__(self,name,cell_dameter,sheartolerance,mu_max,Yxs,Ks):
        """
        Attribute       | Datatype | Information 
        name            | string   | name of the organism               [-]
        cell_diamter    | int      | cell diameter                      [m]
        sheartolerance  | int      | shear stress tolerance             [Pa]
        mu_max          | float    | max specific growth factor         [s-1]
        Yxs             | float    | yield coefficient                  [kg/kg]
        Ks              | float    | affinity constant                  []              [kg/m3]
        qs              | float    | substrate consumtion over time     [kg/s]               
        """
        self.name = name
        self.cell_dameter = cell_dameter
        self.sheartolerance = sheartolerance
        self.mu_max = mu_max
        self.Yxs = Yxs
        self.Ks = Ks
        self.q_max = round(self.mu_max/self.Yxs,3)
        
    def __str__(self):
        return f'''----------------------------------------------
Organism:    {self.name}
----------------------------------------------
Cell Diameter:              {self.cell_dameter} [m]
Substrate Uptakerate:       {self.q_max} [g/h]
Shear Tolerance:            {self.sheartolerance} [Pa]
Max Specific Growth Rate:   {self.mu_max} [h-1]
Ks:                         {self.Ks} [kg/m3]
qs max:                     {self.q_max} [kg/m3]
----------------------------------------------
'''
    
def get_numbered_directories():
    """
    generates a list of numbered directories
    """
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
    """
    Returns latesttime of a subfolder
    """
    return get_numbered_directories()[-1]

def read_locations(path):
    """
    Given a path to an locations file the epsilons are put in a list
    """
    
    locations = []
    with open(path,'r') as file:
        
        # readfile
        a = file.read()
        b = a.split('\n')
        
        # loop until
        loop = int(b[19])+21
        
        for i in range(21,loop):
            c = b[i].replace('(','').replace(')','').split(' ')
            loc = [float(c[0]), float(c[1]), float(c[2])]
            locations.append(loc)
                 
    return locations

def read_epsilons(path):
    """
    Given a path to an epsilon file the epsilons are put in a list
    """
    epsilons = []
    with open(path,'r') as file:
        
        # readfile
        a = file.read()
        b = a.split('\n')
        
        # loop until
        loop = int(b[19])+21
        
        for i in range(21,loop):
            epsilons.append(float(b[i]))
                 
    return epsilons

def read_V(path):
    """
    Works identical tehn read_pressures. Difference is output list is a list of folats
    """
    return read_epsilons(path)

def read_gradU(path):
    """
    Works identical tehn read_locations. Difference is output list is a list of floats
    """
    grad = []
    with open(path,'r') as file:
        
        # readfile
        a = file.read()
        b = a.split('\n')
        
        # loop until
        loop = int(b[19])+21
        for i in range(21,loop):
            c = b[i].replace('(','').replace(')','').split(' ')
            gra = [float(c[0]), float(c[1]), float(c[2]), float(c[3]),float(c[4]),float(c[5]),float(c[6]),float(c[7]),float(c[8])]
            grad.append(gra)
                 
    return grad

def read_U(path):
    """
    Works identical to read_locations. Difference is output list is a list of floats
    """
    return read_locations(path)

def read_tail(File):
    """
    Gets tail of a file
    """
    
    #initialize a counter and tail of a document
    c, tail = 0, ')\n;\n\n'
    
    #try to read the file in windows way
    try:
        with open(f'0\\{File}','r') as f:
            for line in f:
                c += 1
                if c > 21:
                    if 'value' in line:
                        continue
                    elif 'symmetry' in line:
                        tail += '        type symmetry;\n'
                        
                    elif 'type' in line:
                        tail += '        type zeroGradient;\n'
                    else:
                        tail += line
    except:
        with open(f'0/{File}','r') as f:
            for line in f:
                c += 1
                if c > 21:
                    if 'value' in line:
                        continue
                    elif 'symmetry' in line:
                        tail += '        type symmetry;\n'
                        
                    elif 'type' in line:
                        tail += '        type zeroGradient;\n'
                    else:
                        tail += line
    
    return tail

def calculate_local_kolmogorov(epsilons, kinematicViscosity = 8.005e-7, unit = 'µm'): #assuming 30°C
    """
    Calculates local kolmogorov length scale based on the given kinematic viscosity and energydissipation
    rate and stores the values in a list
    """
    kolmogorov = []
    if unit == 'µm':
        print('Kolmogorv length scale in µm')
    else:
        print('Kolmogorv length scale in m')
        
    for i in range(0,len(epsilons)):
        if unit == 'µm':
            kolmogorov.append((kinematicViscosity/(epsilons[i]**3))**(1/4)*1000000)
        else:
            kolmogorov.append((kinematicViscosity/(epsilons[i]**3))**(1/4))
        
    return kolmogorov

def calculate_stresses():
    
    """
    Calculates normal and shear stress values given a speed file U
    and stores values in lists.
    """
    #access latest time
    latest_time = get_latesttime()
    
    #generate needed data
    try:
        os.system('postProcess -func time')
        os.system('postProcess -func writeCellVolumes')
        os.system('simpleFoam -postProcess -func "grad(U)" -field U')
        print()
    except:
        print('Something went worng')
        
    #open needed files for further calculation
    U = read_U(str(latest_time)+'/U')
    try:
        gradU = read_gradU(str(latest_time)+'/grad(U)')
    except:
        gradU = read_gradU(str(latest_time)+'/grad(U)')
    
    #generate dataframe
    grad = pd.DataFrame(gradU)
    grad['normal'] = np.nan
    grad['shear'] = np.nan
    
    #for progressbar
    c, pb = 0,0
    print('\nCalculating Normal- and Shearforces:')
    for i in range(len(U)):
        
        #determine normal stress
        Umag = (U[i][0]**2 + U[i][1]**2 + U[i][2]**2)**(1/2) + 10**(-15)
        Unull = ((U[i][1]-U[i][2])**2 + (U[i][2]-U[i][0])**2 + (U[i][0]-U[i][1])**2)**(1/2) + 10**(-15)
        N1 = (U[i][0]/Umag*gradU[i][0] + U[i][1]/Umag*gradU[i][3] + U[i][2]/Umag*gradU[i][6]) * U[i][0]/Umag 
        N2 = (U[i][0]/Umag*gradU[i][1] + U[i][1]/Umag*gradU[i][4] + U[i][2]/Umag*gradU[i][7]) * U[i][1]/Umag 
        N3 = (U[i][0]/Umag*gradU[i][2] + U[i][1]/Umag*gradU[i][5] + U[i][2]/Umag*gradU[i][8]) * U[i][2]/Umag 
        #N1 = (Ux/Umag*xx + Uy/Umag*yx + Uz/Umag*zx) * Ux/Umag
        #N2 = (Ux/Umag*xy + Uy/Umag*yy + Uz/Umag*zy) * Uy/Umag
        #N3 =(Ux/Umag*xz + Uy/Umag*yz + Uz/Umag*zz) * Uz/Umag
        N4 = (N1 + N2 + N3)**2
        grad['normal'][i] = (2*N4)**(1/2)
        
        #determine shear stress
        dusdys1 = (U[i][0]/Umag*gradU[i][0] + U[i][1]/Umag*gradU[i][3] + U[i][2]/Umag*gradU[i][6]) * (U[i][1]-U[i][2])/Unull             
        dusdys2 = (U[i][0]/Umag*gradU[i][1] + U[i][1]/Umag*gradU[i][4] + U[i][2]/Umag*gradU[i][7]) * (U[i][2]-U[i][0])/Unull             
        dusdys3 = (U[i][0]/Umag*gradU[i][2] + U[i][1]/Umag*gradU[i][5] + U[i][2]/Umag*gradU[i][8]) * (U[i][0]-U[i][1])/Unull             
        #dusdys1 = (Ux/Umag*xx + Uy/Umag*yx + Uz/Umag*zx) *(Uy - Uz)/Unull
        #dusdys2 = (Ux/Umag*xy + Uy/Umag*yy + Uz/Umag*zy) *(Uz - Ux)/Unull
        #dusdys3 = (Ux/Umag*xz + Uy/Umag*yz + Uz/Umag*zz) *(Ux - Uy)/Unull
        
        dvsdxs1 = ((U[i][1] - U[i][2])/Unull*gradU[i][0] + (U[i][2] - U[i][0])/Unull*gradU[i][3] + (U[i][0] - U[i][1])/Unull*gradU[i][6])*U[i][0]/Umag 
        dvsdxs2 = ((U[i][1] - U[i][2])/Unull*gradU[i][1] + (U[i][2] - U[i][0])/Unull*gradU[i][4] + (U[i][0] - U[i][1])/Unull*gradU[i][7])*U[i][1]/Umag 
        dvsdxs3 = ((U[i][1] - U[i][2])/Unull*gradU[i][2] + (U[i][2] - U[i][0])/Unull*gradU[i][5] + (U[i][0] - U[i][1])/Unull*gradU[i][8])*U[i][2]/Umag 
        #dvsdxs1 = ((Uy - Uz)/Unull*xx + (Uz - Ux)/Unull*yx + (Ux - Uy)/Unull*zx)*Ux/Umag
        #dvsdxs2 = ((Uy - Uz)/Unull*xy + (Uz - Ux)/Unull*yy + (Ux - Uy)/Unull*zy)*Uy/Umag
        #dvsdxs3 = ((Uy - Uz)/Unull*xz + (Uz - Ux)/Unull*yz + (Ux - Uy)/Unull*zz)*Uz/Umag
        
        dusdzs1 = (U[i][0]/Umag*gradU[i][0] + U[i][1]/Umag*gradU[i][3] + U[i][2]/Umag*gradU[i][6]) * (U[i][1]*(U[i][0] - U[i][1])/Umag/Unull - U[i][2]*(U[i][2]-U[i][0])/Umag/Unull) 
        dusdzs2 = (U[i][0]/Umag*gradU[i][1] + U[i][1]/Umag*gradU[i][4] + U[i][2]/Umag*gradU[i][7]) * (U[i][2]*(U[i][1] - U[i][2])/Umag/Unull - U[i][0]*(U[i][0]-U[i][1])/Umag/Unull) 
        dusdzs3 = (U[i][0]/Umag*gradU[i][2] + U[i][1]/Umag*gradU[i][5] + U[i][2]/Umag*gradU[i][8]) * (U[i][0]*(U[i][2] - U[i][0])/Umag/Unull - U[i][1]*(U[i][1]-U[i][2])/Umag/Unull) 
        #dusdzs1 = (Ux/Umag*xx + Uy/Umag*yx + Uz/Umag*zx) * (Uy*(Ux - Uy)/Umag/Unull - Uz*(Uz-Ux)/Umag/Unull)
        #dusdzs2 = (Ux/Umag*xy + Uy/Umag*yy + Uz/Umag*zy) * (Uz*(Uy - Uz)/Umag/Unull - Ux*(Ux-Uy)/Umag/Unull)
        #dusdzs3 = (Ux/Umag*xz + Uy/Umag*yz + Uz/Umag*zz) * (Ux*(Uz - Ux)/Umag/Unull - Uy*(Uy-Uz)/Umag/Unull)
        
        dwsdxs11 = ((U[i][0] - U[i][1])*U[i][1]/Umag/Unull - (U[i][2] - U[i][0])*U[i][2]/Umag/Unull)*gradU[i][0]
        dwsdxs12 = ((U[i][1] - U[i][2])*U[i][2]/Umag/Unull - (U[i][0] - U[i][1])*U[i][0]/Umag/Unull)*gradU[i][3]
        dwsdxs13 = ((U[i][2] - U[i][0])*U[i][0]/Umag/Unull - (U[i][1] - U[i][2])*U[i][1]/Umag/Unull)*gradU[i][6]
        # dwsdxs11 = ((Ux - Uy)*Uy/Umag/Unull - (Uz - Ux)*Uz/Umag/Unull)*xx
        # dwsdxs12 = ((Uy - Uz)*Uz/Umag/Unull - (Ux - Uy)*Ux/Umag/Unull)*yx            
        # dwsdxs13 = ((Uz - Ux)*Ux/Umag/Unull - (Uy - Uz)*Uy/Umag/Unull)*zx
        dwsdxs1g = (dwsdxs11 + dwsdxs12 + dwsdxs13)*U[i][0]/Umag;
        
        dwsdxs21 = ((U[i][0] - U[i][1])*U[i][1]/Umag/Unull - (U[i][2] - U[i][0])*U[i][2]/Umag/Unull)*gradU[i][1]
        dwsdxs22 = ((U[i][1] - U[i][2])*U[i][2]/Umag/Unull - (U[i][0] - U[i][1])*U[i][0]/Umag/Unull)*gradU[i][4]
        dwsdxs23 = ((U[i][2] - U[i][0])*U[i][0]/Umag/Unull - (U[i][1] - U[i][2])*U[i][1]/Umag/Unull)*gradU[i][7]
        # dwsdxs21 = ((Ux - Uy)*Uy/Umag/Unull - (Uz - Ux)*Uz/Umag/Unull)*xy
        # dwsdxs22 = ((Uy - Uz)*Uz/Umag/Unull - (Ux - Uy)*Ux/Umag/Unull)*yy
        # dwsdxs23 = ((Uz - Ux)*Ux/Umag/Unull - (Uy - Uz)*Uy/Umag/Unull)*zy
        dwsdxs2g = (dwsdxs21 + dwsdxs22 + dwsdxs23)*U[i][1]/Umag;
        
        dwsdxs31 = ((U[i][0] - U[i][1])*U[i][1]/Umag/Unull - (U[i][2] - U[i][0])*U[i][2]/Umag/Unull)*gradU[i][2]
        dwsdxs32 = ((U[i][1] - U[i][2])*U[i][2]/Umag/Unull - (U[i][0] - U[i][1])*U[i][0]/Umag/Unull)*gradU[i][5]
        dwsdxs33 = ((U[i][2] - U[i][0])*U[i][0]/Umag/Unull - (U[i][1] - U[i][2])*U[i][1]/Umag/Unull)*gradU[i][8]
        # dwsdxs31 = ((Ux - Uy)*Uy/Umag/Unull - (Uz - Ux)*Uz/Umag/Unull)*xz
        # dwsdxs32 = ((Uy - Uz)*Uz/Umag/Unull - (Ux - Uy)*Ux/Umag/Unull)*yz
        # dwsdxs33 = ((Uz - Ux)*Ux/Umag/Unull - (Uy - Uz)*Uy/Umag/Unull)*zz
        dwsdxs3g = (dwsdxs31 + dwsdxs32 + dwsdxs33)*U[i][2]/Umag;
        TermS1 = (dusdys1 + dusdys2 + dusdys3 + dvsdxs1 + dvsdxs2 + dvsdxs3)**2
        TermS2 = (dusdzs1 + dusdzs2 + dusdzs3 + dwsdxs1g + dwsdxs2g + dwsdxs3g)**2
        grad['shear'][i]=(TermS1 + TermS2)**(1/2)
        
        c+= 1
        
        if c == len(U)//50:
                pb += 1
                sys.stdout.write('\r')
                sys.stdout.write("[%-50s] %d%%" % ('='*pb, pb*2))
                sys.stdout.flush()
                c = 0
                
    # print(grad.head(5))
    # print(U)
    # print(gradU)
    # print(latest_time)
    return grad['normal'].values.tolist(), grad['shear'].values.tolist()

def calculate_media_dist(Media_Con):
    """
    Calculates media distribution given the U file by taking the inverse of it
    Media_Con in g/L
    """
    
    #access latest time
    latest_time = get_latesttime()
    
    U = read_U(str(latest_time)+'/U')
    Vol = read_V(str(latest_time)+'/V')  
    
    #calculate individual speed maginitude
    U_inv = []
    for entry in U:
        U_inv.append(1/((entry[0]**2+entry[1]**2+entry[2]**2)**(1/2)))
    U_inv = np.array(U_inv)
    U_inv_tot = sum(U_inv)
    
    #calculate total vloume
    Vol_tot = sum(Vol)                          # m3 Reactor
    media_tot = (Media_Con*1000)*Vol_tot        # kg Substrate
    
    # #determine speed percentage
    U_inv_tot = sum(U_inv)                      # 100%
    speed_percent = np.array(U_inv)/U_inv_tot   # cell%
    
  
    #media dist per cell kg/m3
    media_parts = media_tot*speed_percent
      
    
    return media_parts

def write_kolmogorov(list_of_data,tail):
    """
    Takes list a list of kolmogorov length scale values and writes a OpenFOAM file which is
    interpretable by paraview
    """
    head=f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    {get_latesttime()};
    object      kolmogorov;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 0 0 0 0 0];

internalField   nonuniform List<scalar>""" 
    middle = f"""\n{len(list_of_data)}
(
"""
    try:
        with open(f'{get_latesttime()}/kolmogorov','w') as file:
            file.write(head)
            file.write(middle)
            for element in list_of_data:
                file.write(str(element)+'\n')
            file.write(tail)
            file.close()
    except:
        with open(f'{get_latesttime()}\\kolmogorov','w') as file:
            file.write(head)
            file.write(middle)
            for element in list_of_data:
                file.write(str(element)+'\n')
            file.write(tail)
            file.close()       
    return

def write_normalstress(list_of_data,tail):
    """
    Takes list a list of normalstress values and writes a OpenFOAM file which is
    interpretable by paraview
    """
    head=f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    {get_latesttime()};
    object      normalstress;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [1 -1 -2 0 0 0 0];

internalField   nonuniform List<scalar>""" 
    middle = f"""\n{len(list_of_data)}
(
"""
    try:
        with open(f'{get_latesttime()}/normalstress','w') as file:
            file.write(head)
            file.write(middle)
            for element in list_of_data:
                file.write(str(element)+'\n')
            file.write(tail)
            file.close()
    except:
        with open(f'{get_latesttime()}\\normalstress','w') as file:
            file.write(head)
            file.write(middle)
            for element in list_of_data:
                file.write(str(element)+'\n')
            file.write(tail)
            file.close()       
    return

def write_shearstress(list_of_data,tail):
    """
    Takes list a list of shear stress values and writes a OpenFOAM file which is
    interpretable by paraview
    """
    head=f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    {get_latesttime()};
    object      shearstress;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [1 -1 -2 0 0 0 0];

internalField   nonuniform List<scalar>""" 
    middle = f"""\n{len(list_of_data)}
(
"""
    try:
        with open(f'{get_latesttime()}/shearstress','w') as file:
            file.write(head)
            file.write(middle)
            for element in list_of_data:
                file.write(str(element)+'\n')
            file.write(tail)
            file.close()
    except:
        with open(f'{get_latesttime()}\\shearstress','w') as file:
            file.write(head)
            file.write(middle)
            for element in list_of_data:
                file.write(str(element)+'\n')
            file.write(tail)
            file.close()       
    return

def write_media_dist(list_of_data,Media_Con,tail):
    """
    Takes list a list of kolmogorov length scale values and writes a OpenFOAM file which is
    interpretable by paraview
    """
    head=f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    {get_latesttime()};
    object      media_{Media_Con};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [1 -3 0 0 0 0 0];

internalField   nonuniform List<scalar>""" 
    middle = f"""\n{len(list_of_data)}
(
"""
    try:
        with open(f'{get_latesttime()}/media_{Media_Con}','w') as file:
            file.write(head)
            file.write(middle)
            for element in list_of_data:
                file.write(str(element)+'\n')
            file.write(tail)
            file.close()
    except:
        with open(f'{get_latesttime()}\\media_{Media_Con}','w') as file:
            file.write(head)
            file.write(middle)
            for element in list_of_data:
                file.write(str(element)+'\n')
            file.write(tail)
            file.close()       
    return

def write_all_zones(Org, kolmogorov, normalstresses, shearstresses, media, concentration, tail):
    """
    All DangerZones for each Organism are based on the values given by the 
    individual parameters, they are treaded as sorted lists where least 
    dangerous values greated with number 1 and more risky values are grated 
    numbers up to 3 which imply high risk.
    
    For the individual analized parameters the following is true:
        kolmogorov length scale --> the smaller the value the higher the risk
        normalstress            --> the higher the value the higher the risk
        shearstress             --> the higher the value the higher the risk
        media                   --> the smaller the value the higher the risk 
    
    """
    head0=f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    {get_latesttime()};
    object      DangerZonesAveCombined_con_{concentration}_{Org.name};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar>""" 

    head01=f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    {get_latesttime()};
    object      DangerZonesSumCombined_con_{concentration}_{Org.name};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar>""" 

    head02=f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    {get_latesttime()};
    object      Ave_Growth;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar>""" 

    head03=f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    {get_latesttime()};
    object      Sum_Growth;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar>""" 

    head1=f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    {get_latesttime()};
    object      DangerZonesKolmogorov_{Org.name};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar>""" 

    head2=f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    {get_latesttime()};
    object      DangerZonesNormalstress_{Org.name};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar>""" 

    head3=f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    {get_latesttime()};
    object      DangerZonesShearstress_{Org.name};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar>""" 

    middle = f"""\n{len(kolmogorov)}
(
"""

    head4=f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    {get_latesttime()};
    object      DangerZones_qs_{concentration}_{Org.name};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar>""" 

    #get length of data
    len_data = len(kolmogorov)
    # #determine stepsze for zones
    # step = len_data//3 
    
    #kolmogorov zones
    kolZones = []
    #cellsize for interpretation 
    cell_size = Org.cell_dameter
    
    #zone definitions
    zone1 = 5*cell_size
    zone2 = 2*cell_size
    
    for i in kolmogorov:
        if i >= zone1:
            kolZones.append(1)
        elif i >= zone2:
            kolZones.append(2)
        else:
            kolZones.append(3)
    
    #normalstress zones
    normZones = []
    #maximal normal and shear force endurable
    norm_max = shear_max = Org.sheartolerance

    #zone definitions
    zone1 = norm_max*0.05
    zone2 = norm_max*0.95
    
    for i in normalstresses:
        if i <= zone1:
            normZones.append(1)
        elif i <= zone2:
            normZones.append(2)
        else:
            normZones.append(3)
    
    #shearstress zones
    shearZones = []
    
    #zone definitions
    zone1 = shear_max*0.05
    zone2 = shear_max*0.95
    
    for i in shearstresses:
        if i <= zone1:
            shearZones.append(1)
        elif i <= zone2:
            shearZones.append(2)
        else:
            shearZones.append(3)
    
    #media zones
    mediaZones = []
        
    # evaluation of mediaZones
    q_max = Org.q_max            #see paper malairuang_high_2020     #kg/s
    Ks    = Org.Ks               #from g/L to mol/kg assuming        #kg/m3
                                   
    #zone definitions
    zone1 = 0.95*q_max
    zone2 = 0.05*q_max
    
    for i in media:
        
        m = q_max*((i*1000)/((i*1000)+Ks))
        
        if  m >= zone1:
            mediaZones.append(1)
        elif m >= zone2:
            mediaZones.append(2)
        else:
            mediaZones.append(3)   
    
    
    #combined zones average and summed-up
    combinedAveZones, combinedSumZones = [],[]
    #summantion of all zones divided by 3
    
    for i in range(len_data):
        combinedAveZones.append(round((kolZones[i]+normZones[i]+shearZones[i]+mediaZones[i])/4,0))
        combinedSumZones.append(kolZones[i]+normZones[i]+shearZones[i]+mediaZones[i]-3)
    #Writing the files
    
    #combinedAveZones
    try:
        with open(f'{get_latesttime()}/DangerZonesAveCombined_{concentration}g_{Org.name}','w') as file:
            file.write(head0)
            file.write(middle)
            for element in combinedAveZones:
                    file.write(str(element)+'\n')
            file.write(tail)
            file.close()
        with open('Ave_Growth','w') as file:
            file.write(head02)
            file.write(middle)
            for element in combinedAveZones:
                    file.write(str(element)+'\n')
            file.write(tail)
            file.close()
    except:
        with open(f'{get_latesttime()}\\DangerZonesAveCombined_{concentration}g_{Org.name}','w') as file:
            file.write(head0)
            file.write(middle)
            for element in combinedAveZones:
                    file.write(str(element)+'\n')
            file.write(tail)
            file.close()
        with open('Ave_Growth','w') as file:
            file.write(head02)
            file.write(middle)
            for element in combinedAveZones:
                    file.write(str(element)+'\n')
            file.write(tail)
            file.close()
    #combinedSumZones
    try:
        with open(f'{get_latesttime()}/DangerZonesSumCombined_{concentration}g_{Org.name}','w') as file:
            file.write(head01)
            file.write(middle)
            for element in combinedSumZones:
                    file.write(str(element)+'\n')
            file.write(tail)
            file.close()
        with open('Sum_Growth','w') as file:
            file.write(head03)
            file.write(middle)
            for element in combinedSumZones:
                    file.write(str(element)+'\n')
            file.write(tail)
            file.close()
    except:
        with open(f'{get_latesttime()}\\DangerZonesSumCombined_{concentration}g_{Org.name}','w') as file:
            file.write(head01)
            file.write(middle)
            for element in combinedSumZones:
                    file.write(str(element)+'\n')
            file.write(tail)
            file.close()
        with open('Sum_Growth','w') as file:
            file.write(head03)
            file.write(middle)
            for element in combinedSumZones:
                    file.write(str(element)+'\n')
            file.write(tail)
            file.close()
            
    #kolmogorovZones
    try:
        with open(f'{get_latesttime()}/DangerZonesKolmogorov_{Org.name}','w') as file:
            file.write(head1)
            file.write(middle)
            for element in kolZones:
                    file.write(str(element)+'\n')
            file.write(tail)
            file.close()
    except:
        with open(f'{get_latesttime()}\\DangerZonesKolmogorov_{Org.name}','w') as file:
            file.write(head1)
            file.write(middle)
            for element in kolZones:
                    file.write(str(element)+'\n')
            file.write(tail)
            file.close()    
    
    #normalstressZones
    try:
        with open(f'{get_latesttime()}/DangerZonesNormalstress_{Org.name}','w') as file:
            file.write(head2)
            file.write(middle)
            for element in normZones:
                    file.write(str(element)+'\n')
            file.write(tail)
            file.close()
    except:
        with open(f'{get_latesttime()}\\DangerZonesNormalstress_{Org.name}','w') as file:
            file.write(head2)
            file.write(middle)
            for element in normZones:
                    file.write(str(element)+'\n')
            file.write(tail)
            file.close()
            
    #shearstressZones
    try:
        with open(f'{get_latesttime()}/DangerZonesShearstress_{Org.name}','w') as file:
            file.write(head3)
            file.write(middle)
            for element in shearZones:
                    file.write(str(element)+'\n')
            file.write(tail)
            file.close()
    except:
        with open(f'{get_latesttime()}\\DangerZonesShearstress_{Org.name}','w') as file:
            file.write(head3)
            file.write(middle)
            for element in shearZones:
                    file.write(str(element)+'\n')
            file.write(tail)
            file.close()  
    
    #qs_Zones
    try:
        with open(f'{get_latesttime()}/DangerZones_qs_{concentration}g_{Org.name}','w') as file:
            file.write(head4)
            file.write(middle)
            for element in mediaZones:
                    file.write(str(element)+'\n')
            file.write(tail)
            file.close()
    except:
        with open(f'{get_latesttime()}\\DangerZones_qs_{concentration}g_{Org.name}','w') as file:
            file.write(head4)
            file.write(middle)
            for element in mediaZones:
                    file.write(str(element)+'\n')
            file.write(tail)
            file.close() 
    
    return mediaZones

def overview_plot(kolmogorov, normal, shear):
    
    fig, ax0 = plt.subplots(figsize = (7,7))
    #first axis
    l1 = ax0.plot(sorted(normal),'-k',label='Normalstress')
    l2 = ax0.plot(sorted(shear),'--k',label='Shearstress')
    
    ax1 = ax0.twinx()
    #second axis
    l3 = ax1.plot(sorted(np.array(kolmogorov)*1_000_000),':k',label='Kolmogorov lengthscale')
    
    #set titles for axis
    ax0.set_ylabel('Normal- and Shearstress [ kg$\cdot$m$^{-1}$$\cdot$s$^{-2}$ ]',fontsize=14)
    ax1.set_ylabel('Kolmogorov lengthscale [ $\mu$m ]',fontsize=14)
    ax0.set_xlabel('Sorted values',fontsize=14)
    # ax0.set_xlim(0,100000)
    # ax0.set_ylim(-20)
    ax1.set_ylim()
    ax0.tick_params(axis='both', which='major', labelsize=12)
    #labels for legend
    ls = l1+l2+l3
    labels = [l.get_label() for l in ls]
    ax0.legend(ls,labels,fontsize=12)
    #plt.show()
    
    fig.savefig('overviewplot')
    
    return


def main(list_of_organisms,list_of_concentraitons):
    
    print(f'Running Script: {__file__}\n')
    #select current working directory
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    
    #get tail of a document
    tail = read_tail('epsilon')
    #get data from latest time point
    try:
        epsilons = read_epsilons(f'{get_latesttime()}\\epsilon')
    except:
        epsilons = read_epsilons(f'{get_latesttime()}/epsilon')
    #calculate kolmogorow values
    kolmogorov = calculate_local_kolmogorov(epsilons, unit = 'm')
    print('\nWriting local_kolmogorov_lengthscale file\n')
    write_kolmogorov(kolmogorov, tail)
    
    #calculating normal- and shearstress
    normal, shear = calculate_stresses()
    print('\n\nWriting local normal stress\n')
    write_normalstress(normal,tail)
    print('\nwriting local shear stress\n')
    write_shearstress(shear,tail)
    
    for organism in list_of_organisms:
        print(organism)
        for concentration in list_of_concentraitons:
            media = calculate_media_dist(concentration)
            write_media_dist(media,concentration,tail)
        
            print(f'\nwriting dangerzones files with concentration {concentration}\n')
            zones = write_all_zones(organism, kolmogorov,normal,shear,media,concentration,tail)
    
    # print('\ngenerating overview plot')
    # ploting values on graph
    # overview_plot(kolmogorov, normal, shear)
    
    return

if __name__ == '__main__':
    
    """
    Defining Organisms and Substrate Concentrations to check
    """
    yeast           = Organism('SC' ,10e-6,2700,0.54 ,0.17/3600, 0.034)  # https://www.researchgate.net/publication/6589044_Substrate_inhibition_kinetics_of_Saccharomyces_cerevisiae_in_fed-batch_cultures_operated_at_constant_glucose_and_maltose_concentration_levels
    cho             = Organism('CHO',15e-6, 400,0.038,0.05/3600, 0.286) 
    
    #puting organisms in a list
    organisms       = [yeast,cho]     #Organisms
    
    concentraitons  = [20, 10, 1]     #kg/m3
    
    """
    Runing the main function
    """
    main(organisms,concentraitons)
    
    # calculate_media_dist(10)
    
