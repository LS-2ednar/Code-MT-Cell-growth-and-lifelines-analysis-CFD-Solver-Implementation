# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 21:10:14 2021

@author: wiese
"""

"""
Generate a kinematicCloud File to work with openFoam this
"""
import os

def get_internalmesh():
    
    #try to read the C file in 0
    if os.path.isfile('0\C'):
        with open('0\C','r') as cFile:
            intMesh = cFile.read().split('\n')
            elements = 20+int(intMesh[19])
            
            cFile.close()
    
    else:
        print('No C file found trying to create it')
        try:
            os.system('postProcess -func writeCellCentres -latestTime')
            with open('0\C','r') as cFile:
                intMesh = cFile.read().split('\n')
                elements = 20+int(intMesh[19])
                cFile.close()
        except:
            return 'something went wrong'
        
    
    return intMesh[21:elements]


def get_range_xy(internalMesh):
    
    check = int(round(len(internalMesh)/10,0))
    
    x_s, y_s = [], []
    for element in internalMesh:
        # print(element)
        entrie = element.replace('(','').replace(')','').split(' ')
        # print(entrie)
        x_s.append(float(entrie[0]))
        y_s.append(float(entrie[1]))
    
    x_s, y_s = sorted(x_s), sorted(y_s)
    
    x_min,x_max,y_min,y_max = x_s[0],x_s[-1],y_s[0],y_s[-1]
    
    c, temp = 0,0
    for min_ in x_s:
        if min_ > x_min and temp != min_:
            c += 1
            temp = min_
            if c > check:
                break

    c, temp = 0, 999
    for max_ in x_s[::-1]:
        if max_ < x_max and temp != max_:
            c += 1
            temp = max_
            if c > check:
                break
    #write x_range
    x_range = (min_, max_)
    
    c, temp = 0,0
    for min_ in y_s:
        if min_ > y_min and temp != min_:
            c += 1
            temp = min_
            if c > check:
                break
    c, temp = 0, 999
    for max_ in y_s[::-1]:
        if max_ < y_max and temp != max_:
            c += 1
            temp = max_
            if c > check:
                break
    #write y_range
    y_range = (min_, max_)    
    
    print(x_range,y_range,x_min,y_min)
    print(f'x_range: {x_range} x_min={x_min} x_max={x_max}\ny_range: {y_range} y_min={y_min} y_max:{y_max}')
    return x_range, y_range
    

def set_particles(number_of_particles,internalMesh,x_range,y_range):
    
    #calculations for startposition seleciton
    pos_postions = len(internalMesh)
    step_difference = pos_postions // number_of_particles
    print(step_difference) 
    #setup for the file to be written
    head = """/*--------------------------------*- C++ -*----------------------------------*\\
 =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  9
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       vectorField;
    object      cloudPositions;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

(
"""
    bottom = ')\n\n// ************************************************************************* //'
    
    #writing the file 
    with open('cloudPositions','w') as file:
        file.write(head)
        index, moved_particles = 0,0
        for i in range(number_of_particles):
            check_x = float(internalMesh[index].replace('(','').replace(')','').split(' ')[0])
            check_y = float(internalMesh[index].replace('(','').replace(')','').split(' ')[1])
            if  (check_x > x_range[0]) and (check_x < x_range[1]) and (check_y > y_range[0]) and (check_y < y_range[1]):

                # print(type(internalMesh[index].split(' ')[0]))
                file.write(internalMesh[index]+'\n')
                index += step_difference
            else:
                print('Particle to close to Wall move away')
                step = 0
                while (check_x < x_range[0]) and (check_x > x_range[1]) and (check_y < y_range[0]) and (check_y > y_range[1]):
                    step += 1
                    check_x = float(internalMesh[index+step].replace('(','').replace(')','').split(' ')[0])
                    check_y = float(internalMesh[index+step].replace('(','').replace(')','').split(' ')[1])
                print('Particle moved... placing next particle')
                file.write(internalMesh[index]+'\n')
                index += step_difference
                moved_particles += 1
        file.write(bottom)
        print(f'Moved {moved_particles} Particles')
    return 

if __name__ == '__main__':
    #set current working directory to file location
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)   
    
    mesh = get_internalmesh()
    x_range,y_range = get_range_xy(mesh)
    #set a number of particles in the internal mesh at the cell centers given.
    set_particles(20,mesh, x_range,y_range)
    # print(get_range_xy(get_internalmesh()))
