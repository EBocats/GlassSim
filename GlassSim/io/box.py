import numpy as np
import os
os.environ['OVITO_GUI_MODE'] = '1'
from ovito.io import import_file

def GetBox(fin: str, mult = False):
    pipline = import_file(fin, multiple_frames = mult)
    if mult:
        box = []
        for frame in range(pipline.source.num_frames):
            data = pipline.compute(frame)
            boxi = data.cell.matrix
            box.append(boxi)
    else:
        data = pipline.compute()
        box = data.cell.matrix
    return box
            

def GetLogBoxLen(logfile: str, line):
    linei = open(logfile).readlines()[line]
    L = linei.split()[-1] # cubic box length
    return float(L)


def scale_the_box(xyz = None, L = None, ratio = 1.0, norm = False):
    
    if norm:
        box_o = [L[0,3], L[1,3], L[2,3]]
        box_l = [L[0,0], L[1,1], L[2,2]]
        xyz = xyz - box_o
        xyz = xyz / box_l
        L = np.diag([1.0, 1.0, 1.0])
        L_o = [0.0, 0.0, 0.0]
        L = np.column_stack((L, L_o))
    else:
        xyz = xyz * ratio
        L = L * ratio
    return xyz , L




##############################################Old version##############################################

from GlassSim.io.input import read_stru
from GlassSim.utils.distence import distence_3d

def joint_box(stru1_file, stru2_file, plane = None, type_map = None, gap = 1.0, charge = False, chg_map = None, molecule = False, mol_cuts = None):
    '''
    Joint two structures in a specific plane, and generate a new structure.
    For example, joint a x-z plane of two melten soilds.
    '''

    if plane == None:
        print('No plane specified, Must specify a plane! eg: [1, 0, 1] for x-z plane')
        return
    if type_map == None:
        print('No type map specified, using the same type map!')
    
    type1, xyz1, cell1 = read_stru(stru1_file)
    type2, xyz2, cell2 = read_stru(stru2_file)

    # read cell_matrixs and transform its to a new cell_matrix
    new_cell = np.zeros((3, 2))
    for i in range(3):
        if plane[i] == 1:
            new_cell[i, 0] = cell1[i, -1]
            new_cell[i, 1] = np.max([cell1[i, i], cell2[i, i]]) + new_cell[i, 0]
        else:
            new_cell[i, 0] = cell1[i, -1]
            new_cell[i, 1] = cell1[i, i] + cell2[i, i] + gap + new_cell[i, 0]

    # transform the xyz2 to new xyz2 and joint xyz1 and new xyz2
    new_xyz2 = np.zeros_like(xyz2)
    for i in range(3):
        if plane[i] == 0:
            new_xyz2[:, i] = xyz2[:, i] + cell1[i, i] + new_cell[i, 0] + gap
        else:
            new_xyz2[:, i] = xyz2[:, i] + new_cell[i, 0]
    xyz = np.vstack((xyz1, new_xyz2))
    
    # map the type of atoms
    if type_map == None:
        ts = set(type1)
        type_map = [(i, i) for i in range(1, len(ts)+1)]

    types = [type1, type2]
    new_types = [np.zeros_like(types[0]), np.zeros_like(types[1])]
    for i in range(2):
        map = type_map[i]
        for j in range(len(types[i])):
            new_types[i][j] = map[types[i][j]-1]
    type = np.hstack((new_types[0], new_types[1]))

    if molecule == True and charge == False:
        print('Charge must be specified when molecule is True!')
        return

    # generate charges
    if charge:
        charges = np.zeros(len(type))
        for i in range(len(type)):
            charges[i] = chg_map[type[i]-1]
        
        if molecule != True:
            return type, charges, xyz, new_cell

    # generate molecule id
    if molecule:
        if mol_cuts == None:
            print('No molecule cuts specified, will use the 0.5-1.5 AA as the default!')
            mol_cuts = [(0.5, 1.5)]
        else:
            dis = np.zeros((len(xyz), len(xyz)))
            for i in range(len(xyz)):
                for j in range(i+1, len(xyz)):
                    ls = [new_cell[i, 1] - new_cell[i, 0] for i in range(3)]
                    dis[i, j] = distence_3d(xyz[i], xyz[j], pbc=True, Ls=ls)
            mols = np.zeros(len(xyz))
            bonds = np.zeros((len(xyz), 2))
            mol_id = 1
            for i in range(len(xyz)):
                if mols[i] == 0:
                    mols[i] = mol_id
                    js = []
                    for j in range(i+1, len(xyz)):
                        for mol_cut in mol_cuts:
                            if dis[i, j] < mol_cut[1] and dis[i, j] > mol_cut[0]:
                                js.append(j)
                                bonds[j, 0] = i
                                bonds[j, 1] = j
                    if len(js) > 0:
                        for j in js:
                            mols[j] = mol_id
                        mol_id += 1
                    else:
                        mols[i] = 0

            bonds = bonds[bonds[:, 0] != 0]
            bonds = bonds.astype(int)
            mols = mols.astype(int)
            
            # generate angles using the bonds
            angles = []
            for i in range(len(bonds)):
                for j in range(i+1, len(bonds)):
                    if bonds[i, 1] == bonds[j, 0]:
                        angles.append([bonds[i, 0], bonds[i, 1], bonds[j, 1]])
                    if bonds[i, 0] == bonds[j, 0]:
                        angles.append([bonds[i, 1], bonds[i, 0], bonds[j, 1]])
                    if bonds[i, 0] == bonds[j, 1]:
                        angles.append([bonds[i, 1], bonds[i, 0], bonds[j, 0]])
                    if bonds[i, 1] == bonds[j, 1]:
                        angles.append([bonds[i, 0], bonds[i, 1], bonds[j, 0]])
            angles = np.array(angles)
            
            return mols, type, charges, xyz, new_cell, bonds, angles

    return type, xyz, new_cell