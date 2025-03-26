import numpy as np
import os

def bader_charge(POSCAR = 'POSCAR', ACF = 'ACF.dat', POTCAR = 'POTCAR', output_file='bader.dat'):

    files = [POSCAR, ACF, POTCAR]
    for file in files:
        if not os.path.exists(file):
            print(f'Error: {file} not found')
            return

    os.system("awk '/PAW_PBE/{getline a; print a}' %s > valence.txt" % POTCAR)
    lines = open('valence.txt').readlines()
    types = int(len(lines)/2)
    valence = np.zeros(types)
    for i in range(types):
        valence[i] = float(lines[2*i])

    # Read the POSCAR file
    with open(POSCAR, 'r') as f:
        lines = f.readlines()

    # Read the ACF file
    with open(ACF, 'r') as f:
        lines_acf = f.readlines()
        lines_acf = lines_acf[2:]

    # Get the number of atoms
    num_atoms = [int(x) for x in lines[6].split()]

    # Get the atomic symbols
    symbols = lines[5].split()

    # Get the atomic charges and volumes
    charges = []
    volumes = []
    for t in range(types):
        if t == 0:
            tlines = lines_acf[0:num_atoms[t]]
            charge = float(sum([float(x.split()[-3]) for x in tlines] - valence[t]))/num_atoms[t]
            volume = sum([float(x.split()[-1]) for x in tlines])/num_atoms[t]
        else:
            tlines = lines_acf[sum(num_atoms[0:t]):sum(num_atoms[0:t+1])]
            charge = float(sum([float(x.split()[-3]) for x in tlines] - valence[t]))/num_atoms[t]
            volume = sum([float(x.split()[-1]) for x in tlines])/num_atoms[t]
        charges.append(charge)
        volumes.append(volume)

    # Write the Bader charges to the output file
    with open(output_file, 'w') as f:
        for i in range(len(symbols)):
            f.write(f'{symbols[i]} {charges[i]:.6f} {volumes[i]:.6f}\n')
    os.system('rm valence.txt')
    
    return charges, volumes







