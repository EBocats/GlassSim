# NO spin-polarized calculation!!!
# This script is used to abstract COHP and ICOHP data from LOBSTER output file.  by L.Gao
# Usage: python icohp.py [filename]

import numpy as np
import re
import sys

def export_icohp(fin, pairs=None, export_cohp = True):

    print("This script is used to abstract COHP and ICOHP data from LOBSTER output file. (NO spin-polarized!!!)")
    
    lines = open(fin, "r").readlines()
    nn = int(lines[1].split()[0]) - 1
    print("Number of atom-pairs: ", nn)
    pair_lines = lines[3:3+nn]
    pairs = []
    for line in pair_lines:
        pair = {}
        temp = line.split('->')
        temp1 = temp[0].split(':')[-1]
        temp2 = temp[1].split('(')[0]
        distence = temp[1].split('(')[1].split(')')[0]
        atom1 = re.findall(r'[A-Za-z]+', temp1)[0]
        atom2 = re.findall(r'[A-Za-z]+', temp2)[0]
        id1 = re.findall(r'\d+', temp1)[0]
        id2 = re.findall(r'\d+', temp2)[0]
        pair['pair'] = atom1 + '-' + atom2
        pair['ids'] = id1 + '-' + id2
        pair['distance'] = distence
        pairs.append(pair)

    checks = [pair['ids'] for pair in pairs]
    for check in checks:
        new_check = check.split('-')[1] + '-' + check.split('-')[0]
        if new_check in checks:
            print("Warning: ", new_check, " is in the list!")

    pair_types = [pair['pair'] for pair in pairs]
    pairs_types = list(set(pair_types))
    pairs_types.sort()
    print("Number of atom-pair types: ", len(pairs_types))
    print("Atom-pair types: ", pairs_types)

    cohp_lines = lines[3+nn:]
    ncohp = len(cohp_lines)
    cohp_data = np.loadtxt(cohp_lines)
    if cohp_data.shape[1] != 2*nn+3:
        print("Warning: The number of columns in COHP data is not equal to 2*nn+3!")
        exit()
    zero_id = np.argmin(abs(cohp_data[:, 0]))
    fave = "iCOHP_average.dat"
    fp = open(fave, "w")
    fp.write("# Atom-pair Bonding Antibonding ICOHP\n")
    for type in range(len(pairs_types)):
        n_type = [pair['pair'] for pair in pairs].count(pairs_types[type])
        print("Atom-pair type: ", pairs_types[type], " number ", n_type)
        cohp = np.zeros((ncohp, 2))
        icohp = []
        for i in range(nn):
            if pairs[i]['pair'] == pairs_types[type]:
                cohp[:, 0] = cohp_data[:, 0]
                cohp[:, 1] += cohp_data[:, 2*i+3]
                icohp.append(cohp_data[zero_id, 2*i+4])
        
        icohp = np.array(icohp)
        np.savetxt("ICOHP_" + pairs_types[type] + ".dat", -icohp, fmt='%.4f')
        cohp[:, 1] /= n_type
        if export_cohp:
            fout = "COHP_" + pairs_types[type] + ".dat"
            np.savetxt(fout, cohp, fmt='%.4f')

        cohp = cohp[cohp[:, 0] <= 0, :]
        bonding = cohp[:, 1] < 0
        antibonding = cohp[:, 1] > 0
        integral_bonding = abs(np.trapz(cohp[bonding, 1], cohp[bonding, 0]))
        integral_antibonding = abs(
            np.trapz(cohp[antibonding, 1], cohp[antibonding, 0]))
        print("Average: Bonding = %.4f, Antibonding = %.4f, ICOHP = %.4f" % (
            integral_bonding, integral_antibonding, integral_bonding-integral_antibonding))
        fp.write("%s\t%.4f\t%.4f\t%.4f\n" % (
            pairs_types[type], integral_bonding, integral_antibonding, integral_bonding-integral_antibonding))
    fp.close()

    # if pairs is not None:
    #     pairs = pairs_types
        

