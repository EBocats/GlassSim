import os
os.environ['OVITO_GUI_MODE'] = '1'
from ovito.io import import_file
from ovito.modifiers import VoronoiAnalysisModifier
import re
import numpy as np


def to_voro(filename, radii=[0.5, 0.5, 0.5]):

    if not isinstance(radii, (list, np.ndarray)):
        raise ValueError('radii must be a list or numpy array')

    pipline = import_file(filename)
    data = pipline.compute()
    cell = data.cell
    xl = cell[0, -1]
    xh = cell[0, -1] + cell[0, 0]
    yl = cell[1, -1]
    yh = cell[1, -1] + cell[1, 1]
    zl = cell[2, -1]
    zh = cell[2, -1] + cell[2, 2]

    ids = data.particles.identifiers.array
    types = data.particles.particle_types.array
    positions = data.particles.positions.array
    data = np.column_stack((ids, types, positions))
    data = np.insert(data, 5, values=0, axis=1)  

    for i in range(len(radii)):
        data[data[:, 1] == i+1, 5] = radii[i]
    data_out = np.delete(data, 1, axis=1)
    np.savetxt('nxyzr.dat', data_out, fmt='%d %.6f %.6f %.6f %.2f')
    os.system(
        'voro++ -p -r -c "%i#%A#%n#%v#%c" {} {} {} {} {} {} nxyzr.dat'.format(xl, xh, yl, yh, zl, zh))
    os.system('del nxyzr.dat')  # windows

    atoms = []
    voro_lines = open('nxyzr.dat.vol').readlines()
    for dat in voro_lines:
        dat = dat.split('#')
        atom_id = int(dat[0])
        atom_type = data[data[:, 0] == atom_id, 1]
        atom_radius = data[data[:, 0] == atom_id, 5]
        voro_index = dat[1].split()
        if len(voro_index) <= 3:
            print('Atom %d has no Voronoi index!' % atom_id)
            voro_index = '0 0 0 0'
        elif len(voro_index) == 4:
            voro_index = voro_index[3] + ' 0 0 0'
        elif len(voro_index) == 5:
            voro_index = voro_index[3] + ' ' + voro_index[4] + ' 0 0'
        elif len(voro_index) == 6:
            voro_index = voro_index[3] + ' ' + \
                voro_index[4] + ' ' + voro_index[5] + ' 0'
        else:
            voro_index = voro_index[3] + ' ' + voro_index[4] + \
                ' ' + voro_index[5] + ' ' + voro_index[6]
        neighbor = dat[2].split()
        # nbtype = np.zeros(len(neighbor), dtype=int)
        # for j in range(len(neighbor)):
        #     nbtype[j] = data[data[:, 0] == int(neighbor[j]), 1]
        # nb_type, nb_counts = np.unique(nbtype, return_counts=True)
        volume = float(dat[3])
        centroid = [float(x) for x in dat[4].split()]
        centroid = np.sqrt(centroid[0]**2 + centroid[1]**2 + centroid[2]**2)


        atom = {}
        atom['id'] = atom_id
        atom['type'] = atom_type[0]
        atom['radius'] = atom_radius[0]
        atom['voro_index'] = voro_index
        atom['neighbors'] = neighbor
        # atom['nb_type'] = nb_type
        # atom['nb_counts'] = nb_counts
        atom['volume'] = volume
        atom['centroid'] = centroid
        atoms.append(atom)
    
    atoms = sorted(atoms, key=lambda x: x['id'])
    os.system('del nxyzr.dat.vol')  # windows
    return atoms

def voro_cns(atoms, type = True):
    if type == True:
        types = [atom['type'] for atom in atoms]
        unique_types = list(set(types))
        nns_avg = []
        for i in range(len(unique_types)):
            atoms_i = [atom for atom in atoms if atom['type'] == unique_types[i]]
            natom_neb = [atom['neighbors'] for atom in atoms_i]
            natom_nn = [len(neb) for neb in natom_neb]
            nn_avg = np.mean(natom_nn)
            nn_avg = round(nn_avg, 2)
            nns_avg.append(nn_avg)
        return nns_avg
    else:
        natom_neb = [atom['neighbors'] for atom in atoms]
        natom_nn = [len(neb) for neb in natom_neb]
        nn_avg = np.mean(natom_nn)
        return nn_avg

def _row_histogram(a):
    ca = np.ascontiguousarray(a).view([('', a.dtype)] * a.shape[1])
    unique, indices, inverse = np.unique(ca, return_index=True, return_inverse=True)
    counts = np.bincount(inverse)
    sort_indices = np.argsort(counts)[::-1]
    return (a[indices[sort_indices]], counts[sort_indices])

def vor_list(filename, radii = None, use_radii = False, edge_threshold = 0.0, face_threshold = 0.0):
    
    if use_radii == True:
        if not isinstance(radii, (list, np.ndarray)):
            raise ValueError('radii must be a list or numpy array')
    radii = np.array(radii)

    def _assign_particle_radii(frame, data):
        atom_types = data.particles_.particle_types_
        for i in range(len(radii)):
            atom_types.type_by_id_(i+1).radius = radii[i]
        
    pipeline = import_file(filename)
    data = pipeline.compute()
    if use_radii == True:
        pipeline.modifiers.append(_assign_particle_radii)
    voro = VoronoiAnalysisModifier(
        compute_indices = True,
        use_radii = use_radii,
        edge_threshold = edge_threshold,
        face_threshold = face_threshold
    )
    pipeline.modifiers.append(voro)
    data = pipeline.compute()
    voro_indices = data.particles['Voronoi Index']
    voro_indices = voro_indices[:,2:6]
    unique_indices, counts = _row_histogram(voro_indices)
    
    unique_indices = ['<%d %d %d %d>' % (i[0], i[1], i[2], i[3]) for i in unique_indices]
    voros = []
    for i in range(len(unique_indices)):
        voro = {}
        voro['index'] = unique_indices[i]
        voro['count'] = counts[i]
        voro['percent'] = counts[i] / len(voro_indices) * 100
        voros.append(voro)
    
    return voros


    
    

        
    