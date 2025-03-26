import os
os.environ['OVITO_GUI_MODE'] = '1'
from ovito.data import CutoffNeighborFinder
from ovito.modifiers import ComputePropertyModifier
import numpy as np
from ovito.io import import_file

def LocalEntropy(file, cutoff=6.0, sigma=0.2, use_local_density=True, compute_average=False, average_cutoff=5.0):
    # modified from OVITO's LocalEntropy

    assert (cutoff > 0.0)
    assert (sigma > 0.0 and sigma < cutoff)
    assert (average_cutoff > 0)

    pipline = import_file(file)
    data = pipline.compute()
    global_rho = data.particles.count / data.cell.volume
    finder1 = CutoffNeighborFinder(cutoff, data)
    if compute_average:
        finder2 = CutoffNeighborFinder(average_cutoff, data)

    local_entropy = np.empty(data.particles.count)
    nbins = int(cutoff / sigma) + 1
    r = np.linspace(0.0, cutoff, num=nbins)
    rsq = r**2

    # Precompute normalization factor of g_m(r) function:
    prefactor = rsq * (4 * np.pi * global_rho * np.sqrt(2 * np.pi * sigma**2))
    prefactor[0] = prefactor[1]  # Avoid division by zero at r=0.

    # Iterate over input particles:
    cut_num = []
    if compute_average:
        ave_cut_num = []
    for particle_index in range(data.particles.count):
        # yield particle_index / data.particles.count

        # Get distances r_ij of neighbors within the cutoff range.
        r_ij = finder1.neighbor_distances(particle_index)
        cut_num.append(len(r_ij))
        if compute_average:
            r_ij_ave = finder2.neighbor_distances(particle_index)
            ave_cut_num.append(len(r_ij_ave))

        # Compute differences (r - r_ji) for all {r} and all {r_ij} as a matrix.
        r_diff = np.expand_dims(r, 0) - np.expand_dims(r_ij, 1)

        # Compute g_m(r):
        g_m = np.sum(np.exp(-r_diff**2 / (2.0*sigma**2)), axis=0) / prefactor

        # Estimate local atomic density by counting the number of neighbors within the
        # spherical cutoff region:
        if use_local_density:
            local_volume = 4/3 * np.pi * cutoff**3
            rho = len(r_ij) / local_volume
            g_m *= global_rho / rho
        else:
            rho = global_rho

        # Compute integrand:
        integrand = np.where(
            g_m >= 1e-10, (g_m * np.log(g_m) - g_m + 1.0) * rsq, rsq)

        # Integrate from 0 to cutoff distance:
        local_entropy[particle_index] = -2.0 * \
            np.pi * rho * np.trapz(integrand, r)
        # print(local_entropy[particle_index])
    cut_num = np.array(cut_num)
    cut_num = np.mean(cut_num)
    data.particles_.create_property('Entropy', data=local_entropy)

    if compute_average:
        data.apply(ComputePropertyModifier(
            output_property='Entropy',
            operate_on='particles',
            cutoff_radius=average_cutoff,
            expressions=['Entropy / (NumNeighbors + 1)'],
            neighbor_expressions=['Entropy / (NumNeighbors + 1)']))

        ave_cut_num = np.array(ave_cut_num)
        ave_cut_num = np.mean(ave_cut_num)
        return data.particles_['Entropy'].array, cut_num, ave_cut_num
        
    return data.particles_['Entropy'].array, cut_num