from ase.io.lammpsrun import read_lammps_dump
from dscribe.descriptors import SOAP
from dscribe.kernels import AverageKernel
from dscribe.kernels import REMatchKernel
from sklearn.preprocessing import normalize
import numpy as np


def soap_cmp(fins, rcut = None, nmax = 8, lmax = 6, sigma = 0.1, gamma = None, alpha = 1):

    gamma  = np.where(gamma is None, 1.0/len(fins), gamma)
    if rcut is None:
        rcut = [6.0 for i in range(len(fins))]
        print('Warning: rcut is not specified, using default value 6.0 for all input files!')
    else:
        # 判断rcut是否为列表或者数组，如果是单个值，则转换为列表
        if not isinstance(rcut, (list, tuple, np.ndarray)):
            rcut = [rcut for i in range(len(fins))]
            print('Warning: only one rcut value is specified, using the same value for all input files!')
        elif len(rcut) != len(fins):
            print('Error: the number of rcut values does not match the number of input files!')
            return None
        
    print('using SOAP parameters: rcut=%s, nmax=%d, lmax=%d, sigma=%.2f' % (rcut, nmax, lmax, sigma))
    print('using Gaussian kernel with gamma = %.2f' % gamma)
    print('using REMatch kernel with alpha = %.2f' % alpha)

    soaps = []
    nom_soaps = []
    for i in range(len(fins)):
        dump = read_lammps_dump(fins[i])
        num_species = len(np.unique(dump.get_atomic_numbers()))
        species = list(range(1, num_species + 1))
        soap = SOAP(species=species, r_cut=rcut[i], n_max=nmax, l_max=lmax,
                    sigma=sigma, periodic=True, compression={"mode": "off"}, sparse=False)
        soapc = soap.create(dump)
        soaps.append(soapc)
        nom_soaps.append(normalize(soapc))

    rel = AverageKernel(metric="linear")
    rel_kernel = rel.create(soaps)
    reg = AverageKernel(metric="rbf", gamma=gamma)
    reg_kernel = reg.create(soaps)

    rrel = REMatchKernel(metric="linear", alpha=alpha, threshold=1e-6)
    rrel_kernel = rrel.create(nom_soaps)
    rreg = REMatchKernel(metric="rbf", gamma=gamma, alpha=alpha, threshold=1e-6)
    rreg_kernel = rreg.create(nom_soaps)
    print('Comparing completed!')

    return rel_kernel, reg_kernel, rrel_kernel, rreg_kernel

def direct_cmp(r1, r2):

    if len(r1) != len(r2):
        print('Error: the two vectors have different lengths!')
        return None

    dr = np.zeros(len(r1))
    for i in range(len(r1)):
        dr[i] = min(r1[i], r2[i])

    A1 = np.trapz(r1)
    A2 = np.trapz(r2)
    A = np.trapz(dr)
    AA = max(A1, A2)
    return A/AA

def eigenvalue_cmp(r1, r2):
    
    if len(r1) != len(r2):
        print('Error: the two vectors have different lengths!')
        return None

    r1 = r1/np.sqrt(np.sum(r1**2))
    r2 = r2/np.sqrt(np.sum(r2**2))
    return np.dot(r1, r2)