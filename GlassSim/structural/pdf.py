import numpy as np
import os
os.environ['OVITO_GUI_MODE'] = '1'
from ovito.io import import_file
from ovito.modifiers import CoordinationAnalysisModifier

# def pdf_all(xyz, L):


def pdf_all(file, from_file = True, cutoff = 10.0, bins = 100, xyz = None, L = None):

    if from_file:
        pipeline = import_file(file)
        modify = CoordinationAnalysisModifier(cutoff = cutoff, number_of_bins = bins)
        pipeline.modifiers.append(modify)
        data = pipeline.compute()
        pdf = data.tables['coordination-rdf'].xy()
        
    else:
        dis = np.zeros((len(xyz), len(xyz)))
        dis = np.sqrt(np.sum((xyz[:, np.newaxis] - xyz[np.newaxis, :])**2, axis=-1))

        rh = 10.0
        bins = 100
        pdf = np.zeros(bins)
        for i in range(len(xyz)):
            for j in range(i+1, len(xyz)):
                if dis[i,j] < rh:
                    pdf[int(dis[i,j]/rh*bins)] += 1
        
        N = len(xyz)
        rho = N / (L[0,0]*L[1,1]*L[2,2])
        for i in range(bins):
            pdf[i] /= N
            pdf[i] /= 4*np.pi*rh**2*rho
        ri = np.linspace(0, rh, bins)
        pdf = np.column_stack((ri, pdf))
        return pdf
    return pdf
    
    

