import os
os.environ['OVITO_GUI_MODE'] = '1'
from ovito.io import import_file








def read_stru(fn):
    '''
    Old version! Not recommended.
    '''
    pipline = import_file(fn)
    data = pipline.compute()
    return data.particles['Particle Type'].array, data.particles['Position'].array, data.cell.matrix
