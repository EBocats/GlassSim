import os
os.environ['OVITO_GUI_MODE'] = '1'
from ovito.io import import_file, export_file
import numpy as np


def ovito_to_dump(fname, muti_frames = False):
    pipline = import_file(fname, multiple_frames = muti_frames)
    name = fname.split('.')[0]
    export_file(pipline, name + '.dump', 'lammps/dump', columns = ['Particle Identifier', 'Particle Type', 'Position.X', 'Position.Y', 'Position.Z'], multiple_frames = muti_frames)
    return name + '.dump'


def add_property_todump(file_name, property_name, property_value, multiple_propertys=False, frame=0, multiple_frames=False):

    name = file_name.split('.')[0] + '_%s_modified' % frame
    pipeline = import_file(file_name, multiple_frames=multiple_frames)
    data = pipeline.compute(frame)
    
    if multiple_propertys:
        head = ['Particle Identifier', 'Particle Type',
                'Position.X', 'Position.Y', 'Position.Z']
        for i in range(len(property_name)):
            property_valuei = np.array(property_value[i])
            data.particles_.create_property(property_name[i], data=property_valuei)
            head.append(property_name[i])
        export_file(data, name + '.dump', 'lammps/dump',
                    columns=head, frame=frame)
    else:
        property_value = np.array(property_value)
        data.particles_.create_property(property_name, data=property_value)
        export_file(data, name + '.dump', 'lammps/dump', columns=['Particle Identifier', 'Particle Type',
                                                                  'Position.X', 'Position.Y', 'Position.Z', property_name], frame=frame)

    return name + '.dump'