import os
os.environ['OVITO_GUI_MODE'] = '1'
from ovito.io import import_file
from ovito.modifiers import CalculateDisplacementsModifier
import numpy as np
from mgs.tools.pbc import apply_PBC


def get_dislpacement(file_path, frame = 1, frame_offset = None, reference_frame = None, particle_type=None):

    if not os.path.exists(file_path):
        raise FileNotFoundError("The file does not exist")
    try:
        pipeline = import_file(file_path, multiple_frames=True)
    except:
        raise ValueError("The file is not a valid trajectory file")

    if frame_offset is None and reference_frame is None:
        modifier = CalculateDisplacementsModifier()
    elif frame_offset is None and reference_frame is not None:
        modifier = CalculateDisplacementsModifier(reference_frame=reference_frame)
    elif frame_offset is not None and reference_frame is None:
        modifier = CalculateDisplacementsModifier(
            use_frame_offset=True, frame_offset=frame_offset)
    
    pipeline.modifiers.append(modifier)
    data = pipeline.compute(frame)

    if particle_type is None:
        displacements = data.particles['Displacement Magnitude']
    else:
        types = data.particles['Particle Type']
        if particle_type not in types:
            raise ValueError("The particle type does not exist")
        displacements = data.particles['Displacement Magnitude'][data.particles['Particle Type'] == particle_type]

    return displacements.array


def pu(displacement=None, bin=200, displacement_max=10.0):

    try:
        if not isinstance(displacement, np.ndarray):
            displacement = np.array(displacement)
    except:
        raise ValueError(
            "The input is not a valid displacement data, must be a 1-D numpy array or a list")

    displacement = displacement[displacement <= displacement_max]
    u_max = displacement_max
    detu = u_max/bin
    up = 0
    if u_max >= detu:
        hist, edges = np.histogram(displacement, bin, density=True)
        up = edges[np.argmax(hist)]
        sts = np.column_stack((edges[1:], hist))
    return up, sts


def pus(displacements=None, dump_file_path=None, frame_offset=0, bin=200, particle_type=None):

    if dump_file_path is not None:
        displacements = get_dislpacement(
            dump_file_path, frame_offset, particle_type)
        ups = np.zeros(len(displacements))
        stss = []
        for i in range(len(displacements)):
            up, sts = pu(displacements[i], bin)
            ups[i] = up
            stss.append(sts)

    elif displacements is not None:
        ups = np.zeros(len(displacements))
        stss = []
        for displacement in displacements:
            up, sts = pu(displacement, bin)
            ups[i] = up
            stss.append(sts)

    ups = np.array(ups).mean()
    out = stss[0]
    for i in range(1, len(stss)):
        out = np.row_stack((out, stss[i]))
    out = out[out[:, 0].argsort()]
    return ups, out


def pu_auto(displacement = None, bin = 200):
    
    try:
        if not isinstance(displacement, np.ndarray):
            displacement = np.array(displacement)
    except:
        raise ValueError("The input is not a valid displacement data, must be a 1-D numpy array or a list")

    u_max = np.max(displacement)
    detu = u_max/bin 
    up=0
    if u_max>=detu:
        hist,edges = np.histogram(displacement, bin, density = True)
        up = edges[np.argmax(hist)]
        sts = np.column_stack((edges[1:],hist))    
    return up, sts


def displacement_xyz(A, B, Ls, PBC = True):

    if not isinstance(A, np.ndarray):
        A = np.array(A)
    if not isinstance(B, np.ndarray):
        B = np.array(B)
    
    dis = A - B
    if PBC:
        dis = apply_PBC(dis, Ls)

    return dis


        
            

            
    
    
    

    

            
