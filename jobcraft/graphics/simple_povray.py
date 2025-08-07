import sys
from ase.io import read,write
import numpy as np
from matplotlib import cm

def povray_struc(atoms,save_name,show_unit_cell,rot = '0x, 0y,0z'):
    '''
    This is a function which takes list of ase.atoms objects, read charges
    from them and povray the picked one
    Charges will be normalized for each element
    '''
    if show_unit_cell is True:
        _tmp_show = 2
    else:
        _tmp_show = 0
    generic_projection_settings = {
            'rotation': rot,
            'show_unit_cell': show_unit_cell}

    povray_settings = {
        'display': False,  # Display while rendering
        'pause': False,  # Pause when done rendering (only if display)
        'canvas_width': 512,  # Width of canvas in pixels
        'camera_dist' : 50,
        'image_plane': None,
        'transparent': False,
        'camera_type': 'perspective',  # perspective, ultra_wide_angle
        'textures': ["simple" for a in atoms],
        'transmittances': [0.0 for a in atoms],
        # 'colors' : ,
        'background':'White'
    }
    rendere = write(f'{save_name}.pov',
                    atoms,
                    **generic_projection_settings,
                    povray_settings=povray_settings)


n_snap = 13
mols = read(f"{sys.argv[1]}@:",format="extxyz")

for n_snap,mol in enumerate(mols):
    povray_struc(mol,f"mol{n_snap}",True)
