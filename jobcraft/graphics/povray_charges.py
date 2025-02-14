from ase.io import read,write
import numpy as np
from matplotlib import cm

def normData(arr,allC):
    minVal = np.min(allC)
    maxVal = np.max(allC)
    newArr = []
    for xi in arr:
        newArr.append((xi-minVal)/(maxVal-minVal))
    return newArr


def povray_struc_charges(list_ase,charges_key,n_snap):
    '''
    This is a function which takes list of ase.atoms objects, read charges
    from them and povray the picked one
    Charges will be normalized for each element
    '''
    charges_from_everywhere = []
    symbols_from_everywhere = []
    for id_mol in range(len(list_ase)):
        charges_from_everywhere.extend(list_ase[id_mol].arrays[q_name])
        symbols_from_everywhere.extend(list_ase[id_mol].symbols)
    charges_from_everywhere_dict = {"H":[],"O":[]}

    for q,el in zip(charges_from_everywhere,symbols_from_everywhere):
        charges_from_everywhere_dict[el].append(abs(q))

    
    struc_charges = list_ase[n_snap].arrays[charges_key]
    struc_symbols = list_ase[n_snap].symbols
    qs_dict  = {"H":[],"O":[]}
    for q,el in zip(struc_charges,struc_symbols):
        qs_dict[el].append(abs(q))

    normed_POS = normData(qs_dict["H"],charges_from_everywhere_dict["H"])
    normed_NEG = normData(qs_dict["O"],charges_from_everywhere_dict["O"])
    colorsPOS = [ cm.Reds(x) for x in normed_POS ]
    colorsNEG = [ cm.Blues(x) for x in normed_NEG ]

    cols_struc = []
    countPOS = 0
    countNEG = 0
    for atom in list_ase[n_snap]:
        if atom.symbol == "O":
            cols_struc.append(colorsNEG[countNEG])
            countNEG += 1
        elif atom.symbol == "H":
            print("count",countPOS+countNEG)
            cols_struc.append(colorsPOS[countPOS])
            countPOS +=1
        generic_projection_settings = {
            'rotation': '0x, 90y,0z',
            'show_unit_cell': 0}

        povray_settings = {
            'display': False,  # Display while rendering
            'pause': False,  # Pause when done rendering (only if display)
            'canvas_width': 512,  # Width of canvas in pixels
            'camera_dist' : 50,
            'image_plane': None,
            'transparent': False,
            'camera_type': 'perspective',  # perspective, ultra_wide_angle
            'textures': ["simple" for a in list_ase[n_snap]],
            'transmittances': [0.0 for a in list_ase[n_snap]],
            'colors' : cols_struc,
            'background':'White'
        }
        rendere = write(f'{q_name}pic{n_snap}.pov',
                        list_ase[n_snap],
                        **generic_projection_settings,
                        povray_settings=povray_settings)


n_snap = 13
# n_snap = 0
neut = read("neut_with_charges.xyz@:",format="extxyz")
neg = read("neg_with_charges.xyz@:",format="extxyz")
pos = read("pos_with_charges.xyz@:",format="extxyz")
# neg = read("testik.xyz@:",format="extxyz")
print("neut",len(neut))
print("neg",len(neg))
print("pos",len(pos))
all_strucs =[]
all_strucs.extend(neut)
all_strucs.extend(neg)
all_strucs.extend(pos)
write("test.xyz",all_strucs)
q_name = "kqeq_charges"
# q_name = "dft_hirshfeld"
for n_snap in range(len(all_strucs)):
    print(f"start {n_snap}")
    povray_struc_charges(all_strucs,q_name,n_snap)
