import numpy as np


# This was originally written only for ethanol - benzene mixed system, if you reading this I was just lazy to ever do it for a general system
N_at_ethanol = 9
benzene_start = 1305+9
N_at_benzene = 12

def get_separation(mol,supercell_size=(3, 3, 3)):
    mults = supercell_size[0]*supercell_size[1]*supercell_size[2]
    nAts = int(len(mol)/mults)
    ethanol_nAts = benzene_start
    benzene_nAts = nAts - ethanol_nAts
    ethanol_part = []
    benzene_part = []
    for mult in range(mults):
        ethanol_part.extend(range(mult*(nAts),mult*(nAts)+ethanol_nAts))
        benzene_part.extend(range(mult*(nAts)+ethanol_nAts,mult*(nAts)+nAts))
    ethanol = mol[ethanol_part]
    benzene = mol[benzene_part]
    return ethanol,benzene
    

def get_molecule_list(mol,supercell_size=[3,3,3]):
    mults = supercell_size[0]*supercell_size[1]*supercell_size[2]
    nAts = int(len(mol)/mults)
    moldy_ethen, benzene = get_separation(mol,supercell_size)
    ethanol_nAts = int(benzene_start*mults)
    benzene_nAts = int((nAts - ethanol_nAts)*mults)
    ethanol_list = []
    benzene_list = []
    for id_ethanol in range(0,ethanol_nAts,9):
        ethanol_list.append(moldy_ethen[id_ethanol:id_ethanol+9])
    for id_benzene in range(0,benzene_nAts,12):
        benzene_list.append(benzene[id_benzene:id_benzene+12])
        # ethanol_list[-1].cell = [False,False,False]
        # ethanol_list[-1].center()
    return benzene_list,ethanol_list

def extract_clusters(mol,r_cut = 8.0):
    print("extracting clusters")
    assert mol.cell, "Works only for PBC systems"
    # cluster_list = []
    center = [mol.cell[0][0]/2,mol.cell[1][1]/2,mol.cell[2][2]/2]
    benzene_list,ethanol_list = get_molecule_list(mol,[1,1,1])
    
    for id_benzene,benzene in enumerate(benzene_list):
        positions = benzene.get_positions()
        distances = np.linalg.norm(positions - center,axis=1)
        is_in_question_mark = distances <= r_cut
        if True in is_in_question_mark:
            try:
                final_cluster += benzene
            except:
                final_cluster = benzene
    for id_ethanol,ethanol in enumerate(ethanol_list):
        positions = ethanol.get_positions()
        distances = np.linalg.norm(positions - center,axis=1)
        is_in_question_mark = distances <= r_cut
        if True in is_in_question_mark:
            try:
                final_cluster += ethanol
            except:
                final_cluster = ethanol
    return final_cluster
