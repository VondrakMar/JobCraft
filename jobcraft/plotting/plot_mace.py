import numpy as np
from ase.io import read,write
import matplotlib.pyplot as plt


# If you are not me and reading this, I know these are wrong values, don't judge me
atomic_energies = {"H": -13.59803017, "O": -2043.567039796}

def preprocess_energy_for_per_atom(atoms_list,dft_energy,ml_energy,plot_per_atom,atomic_energies):
    dft_energy_list = []
    ml_energy_list = []
    energy_label = "something went wrong, you would wish this is Rust now eh?!"
    if atomic_energies is None and not plot_per_atom:
        energy_label = "energy eV"
        for atoms in atoms_list:
            dft_energy_list.append(atoms.info[dft_energy])
            ml_energy_list.append(atoms.info[ml_energy])
    elif atomic_energies is not None and not plot_per_atom:
        energy_label = "energy-atomization eV"
        for atoms in atoms_list:
            cur_atE = sum([atomic_energies[sym] for sym in atoms.symbols])
            dft_energy_list.append(atoms.info[dft_energy]-cur_atE)
            ml_energy_list.append(atoms.info[ml_energy]-cur_atE)
    elif atomic_energies is not None and plot_per_atom:
        energy_label = "energy-atomization eV/atom"
        for atoms in atoms_list:
            cur_atE = sum([atomic_energies[sym] for sym in atoms.symbols])
            dftE = (atoms.info[dft_energy] - cur_atE)/len(atoms)
            mlE = (atoms.info[ml_energy] - cur_atE)/len(atoms)
            dft_energy_list.append(dftE)
            ml_energy_list.append(mlE)
    elif atomic_energies is None and plot_per_atom:
        energy_label = "energy eV/atom"
        for atoms in atoms_list:
            dft_energy_list.append(atoms.info[dft_energy]/len(atoms))
            ml_energy_list.append(atoms.info[ml_energy]/len(atoms))
    return np.array(ml_energy_list),np.array(dft_energy_list),energy_label
            
def plot_correlation(atoms_list,
                     ml_energy = "MACE_energy",
                     dft_energy = "dft_energy",
                     dft_info_list=["dft_energy"],
                     dft_arrays_list=["dft_forces"],
                     ml_info_list=["MACE_energy"],
                     ml_arrays_list=["MACE_forces"],
                     dft_pre="dft_",
                     ml_pre="ml_",
                     info_names=["energy"],
                     arrays_names=["forces"],
                     atom_Es = None,
                     plot_per_atom=False):#["ml_energy","dft_energy"]):
    # if plot_per_atom and atom_Es is None:
        # print("If you want to plot atomization energies, you have to provide atomic energies, I will plot whole energy")
        # plot_per_atom = False
    assert all([ml_name in atoms_list[0].arrays for ml_name in ml_arrays_list]), "something wrong with ml arrays labels"
    assert all([dft_name in atoms_list[0].arrays for dft_name in dft_arrays_list]), "something wrong with dft arrays labels"
    assert all([ml_name in atoms_list[0].info for ml_name in ml_info_list]), "something wrong with ml info labels"
    assert all([dft_name in atoms_list[0].info for dft_name in dft_info_list]), "something wrong with dft info labels"
    assert len(dft_info_list) == len(ml_info_list)
    assert len(dft_arrays_list) == len(ml_arrays_list)
    assert len(info_names) == len(dft_info_list)
    assert len(arrays_names) == len(ml_arrays_list)
    if ml_energy in dft_info_list:
        print("removing ml energy tag from the list you provided, so it is not plotted twice")
        ml_info_list.remove(ml_energy)
    if dft_energy in dft_info_list:
        print("removing dft energy tag from the list you provided, so it is not plotted twice")
        dft_info_list.remove(dft_energy)
    
    ml_energy,dft_energy,energy_label = preprocess_energy_for_per_atom(atoms_list,dft_energy,ml_energy,plot_per_atom,atom_Es)
    dft_values = {}
    ml_values = {}
    config_types = []
    atoms_lens = []
    for atoms in atoms_list:
        config_types.append(atoms.info.get("config_type", "unknown"))
        atoms_lens.append(len(atoms))
    config_types = np.array(config_types)
    unique_types = np.unique(config_types)


    for info_dft,info_ml in zip(dft_info_list,ml_info_list):
        dft_values[info_dft] = []
        ml_values[info_ml] = []
        for atoms in atoms_list:
            dft_values[info_dft].append(atoms.info[info_dft])
            ml_values[info_ml].append(atoms.info[info_ml])
        dft_values[info_dft] = np.array(dft_values[info_dft])
        ml_values[info_ml] = np.array(ml_values[info_ml])

    for arrays_dft,arrays_ml in zip(dft_arrays_list,ml_arrays_list):
        dft_values[arrays_dft] = []
        ml_values[arrays_ml] = []
        for atoms in atoms_list:
            dft_values[arrays_dft].extend(atoms.arrays[arrays_dft].flatten())
            ml_values[arrays_ml].extend(atoms.arrays[arrays_ml].flatten())
        dft_values[arrays_dft] = np.array(dft_values[arrays_dft])
        ml_values[arrays_ml] = np.array(ml_values[arrays_ml])


    
    plt.figure(figsize=(7,7))
    for ctype in unique_types:
        mask = config_types == ctype
        plt.scatter(dft_energy[mask],
                    ml_energy[mask],
                    alpha=0.7,
                    edgecolor="k",
                    label=ctype)
    
    min_e = min(dft_energy.min(), ml_energy.min())
    max_e = max(dft_energy.max(), ml_energy.max())
    plt.plot([min_e, max_e], [min_e, max_e], "r--", label="y = x")
    plt.xlabel(f"{dft_pre}{energy_label}")
    plt.ylabel(f"{ml_pre}{energy_label}")
    plt.legend()
    plt.tight_layout()
    plt.savefig("ftE.png", dpi=300)
    plt.show()



    
    for info_dft,info_ml in zip(dft_info_list,ml_info_list):
        for ctype in unique_types:
            mask = config_types == ctype
            plt.scatter(dft_values[info_dft][mask],
                        ml_values[info_ml][mask],
                        alpha=0.7,
                        edgecolor="k",
                        label=ctype)

        min_e = min(dft_values[info_dft].min(), ml_values[info_ml].min())
        max_e = max(dft_values[info_dft].max(), ml_values[info_ml].max())
        plt.plot([min_e, max_e], [min_e, max_e], "r--", label="y = x")


        plt.xlabel("DFT Energy")
        plt.ylabel("ML Energy")
        plt.legend()
        plt.tight_layout()
        plt.savefig("ftE.png", dpi=300)
        plt.show()

    expanded_config_types = np.repeat(config_types, atoms_lens)
    expanded_config_types_forces = np.repeat(config_types, np.array(atoms_lens) * 3)
    # Now you can make masks directly
    unique_types = np.unique(expanded_config_types)

    for arrays_dft, arrays_ml,array_name in zip(dft_arrays_list, ml_arrays_list,arrays_names):
        for ctype in unique_types:
            if "forces" in arrays_dft:
                mask = expanded_config_types_forces == ctype
            else:
                mask = expanded_config_types == ctype
            plt.scatter(dft_values[arrays_dft][mask],
                        ml_values[arrays_ml][mask],
                        alpha=0.7,
                        edgecolor="k",
                        label=ctype)

                
        min_e = min(dft_values[arrays_dft].min(), ml_values[arrays_ml].min())
        max_e = max(dft_values[arrays_dft].max(), ml_values[arrays_ml].max())
        plt.plot([min_e, max_e], [min_e, max_e], "r--", label="y = x")

        plt.xlabel(f"{dft_pre}{array_name}")
        plt.ylabel(f"{ml_pre}{array_name}")
        # plt.title(f"Correlation between DFT and ML energies (R² = {r2:.3f})")
        plt.legend()
        plt.tight_layout()
        plt.savefig("saved_pic.png", dpi=300)
        plt.show()


# mols = read("test_res_positive_bulk_scratch.xyz@:",format="extxyz")
mols = read("ft.xyz@:20",format="extxyz")
# print(sum(mols[0].arrays["MACE_charges"]))
print([k for k in mols[10].arrays.keys()])
print([k for k in mols[10].info.keys()])
plot_correlation(mols,
                 ml_energy = "MACE_energy",
                 dft_energy = "dft_energy",
                 atom_Es = atomic_energies,
                 plot_per_atom = True,
                 dft_info_list=["dft_energy"],
                 ml_info_list=["MACE_energy"],
                 ml_arrays_list=["MACE_forces","MACE_charges"],
                 dft_arrays_list=["dft_forces","dft_hirshfeld_q"],
                 info_names=["energy (eV/atom)"],
                 arrays_names=["forces (eV/A)", "charges (e)"],
                 )
