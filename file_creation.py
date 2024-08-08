def prep_aims_ase_file(base_file="aims_ase_template.py",final_name="runAims.py",aims_command="raven",aims_species="raven/species"):
    with open(base_file, 'r') as base_file:
        base_lines = base_file.read()
    with open(f'{final_name}', 'w') as final_file:
        final_file.write(f'aims_command="{aims_command}"\n')
        final_file.write(f'aims_species="{aims_species}"\n')
        final_file.write("\n"+base_lines)

prep_aims_ase_file()
