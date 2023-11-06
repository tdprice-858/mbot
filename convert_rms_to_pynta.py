import yaml
import sys
sys.path.insert(0, "/Users/tdprice/Desktop/Desktop/RMG-Py/")
from rmgpy.molecule import Molecule
from rmgpy.molecule.draw import MoleculeDrawer
#import rmgpy
from rmgpy.molecule.adjlist import to_adjacency_list
from ase import Atoms
from ase.data.pubchem import pubchem_atoms_search
from ase.visualize import view
import numpy as np
import os

print('hello')

file = '/Users/tdprice/Desktop/Desktop/RMG-Py/examples/rmg/catalysis/\
MeOH_ox_Ag_no_add_species/rms/chem91.rms'
chemkin_species = '/Users/tdprice/Desktop/Desktop/RMG-Py/examples/\
rmg/catalysis/MeOH_ox_Ag_no_add_species/chemkin/species_dictionary_updated.txt'

with open(file, 'r') as file:
    rms_file = yaml.safe_load(file)

with open(chemkin_species, 'r') as chem_file:
    #print(chem_file.readline())
    i = 0
    species_dictionary = {}
    for line in chem_file.readlines():
        if i==0:
            #print(line)
            if '(' in line:
                #print('('.join(line.split('(')[:-1]))
                species = '('.join(line.split('(')[:-1])
                species_dictionary[species] = ''
            else:
                continue
        if 'multiplicity' in line:
            continue
        if i>=1 and 'multiplicity' not in line and '\n'!=line:
            species_dictionary[species] += line
            print(line)
        #print(line)
        i+=1
        if line=='\n':
            print(line)
            i=0
    #chemkin_file = yaml.safe_load(chem_file)
print(species_dictionary)
#print(chemkin_file.split('\n'))
#for line in chemkin_file.readlines():
#    print(line)


#print(type(rms_file['Phases']))

#for key, values in rms_file.items():
#    print(key)
    #print(values)
    #if key=='Phases':
        #print(key)
     #   for k, values in rms_file[f"{key}"]:
            #print(k)
            #print(rms_file[f"{key}"])'''
#print(rms_file['Phases'])
#print(len(rms_file['Phases']))
#print(rms_file['Phases'][0])
#print(rms_file['Phases'][1])

#print(f"length of reactions {len(rms_file['Reactions'])}")
#print(f"type of reactions == {type(rms_file['Reactions'])}")
#for count, item in enumerate(rms_file['Reactions']):
#    print(count)
#    print(item)

#for key, item in rms_file['Phases'][0].items():
#    print(f"key... {key}")
#    print(f"item .... {item}")

##for key, item in rms_file['Phases'][1].items():
 #   print(f"key... {key}")
    #print(f"item .... {item}")
    #print(f"item length.... {len(item)}")

#Surface species
#print(rms_file['Phases'][1]['name'])
#print(rms_file['Phases'][1]['Species'])
species_dict = {}
for species in rms_file['Phases'][1]['Species']:
    #print(species['smiles'])
    if species['smiles']=='[Pt]' and species['name']=='vacantX':
        species_dict['X'] = species['adjlist']
    else:
        species_dict[f"{species['smiles'].strip()}"] = species['adjlist']


    #print(f"name {species['name']}")
    #print(f"smiles {species['smiles']}")
    #print(species['adjlist'])
for species in rms_file['Phases'][0]['Species']:
    #try:
    #    atoms = pubchem_atoms_search(smiles=f"{species['smiles'].strip()}")
    #except:
    #    continue
    #view(atoms)
    adjlist = Molecule().from_smiles(species['smiles'].strip()).to_adjacency_list()
    #adj_list = to_adjacency_list(species['smiles'].strip(), multiplicity=species['radicalelectrons'])
    #print(adjlist)
    #print('hhhihihihihihihiuh')
    species_dict[f"{species['smiles'].strip()}"] = adjlist

#print(species_dict)
# For many of the mondentate species, we could really have a bidentate
rxn_yaml = {'rxn_1' : 'CH3OH + vacantX + OX <-> CO[Pt] + O[Pt]',
'rxn_2' : 'CH3OH + vacantX + O[Pt] <-> CO[Pt] + O + vacantX',
'rxn_3' : 'CO[Pt] + OX <-> CH2O + vacantX + O[Pt]',
'rxn_4' : 'CO[Pt] + O[Pt] <-> CH2O + vacantX + O + vacantX',
'rxn_5' : 'CH3OH + vacantX + OX <-> OC[Pt] + O[Pt]',
'rxn_6' : 'CH3OH + vacantX + O[Pt] <-> OC[Pt] + O + vacantX',
'rxn_7' : 'OC[Pt] + OX <-> CH2O + vacantX + O[Pt]',
'rxn_8' : 'OC[Pt] + O[Pt] <-> CH2O + vacantX + O + vacantX',
'rxn_9' : 'OC[Pt] + OX <-> OCX + O[Pt]',
'rxn_10' : 'OC[Pt] + O[Pt] <-> OCX + O + vacantX',
'rxn_11' : 'OCX + OX <-> CHOX + O[Pt]',
'rxn_12' : 'OCX + O[Pt] <-> CHOX + O + vacantX',
'rxn_13' : 'OCX + OX <-> HOCX + O[Pt]',
'rxn_14' : 'OCX + O[Pt] <-> HOCX + O + vacantX',
'rxn_15' : 'HOCX + OX <-> COX + O[Pt]',
'rxn_16' : 'HOCX + O[Pt] <-> COX + O + vacantX',
'rxn_17' : 'CH2O + vacantX + OX <-> CHOX + O[Pt]',
'rxn_18' : 'CH2O + vacantX + O[Pt] <-> CHOX + O + vacantX',
'rxn_19' : 'CHOX + OX <-> COX + O[Pt]',
'rxn_20' : 'CHOX + O[Pt] <-> COX + O + vacantX',
'rxn_21' : 'CH3OH + vacantX + vacantX <-> CO[Pt] + [Pt]',
'rxn_22' : 'CH3OH + vacantX + vacantX <-> OC[Pt] + [Pt]',
'rxn_23' : 'CH3OH + vacantX + vacantX <-> O[Pt] + C[Pt] ',
'rxn_24' : 'CO[Pt] + vacantX <-> CH2O + vacantX + [Pt]',
'rxn_25' : 'CO[Pt] + vacantX <-> C[Pt] + OX',
'rxn_26' : 'Formaldehyde desorption',
'rxn_27' : 'CH2O + vacantX + vacantX <-> CHOX + [Pt]',
'rxn_28' : 'CH2O + vacantX + vacantX <-> CH2X + OX',
'rxn_29' : 'CHOX + vacantX <-> COX + [Pt]',
'rxn_30' : 'CHOX + vacantX <-> CHX + OX',
'rxn_31' : 'OC[Pt] + vacantX <-> CH2O + vacantX + [Pt]',
'rxn_32' : 'OC[Pt] + vacantX <-> OCX + [Pt]',
'rxn_33' : 'OC[Pt] + vacantX <-> CH2X + O[Pt]',
'rxn_34' : 'OCX + vacantX <-> HOCX + [Pt]',
'rxn_35' : 'OCX + vacantX <-> CHOX + [Pt]',
'rxn_36' : 'OCX + vacantX <-> CHX + O[Pt]',
'rxn_37' : 'HOCX + vacantX <-> COX + [Pt]',
'rxn_38' : 'HOCX + vacantX <-> CX + O[Pt]',
'rxn_39' : 'C[Pt] + vacantX <-> CH2X + [Pt]',
'rxn_40' : 'CH2X + vacantX <-> CHX + [Pt]',
'rxn_41' : 'CHX + vacantX <-> CX + [Pt]',
'rxn_42' : 'C[Pt] + OX <-> CH2X + O[Pt]',
'rxn_43' : 'CH2X + OX <-> CHX + O[Pt]',
'rxn_44' : 'CHX + OX <-> CX + O[Pt]',
'rxn_45' : 'C[Pt] + O[Pt] <-> CH2X + O + vacantX',
'rxn_46' : 'CH2X + O[Pt] <-> CHX + O + vacantX',
'rxn_47' : 'CHX + O[Pt] <-> CX + O + vacantX',
'rxn_48' : 'CO + vacantX <-> COX',
'rxn_49' : 'CO2 + vacantX <-> COX + OX',
'rxn_50' : 'COX + vacantX <-> CX + OX',
'rxn_51' : 'O + vacantX + vacantX <-> O[Pt] + [Pt]',
'rxn_52' : 'O[Pt] + vacantX <-> OX + [Pt]',
'rxn_53' : 'O2 + vacantX + vacantX <-> OX + OX',
'rxn_54' : '[H][H] + vacantX + vacantX <-> [Pt] + [Pt]',
'rxn_55' : 'O[Pt] + O[Pt] <-> O + vacantX + OX'}

'''#What is CH2OH in smiles?
'''
#rxns = ['rxn_23','rxn_26','rxn_33','rxn_36']
rxns = ['rxn_1']
for ind, rxn in enumerate(rxns):
    reactants = rxn_yaml[rxn].split('<->')[0]
    products = rxn_yaml[rxn].split('<->')[1]
    react = [reactant.strip() for reactant in reactants.split('+')]
    prods = [product.strip() for product in products.split('+')]
    react_adj_list = [species_dictionary[f"{r}"] for r in react]
    prod_adj_list = [species_dictionary[f"{r}"] for r in prods]
    new_react_adj =''
    for r in react_adj_list:
        new_react_adj += r
    new_prod_adj = ''
    for r in prod_adj_list:
        new_prod_adj += r
    #print(new_react_adj)

    new_react_adj = new_react_adj.split('\n')[:-1]

    new_prod_adj = new_prod_adj.split('\n')[:-1]
    # Renumbering individual reactant adj list for total adj list
    #print(test)

    def reorder_adj_list(test):
        for count,t in enumerate(test):
            if not str(count+1)==t.split()[0]:
                line = test[count]
                split_line = line.split()
                old_index = split_line[0]
                split_line[0] = str(count+1)
                new_line = ' '.join(split_line)
                test[count] = new_line
                if '{' in new_line:
                    connections = new_line.split('c0 ')[1:]
                    conns = connections[0].split()
                    # This for loop is for replacing the old connections so the final adjacency list is correct
                    for new_count,c in enumerate(connections[0].split()):
                        new_connection = int(c.split('{')[1].split(',')[0]) - int(old_index) + int(new_line.split()[0])
                        new_conn = '{' + str(new_connection) + ',' + c.split('{')[1].split(',')[1]
                        conns[new_count] = new_conn
                        newer_line = new_line.split('c0 ')[0] + 'c0 ' + ' '.join(conns)
                        test[count] = newer_line
        return test
    old_r_adj = '\n'.join(new_react_adj)
    print(f"Old  reactant adj {old_r_adj}")
    ordered_react_adj = reorder_adj_list(new_react_adj)
    ord_r_adj = '\n'.join(ordered_react_adj)
    print(f"Ordered reactant adj {ord_r_adj}")
    ordered_prod_adj = reorder_adj_list(new_prod_adj)
    old_p_adj = '\n'.join(new_prod_adj)
    print(f"Old product adj {old_p_adj}")
    ord_p_adj = '\n'.join(ordered_prod_adj)
    print(f"Ordered product adj {ord_p_adj}")

    def get_atom_connect_dict(adj):
        '''
        Function to extract the atoms that each atom is bonded to.
        Returns a dictionary with the format...
        {'1 O': ['{C,S}', '{H,S}'], ...}
        {'index element': ['{bonded atom, type of bond}, ...]
        '''
        ordered_react_connections = {}
        for line in adj:
            connected_atoms = []
            for connection in line.split()[5:]:
                try:
                    connect_atom = (adj[int(connection.split(',')[0][1]) - 1]).split()[1]
                    connected_atoms.append(connect_atom)
                except:
                    continue
            ordered_react_connections[f"{' '.join(line.split()[:2])}"] = connected_atoms
        return ordered_react_connections

    react_atom_connects = get_atom_connect_dict(ordered_react_adj)
    prod_atom_connects = get_atom_connect_dict(ordered_prod_adj)

    mismatched_lines = []
    for react in react_atom_connects:
        try:
            if prod_atom_connects[react]:
                match = True
        except:
            mismatched_lines.append(react.split()[0])
            continue

    prod_line_map = []
    for line in mismatched_lines:
        index_of_react_adj = int(line)-1
        for new_line in mismatched_lines:
            prod_index_of_adj = int(new_line) - 1
            react_element = ordered_react_adj[index_of_react_adj].split()[1]
            prod_element = ordered_prod_adj[prod_index_of_adj].split()[1]
            if react_element==prod_element:
                prod_line_map.append([prod_index_of_adj, index_of_react_adj])
    mismatched_prod_lines = {int(line) - 1:ordered_prod_adj[int(line) - 1] for line in mismatched_lines}
    for map in prod_line_map:
        old_line = mismatched_prod_lines[map[0]].split()
        old_line[0] = str(map[1] + 1)
        new_line = ' '.join(old_line)
        ordered_prod_adj[map[1]] = new_line

    for count, line in enumerate(ordered_prod_adj):
        try:
            connections = line.split()[5:]
            for new_count,con in enumerate(connections):
                changed = False
                for map in prod_line_map:
                    if str(map[0] + 1) in con and not changed:
                        con = con.replace(str(map[0] + 1), str(map[1] + 1))
                        connections[new_count] = con
                        changed = True
            new_line = ' '.join(line.split()[:5] + connections)
            ordered_prod_adj[count] = new_line

        except:
            continue

    bond_break_pos = []

    for react in react_atom_connects:
        bond_break_pos.append(react.split()[0])
        sorted_react = sorted(react_atom_connects[f"{react}"])
        for prod in prod_atom_connects:
            sorted_prod = sorted(prod_atom_connects[f"{prod}"])
            if sorted_react == sorted_prod and prod==react:
                bond_break_pos.pop()

    for count, pos in enumerate(bond_break_pos):
        old_react_line = ordered_react_adj[int(pos) - 1].split()
        old_prod_line = ordered_prod_adj[int(pos) - 1].split()

        new_react_line = old_react_line[0] + f" *{count + 1} " + ' '.join(old_react_line[1:])
        new_prod_line = old_prod_line[0] + f" *{count + 1} " + ' '.join(old_prod_line[1:])

        ordered_react_adj[int(pos) - 1] = new_react_line
        ordered_prod_adj[int(pos) - 1] = new_prod_line

    # Need to get multiplicity
    react_unpaired_electrons = 0
    prod_unpaired_electrons = 0
    for index in range(len(ordered_react_adj)):
        react_unpaired_electrons += int(ordered_react_adj[index].split('u')[1].split()[0])
        prod_unpaired_electrons += int(ordered_prod_adj[index].split('u')[1].split()[0])

    reactant_mult = 1 + react_unpaired_electrons
    product_mult = 1 + prod_unpaired_electrons
    react_adj = '\n\n    '.join(ordered_react_adj)
    prod_adj = '\n\n    '.join(ordered_prod_adj)
    '''os.chdir('/Users/tdprice/Desktop/ML_training_data/')
    os.mkdir(f"{rxn}")
    os.chdir(f"/Users/tdprice/Desktop/ML_training_data/{rxn}")
    with open(f"rxn.yaml", 'w') as l:
        l.write(f'- index: 0\n')
        l.write(f"  reactant: 'multiplicity {reactant_mult}\n\n")
        l.write(f"    {react_adj} \n\n    '\n")
        l.write(f"  product: 'multiplicity {product_mult}\n\n")
        l.write(f"    {prod_adj} \n\n    '\n")
        l.write(f"  reaction: {rxn_yaml[rxn]}\n")
        l.write(f"  reaction_family: Eley-Rideal")
    os.chdir('/Users/tdprice/Desktop/')'''
    print(f'- index: {ind}')
    print(f"  reactant: 'multiplicity {reactant_mult}\n")
    react_adj = '\n\n    '.join(ordered_react_adj)
    print(f"    {react_adj} \n\n    '")
    print(f"  product: 'multiplicity {product_mult}\n")
    prod_adj = '\n\n    '.join(ordered_prod_adj)
    print(f"    {prod_adj} \n\n    '")
    print(f"  reaction: {rxn_yaml[rxn]}")
    print(f"  reaction_family: Eley-Rideal")

    react_molecule = Molecule().from_adjacency_list(react_adj)
    MoleculeDrawer().draw(molecule=react_molecule, file_format='png', target=f'react.png')
    prod_molecule = Molecule().from_adjacency_list(prod_adj)
    MoleculeDrawer().draw(molecule=prod_molecule, file_format='png', target=f'prod.png')
    from PIL import Image
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    #mpl.rcParams['figure.dpi'] = 800
    import matplotlib as mpl
    f, axs = plt.subplots(1, 3, figsize=(3.5, 1.5))
    f.suptitle(f"{rxn}:{rxn_yaml[rxn]}", fontsize=8)
    react_im = Image.open(f"react.png")
    #width, height = react_im.size
    #react_im.thumbnail((width * 0.5,height * 0.5))
    axs[0].imshow(react_im)
    #plt.imshow(react_im)
    #plt.subplot(1,2,2)
    arrow_im = Image.open("equil_arrows.png")
    width, height = arrow_im.size
    arrow_im.resize((int(width * 0.25),int(height * 0.25)))
    axs[1].imshow(arrow_im)
    prod_im = Image.open(f"prod.png")
    #width, height = prod_im.size
    #prod_im.thumbnail((width * 0.5,height * 0.5))
    axs[2].imshow(prod_im)
    for count, a in enumerate(axs.flat):
            a.spines['top'].set_visible(False)
            a.spines['right'].set_visible(False)
            a.spines['bottom'].set_visible(False)
            a.spines['left'].set_visible(False)
            a.set_xticks([])
            a.set_yticks([])
    plt.savefig(f"{rxn}.png")
    plt.show()


        #elif prod == react:
            #print(react)
        #    print('Have same atom connections')
        #    print(prod)
        #    print(react)
        #elif react==prod:
        #    print('reactant and product are the same')
        #    print(prod)
        #    print(react)

    #for prod in prod_atom_connects:


#for item in test:
#    print(item)

#for count, line in enumerate(react_adj):
    #print(line)
#print(products)
#for species in rxn_1.split('<->'):
#    print('Hi')
#    print(species)
    #print(species.strip('<->'))
#'\n'.join(reactants
'''for species in reactants.split():
    #print(species.strip())
    if species not in ['+','<', '-', '>', '<->']:
        print(species.strip())
        #print(species.strip())
        print(species_dict[f"{species.strip()}"])
        #if species.strip() == 'CO[Pt]':
        #    print(species.strip())
        #    print(species_dict[f"{species.strip()}"])
        #    print(len(species_dict[f"{species.strip()}"]))
        #    print(type(species_dict[f"{species.strip()}"]))
            #for character in species_dict[f"{species.strip()}"]:
            #    print(character)
            #    if character== '\n':
            #        print('I am a new paragraph')
'''

'''string = '1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}'
print(string)
connections = string.split('c0 ')[1:]
print(connections[0].split())
conns = connections[0].split()
for count, c in enumerate(connections[0].split()):
    print(c)
    print(c.split('{'))
    print(c.split('{')[1].split(','))
    print(c.split('{')[1].split(',')[0])
    new_conn = '{' + str(int(c.split('{')[1].split(',')[0]) + 1) +','+ c.split('{')[1].split(',')[1]
    conns[count] = new_conn

print(conns)
print(string.split('c0 ')[0] + 'c0 ' + ' '.join(conns))
'''