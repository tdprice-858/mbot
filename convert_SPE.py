import os, sys
import json
import subprocess



rootdir = '/global/cfs/cdirs/m3548/tdprice/10_ML_training_loop/iter1'
number_of_spes = []
paths_of_completed_rxns = []
unfinished = []
os.chdir(rootdir)
path_to_training_data = '/global/cfs/cdirs/m3548/tdprice/10_ML_training_loop/iter1/Ag_data_for_trevor'
for dir in os.listdir(rootdir):
    if dir.split('_')[0] == 'rxn':
        paths_of_completed_rxns = []
        # print(dir)
        if os.path.exists(os.path.join(rootdir, dir, 'TS0/structures_for_ML.json')):
            with open(os.path.join(rootdir, dir, 'TS0/structures_for_ML.json')) as f:
                structs = json.load(f)
                for struct in structs:
                    path = '/'.join(struct.split('/')[:-1])

                    number_of_spes.append(path)
                    if os.path.exists(os.path.join(path, 'sella_fs/spe_RPBE_500_441/vasprun.xml')):
                        paths_of_completed_rxns.append(os.path.join(path, 'sella_fs/spe_RPBE_500_441/vasprun.xml'))
                    else:
                        unfinished.append(path)
            paths = ' '.join(paths_of_completed_rxns)

            #os.system(f"dptools convert {paths} {dir}_iter1.traj")
            #os.system(f"dptools input -a -p {path_to_training_data} {rootdir}/{dir}_iter1.traj")
print(len(unfinished))
#print(' '.join(paths_of_completed_rxns))
#paths= ' '.join(paths_of_completed_rxns)
#os.system(f"dptools convert {paths} all_spes_for_training.traj")