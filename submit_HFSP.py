import os, sys
import json
rootdir = '/global/cfs/cdirs/m3548/tdprice/10_ML_training_loop/iter1'
number_of_spes = []
number_of_completed_rxns = []
os.system('export PATH="/global/cfs/cdirs/m3548/tdprice/10_ML_training_loop/ml_pynta/pynta:$PATH"')
os.system('export PYTHONPATH="/global/cfs/cdirs/m3548/tdprice/10_ML_training_loop/ml_pynta/pynta:$PYTHONPATH"')
for dir in os.listdir(rootdir):
    if dir.split('_')[0]=='rxn':
        #print(dir)
        if os.path.exists(os.path.join(rootdir,dir,'TS0/structures_for_ML.json')):
            print(dir)
            print('HFSP finished')
            '''number_of_completed_rxns.append(dir)
            with open(os.path.join(rootdir,dir,'TS0/structures_for_ML.json')) as f:
                structs = json.load(f)
                for struct in structs:
                    path = '/'.join(struct.split('/')[:-1])
                    #print(path)
                    os.chdir(path)
                    #print(os.getcwd())
                    number_of_spes.append(path)
                    #if os.path.exists(os.path.join(path,'spe_RPBE_500_441')):

                    if not os.path.exists(os.path.join(path,'spe_RPBE_500_441')):
                        print('spe does not exist')
                        print(dir)
                        #os.system('cp ../../../submit_spe.py .')
                        #os.system('sbatch submit_spe.py')'''
        else:
            if os.path.exists(os.path.join(rootdir, dir)):
                print(os.path.join(rootdir, dir))
                os.chdir(os.path.join(rootdir, dir))
                os.system('sbatch run.script')

            #print(dir)
            #       os.system('sbatch submit_spe.py')
            #print('TS opt did not finish ^^^^')




print(len(number_of_spes))
print(len(number_of_completed_rxns))