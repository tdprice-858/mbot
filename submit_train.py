#!/usr/bin/env python
# -*- coding: utf-8 -*-

#SBATCH --job-name=tm_mgo_h2
#SBATCH --nodes=1
#SBATCH --qos=regular --time=00-24:00:00
#SBATCH --constraint=cpu
#SBATCH --account=m4126
#SBATCH --output=job.out
#SBATCH --error=job.err

export OMP_NUM_THREADS=16
export TF_INTRA_OP_PARALLELISM_THREADS=16
export TF_INTER_OP_PARALLELISM_THREADS=4

if [[ -e 'checkpoint' ]]; then
        dp train --restart model.ckpt in.json
else
        dp train in.json
fi

dp freeze -o graph.pb
