import os
import time


def submit(jobid,cmd, ncores=1):
    id = str(jobid)
    jobname = 'job_' + id
    memcore = 10000
    email = 'amaal@dtu.dk'
    block_size = ncores if ncores <=24 else 24
 
    # begin str for jobscript
    strcmd = '#!/bin/sh\n'
    strcmd += '#BSUB -J ' + jobname + '\n'
    strcmd += '#BSUB -q compute\n'
    strcmd += '#BSUB -n ' + str(ncores) + '\n'
    strcmd += '#BSUB -R "span[block='+str(block_size)+']"\n'
    strcmd += '#BSUB -R "rusage[mem=' + str(memcore) + 'MB]"\n'
    strcmd += '#BSUB -W 60:00\n'
    strcmd += '#BSUB -u ' + email + '\n'
    strcmd += '#BSUB -N \n'
    strcmd += '#BSUB -B \n'
    strcmd += '#BSUB -o hpc/output/output_' + id + '.out\n'
    strcmd += '#BSUB -e hpc/error/error_' + id + '.err\n'
    strcmd += "source /dtu/3d-imaging-center/courses/conda/conda_init.sh \n"
    strcmd += 'conda activate /zhome/0f/0/161811/.conda/envs/fenicsproject \n'  
    strcmd += cmd
 
    jobscript = 'hpc/submit_'+ jobname + '.sh'
    f = open(jobscript, 'w')
    f.write(strcmd)
    f.close()
    os.system('bsub < ' + jobscript)

if __name__ == "__main__":
    smoothing_factor_list = [0.0001, 0.001, 0.01, 0.1, 1, 10]
    regularization_strength_list = [1, 5, 7, 10, 20, 30]
    par_dim_list = [int(32**2)]

    idx = 0
    for smoothing_factor in smoothing_factor_list:
        for regularization_strength in regularization_strength_list:
            for par_dim in par_dim_list:
                idx += 1
                cmd = f"NB_ARGS=' --smoothing_factor {smoothing_factor} --regularization_strength {regularization_strength} --par_dim {par_dim}' jupyter nbconvert --execute --to notebook Poisson_2D_MYULA.ipynb --output Poisson_2D_MYULA{smoothing_factor}_{regularization_strength}_{par_dim}.ipynb"
                tag = "Poisson_2D_MYULA"+str(idx)+"_" + str(smoothing_factor) + "_" + str(regularization_strength) + "_" + str(par_dim)
                print(cmd)
                print(tag)
                submit(tag, cmd)
