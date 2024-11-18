import os
import time

  
def tag(experiment):
    tag_string = f"{experiment[0]}_domain_N{experiment[1]}_O{experiment[3]}"
    return tag_string

def submit(jobid,cmd, ncores=1):
    id = str(jobid)
    jobname = 'job_' + id
    memcore = 10000
    maxmem = 8000
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
    cmd = 'jupyter nbconvert --execute --to notebook intro_example_Poisson_2D_long.ipynb --output intro_example_Poisson_2D_long_cpy2_2.ipynb'
    timestr = time.strftime("%Y%m%d_%H%M%S")
    print(timestr)
    tag = 'myula_'+timestr
    print(cmd)
    print(tag)
    submit(tag, cmd)