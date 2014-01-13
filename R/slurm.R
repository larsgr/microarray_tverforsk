# generate header for SLURM batch script

generateBatchScript <- function(jobName,command,cores=1){
  return( sprintf("#!/bin/sh\n#SBATCH -n %-12i # reserved cores(CPUs)\n#SBATCH -J %-12s # name for the job\n%s\n",
                  cores,jobName,command))
}
