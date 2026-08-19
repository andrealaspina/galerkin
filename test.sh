#!/bin/sh
#SBATCH --job-name=test                          # Job name
#SBATCH -c 2                                     # Number of cores
#SBATCH --mem=6gb                                # Memory used
#SBATCH -t 0-00:20:00                            # Time to run (d-hh:mm:ss)
#SBATCH -A free                                  # Account name
#SBATCH --nodes=1                                # Number of nodes
#SBATCH --ntasks=1                               # Number of tasks
#SBATCH -e %x.err                                # Error file
#SBATCH -o %x.out                                # Output file
#SBATCH --mail-type=ALL                          # Email type (ALL/BEGIN/END/FAIL/NONE)
#SBATCH --mail-user=username@domain.com          # Email address

# Root of the galerkin installation
GALERKIN=/galerkin

module load matlab

matlab -nodisplay -nosplash -r                                                            \
  "cd('${GALERKIN}'); init;                                                               \
   distcomp.feature('LocalUseMpiexec',false);                                             \
   delete(gcp);                                                                           \
   parpool(parcluster('local'),str2num(getenv('SLURM_CPUS_ON_NODE')),'IdleTimeout',6*60); \
   run(sprintf('./%s.m',getenv('SLURM_JOB_NAME'))); quit"
