#!/bin/bash                                                                                                                                                

#SBATCH -N 1                          ##  Use N nodes                                                                                                      
#SBATCH -p debug                                                                                                                                          
#SBATCH -t 00:30:00                   ##  xx mins                                                                                                          
#SBATCH -J Rafelski                                                                                                                                    
#SBATCH -C haswell                                                                                                                                          

#SBATCH -o "/global/homes/m/mjwilson/Rafelski15/eazy-photoz/out.txt"                                                                 
#SBATCH -e "/global/homes/m/mjwilson/Rafelski15/eazy-photozxs/err.txt"                                                                

cd $SLURM_SUBMIT_DIR                                                                                            

source env.sh

srun -n 1 python py/genmags.py
