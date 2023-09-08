#!/bin/bash                                                                                                                                                
#                                                                                                                                                          
#SBATCH --job-name=wheat_cf_selhausen
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=cpu256
#SBATCH --time=1000:00:00
#SBATCH --mem=200G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END
#SBATCH --mail-user=s.ullah@fz-juelich.de
 
cd $HOME/Documents/Dumux/CPlantBox/experimental/selhausen_wheat_CF/

python3 multiruns.py