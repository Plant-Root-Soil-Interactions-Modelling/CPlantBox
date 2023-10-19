#!/bin/bash                                                                                                                                                
#                                                                                                                                                          
#SBATCH --job-name=Zea_mays_1_Leitner_2010_cf_selhausen
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=cpu256
#SBATCH --time=1000:00:00
#SBATCH --mem=200G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END
#SBATCH --mail-user=s.ullah@fz-juelich.de
 
cd $HOME/Documents/Dumux/CPlantBox/experimental/selhausen_maize_CF
python3 multiruns.py