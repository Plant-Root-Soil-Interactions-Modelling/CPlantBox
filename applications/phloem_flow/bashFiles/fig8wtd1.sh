#!/bin/bash
#
#SBATCH --job-name=fig8wtd_1
#SBATCH --ntasks=256
#SBATCH --nodes=1
#SBATCH --partition=cpu256
#SBATCH --time=20-00:00:00
#SBATCH --mem=450G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END,FAIL,ALL
#SBATCH --mail-user=m.giraud@fz-juelich.de


cd $HOME/auxin_sizeBuds/CPlantBox/applications/phloem_flow

python3 ExpCalib_WTD06062023.py 1
#zip AllAuxC1 AllAuxC1/*.txt

#sbatch --nodelist=node03 All_C1.sh
