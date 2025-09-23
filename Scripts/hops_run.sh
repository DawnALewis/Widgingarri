#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH --time=72:00:00
#SBATCH --mem=199GB
#SBATCH -c 32

ml Singularity
singularity exec -B /gpfs/ /hpcfs/groups/acad_users/shyrav/resources/singularity/nf-core-eager_2.4.5-sharding.sif hops -i *.rma6 -output HOPS -m me_po -c HOPS_configFile_Melioidosis.txt
