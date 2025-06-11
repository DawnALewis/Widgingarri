#!/bin/bash

ml Singularity
singularity exec -B /gpfs/ /hpcfs/groups/acad_users/containers/nf-core-eager_2.4.5-sharding.sif hops \
  -m me_po \
  -i /gpfs/users/a1867445/Widgingarri/malt_conda/*.rma6 \
  -o /gpfs/users/a1867445/Widgingarri/malt_conda/hops/
