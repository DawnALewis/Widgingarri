#!/bin/bash
ml Anaconda3
source activate malt_0.61
malt-run \
        -i ${INPUT} \
        --index /hpcfs/groups/acad_users/Sediments_shared/nt_2019Nov19_step3 \
        -o ${OUTPUT}.rma6 \
        -at SemiGlobal \
	-mem load \
	--mode BlastN \
	--alignments /gpfs/users/a1867445/Widgingarri/malt_conda/alignments/ \
	--format SAM \
	-t 72 -v -mpi 85 -supp 0.001
