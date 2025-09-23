To run MALT in the eager pipeline, you need a taxon list (likely due to primary use in pathogen studies). Widgingarri did not have any target taxa so I ran Eager per Eager_pipeline.sh and ran malt the unmapped, post-filtering (low complexity removed) libraries through MALT alignment in the following way:

##### Metagenomic alignment after deduplication through eager pipeline
```
Libraries="Seda24WidgLib_10 Seda24WidgLib_11 Seda24WidgLib_12 Seda24WidgLib_13 Seda24WidgLib_14 Seda24WidgLib_15 Seda24WidgLib_16 Seda24WidgLib_17 Seda24WidgLib_18 Seda24WidgLib_19 Seda24WidgLib_1 Seda24WidgLib_20 Seda24WidgLib_21 Seda24WidgLib_22 Seda24WidgLib_23 Seda24WidgLib_24 Seda24WidgLib_2 Seda24WidgLib_3 Seda24WidgLib_4 Seda24WidgLib_5 Seda24WidgLib_6 Seda24WidgLib_7 Seda24WidgLib_8 Seda24WidgLib_9 Seda24WidgLibs_25 Seda24WidgLibs_26 Seda24WidgLibs_27 Seda24WidgLibs_28 Seda24WidgLibs_29 Seda24WidgLibs_30 Seda24WidgLibs_31 Seda24WidgLibs_32 Seda24WidgLibs_3"

for l in ${Libraries}; do  sbatch -J malt_${l} -o $PWD/malt_${l}.log -D $PWD -N 1 -A strategic -p highmem  --mem=1200G -c 72 --time=72:00:00 --export INPUT=/gpfs/users/a1867445/Widgingarri/results/metagenomic_complexity_filter/${l}.unmapped.fastq.gz_lowcomplexityremoved.fq.gz,OUTPUT=/gpfs/users/a1867445/Widgingarri/malt_conda/${l}_malt malt.sh ; done
```
malt.sh
```
#!/bin/bash

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=dawn.lewis@adelaide.edu.au 

ml Singularity

singularity exec -B malt-run.vmoptions:/usr/local/opt/malt-0.61/malt-run.vmoptions /hpcfs/groups/acad_users/containers/malt_0.61--hdfd78af_0.sif malt-run -i /gpfs/users/a1867445/Widgingarri/results/metagenomic_classification/malt/$l.unmapped.fastq.gz_lowcomplexityremoved.rma6 --index /hpcfs/groups/acad_users/Sediments_shared/nt_2019Nov19_step3 -o /gpfs/users/a1867445/Widgingarri/results/metagenomic_classification/malt/test/ -at SemiGlobal -mem load --mode BlastN --alignments /gpfs/users/a1867445/Widgingarri/results/metagenomic_classification/malt/test/alignments/ --format SAM -t 72 -v -mpi 85 -wlca -supp 0.001 -lcp 80
```


### meta comparison after MALT run
```
singularity exec -B /gpfs/ /hpcfs/groups/acad_users/containers/megan_6.24.20--h9ee0642_0.sif compute-comparison -i /gpfs/users/a1867445/Widgingarri/results/metagenomic_classification/malt/*.rma6 -o /gpfs/users/a1867445/Widgingarri/results/metagenomic_classification/malt/Widg85samples.megan -n false
```
### To make an extensive taxa list for HOPS
find_unique_taxa.sh 
```
#!/bin/bash

# Define input directory and output file
INPUT_DIR="/gpfs/users/a1867445/Widgingarri/results/metagenomic_classification/malt"
OUTPUT_FILE="widgingarri_unique_taxa.txt"
TEMP_FILE="all_taxa_temp.txt"

# Remove previous output files if they exist
rm -f $OUTPUT_FILE $TEMP_FILE

# Loop through all .rma6 files in the directory
for file in "$INPUT_DIR"/*.rma6; do
    echo "Processing: $file"
    
    # Extract sample name (removing the path and extension)
    sample_name=$(basename "$file" .unmapped.fastq.gz_lowcomplexityremoved.rma6)

    # Run MEGAN's rma2info to extract taxonomy classifications
    singularity exec -B /gpfs/ /hpcfs/groups/acad_users/containers/megan_6.24.20--h9ee0642_0.sif \
    rma2info -i "$file" -r2c Taxonomy -n >> "$TEMP_FILE"

done

# Extract only unique taxa
sort "$TEMP_FILE" | uniq > "$OUTPUT_FILE"

# Remove temporary file
rm -f "$TEMP_FILE"

echo "Unique taxa saved in: $OUTPUT_FILE"
```
### run HOPS
```
sbatch -A strategic -p highmem --mem 600G run_hops.sh
```

```
#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH --time=72:00:00
#SBATCH --mem=199GB
#SBATCH -c 32

ml Singularity
singularity exec -B /gpfs/ /hpcfs/groups/acad_users/shyrav/resources/singularity/nf-core-eager_2.4.5-sharding.sif hops -i *.rma6 -output HOPS -m me_po -c HOPS_configFile_Melioidosis.txt



