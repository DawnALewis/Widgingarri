To run MALT in the eager pipeline, you need a taxon list (likely due to primary use in pathogen studies). Widgingarri did not have any target taxa so I ran Eager per Eager_pipeline.sh and ran malt the unmapped, post-filtering (low complexity removed) libraries through MALT alignment in the following way:

##### Metagenomic alignment after deduplication through eager pipeline
run_malt.sh
```
for l in /gpfs/users/a1867445/Widgingarri/results/metagenomic_complexity_filter/*.gz; do
sbatch -J malt_${l} -o $PWD/malt_${l}.log -D $PWD -N 1 -A strategic -p highmem  --mem=1200G -c 72 --time=72:00:00 malt.sh
done
```
malt.sh
```
#!/bin/bash

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=dawn.lewis@adelaide.edu.au 

ml Singularity

singularity exec -B /gpfs/ /hpcfs/groups/acad_users/containers/nf-core-eager_2.4.5-sharding.sif malt-run -i ${l} --index /hpcfs/groups/acad_users/Metagenomic_screening_db/malt_nt_2019Nov19_step3/ -o /gpfs/users/a1867445/Widgingarri/results/metagenomic_classification/malt/ -at SemiGlobal -mem load --mode BlastN --alignments /gpfs/users/a1867445/Widgingarri/results/metagenomic_classification/malt/test/alignments/ --format SAM -t 72 -v -mpi 85 -wlca -supp 0.001 -lcp 80

```
change the LCA parameters with malt alignments
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



