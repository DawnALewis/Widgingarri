### run_krakenuniq.sh
```
#!/bin/bash

# Directory containing the FASTQ files
FASTQ_DIR="/gpfs/users/a1867445/Cloggs_malt/lane_merged/"

# Find all FASTQ files in the specified directory
FASTQ_FILES=$(find "$FASTQ_DIR" -type f -name "*.gz")

# Loop over each FASTQ file
for FASTQ in $FASTQ_FILES; do
    PREFIX=$(basename "$FASTQ" .gz) # Strip .gz extension to use as prefix

    sbatch -J krakenUniq_${PREFIX} -D /hpcfs/groups/acad_users/dawn/cloggs_lane_merged/krakenuni_nt -o /hpcfs/groups/acad_users/dawn/cloggs_lane_merged/krakenuni_nt/${PREFIX}_krakenUniq.out -N 1 -c 64 -p icelake \
        --mem=150GB --time=12:00:00 \
        --export DB=/hpcfs/groups/acad_users/Metagenomic_screening_db/KRAKENUNIQ_databases/KrakenUniq_database_based_on_full_NCBI_NT_from_December_2020/,fastq=${FASTQ},prefix=${PREFIX},SIZE=100,CPU=64 \
        kraken.sh
done
```

### kraken.sh
```
#!/bin/bash

ml Singularity

singularity exec -B /gpfs/ /hpcfs/groups/acad_users/containers/krakenuniq_1.0.4--pl5321h6dccd9a_1.sif krakenuniq --preload-size 100G --db /hpcfs/groups/acad_users/Metagenomic_screening_db/KRAKENUNIQ_databases/KrakenUniq_database_based_on_full_NCBI_NT_from_December_2020/ \
    --fastq-input ${FASTQ} \
    --threads 64 \
    --output /gpfs/users/a1867445/Cloggs_May_2025/eager/krakenUniq/${PREFIX}.kraKenOut \
    --report-file /gpfs/users/a1867445/Cloggs_May_2025/eager/krakenUniq/${PREFIX}.report \
    --gzip-compressed --only-classified-out
```


### run_krakenuniq_processing.sh 
```
#!/bin/bash


indir=$1
module purge
module load Singularity

mkdir -p output/{formatted,abundance,plots}
for report in $(ls ${indir}/*report); do
    new_report=$(basename $report | sed 's/.report/.formatted.report/')
    sed -e "s/\r//g" ${report} > output/formatted/${new_report}

    singularity exec -B /gpfs/ /hpcfs/groups/acad_users/containers/ngspy_0.3.sif \
        python3 /hpcfs/groups/acad_users/shyrav/projects/metagenomic_screening/krakenuniq/krakenuniq_filter.py \
            output/formatted/${new_report} 1000 100
done

module purge
module load R/4.3.1-foss-2021b

Rscript /hpcfs/groups/acad_users/shyrav/projects/metagenomic_screening/krakenuniq/krakenuniq_abundances.R output/formatted ./output/abundance
Rscript /hpcfs/groups/acad_users/shyrav/projects/metagenomic_screening/krakenuniq/plot_abundances.R ./output/abundance ./output/plots
```
