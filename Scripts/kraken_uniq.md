
```
bash krakenuniq.sh
```

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

```
bash run_krakenuniq_processing.sh /gpfs/users/a1867445/Cloggs_May_2025/eager/krakenUniq/
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

### krakenuniq_abundances.R
```
#This is a script for plotting KrakenUniq abundance matrix.
#Run this script as:
#Rscript plot_krakenuniq_abundance_matrix.R in_dir out_dir

args = commandArgs(trailingOnly=TRUE)
in_dir<-as.character(args[1])
out_dir<-as.character(args[2])

library("pheatmap")
ku_abundance<-read.delim(paste0(in_dir,"/krakenuniq_abundance_matrix.txt"),header=TRUE,row.names=1,check.names=FALSE,sep="\t")


#FUNCTION FOR TUNING FONTSIZE ON MICROBIAL ABUNDANCE HEATMAP
my_fontsize<-function(ku_abundance)
{
  if(dim(ku_abundance)[1]<50){return(12)}
  else if(dim(ku_abundance)[1]>=50 & dim(ku_abundance)[1]<100){return(10)}
  else if(dim(ku_abundance)[1]>=100 & dim(ku_abundance)[1]<150){return(8)}
  else{return(6)}
}


#ABSOLUTE ABUNDANCE HEATMAP
pdf(paste0(out_dir,"/krakenuniq_absolute_abundance_heatmap.pdf"),paper="a4r",width=297,height=210)
if(dim(ku_abundance)[1]>1 & dim(ku_abundance)[2]>1)
{
  pheatmap(ku_abundance, display_numbers=TRUE,fontsize=my_fontsize(ku_abundance),
           main="KrakenUniq Absolute Abundance",cluster_rows=FALSE,cluster_cols=FALSE,number_format="%i")
}else
{
  pheatmap(ku_abundance, display_numbers=TRUE,fontsize=8,
           main="KrakenUniq Absolute Abundance",cluster_rows=FALSE,cluster_cols=FALSE,number_format="%i",breaks=c(0,1))
}
dev.off()


#NORMALIZE BY SEQUENCING SEPTH
for(i in 1:dim(ku_abundance)[2])
{
  ku_abundance[,i]<-ku_abundance[,i]/sum(ku_abundance[,i])
}


#ABSOLUTE ABUNDANCE HEATMAP
pdf(paste0(out_dir,"/krakenuniq_normalized_abundance_heatmap.pdf"),paper="a4r",width=297,height=210)
if(dim(ku_abundance)[1]>1 & dim(ku_abundance)[2]>1)
{
  pheatmap(ku_abundance, display_numbers=TRUE,fontsize=my_fontsize(ku_abundance),
           main="KrakenUniq Normalized Abundance",cluster_rows=FALSE,cluster_cols=FALSE,number_format="%.3f")
}else
{
  pheatmap(ku_abundance, display_numbers=TRUE,fontsize=8,
           main="KrakenUniq Normalized Abundance",cluster_rows=FALSE,cluster_cols=FALSE,number_format="%.3f",breaks=c(0,1))
}
```

### Get full lineage from uniq IDs
```
for tax in $(cat /gpfs/users/a1867445/Cloggs_May_2025/eager/krakenUniq/output/abundance/unique_species_taxid_list.txt); do grep -P "^$tax\s" /hpcfs/groups/acad_users/dawn/scripts/krakenuniq/fullnamelineage.dmp ; done | tee full_names_for_my_species.txt
```
### get into order
```
for species in $(tail -n+1 output/abundance/krakenuniq_abundance_matrix.txt | awk '{print $1;}'); do grep -P "^$species\s" full_names_for_my_species.txt ; done | awk -F'|' '{print $1","$2","$3}'  | sed 's/\t//g' > fullnames_same_order_as_matrix.csv
```

##### fullnames_same_order_as_matrix.csv and krakenuniq_abundance_matrix can be merged manually if you want.
