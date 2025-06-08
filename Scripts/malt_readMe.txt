### Metagenomic alignment after deduplication through eager pipeline

Libraries="Seda24WidgLib_10 Seda24WidgLib_11 Seda24WidgLib_12 Seda24WidgLib_13 Seda24WidgLib_14 Seda24WidgLib_15 Seda24WidgLib_16 Seda24WidgLib_17 Seda24WidgLib_18 Seda24WidgLib_19 Seda24WidgLib_1 Seda24WidgLib_20 Seda24WidgLib_21 Seda24WidgLib_22 Seda24WidgLib_23 Seda24WidgLib_24 Seda24WidgLib_2 Seda24WidgLib_3 Seda24WidgLib_4 Seda24WidgLib_5 Seda24WidgLib_6 Seda24WidgLib_7 Seda24WidgLib_8 Seda24WidgLib_9 Seda24WidgLibs_25 Seda24WidgLibs_26 Seda24WidgLibs_27 Seda24WidgLibs_28 Seda24WidgLibs_29 Seda24WidgLibs_30 Seda24WidgLibs_31 Seda24WidgLibs_32 Seda24WidgLibs_3

for l in ${Libraries}; do  sbatch -J malt_${l} -o $PWD/malt_${l}.log -D $PWD -N 1 -A strategic -p highmem  --mem=1200G -c 72 --time=72:00:00 --export INPUT=/gpfs/users/a1867445/Widgingarri/results/metagenomic_complexity_filter/${l}.unmapped.fastq.gz_lowcomplexityremoved.fq.gz,OUTPUT=/gpfs/users/a1867445/Widgingarri/malt_conda/${l}_malt malt_conda.sh ; done

### meta comparison after MALT run
singularity exec -B /gpfs/ /hpcfs/groups/acad_users/containers/megan_6.24.20--h9ee0642_0.sif compute-comparison -i /gpfs/users/a1867445/Widgingarri/malt_conda/*.rma6 -o /gpfs/users/a1867445/Widgingarri/malt_conda/Widg85samples.megan -n false
