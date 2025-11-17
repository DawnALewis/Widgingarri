# BIOM data 85% ID and 0.001% abundance

These chunks are set up to rely on previous chunks. May need to re-read files or load libs if only using a single chunk.

## Decontam

Use megan to extract BIOM1 files from for:

ALL EBCs and Samples (Widg85and001samples-Taxonomy.biom)

### otu_plot

```{r}
setwd("/Users/dawnlewis/Library/CloudStorage/Box-Box/Widgingarri/2025/R_stuff")

#PIPELINE CON BIOM OFICIAL
library(phyloseq) 
library(tidyr)
library(decontam)
library(ggplot2)

all_data <- import_biom("/Users/dawnlewis/Library/CloudStorage/Box-Box/Widgingarri/Widg_OctResults/85ID/Widg85and001samples-Taxonomy.biom")
data_t <- t(all_data)

metadata <- import_qiime_sample_data("/Users/dawnlewis/Library/CloudStorage/Box-Box/Widgingarri/2025/R_stuff/decontamMetadata_Widg.txt")

data_metadata <- merge_phyloseq(data_t, metadata)

df <- as.data.frame(sample_data(data_metadata)) 
df$LibrarySize <- sample_sums(data_metadata)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

data_metadata_RA <- transform_sample_counts(data_metadata, function(OTU)(OTU/sum(OTU)*100))

OT <- otu_table(data_metadata_RA)
mat <- as(OT, "matrix")
mat[!is.finite(mat)] <- 0  # handles NA, NaN, Inf
otu_table(data_metadata_RA) <- otu_table(mat, taxa_are_rows = taxa_are_rows(data_metadata_RA))

otu_table(data_metadata_RA) <- replace_na(otu_table(data_metadata_RA), 0)
```

### Prevalence

False means not a contaminant, True means it is a contaminant

```{r}
#cont. from above
#prevalence (default threshold=0.1)
sample_data(data_metadata_RA)$is.neg <- sample_data(data_metadata_RA)$Sample_or_Control == "Control"
contamdf.prev <- isContaminant(data_metadata_RA, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)


#prevalence threshold=0.5
sample_data(data_metadata_RA)$is.neg <- sample_data(data_metadata_RA)$Sample_or_Control == "Control"
contamdf.prev05 <- isContaminant(data_metadata_RA, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)

#change contamdf.prev if want to change threshold
write.table(contamdf.prev05, file = "contam_Widg_05.txt", sep = "\t", quote = FALSE)
```

### Visualising prevalence

```{r}
#Look at the number of times several of these taxa were observed in negative controls and positive samples.
data_metadata_RA_samples <- prune_samples(sample_data(data_metadata_RA)$Sample_or_Control == "Sample", data_metadata_RA)
data_metadata_RA_control <- prune_samples(sample_data(data_metadata_RA)$Sample_or_Control == "Control", data_metadata_RA)

# Make data.frame of prevalence in positive and negative samples
df.RA <- data.frame(RA.samples=taxa_sums(data_metadata_RA_samples), RA.control=taxa_sums(data_metadata_RA_control),
contaminant=contamdf.prev05$contaminant)
ggplot(data=df.RA, aes(x=RA.control, y=RA.samples, color=contaminant)) + geom_point() +
                      xlab("Prevalence (Controls)") + ylab("Prevalence (Samples)")
                    
ggplot(data = contamdf.prev05, aes(x=p)) +
                      geom_histogram(binwidth = 0.01) +
                      labs(x = 'decontam Score', y='Number of species')
                    
```

## Data Analysis

Use megan to extract BIOM1 files from for:

\>ALL EBCs and Samples (Widg85and001samples-Taxonomy.biom)

\>EBCs only (WidgEBC85and001samples-Taxonomy.biom)

\>BONE samples only (WidgBONES85and001samples-Taxonomy.biom)

Unhash #all_data and #data_t if not already loaded.

This uses a different metadata file with EVERYTHING you have. Row names need to be unique and match the biom files.

```{r}
library(phyloseq)
library(readr)
library(dplyr)
library(vegan)
library(colorspace)

EBCs <-import_biom("/Users/dawnlewis/Library/CloudStorage/Box-Box/Widgingarri/Widg_OctResults/85ID/WidgEBC85and001samples-Taxonomy.biom")
#all_data <- import_biom("Widg85and001samples-Taxonomy.biom")
#data_t <- t(all_data)
metadata<-import_qiime_sample_data("/Users/dawnlewis/Library/CloudStorage/Box-Box/Widgingarri/2025/R_stuff/Sample_metadata.txt")
data_metadata <- merge_phyloseq(data_t, metadata)

```

Check the contamination file "contam_Widg.txt" is the same as the output form Decontam above. And it is located in the working directory

```{r}
## 1) Make sure taxa orientation matches across objects & flip if needed (unhash the first two and hash the second two #if)
#if (!taxa_are_rows(EBCs))  EBCs  <- tax_glom(t(EBCs), "none") %>% t()  
#if (!taxa_are_rows(data_t)) data_t <- tax_glom(t(data_t), "none") %>% t()

if (!taxa_are_rows(EBCs))  EBCs  <- t(EBCs)
if (!taxa_are_rows(data_t)) data_t <- t(data_t)

## 2) Get taxa that appear (non-zero) in EBCs
ebc_taxa <- taxa_names(EBCs)

## 3) Remove all taxa present in EBCs from your main object
data_noEBC <- prune_taxa(!(taxa_names(data_metadata) %in% ebc_taxa), data_metadata)

## 4) Read decontam results and pull taxa flagged TRUE
# Expecting a column named 'contaminant' and a taxa/feature ID that matches phyloseq taxa_names
decontam <- read_tsv("contam_Widg_05.txt", show_col_types = FALSE)
colnames(decontam)[1] <- "TaxaID" #add TaxaID in the first column name
decontam$contaminant <- as.logical(decontam$contaminant)

contam_ids <- decontam %>%
  filter(contaminant) %>%
  pull(TaxaID) %>%
  as.character()

## 5) Remove taxa flagged as contaminants
data_clean <- prune_taxa(!(taxa_names(data_noEBC) %in% contam_ids), data_noEBC)

####Remove samples####
## Drop zero-abundance taxa/samples just in case
data_clean <- prune_taxa(taxa_sums(data_clean) > 0, data_clean)
data_clean <- prune_samples(sample_sums(data_clean) > 0, data_clean)

# Remove specific samples by name (unhash as needed)
#data_clean <- prune_samples(!(sample_names(data_clean) %in% c("440cm_VP", "25cm_VP")), data_clean)
#data_clean <- prune_taxa(taxa_sums(data_clean) >= 2, data_clean)


## Quick before/after summary
cat("Before: ", nsamples(data_metadata), "samples,", ntaxa(data_metadata), "taxa\n")

cat("After : ", nsamples(data_clean),    "samples,", ntaxa(data_clean),    "taxa\n")


```

# Rarefaction

Set the target based on what makes sense for your data.

```{r}
#cont from above.
mat <- as(otu_table(data_clean), "matrix")
# Ensure samples are rows
if (taxa_are_rows(data_clean)) {
  mat <- t(mat)
}

depths <- unique(pmin(rowSums(mat), round(seq(1, max(rowSums(mat)), length.out = 50))))

df_list <- lapply(depths, function(d){
  rd <- rarefy(mat, sample = d, se = FALSE)
  data.frame(Sample = rownames(mat), Depth = d, Richness = as.numeric(rd))
})

df <- dplyr::bind_rows(df_list)

ggplot(df, aes(x = Depth, y = Richness, color = Sample)) +
  geom_line(alpha = 0.7) +
  labs(x = "Sequencing depth (reads)", y = "Expected richness (rarefied)") +
  theme_minimal()


###check how many reads in each library before rarefying
# Number of reads per sample
reads_per_sample <- sample_sums(data_clean)

# Show summary
summary(reads_per_sample)

# Show all samples with counts
reads_per_sample

##Rarefy
# 1) Choose target = min library size so no sample is dropped
#target <- min(sample_sums(data_clean))
#target <- mean(sample_sums(data_clean))
target <- median(sample_sums(data_clean))

# 2) Rarefy (without replacement; reproducible rngseed)
set.seed(123)
data_clean_rare <- rarefy_even_depth(
  data_clean,
  sample.size = target,
  rngseed = 123,
  replace = FALSE,
  trimOTUs = TRUE,   # drop taxa that end up 0 after rarefaction
  verbose = TRUE
)

# 3) Sanity check
summary(sample_sums(data_clean_rare))

```

## RelAbund and Taxonomy Fix

shouldn't need to change anything here

```{r}

#Relative abundance of data_clean_rare
data_clean_rare_RA <- transform_sample_counts(data_clean_rare, function(OTU)(OTU/sum(OTU)*100))

###Fix taxonomy####
library(phyloseq)
# ps = your phyloseq object (here called data_clean)
ps <- data_clean_rare_RA
tx <- as.matrix(tax_table(ps))
# Helper: extract the first non-NA cell in a row that starts with a given prefix,
# and strip the prefix (e.g., "p__" -> "Bacillota")
extract_first_pref <- function(row, pref) {
  v <- row[!is.na(row) & startsWith(row, pref)]
  if (length(v)) sub(paste0("^", pref), "", v[1]) else NA_character_
}
# Build a new, tidy taxonomy by prefix (domain/kingdom/phylum/.../species)
wanted <- c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species")

new_tx <- t(apply(tx, 1, function(row) c(
  Domain  = extract_first_pref(row, "d__"),
  Kingdom = extract_first_pref(row, "k__"),
  Phylum  = extract_first_pref(row, "p__"),
  Class   = extract_first_pref(row, "c__"),
  Order   = extract_first_pref(row, "o__"),
  Family  = extract_first_pref(row, "f__"),
  Genus   = extract_first_pref(row, "g__"),
  Species = extract_first_pref(row, "s__")
)))

new_tx <- as.data.frame(new_tx, stringsAsFactors = FALSE, check.names = FALSE)
colnames(new_tx) <- wanted

valid_domain <- c("Bacteria", "Archaea", "Eukaryota")
#Remove NA (viruses, uncultured, etc)
new_tx <- new_tx[!is.na(new_tx$Domain) & new_tx$Domain %in% valid_domain, ]

fill_down_taxonomy <- function(df) {
  rank_order <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  df <- df[, rank_order, drop = FALSE]
  for (i in seq_len(nrow(df))) {
    # get first non-NA rank going top-down and fill below it if NA
    row <- df[i, ]
    for (r in seq_along(rank_order)) {
      if (r > 1 && is.na(row[[rank_order[r]]])) {
        row[[rank_order[r]]] <- row[[rank_order[r - 1]]]
      }
    }
    df[i, ] <- row
  }
  df
}

new_tx_filled <- fill_down_taxonomy(new_tx)
# Write back to the phyloseq object with clean rank names
tax_table(ps) <- as.matrix(new_tx_filled[, c("Domain", "Kingdom","Phylum","Class","Order","Family","Genus","Species")])
# Put standard ranks first
tax_table(ps) <- tax_table(ps)[, c("Domain", "Kingdom","Phylum","Class","Order","Family","Genus","Species")]
rank_names(ps)
```

# Bar Plot ALL

```{r}
#| label: Data
library(microbiome)
library(colorspace)
library(dplyr)
library(ggplot2)


ps.family <- tax_glom(ps, taxrank = "Family")


# --- 1. Make sure we're working at Family level ---
ps.family <- tax_glom(ps, taxrank = "Family")

# --- 2. Identify top 5 families by total abundance ---
tax_sums <- taxa_sums(ps.family)
top5_otus <- names(sort(tax_sums, decreasing = TRUE)[1:5])

# Get the Family names of these OTUs
top5_families <- tax_table(ps.family)[top5_otus, "Family"] %>% as.character() %>% na.omit()

# --- 3. Specify non-specific families ---
non_specific <- c("Bacteria", "Archaea")

# --- 4. Combine all families of interest ---
families_to_keep <- c(top5_families, non_specific)

# --- 5. Identify taxa IDs that belong to these families ---
keep_taxa <- rownames(tax_table(ps.family))[tax_table(ps.family)[, "Family"] %in% families_to_keep]

# --- 6. Prune phyloseq object to keep only these taxa ---
ps.keep <- prune_taxa(keep_taxa, ps.family)

# --- 7. Replace OTU IDs with Family names for clarity ---
taxa_names(ps.keep) <- tax_table(ps.keep)[, "Family"]

# --- 8. Extract counts (taxa x samples) ---
counts <- otu_table(ps.keep)
if (!taxa_are_rows(ps.keep)) counts <- t(counts)  # ensure taxa are rows

# --- 9. Compute proportions relative to total reads per sample ---
total_reads <- sample_sums(ps.family)   # total reads per sample, including all families
prop_matrix <- sweep(counts, 2, total_reads, "/")  # divide each column by total reads

# 1. Start from the prop_matrix (taxa x samples)
prop_wide <- as.data.frame(t(prop_matrix))  # samples as rows, families as columns
colnames(prop_wide) <- taxa_names(ps.keep)  # set family names as column headers

# 2. Add Library/sample IDs as a column
prop_wide$Library <- rownames(prop_wide)

# 3. Join metadata
prop_wide <- prop_wide %>%
  left_join(metadata_df %>% mutate(Library = rownames(metadata_df)), by = "Library")

# 4. Optional: reorder columns so Library + metadata come first
metadata_cols <- setdiff(colnames(metadata_df), "Library")
prop_wide <- prop_wide %>%
  select(Library, all_of(metadata_cols), everything())

# 5. Check
head(prop_wide)



## Top 40 
ps.family.40 <- prune_taxa(names(sort(taxa_sums(ps.family),TRUE)[1:40]), ps.family)
#What % of reads do the top 40 Family have? ## 0.6752479, 0.03433144
print("mean:")
mean(sort(sample_sums(ps.family.40)/sample_sums(ps.family))) %>% print()
print("SD:")
sd(sort(sample_sums(ps.family.40)/sample_sums(ps.family))) %>% print() 
#Looks like the top 40 phyla account for the vast majority of sequences 


#Convert to compositional (relative abundance)
ps.family_RA40 <- microbiome::transform(ps.family.40, "compositional")

#Melt into a dataframe
pdfamily_RA40 <- psmelt(ps.family_RA40)

pd.ordered40 <- pdfamily_RA40 %>%
  group_by(Depth_m) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%  # normalize per depth
  ungroup()

#Select depth column
pd.ordered40$Depth_m <- factor(pd.ordered40$Depth_m, 
                                levels = sort(unique(pd.ordered40$Depth_m)),
                                ordered = TRUE)

```
#### Colour Pallette - all plots
#| label: ColourPalette
family_cols <- c(
  "Burkholderiales"="#CAB2D6",
  "Nitrososphaeraceae"="#FB8072",
  "Streptomycetaceae"="#D95F02",
  "Mycobacteriaceae"="#66A61E",
  "Pseudonocardiaceae"="#8B8000",
  "Gemmatimonadaceae"="#E69F00",
  "Bacteria"="#CC79A7",
  "Malvaceae"="#00FF7F",
  "Bradyrhizobiaceae"="#B8860B",
  "Burkholderiaceae"="#FF1493",
  "Comamonadaceae"="#A540AB",
  "Gemmataceae"="#D55E00",
  "Streptosporangiaceae"="#9467BD",
  "Hyphomicrobiaceae"="#8C564B",
  "Acidobacteriaceae"="#FFD92F",
  "Alcaligenaceae"="#E6AB02",
  "Micromonosporaceae"="#191970",
  "Conexibacteraceae"="#228B22",
  "Myxococcaceae"="#56B4E9",
  "Thermomonosporaceae"="#009E73",
  "Archaea"="#80B1D3",
  "Vicinamibacteraceae"="#A6761D",
  "Baekduiaceae"="#FDB462",
  "Isosphaeraceae"="#BC80BD",
  "Rhizobiaceae"="#1B9E77",
  "Microbacteriaceae"="#E41A1C",
  "Rhodospirillaceae"="#FF6347",
  "Phyllobacteriaceae"="#FB9A99",
  "Nitrospiraceae"="#6B83D5",
  "Nocardioidaceae"="#8DD3C7",
  "Nocardiaceae"="#B3B3B3",
  "Pseudomonadaceae"="#B2DF8A",
  "Sphingomonadaceae"="#B3DE69",
  "Stellaceae"="#999999",
  "Xanthobacteraceae"="#FF7F00",
  "Archangiaceae"="#7FCDBB",
  "Anaeromyxobacteraceae"="#CAB2D6",
  "Oxalobacteraceae"="#1F78B4",
  "Xanthomonadaceae"="#41AE76",
  "Planctomycetaceae"="#A6D96A",
  "Acidiferrobacteraceae"="#66C2A5",
  "Micrococcaceae"="#3288BD",
  "Opitutaceae"="#2CA25F",
  "Thaumarchaeota"="#006D2C",
  "Euryarchaeota"="#55AB98",
  "Natrialbaceae"="#CCEBC5",
  "Halobacteria"="#8DA0CB",
  "Halobacteriaceae"="#FBB4AE",
  "Halorubraceae"="#B7AAE8",
  "Thermococcaceae"="#E5C494",
  "Haloarculaceae"="#F2B447",
  "Haloferacaceae"="#FED9A6",
  "Crenarchaeota"="#E78AC3",
  "Nitrosopumilaceae"="#FDD0A2",
  "Methanomicrobiaceae"="#FF7F00",
  "Methanocellaceae"="#A60026",
  "Methanotrichaceae"="#B62E2E",
  "Methanomassiliicoccaceae"="#C4543C",
  "Methanosarcinaceae"="#D0794A",
  "Methanoregulaceae"="#D99F5A",
  "Methanomicrobiales"="#E2C56B",
  "Haloferacales"="#E8D97D",
  "Candidatus Methanomethylophilaceae"="#D1E289",
  "Picrophilaceae"="#B4E699",
  "Cenarchaeaceae"="#96E9A8",
  "Thermoplasmatales"="#78EBB7",
  "Thermoproteaceae"="#5CEAC5",
  "Thermoplasmata"="#41E7D2",
  "Methanobacteriaceae"="#2DE1DD",
  "Candidatus Nitrosocaldaceae"="#26D6E3",
  "Halobacteriales"="#2EC7E2",
  "Sulfolobaceae"="#3BB5DA",
  "Methanocaldococcaceae"="#4A9FCF",
  "Archaeoglobaceae"="#5988C2",
  "Candidatus Thermoplasmatota"="#676EB3",
  "Thermoprotei"="#7454A2",
  "Conexivisphaeraceae"="#803A90",
  "Methanomicrobia"="#87257E",
  "Desulfurococcaceae"="#8B176C",
  "Acidilobaceae"="#8C125B",
  "Methanopyraceae"="#8A154B",
  ## Eukaryotes
  "Glomeraceae"="#7D2C3E",
  "Moraceae"="#6C3A33",
  "Phaffomycetaceae"="#58462A",
  "Poaceae"="#425024",
  "Rubiaceae"="#2C571F",
  "Mamiellaceae"="#1A5C1F",
  "Fabaceae"="#146020",
  "Pelagomonadales"="#1E632C",
  "Chaetomiaceae"="#2E653B",
  "Aspergillaceae"="#41664A",
  "Sporidiobolaceae"="#546459",
  "Noelaerhabdaceae"="#675F67",
  "Chlorellaceae"="#795877",
  "Rosales"="#8A4F85",
  "Herpotrichiellaceae"="#9A4591",
  "Salicaceae"="#A73A9C",
  "Culicidae"="#B22CA4",
  "Muridae"="#BB1BAB",
  "Clavicipitaceae"="#C300B0",
  "Euphorbiaceae"="#C312A0",
  "Solanaceae"="#C2278E",
  "Acanthamoebidae"="#BF3D7C",
  "Saprolegniaceae"="#BA526B",
  "Ustilaginaceae"="#B2675B",
  "Brassicaceae"="#A87A4D",
  "Formicidae"="#9C8B42",
  "Hyalellidae"="#8D993B",
  "Vitaceae"="#7CA437",
  "Oleaceae"="#68AC36",
  "Gadidae"="#51B137",
  "Rosaceae"="#39B33B",
  "Fagaceae"="#20B343",
  "Plectosphaerellaceae"="#00B04E",
  "Rhamnaceae"="#00AB5D",
  "Bathycoccaceae"="#00A36E",
  "Cucurbitaceae"="#009A80",
  "Convolvulaceae"="#008F92",
  "Cannabaceae"="#0083A4",
  "Rutaceae"="#0075B4"
)
