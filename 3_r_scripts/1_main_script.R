###
# Written by: Conor Taff, Jenny Uehling, and Paige Becker
# Last updated: 3/13/2021
# Run under R Studio XX on XX

# This is the main analysis script for processing the COI data after running through
# AMPtk. Just a few files are saved from AMPtk to use here.

# The files are then used for analyses related to Paige Becker's honors thesis
# on tree swallow parental investment and diet quality.

################################################################################
# Load libraries ----

pacman::p_load("tidyverse", "phyloseq", "plyr", "vegan", "here", "ggpubr", "igraph", "data.table", "dplyr", "car")
# tidyverse & plyr for data wrangling
# phyloseq & vegan for community analyses and plotting
# here for file reference by relative paths
# ggpubr for plotting
# igraph for network plots

# Will need to install these libraries first if not yet installed
# To install phyloseq, must install BiocManager and then use the following command: BiocManager::install(c("phyloseq"))

################################################################################
# Load and wrangle tree swallow data ----

s_info <- read.delim(here("1_raw_data", "tres_sample_info.txt"))    # prepared outside of AMPtk

# Create column for unit_box_year to use to match up treatments with nestlings/nests
s_info$site_box_year <- paste(s_info$site, s_info$nest, s_info$year, sep="_")

# Import nest information, which has treatment categories
nest_info <- read.csv(here("1_raw_data", "Nest_Records_11.18.2020.csv"))
# Note: 4.7 from 2019 is not in the database. This error is going to be corrected.
# For now, I have put a row in this raw data file that contains the bare minimum information necessary.

# Create a column in nest_info that has the same information as site_box_year
nest_info$site_box_year <- paste(nest_info$Site, nest_info$Nest, nest_info$Exp_Year, sep="_")

# We're going to just delete the column that's already in "s_info" that has experiment and treatment.
# This is because it does not have experiment/treatment information for nestlings, only for adults.
# We're just going to re-import all the experiment info from the nest spreadsheet.
s_info <- subset(s_info, select=-c(experiment, treatment))

# Now, we'll extract just the columns we need from nest_info
nest_info <- subset(nest_info, select=c(Nest_Experiment, Nest_Treatment, site_box_year))

# Put nest experiment info in "s_info"
s_info <- left_join(s_info, nest_info, by = "site_box_year")

# Create a column that lists both the experiment AND the treatment
s_info$exp_treat <- paste(s_info$Nest_Experiment, s_info$Nest_Treatment, sep="_")

################################################################################
# Load and wrangle sequencing data ----

# the prefix used for AMPtk processing
# this is set in amptk and all files produced there have this prefix
# when adapting this to a different set of data, change prefix
amptk_prefix <- "trescoi"

# For now, we're going to need to process the sequencing runs from November and 
# December separately.

############################# November Run (11)

# Load the number of reads by taxa per sample table. Format for phyloseq.
otu_ab_11 <- read.delim(here("1_raw_data", paste0(amptk_prefix, "_2020_11.cluster.otu_table.txt")))
otu_ab_11$X.OTU.ID <- paste0(otu_ab_11$X.OTU.ID, "_11") # delineate as Nov seq run
rownames(otu_ab_11) <- otu_ab_11$X.OTU.ID   # give rownames from sample names
otu_ab_11 <- otu_ab_11[, 2:ncol(otu_ab_11)]    # remove the column of sample names

# Read the mapping table
# this 'mapping' table from AMPtk is mostly blank but has all the sample
# names so it can be joined to actual sample metadata. It's also possible
# to merge in the sample metadata within the AMPtk pipeline.
map_11 <- read.delim(here("1_raw_data", paste0(amptk_prefix, "_2020_11.mapping_file.txt")))

for(i in 1:nrow(map_11)){
  map_11$sampleID[i] <- strsplit(map_11$X.SampleID[i], "x")[[1]][2]
}

# write.table(map_11, "map_11.txt", sep = "\t") # if you want to save a copy of the mapping file

# Read the otu taxonomy table
otu_tax_11 <- read.delim(here("1_raw_data", paste0(amptk_prefix, "_2020_11.cluster.taxonomy.txt")))
otu_tax_11$X.OTUID <- paste0(otu_tax_11$X.OTUID, "_11") # delineate as Nov seq run
rownames(otu_tax_11) <- otu_tax_11$X.OTUID

# The taxonomy result from AMPtk is in one long string of text. This is splitting up the string
# and filling in a bunch of different columns. 
for(i in 1:nrow(otu_tax_11)){
  temp <- otu_tax_11$taxonomy[i]
  temp2 <- strsplit(temp, "\\|")[[1]][3]
  temp3 <- strsplit(temp2, ":")[[1]][2]
  otu_tax_11$search_hit[i] <- strsplit(temp, "\\|")[[1]][1]
  otu_tax_11$hit_score[i] <- strsplit(temp, "\\|")[[1]][2]
  otu_tax_11$database[i] <- strsplit(temp2, ":")[[1]][1]
  otu_tax_11$accession[i] <- strsplit(temp3, ";")[[1]][1]
  otu_tax_11$kingdom[i] <- strsplit(strsplit(temp2, "k:")[[1]][2], ",")[[1]][1]
  otu_tax_11$phylum[i] <- strsplit(strsplit(temp2, "p:")[[1]][2], ",")[[1]][1]
  otu_tax_11$class[i] <- strsplit(strsplit(temp2, "c:")[[1]][2], ",")[[1]][1]
  otu_tax_11$order[i] <- strsplit(strsplit(temp2, "o:")[[1]][2], ",")[[1]][1]
  otu_tax_11$family[i] <- strsplit(strsplit(temp2, "f:")[[1]][2], ",")[[1]][1]
  otu_tax_11$genus[i] <- strsplit(strsplit(temp2, "g:")[[1]][2], ",")[[1]][1]
  otu_tax_11$species[i] <- strsplit(strsplit(temp2, "s:")[[1]][2], ",")[[1]][1]		
}

# Replace database mismatches caused by matches that aren't from BOLD records
otu_tax_11$database <- gsub("None;k", "None", otu_tax_11$database)
otu_tax_11$accession <- gsub("Animalia,p", "None", otu_tax_11$accession)

# For phyloseq this has to be added as a matrix rather than data frame    
otu_tax_11 <- as.matrix(otu_tax_11)
# This is saving just the taxonomic ranks rather than the database info.
otu_tax_11_2 <- otu_tax_11[, 10:16]

# Add in the sample information to each sample
map_11_2 <- join(map_11, s_info, "sampleID", "left", "first")
rownames(map_11_2) <- map_11_2$X.SampleID

# Print out some summary information about this dataset
summary_table <- map_11_2 %>% dplyr::count(site, age.1, cap_num)

# This will identify how many samples there are in each category of site,
# age, and capture number.

# This will also identify samples that don't match the metadata file and write them as a separate output
# missing <- subset(map2, is.na(map2$band) == TRUE)
# write.table(missing, "missing_info.txt", sep = "\t") 

############################# December Run (12)

# Load the number of reads by taxa per sample table. Format for phyloseq.
otu_ab_12 <- read.delim(here("1_raw_data", paste0(amptk_prefix, "_2020_12.cluster.otu_table.txt")))
otu_ab_12$X.OTU.ID <- paste0(otu_ab_12$X.OTU.ID, "_12") # delineate as Dec seq run
rownames(otu_ab_12) <- otu_ab_12$X.OTU.ID   # give rownames from sample names
otu_ab_12 <- otu_ab_12[, 2:ncol(otu_ab_12)]    # remove the column of sample names

# Read the mapping table
# this 'mapping' table from AMPtk is mostly blank but has all the sample
# names so it can be joined to actual sample metadata. It's also possible
# to merge in the sample metadata within the AMPtk pipeline.
map_12 <- read.delim(here("1_raw_data", paste0(amptk_prefix, "_2020_12.mapping_file.txt")))

for(i in 1:nrow(map_12)){
  map_12$sampleID[i] <- strsplit(map_12$X.SampleID[i], "x")[[1]][2]
}

# write.table(map_12, "map_12.txt", sep = "\t") # if you want to save a copy of the mapping file

# Read the otu taxonomy table
otu_tax_12 <- read.delim(here("1_raw_data", paste0(amptk_prefix, "_2020_12.cluster.taxonomy.txt")))
otu_tax_12$X.OTUID <- paste0(otu_tax_12$X.OTUID, "_12") # delineate as Dec seq run
rownames(otu_tax_12) <- otu_tax_12$X.OTUID

# The taxonomy result from AMPtk is in one long string of text. This is splitting up the string
# and filling in a bunch of different columns. 
for(i in 1:nrow(otu_tax_12)){
  temp <- otu_tax_12$taxonomy[i]
  temp2 <- strsplit(temp, "\\|")[[1]][3]
  temp3 <- strsplit(temp2, ":")[[1]][2]
  otu_tax_12$search_hit[i] <- strsplit(temp, "\\|")[[1]][1]
  otu_tax_12$hit_score[i] <- strsplit(temp, "\\|")[[1]][2]
  otu_tax_12$database[i] <- strsplit(temp2, ":")[[1]][1]
  otu_tax_12$accession[i] <- strsplit(temp3, ";")[[1]][1]
  otu_tax_12$kingdom[i] <- strsplit(strsplit(temp2, "k:")[[1]][2], ",")[[1]][1]
  otu_tax_12$phylum[i] <- strsplit(strsplit(temp2, "p:")[[1]][2], ",")[[1]][1]
  otu_tax_12$class[i] <- strsplit(strsplit(temp2, "c:")[[1]][2], ",")[[1]][1]
  otu_tax_12$order[i] <- strsplit(strsplit(temp2, "o:")[[1]][2], ",")[[1]][1]
  otu_tax_12$family[i] <- strsplit(strsplit(temp2, "f:")[[1]][2], ",")[[1]][1]
  otu_tax_12$genus[i] <- strsplit(strsplit(temp2, "g:")[[1]][2], ",")[[1]][1]
  otu_tax_12$species[i] <- strsplit(strsplit(temp2, "s:")[[1]][2], ",")[[1]][1]		
}

# Replace database mismatches caused by matches that aren't from BOLD records
otu_tax_12$database <- gsub("None;k", "None", otu_tax_12$database)
otu_tax_12$accession <- gsub("Animalia,p", "None", otu_tax_12$accession)

# For phyloseq this has to be added as a matrix rather than data frame    
otu_tax_12 <- as.matrix(otu_tax_12)
# This is saving just the taxonomic ranks rather than the database info.
otu_tax_12_2 <- otu_tax_12[, 10:16]

# Add in the sample information to each sample
map_12_2 <- join(map_12, s_info, "sampleID", "left", "first")
rownames(map_12_2) <- map_12_2$X.SampleID

# Print out some summary information about this dataset
summary_table <- map_12_2 %>% dplyr::count(site, age.1, cap_num)

# This will identify how many samples there are in each category of site,
# age, and capture number.

# This will also identify samples that don't match the metadata file and write them as a separate output
# missing <- subset(map_12_2, is.na(map_12_2$band) == TRUE)
# write.table(missing, "missing_info.txt", sep = "\t")

############################# March Run (03)

# Load the number of reads by taxa per sample table. Format for phyloseq.
otu_ab_03 <- read.delim(here("1_raw_data", paste0(amptk_prefix, "_2021_03.cluster.otu_table.txt")))
otu_ab_03$X.OTU.ID <- paste0(otu_ab_03$X.OTU.ID, "_03") # delineate as March seq run
rownames(otu_ab_03) <- otu_ab_03$X.OTU.ID   # give rownames from sample names
otu_ab_03 <- otu_ab_03[, 2:ncol(otu_ab_03)]    # remove the column of sample names

# Read the mapping table
# this 'mapping' table from AMPtk is mostly blank but has all the sample
# names so it can be joined to actual sample metadata. It's also possible
# to merge in the sample metadata within the AMPtk pipeline.
map_03 <- read.delim(here("1_raw_data", paste0(amptk_prefix, "_2021_03.mapping_file.txt")))

for(i in 1:nrow(map_03)){
  map_03$sampleID[i] <- strsplit(map_03$X.SampleID[i], "x")[[1]][2]
}

# write.table(map_03, "map_03.txt", sep = "\t") # if you want to save a copy of the mapping file

# Read the otu taxonomy table
otu_tax_03 <- read.delim(here("1_raw_data", paste0(amptk_prefix, "_2021_03.cluster.taxonomy.txt")))
otu_tax_03$X.OTUID <- paste0(otu_tax_03$X.OTUID, "_03") # delineate as March seq run
rownames(otu_tax_03) <- otu_tax_03$X.OTUID

# The taxonomy result from AMPtk is in one long string of text. This is splitting up the string
# and filling in a bunch of different columns. 
for(i in 1:nrow(otu_tax_03)){
  temp <- otu_tax_03$taxonomy[i]
  temp2 <- strsplit(temp, "\\|")[[1]][3]
  temp3 <- strsplit(temp2, ":")[[1]][2]
  otu_tax_03$search_hit[i] <- strsplit(temp, "\\|")[[1]][1]
  otu_tax_03$hit_score[i] <- strsplit(temp, "\\|")[[1]][2]
  otu_tax_03$database[i] <- strsplit(temp2, ":")[[1]][1]
  otu_tax_03$accession[i] <- strsplit(temp3, ";")[[1]][1]
  otu_tax_03$kingdom[i] <- strsplit(strsplit(temp2, "k:")[[1]][2], ",")[[1]][1]
  otu_tax_03$phylum[i] <- strsplit(strsplit(temp2, "p:")[[1]][2], ",")[[1]][1]
  otu_tax_03$class[i] <- strsplit(strsplit(temp2, "c:")[[1]][2], ",")[[1]][1]
  otu_tax_03$order[i] <- strsplit(strsplit(temp2, "o:")[[1]][2], ",")[[1]][1]
  otu_tax_03$family[i] <- strsplit(strsplit(temp2, "f:")[[1]][2], ",")[[1]][1]
  otu_tax_03$genus[i] <- strsplit(strsplit(temp2, "g:")[[1]][2], ",")[[1]][1]
  otu_tax_03$species[i] <- strsplit(strsplit(temp2, "s:")[[1]][2], ",")[[1]][1]		
}

# Replace database mismatches caused by matches that aren't from BOLD records
otu_tax_03$database <- gsub("None;k", "None", otu_tax_03$database)
otu_tax_03$accession <- gsub("Animalia,p", "None", otu_tax_03$accession)

# For phyloseq this has to be added as a matrix rather than data frame    
otu_tax_03 <- as.matrix(otu_tax_03)
# This is saving just the taxonomic ranks rather than the database info.
otu_tax_03_2 <- otu_tax_03[, 10:16]

# Add in the sample information to each sample
map_03_2 <- join(map_03, s_info, "sampleID", "left", "first")
rownames(map_03_2) <- map_03_2$X.SampleID

# Print out some summary information about this dataset
summary_table <- map_03_2 %>% dplyr::count(site, age.1, cap_num)

# This will identify how many samples there are in each category of site,
# age, and capture number.

# This will also identify samples that don't match the metadata file and write them as a separate output
# missing <- subset(map2, is.na(map2$band) == TRUE)
# write.table(missing, "missing_info.txt", sep = "\t")

################################################################################
# Build initial phyloseq object ----

## Combine tax tables from both sequencing runs
otu_tax2 <- rbind(otu_tax_11_2, otu_tax_12_2, otu_tax_03_2)
TAX = tax_table(otu_tax2)

## Combine map files from both sequencing runs
map2 <- rbind(map_11_2, map_12_2, map_03_2)
SAM = sample_data(map2)

## Combine otu tables from both sequencing runs
# To do this, we'll have to create new columns that match
nov_names <- colnames(otu_ab_11)
dec_names <- colnames(otu_ab_12)
mar_names <- colnames(otu_ab_03)
otu_ab_11[dec_names] <- 0 # fill new columns with 0
otu_ab_11[mar_names] <- 0 # fill new columns with 0
otu_ab_12[nov_names] <- 0 # fill new columns with 0
otu_ab_12[mar_names] <- 0 # fill new columns with 0
otu_ab_03[nov_names] <- 0 # fill new columns with 0
otu_ab_03[dec_names] <- 0 # fill new columns with 0
otu_ab <- rbind(otu_ab_11, otu_ab_12, otu_ab_03)
otu_ab <- as.matrix(otu_ab)
OTU = otu_table(otu_ab, taxa_are_rows = TRUE)

coi_ps <- phyloseq(OTU, TAX, SAM)

# Subset just to arthropods
coi_ps2 <- subset_taxa(coi_ps, phylum == "Arthropoda")

# Extract the unique arthropod families found in the dataset
unique_families <- get_taxa_unique(coi_ps2, taxonomic.rank = "family")
# Save file with list of unique families to use to research aquatic vs.
# terrestrial families.
write.csv(unique_families, here("5_other_outputs/unique_families.csv"))
# This file was then taken out of R to research aquatic vs. terrestrial
# families, and will be re-imported later.

# Subset just to the data we want for this project
coi_ps2 <- subset_samples(coi_ps2, parental_investment_proj == "yes")

################################################################################
# Check sample effort ----

# Histogram of sample reads. This is counting a sum of arthropod reads for each sample.
coi_ps2_wild <- subset_samples(coi_ps2, site != "Gut_passage") # take out gut passage time birds
depth <- data.frame(as(sample_data(coi_ps2_wild), "data.frame"),
                    TotalReads = sample_sums(coi_ps2_wild), keep.rownames = TRUE)
p <- ggplot(depth, aes(log(TotalReads))) + geom_histogram(fill = "slateblue") + 
  ylab("Count of Samples") + xlab("log(Reads)") +
  theme_classic() + geom_vline(xintercept = log(150), linetype = "dashed", col = "coral3", size = 1) + 
  geom_text(x = log(150) - 0.2, y = 40, label = "150 Reads", angle = 90)
p2 <- p + facet_grid(~ age)     # same but splitting out adult/nestling/negative_control

# Save the histograms to file to be added to the markdown
ggsave(here("3_r_scripts/total_reads.png"), plot = p, width = 8, height = 4.5, device = "png")
ggsave(here("3_r_scripts/total_reads_split.png"), plot = p2, width = 8.2, height = 4, device = "png")

# Histogram of sample reads across captive nestling days
coi_ps2_captive <- subset_samples(coi_ps2, site == "Gut_passage") # extract just gut passage time birds
depth <- data.frame(as(sample_data(coi_ps2_captive), "data.frame"),
                    TotalReads = sample_sums(coi_ps2_captive), keep.rownames = TRUE)
p_cap <- ggplot(depth, aes(log(TotalReads))) + geom_histogram(fill = "slateblue") + 
  ylab("Count of Samples") + xlab("log(Reads)") +
  theme_classic() + geom_vline(xintercept = log(150), linetype = "dashed", col = "coral3", size = 1)
p2_cap <- p_cap + facet_wrap(~ age.1, ncol = 1) # same but splitting out age (by day)
# set order for ages for neatness in plots 
age_order = c("6", "7", "8", "9", "10", "11", "12")
p2_cap$data$age.1 <- as.character(p2_cap$data$age.1)
p2_cap$data$age.1 <- factor(p2_cap$data$age.1, levels=age_order)

# Save the histograms to file to be added to the markdown
ggsave(here("3_r_scripts/total_reads_cap.png"), plot = p_cap, width = 8, height = 4.5, device = "png")
ggsave(here("3_r_scripts/total_reads_split_cap.png"), plot = p2_cap, width = 4, height = 6, device = "png")

# Rarefy to even depth of 150
# running this would rarefy to an even depth moving forward
# coi_ps2 <- rarefy_even_depth(coi_ps2, sample.size = 150, rngseed = 92)

################################################################################
# Agglomerate taxa ----

# Depending on the analyses, we may want to agglomerate to different taxonomic ranks.
# Many of the sequence IDs do not go all the way to species, so in those cases analyses
# at the species level wouldn't include those sequences
coi_genus <- tax_glom(coi_ps2, taxrank = "genus")
coi_fam <- tax_glom(coi_ps2, taxrank = "family")
coi_ord <- tax_glom(coi_ps2, taxrank = "order")

# Merge family to life history
life_history <- read.csv(here("5_other_outputs", "unique_families_aquatic_terrestrial_IDs.csv"))    # prepared outside of R
otu_lh <- plyr::join(as.data.frame(tax_table(coi_fam)), life_history, "family", "left", "first")
coi_fam2 <- phyloseq(
  otu_table(coi_fam),
  tax_table(as.matrix(otu_lh)),
  sample_data(coi_fam)
)

glom_ps <- coi_fam2    # change here which agglomeration you want to use for plots below

################################################################################
# Filtering criteria ----

# Remove negative controls
glom_ps <- subset_samples(glom_ps, age != "neg_control")

# Remove OTUs with less than 5 reads in a sample aka 5-tons
# (could change to singletons, 50-tons, or whatever)
coi_ps2 <- prune_taxa(taxa_sums(glom_ps) > 5, glom_ps)

# Create a record of the seqeuncing depth of each sample before transforming to relative abundance
depth_postprune <- data.frame(as(sample_data(coi_ps2), "data.frame"),
                    TotalReads = sample_sums(coi_ps2), keep.rownames = TRUE)

# Transform to relative abundance
coi_ra <- transform_sample_counts(coi_ps2, function(x) x / sum(x))

## Filter out taxa with relative abundance values below some threshold
coi_ra2 <- filter_taxa(coi_ra, function(x) mean(x) > 1e-5, TRUE)

# Transform to presence absence
coi_pa <- transform_sample_counts(coi_ra2, function(x) ceiling(x))

# Limit to genera in 20% of samples
coi_20 <- prune_taxa(genefilter_sample(coi_pa, filterfun_sample(function(x) x > 0.1), A = 0.2 * nsamples(coi_pa)), coi_pa)

################################################################################
# Plot Patterns in Taxonomic Groups ----

# We'll start with some general exploration of the data. What types of families
# are we seeing, and how are they represented? Which families are most common?

# Note that we will use coi_ra2 for plots for relative abundance, and 
# coi_pa for plots for presence/absence.

######### Examine all families

# First, let's only use samples that were collected from birds in the wild for this figure
coi_pa_wild <- subset_samples(coi_pa, site != "Gut_passage")

# All genera
p <- plot_bar(coi_pa_wild, "family") + theme_classic() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  geom_hline(yintercept = 393 * 0.1, linetype = "dotted", col = "coral3") + 
  geom_hline(yintercept = 393 * 0.2, linetype = "dotted", col = "coral3") + 
  geom_hline(yintercept = 393 * 0.3, linetype = "dotted", col = "coral3")
ggsave(here("3_r_scripts/family_bar.png"), width = 10, height = 4.5, device = "png")

# This is pretty overwhelming and not super useful to look at.

######### Examine genera over 20% split by age

# Again, only use samples that were collected from birds in the wild
coi_20_wild <- subset_samples(coi_20, site != "Gut_passage")

p <- plot_bar(coi_20_wild, "family", fill="life_history") + theme_classic() + theme(axis.text.x = element_text(angle = 90)) + 
  facet_wrap(~ age, ncol = 1)

# Create a figure that has families color coded by life history
# phyloseq automatically draws black borders around every sample which makes it impossible to see the colors of 
# aquatic vs. terrestrial samples.
# We will have to redefine the function to get it to stop doing this (I found this code on stack overflow)
plot_bar_2 <-  function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, facet_grid = NULL, border_color = NA) 
{
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack",  color = border_color)
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

p <- plot_bar_2(coi_20, x="family", fill="life_history") + 
  xlab("Family") +
  ylab("Number of samples detected in") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, size = 12)) + 
  theme(axis.title = element_text(size = 14)) +
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12)) +
  theme(strip.text = element_text(size = 12)) +
  facet_wrap(~ age, ncol = 1) +
  scale_fill_manual("Life history", values = c("aquatic" = "skyblue1", "terrestrial" = "palegreen3", "both" = "tan4", 
                                               "unknown" = "gray"))
ggsave(here("3_r_scripts/common_families.png"), width = 10, height = 6, device = "png")

################################################################################
# Convert phyloseq object to data frame  ----

# At this point, we've done all of the filtering/organizing/plotting that we want with phyloseq.
# We will now convert the phyloseq object into a data frame so that we can more
# easily work with it.

plot_ra <- psmelt(coi_ra2) # psmelt makes a phyloseq object into a data frame

plot_pa <- psmelt(coi_pa) # psmelt makes a phyloseq object into a data frame

################################################################################
# Calculate percent aquatic & number of families in each sample ----

######## Percent aquatic using relative abundance
# For this first calculation, we will also save the number of families per sample

# Identify the unique samples
sample <- unique(plot_ra$sampleID)

# Create an empty data frame to store relative abundance and number of families information
aquatic <- data.frame(matrix(NA, ncol = 3, nrow = length(sample)))
aquatic <- data.frame()

for (i in 1:length(sample)){
  sam <- sample[i] # identify sample
  list <- plot_ra[plot_ra$sampleID == sam ,] # pull out all records from plot_ra that have that sample ID
  unique_fam <- unique(list$family[list$Abundance != "0"])
  num_fam <- length(unique_fam)
  list_aq <- list[list$life_history == "aquatic" ,] # of those records, pull out only those that have aquatic life histories
  rel_ab <- sum(list_aq$Abundance) # sum up all of the percentages of aquatic insects for that sample
  aquatic[i,1] <- sam # save the sample name
  aquatic[i,2] <- rel_ab # save the total aquatic relative abundance
  aquatic[i,3] <- num_fam # save the number of families in the sample
}

# Name columns
names(aquatic)[names(aquatic) == "V1"] <- "sampleID"
names(aquatic)[names(aquatic) == "V2"] <- "percent_aquatic_ra"
names(aquatic)[names(aquatic) == "V3"] <- "num_fam"

# Add in sample information
aquatic <- merge(aquatic, s_info, by = "sampleID", all.x = TRUE, all.y = FALSE)

######## Percent aquatic using presence/absence

# Identify the unique samples
sample <- unique(plot_pa$sampleID)

# Create an empty data frame to store presence/absence information
pa_aq <- data.frame(matrix(NA, ncol = 2, nrow = length(sample)))
pa_aq <- data.frame()

for (i in 1:length(sample)){
  sam <- sample[i] # identify sample
  list <- plot_pa[plot_pa$sampleID == sam ,] # pull out all records from plot_ra that have that sample ID
  list <- list[list$Abundance != 0 ,] # Only include OTUs in this calculation if they're actually present in the sample
  list_aq <- list[list$life_history == "aquatic" ,] # of those records, pull out only those that have aquatic life histories
  aq_count <- sum(list_aq$Abundance) # sum up all of the percentages of aquatic insects for that sample
  all_count <- (sum(list$Abundance))
  per_aq <- (aq_count / all_count)
  pa_aq[i,1] <- sam # save the sample name
  pa_aq[i,2] <- per_aq # save the total aquatic relative abundance
}

# Name columns
names(pa_aq)[names(pa_aq) == "V1"] <- "sampleID"
names(pa_aq)[names(pa_aq) == "V2"] <- "percent_aquatic_pa"

# Add in sample information
aquatic <- merge(pa_aq, aquatic, by = "sampleID")

# Check to make sure that the data are classified correctly
# str(aquatic) # This is important to check for modeling purposes, but I'm going to comment it out here so it doesn't spit out a long output in the Markdown document.

# Check on the relationship between percent aquatic when it is calculated 
# with relative abundance vs. presence/absence

p <- ggplot(aquatic, aes(x=percent_aquatic_ra, y = percent_aquatic_pa)) + geom_point() +
  geom_smooth(method=lm) + theme_classic() + xlab("Percent aquatic relative abundance") +
  ylab("Percent aquatic presence/absence")
ggsave(here("3_r_scripts/percent_aquatic_compare.png"), width = 5, height = 4, device = "png")

cor <- cor.test(aquatic$percent_aquatic_ra, aquatic$percent_aquatic_pa)
# There is a low correlation between percent aquatic when calculated with relative
# abundance vs. presence/absence.

################################################################################
# Add in the number of reads per sample (sequencing depth) ----

# Only keep the columns we need
depth_postprune <- depth_postprune[,c("sampleID", "TotalReads")]

# Merge with larger dataset
aquatic <- merge(aquatic, depth_postprune, by = "sampleID")

################################################################################
# Calculate percent aquatic pooled by nest, for use in direct comparisons of adult vs. nestling diet ----

# Subset data to exclude captive nestling samples
aquatic_wild <- aquatic[aquatic$site == "Unit_4" | aquatic$site == "Turkey_Hill" ,]

# Filter data just to adult and nestling captures during provisioning
# (i.e. take out first adult captures during provisioning)
aquatic_wild_ad <- aquatic_wild[aquatic_wild$cap_num != "1" & aquatic_wild$age == "Adult" ,]
aquatic_wild_nest <- aquatic_wild[aquatic_wild$age == "Nestling" ,]
aquatic_wild_prov <- rbind(aquatic_wild_ad, aquatic_wild_nest)

# Pull out nestlings to average across nest
aquatic_wild_prov_n <- aquatic_wild_prov[aquatic_wild_prov$age == "Nestling" ,]

# Pull out adults
aquatic_wild_prov_a <- aquatic_wild_prov[aquatic_wild_prov$age == "Adult" ,]

######## Relative abundance calculations

# Create an empty dataframe to store relative abundance information for each nest
nests <- unique(aquatic_wild_prov_n$site_box_year)
aquatic_wild_prov_npooled <- data.frame(matrix(NA, ncol = 2, nrow = length(nests)))
aquatic_wild_prov_npooled <- data.frame()

# Calculate average aquatic percentage with relative abundance for the nestlings in each nest
for (i in 1:length(nests)){
  nest <- nests[i] # identify nest
  list <- aquatic_wild_prov_n[aquatic_wild_prov_n$site_box_year == nest ,] # pull out all nestlings from that nest
  mean_aq <- mean(list$percent_aquatic_ra) # average % aquatic insects
  aquatic_wild_prov_npooled[i,1] <- nest # save the nest name
  aquatic_wild_prov_npooled[i,2] <- mean_aq # save the total aquatic relative abundance
}

# Name columns
names(aquatic_wild_prov_npooled)[names(aquatic_wild_prov_npooled) == "V1"] <- "site_box_year"
names(aquatic_wild_prov_npooled)[names(aquatic_wild_prov_npooled) == "V2"] <- "percent_aquatic_ra_nestlings_mean"

# Combine adults with pooled nestling percentages
aquatic_wild_prov_a <- merge(aquatic_wild_prov_a, aquatic_wild_prov_npooled, by = "site_box_year")

######## Presence/absence calculations

# Create an empty dataframe to store presence/absence information for each nest
nests <- unique(aquatic_wild_prov_n$site_box_year)
aquatic_wild_prov_npooled <- data.frame(matrix(NA, ncol = 2, nrow = length(nests)))
aquatic_wild_prov_npooled <- data.frame()

# Calculate average aquatic percentage with presence/absence for the nestlings in each nest
for (i in 1:length(nests)){
  nest <- nests[i] # identify nest
  list <- aquatic_wild_prov_n[aquatic_wild_prov_n$site_box_year == nest ,] # pull out all nestlings from that nest
  mean_aq <- mean(list$percent_aquatic_pa) # average % aquatic insects using presence/absence
  aquatic_wild_prov_npooled[i,1] <- nest # save the nest name
  aquatic_wild_prov_npooled[i,2] <- mean_aq # save the total aquatic percent
}

# Name columns
names(aquatic_wild_prov_npooled)[names(aquatic_wild_prov_npooled) == "V1"] <- "site_box_year"
names(aquatic_wild_prov_npooled)[names(aquatic_wild_prov_npooled) == "V2"] <- "percent_aquatic_pa_nestlings_mean"

# Combine adults with pooled nestling percentages
# This time, rename the file to something more 
aquatic_prov_comp <- merge(aquatic_wild_prov_a, aquatic_wild_prov_npooled, by = "site_box_year")

################################################################################
# Check to see if same birds sampled across multiple years  ----

# If there are birds that were sampled across years, then we'd need to probably
# randomly choose one year and throw out the other year for that bird.
# (We don't have enough repeat measures to make it worth using band as a random effect.)

table(aquatic$band[aquatic$age == "Adult" & aquatic$cap_num == "1"]) 
# None of the same birds are in this dataset for incubation captures in both 2019 and 2020

table(aquatic$band[aquatic$age == "Adult" & aquatic$cap_num != "1"])
# None of the same birds are in this dataset for provisioning captures in both 2019 and 2020

table(aquatic_prov_comp$band) # There are none of the same birds caught between years in this dataset

################################################################################
# Export files for analyses for Paige's thesis (performed in R Markdown)  ----

write.csv(aquatic, here("3_r_scripts", "aquatic.csv "))

write.csv(aquatic_prov_comp, here("3_r_scripts", "aquatic_prov_comp.csv"))
