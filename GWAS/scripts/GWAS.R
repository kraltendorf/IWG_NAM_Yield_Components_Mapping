# Project: IWG_NAM_Yield_Components_Mapping
# Analysis - GWAS
# Author: Kayla R. Altendorf
# Date: 8/21/21

# load packages

library("sommer")
library("vcfR")
library("dplyr")
library("tibble")
library("multtest")
library("gplots")
library("LDheatmap")
library("genetics")
library("ape")
library("EMMREML")
library("compiler") #this library is already installed in R
library("scatterplot3d")
library("cowplot")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")
BiocManager::install(c("multtest"))

# download gapit
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

# location of the github directory
dir <- "/users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Yield Components Mapping/Scripts for Github/"
# folder we're on
folder <- "/GWAS"

#### Step 1: Load in Emmeans ####
# create a vector of traits

traits <- c("flag_leaf_area", "floret_site_utilization", "florets_per_spikelet", "height", "reproductive_tiller_ct", 
           "seed_area", "seeds_per_spike", "spikelet_density", "spikelets_per_spike", "stem_diameter", 
           "thousand_grain_weight", "yield_per_spike")
year <- c("2017", "2018", "2017", "2018")
loc <- c("STP", "STP", "TLI", "TLI")
env <- c("stp17", "stp18", "tli17", "tli18")


# read in backbone
backbone <- read.csv(paste(dir, "Phenotypic Data Analysis/Data/", "backbone.csv", sep = ""), header = T) %>% 
  dplyr::select(famID, parent, plantID, plantID3, longID) %>% 
  mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% 
  dplyr::select(-plantID3, -famID) %>% 
  distinct()

# read in emmeans files for flag leaf area
files <- list.files(paste(dir, "/Phenotypic Data Analysis/data/", traits[1], sep = ""), pattern = "emmeans_genet", full.names = TRUE)
flag_leaf_area <- list()
for (j in 1:length(files)) {
  flag_leaf_area[[j]] <- read.table(files[j], header = T) %>% 
    mutate(year = year[j], loc = loc[j]) %>% 
    mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% # create a column for joining
    left_join(backbone, ., by = "famID_plantID3") %>%
    dplyr::select(longID, emmean) %>%
    filter(! is.na(longID), 
           ! is.na(emmean))
}

# read in emmeans files for floret_site_utilization
files <- list.files(paste(dir, "/Phenotypic Data Analysis/data/", traits[2], sep = ""), pattern = "emmeans_genet", full.names = TRUE)
floret_site_utilization <- list()
for (j in 1:length(files)) {
  floret_site_utilization[[j]] <- read.table(files[j], header = T) %>% 
    mutate(year = year[j], loc = loc[j]) %>% 
    mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% # create a column for joining
    left_join(backbone, ., by = "famID_plantID3") %>%
    dplyr::select(longID, emmean) %>%
    filter(! is.na(longID), 
           ! is.na(emmean))
}

# read in emmeans files for florets_per_spikelet
files <- list.files(paste(dir, "/Phenotypic Data Analysis/data/", traits[3], sep = ""), pattern = "emmeans_genet", full.names = TRUE)
florets_per_spikelet <- list()
for (j in 1:length(files)) {
  florets_per_spikelet[[j]] <- read.table(files[j], header = T) %>% 
    mutate(year = year[j], loc = loc[j]) %>% 
    mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% # create a column for joining
    left_join(backbone, ., by = "famID_plantID3") %>%
    dplyr::select(longID, emmean) %>%
    filter(! is.na(longID), 
           ! is.na(emmean))
}

# read in emmeans files for height
files <- list.files(paste(dir, "/Phenotypic Data Analysis/data/", traits[4], sep = ""), pattern = "emmeans_genet", full.names = TRUE)
height <- list()
for (j in 1:length(files)) {
  height[[j]] <- read.table(files[j], header = T) %>% 
    mutate(year = year[j], loc = loc[j]) %>% 
    mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% # create a column for joining
    left_join(backbone, ., by = "famID_plantID3") %>%
    dplyr::select(longID, emmean) %>%
    filter(! is.na(longID), 
           ! is.na(emmean))
}

# read in emmeans files for reproductive_tiller_ct
files <- list.files(paste(dir, "/Phenotypic Data Analysis/data/", traits[5], sep = ""), pattern = "emmeans_genet_transformed_scale", full.names = TRUE)
reproductive_tiller_ct <- list()
for (j in 1:length(files)) {
  reproductive_tiller_ct[[j]] <- read.table(files[j], header = T) %>% 
    mutate(year = year[j], loc = loc[j]) %>% 
    mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% # create a column for joining
    left_join(backbone, ., by = "famID_plantID3") %>%
    dplyr::select(longID, emmean) %>%
    filter(! is.na(longID), 
           ! is.na(emmean))
}

# read in emmeans files for seed_area
files <- list.files(paste(dir, "/Phenotypic Data Analysis/data/", traits[6], sep = ""), pattern = "emmeans_genet", full.names = TRUE)
seed_area <- list()
for (j in 1:length(files)) {
  seed_area[[j]] <- read.table(files[j], header = T) %>% 
    mutate(year = year[j], loc = loc[j]) %>% 
    mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% # create a column for joining
    left_join(backbone, ., by = "famID_plantID3") %>%
    dplyr::select(longID, emmean) %>%
    filter(! is.na(longID), 
           ! is.na(emmean))
}

# read in emmeans files for seeds_per_spike
files <- list.files(paste(dir, "/Phenotypic Data Analysis/data/", traits[7], sep = ""), pattern = "emmeans_genet", full.names = TRUE)
seeds_per_spike <- list()
for (j in 1:length(files)) {
  seeds_per_spike[[j]] <- read.table(files[j], header = T) %>% 
    mutate(year = year[j], loc = loc[j]) %>% 
    mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% # create a column for joining
    left_join(backbone, ., by = "famID_plantID3") %>%
    dplyr::select(longID, emmean) %>%
    filter(! is.na(longID), 
           ! is.na(emmean))
}

# read in emmeans files for spikelet_density
files <- list.files(paste(dir, "/Phenotypic Data Analysis/data/", traits[8], sep = ""), pattern = "emmeans_genet", full.names = TRUE)
spikelet_density <- list()
for (j in 1:length(files)) {
  spikelet_density[[j]] <- read.table(files[j], header = T) %>% 
    mutate(year = year[j], loc = loc[j]) %>% 
    mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% # create a column for joining
    left_join(backbone, ., by = "famID_plantID3") %>%
    dplyr::select(longID, emmean) %>%
    filter(! is.na(longID), 
           ! is.na(emmean))
}

# read in emmeans files for spikelets_per_spike
files <- list.files(paste(dir, "/Phenotypic Data Analysis/data/", traits[9], sep = ""), pattern = "emmeans_genet", full.names = TRUE)
spikelets_per_spike <- list()
for (j in 1:length(files)) {
  spikelets_per_spike[[j]] <- read.table(files[j], header = T) %>% 
    mutate(year = year[j], loc = loc[j]) %>% 
    mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% # create a column for joining
    left_join(backbone, ., by = "famID_plantID3") %>%
    dplyr::select(longID, emmean) %>%
    filter(! is.na(longID), 
           ! is.na(emmean))
}

# read in emmeans files for stem_diameter
files <- list.files(paste(dir, "/Phenotypic Data Analysis/data/", traits[10], sep = ""), pattern = "emmeans_genet", full.names = TRUE)
stem_diameter <- list()
for (j in 1:length(files)) {
  stem_diameter[[j]] <- read.table(files[j], header = T) %>% 
    mutate(year = year[j], loc = loc[j]) %>% 
    mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% # create a column for joining
    left_join(backbone, ., by = "famID_plantID3") %>%
    dplyr::select(longID, emmean) %>%
    filter(! is.na(longID), 
           ! is.na(emmean))
}

# read in emmeans files for thousand_grain_weight
files <- list.files(paste(dir, "/Phenotypic Data Analysis/data/", traits[11], sep = ""), pattern = "emmeans_genet", full.names = TRUE)
thousand_grain_weight <- list()
for (j in 1:length(files)) {
  thousand_grain_weight[[j]] <- read.table(files[j], header = T) %>% 
    mutate(year = year[j], loc = loc[j]) %>% 
    mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% # create a column for joining
    left_join(backbone, ., by = "famID_plantID3") %>%
    dplyr::select(longID, emmean) %>%
    filter(! is.na(longID), 
           ! is.na(emmean))
}

# read in emmeans files for yield_per_spike
files <- list.files(paste(dir, "/Phenotypic Data Analysis/data/", traits[12], sep = ""), pattern = "emmeans_genet", full.names = TRUE)
yield_per_spike <- list()
for (j in 1:length(files)) {
  yield_per_spike[[j]] <- read.table(files[j], header = T) %>% 
    mutate(year = year[j], loc = loc[j]) %>% 
    mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% # create a column for joining
    left_join(backbone, ., by = "famID_plantID3") %>%
    dplyr::select(longID, emmean) %>%
    filter(! is.na(longID), 
           ! is.na(emmean))
}



##### Step 2: Prepare SNP Data in Numeric Format #####
# prepare imputed genotype file in numeric format
vcf <- read.vcfR(paste(dir, folder, "/data/NAM_GATK_imputed.vcf", sep = ""), convertNA = TRUE)

# extract genotypes and id_frame
genotypes <- extract.gt(vcf, convertNA = FALSE)
id_frame <- as.data.frame(vcf@fix[,1:5]) # extract a dataframe with SNP ids -- this will come in handy later
id_frame <- id_frame %>% mutate(ID = paste(CHROM, "_", POS, sep = ""))


# update sample names from their ID in variant calling (e.g. flowcell, lane, barcode)
# to their sample names
# read in key
key <- read.table(paste(dir, folder, "/data/new_key.txt", sep = ""), header = T) %>% 
  mutate(GATK_Sample = paste(Flowcell, "_", Lane, "_", Barcode_ID, "_", sep = "")) %>%
  dplyr::select(GATK_Sample, Sample)

# prepare genotypes
genotypes1 <- t(genotypes) # transpose
genotypes2 <- as.data.frame(genotypes1) %>% rownames_to_column(var = "GATK_Sample")

genotypes3 <- left_join(genotypes2, key, by = "GATK_Sample") %>%
  column_to_rownames("Sample") %>%
  dplyr::select(-GATK_Sample) %>%
  t()

# change all | to / to remove phasing information, if any
genotypes3[genotypes3=="0|1"] <- "0/1"
genotypes3[genotypes3=="1|0"] <- "0/1"
genotypes3[genotypes3=="1|1"] <- "1/1"
genotypes3[genotypes3=="0|0"] <- "0/0"

genotypes_numeric <- genotypes3

for (i in 1:nrow(genotypes3)) {
  gt <- unlist(genotypes3[i,])
  gt1 <- gt
  gt1[gt == "0/0"] <- 0
  gt1[gt == "0/1"] <- 1
  gt1[gt == "1/1"] <- 2
  gt1[gt == "./."] <- NA
  genotypes_numeric[i,] <- gt1 
}

genotypes_numeric <- as.data.frame(t(genotypes_numeric)) %>% rownames_to_column(var = "taxa")

# extract snp information as a separate file
mdp_SNP_information <- id_frame %>% 
  dplyr::select(ID, CHROM, POS) %>%
  dplyr::rename(Name = ID, 
                Chromosome = CHROM, 
                Position = POS) %>%
  mutate(Chromosome = as.numeric(substr(Chromosome, 4, 5)))

# write out result
write.table(genotypes_numeric, paste(dir, folder, "/output/GATK_NAM_snp_matrix_imputed.table.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)
write.table(mdp_SNP_information, paste(dir, folder, "/output/mdp_SNP_information.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)

# read it back in 
myGD <- read.table(paste(dir, folder, "/output/GATK_NAM_snp_matrix_imputed.table.txt", sep = ""), head = T) 
myGM <- read.table(paste(dir, folder, "/output/mdp_SNP_information.txt", sep = ""), head = T) 



#### Step 3: Run GAPIT ####
# edit here only
traits # to reference trait order
trait <- spikelets_per_spike
trait_name <- c("spikelets_per_spike")

# set directory
dir.create(paste(dir, folder, "/output/GAPIT/", trait_name, sep = ""))
setwd(paste(dir, folder, "/output/GAPIT/", trait_name, sep = ""))

# run

for (i in 1:length(trait)) {
  phenotype <- trait[[i]]
  
  # choose only samples that are in common
  P <- phenotype$longID
  G <- myGD$taxa
  common <- Reduce(intersect, list(G, P))
  
  # filter phenotype data
  P2 <- filter(phenotype, longID %in% common) %>% droplevels() %>% arrange(longID)
  
  myGD2 <- myGD %>% 
    filter(taxa %in% common) %>%
    arrange(taxa)
  
  # make sure sample names align
  print(summary(myGD2$taxa == P2$longID))
  
  # rename phenotype header
  colnames(P2)[1] <- "Taxa"
  
  # rename phenotype to correct trait
  colnames(P2)[2] <- paste(trait_name, env[i], sep = "_")
  
  # gapit command
  myGAPIT <- GAPIT(
    Y=P2,
    GD=myGD2,
    GM=myGM,
    SNP.MAF=0.005, # inserting a MAF here since imputation introduced a few low freq alleles post vcftools filtering
    PCA.total = 0)
}



