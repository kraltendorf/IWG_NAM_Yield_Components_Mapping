# Project: IWG_NAM_Yield_Components_Mapping
# Script 2 - Narrowing GWAS Peaks
# Author: Kayla R. Altendorf
# Date: 8/21/21

# load required packages
library("sommer")
library("vcfR")
library("dplyr")
library("ggplot2")
library("tibble")
library("rrBLUP")
library("tidyr")

# location of github directory
dir <- "/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Yield Components Mapping/Scripts for Github/"

# script we're on
folder <- c("/GWAS")

# create a vector of traits
traits <- c("flag_leaf_area", "floret_site_utilization", "florets_per_spikelet", "height", "reproductive_tiller_ct", 
            "seed_area", "seeds_per_spike", "spikelet_density", "spikelets_per_spike", "stem_diameter", 
            "thousand_grain_weight", "yield_per_spike")

year <- c("2017", "2018", "2017", "2018")
loc <- c("STP", "STP", "TLI", "TLI")
env <- c("stp17", "stp18", "tli17", "tli18")

# set the base pair threshold - this was identified in the LD analysis
# LD decays below 0.2 at a median rate of 21 cM. Average distance bp distance per cM across chromosomes x 21 = 
bp_threshold <-  63093157


#### Step 1: Read in and Format Genotypic Data ####
# this is the imputed dataset which was used in GWAS
vcf <- read.vcfR(paste(dir, folder, "/data/NAM_GATK_imputed.vcf", sep = ""), convertNA = TRUE)

# extract genotypes and id_frame
genotypes <- extract.gt(vcf, convertNA = FALSE)
id_frame <- as.data.frame(vcf@fix[,1:5]) %>% mutate(ID = paste(CHROM, "_", POS, sep = ""))

# update sample names from their ID in variant calling (e.g. flowcell, lane, barcode) to their sample names
# read in key
key <- read.table(paste(dir, folder, "/data/new_key.txt", sep = ""), header = T) %>% 
  mutate(GATK_Sample = paste(Flowcell, "_", Lane, "_", Barcode_ID, "_", sep = "")) %>%
  dplyr::select(GATK_Sample, Sample)

# prepare genotypes
genotypes <- as.data.frame(t(genotypes)) %>% 
  rownames_to_column(var = "GATK_Sample") %>%
  left_join(., key, by = "GATK_Sample") %>%
  arrange(Sample) %>%
  column_to_rownames("Sample") %>%
  dplyr::select(-GATK_Sample) %>%
  t()

# change all | to / to remove any phasing information that might be present
genotypes[genotypes=="0|1"] <- "0/1"
genotypes[genotypes=="1|0"] <- "0/1"
genotypes[genotypes=="1|1"] <- "1/1"
genotypes[genotypes=="0|0"] <- "0/0"

# convert to marker matrix
genotypes_matrix <- genotypes # create a new dataframe for the output

for (i in 1:nrow(genotypes)) {
  gt <- unlist(genotypes[i,])
  gt1 <- gt
  gt1[gt == "0/0"] <- 1
  gt1[gt == "0/1"] <- 0
  gt1[gt == "1/1"] <- -1
  gt1[gt == "./."] <- NA
  genotypes_matrix[i,] <- gt1 
}

genotype <- cbind(id_frame, genotypes_matrix)  # add the id_frame back
rownames(genotype) <- NULL # remove rownames


#### Step 2: Read in and Format Phenotypic Data ####
# read in backbone
backbone <- read.csv(paste(dir, "Phenotypic Data Analysis/Data/", "backbone.csv", sep = ""), header = T) %>% 
  dplyr::select(famID, parent, plantID, plantID3, longID) %>% 
  mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% 
  dplyr::select(-plantID3, -famID) %>% 
  distinct()


# read in emmeans files for flag_leaf_area
files <- list.files(paste(dir, "/Phenotypic Data Analysis/data/", traits[1], sep = ""), pattern = "emmeans_genet", full.names = TRUE)
flag_leaf_area <- list()
for (j in 1:length(files)) {
  flag_leaf_area[[j]] <- read.table(files[j], header = T) %>% 
    mutate(year = year[j], loc = loc[j]) %>% 
    mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% # create a column for joining
    left_join(backbone, ., by = "famID_plantID3") %>%
    dplyr::select(longID, emmean) %>%
    filter(! is.na(longID), 
           ! is.na(emmean)) %>%
    arrange(longID)
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
           ! is.na(emmean)) %>%
    arrange(longID)
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
           ! is.na(emmean)) %>%
    arrange(longID)
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
           ! is.na(emmean)) %>%
    arrange(longID)
}

# read in emmeans files for reproductive_tiller_ct
files <- list.files(paste(dir, "/Phenotypic Data Analysis/data/", traits[5], sep = ""), pattern = "emmeans_genet_transformed", full.names = TRUE)
reproductive_tiller_ct <- list()
for (j in 1:length(files)) {
  reproductive_tiller_ct[[j]] <- read.table(files[j], header = T) %>% 
    mutate(year = year[j], loc = loc[j]) %>% 
    mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% # create a column for joining
    left_join(backbone, ., by = "famID_plantID3") %>%
    dplyr::select(longID, emmean) %>%
    filter(! is.na(longID), 
           ! is.na(emmean)) %>%
    arrange(longID)
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
           ! is.na(emmean)) %>%
    arrange(longID)
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
           ! is.na(emmean)) %>%
    arrange(longID)
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
           ! is.na(emmean)) %>%
    arrange(longID)
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
           ! is.na(emmean)) %>%
    arrange(longID)
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
           ! is.na(emmean)) %>%
    arrange(longID)
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
           ! is.na(emmean)) %>%
    arrange(longID)
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
           ! is.na(emmean)) %>%
    arrange(longID)
}


traits

#### Step 3: Read in GWAS Results from GAPIT ####
flag_leaf_area_files <- list.files(paste(dir, folder, "/output/GAPIT/", traits[1], sep = ""), pattern = "GWAS.Results.csv", full.names = TRUE)
floret_site_utilization_files <- list.files(paste(dir, folder, "/output/GAPIT/", traits[2], sep = ""), pattern = "GWAS.Results.csv", full.names = TRUE)
florets_per_spikelet_files <- list.files(paste(dir, folder, "/output/GAPIT/", traits[3], sep = ""), pattern = "GWAS.Results.csv", full.names = TRUE)
height_files <- list.files(paste(dir, folder, "/output/GAPIT/", traits[4], sep = ""), pattern = "GWAS.Results.csv", full.names = TRUE)
reproductive_tiller_ct_files <- list.files(paste(dir, folder, "/output/GAPIT/", traits[5], sep = ""), pattern = "GWAS.Results.csv", full.names = TRUE)
seed_area_files <- list.files(paste(dir, folder, "/output/GAPIT/", traits[6], sep = ""), pattern = "GWAS.Results.csv", full.names = TRUE)
seeds_per_spike_files <- list.files(paste(dir, folder, "/output/GAPIT/", traits[7], sep = ""), pattern = "GWAS.Results.csv", full.names = TRUE)
spikelet_density_files <- list.files(paste(dir, folder, "/output/GAPIT/", traits[8], sep = ""), pattern = "GWAS.Results.csv", full.names = TRUE)
spikelets_per_spike_files <- list.files(paste(dir, folder, "/output/GAPIT/", traits[9], sep = ""), pattern = "GWAS.Results.csv", full.names = TRUE)
stem_diameter_files <- list.files(paste(dir, folder, "/output/GAPIT/", traits[10], sep = ""), pattern = "GWAS.Results.csv", full.names = TRUE)
thousand_grain_weight_files <- list.files(paste(dir, folder, "/output/GAPIT/", traits[11], sep = ""), pattern = "GWAS.Results.csv", full.names = TRUE)
yield_per_spike_files <- list.files(paste(dir, folder, "/output/GAPIT/", traits[12], sep = ""), pattern = "GWAS.Results.csv", full.names = TRUE)


# now iterate through each trait
# filter significant SNPs with -log10(p) 3.6 or p.value < 0.00025

# flag_leaf_area
flag_leaf_area_gwas <- list()
for (i in 1:length(flag_leaf_area_files)) {
  flag_leaf_area_gwas[[i]] <- read.csv(flag_leaf_area_files[i], header = T) %>% filter(P.value < 0.00025)
}

# floret_site_utlization
floret_site_utilization_gwas <- list()
for (i in 1:length(floret_site_utilization_files)) {
  floret_site_utilization_gwas[[i]] <- read.csv(floret_site_utilization_files[i], header = T) %>% filter(P.value < 0.00025)
}

# florets_per_spikelet
florets_per_spikelet_gwas <- list()
for (i in 1:length(florets_per_spikelet_files)) {
  florets_per_spikelet_gwas[[i]] <- read.csv(florets_per_spikelet_files[i], header = T) %>% filter(P.value < 0.00025)
}

# height
height_gwas <- list()
for (i in 1:length(height_files)) {
  height_gwas[[i]] <- read.csv(height_files[i], header = T) %>% filter(P.value < 0.00025)
}

# reproductive_tiller_ct
reproductive_tiller_ct_gwas <- list()
for (i in 1:length(reproductive_tiller_ct_files)) {
  reproductive_tiller_ct_gwas[[i]] <- read.csv(reproductive_tiller_ct_files[i], header = T) %>% filter(P.value < 0.00025)
}

# seed_area
seed_area_gwas <- list()
for (i in 1:length(seed_area_files)) {
  seed_area_gwas[[i]] <- read.csv(seed_area_files[i], header = T) %>% filter(P.value < 0.00025)
}

# seeds_per_spike
seeds_per_spike_gwas <- list()
for (i in 1:length(seeds_per_spike_files)) {
  seeds_per_spike_gwas[[i]] <- read.csv(seeds_per_spike_files[i], header = T) %>% filter(P.value < 0.00025)
}

# spikelet_density
spikelet_density_gwas <- list()
for (i in 1:length(spikelet_density_files)) {
  spikelet_density_gwas[[i]] <- read.csv(spikelet_density_files[i], header = T) %>% filter(P.value < 0.00025)
}

# spikelets_per_spike
spikelets_per_spike_gwas <- list()
for (i in 1:length(spikelets_per_spike_files)) {
  spikelets_per_spike_gwas[[i]] <- read.csv(spikelets_per_spike_files[i], header = T) %>% filter(P.value < 0.00025)
}

# stem_diameter
stem_diameter_gwas <- list()
for (i in 1:length(stem_diameter_files)) {
  stem_diameter_gwas[[i]] <- read.csv(stem_diameter_files[i], header = T) %>% filter(P.value < 0.00025)
}

# thousand_grain_weight
thousand_grain_weight_gwas <- list()
for (i in 1:length(thousand_grain_weight_files)) {
  thousand_grain_weight_gwas[[i]] <- read.csv(thousand_grain_weight_files[i], header = T) %>% filter(P.value < 0.00025)
}

# yield_per_spike
yield_per_spike_gwas <- list()
for (i in 1:length(yield_per_spike_files)) {
  yield_per_spike_gwas[[i]] <- read.csv(yield_per_spike_files[i], header = T) %>% filter(P.value < 0.00025)
}


#### Step 3: Identify Priority SNPs ####
# identify SNPs that are commonly significant across traits and environments 
# sort by most significant pvalues
traits
gwas <- stem_diameter_gwas

priority_snps_pval <- rbind(do.call("rbind", gwas)) %>%
  group_by(SNP) %>% 
  summarise(p_val_sum = sum(P.value)) %>%
  arrange(p_val_sum)

# sort by most frequent
priority_snps_n <- rbind(do.call("rbind", gwas)) %>%
  group_by(SNP) %>% 
  tally() %>% 
  arrange(-n) %>%
  filter(n >= 2) 



# merge togther
priority_snps <- left_join(priority_snps_n, priority_snps_pval, by = "SNP") %>%
  mutate(chrom = substr(SNP, 4, 5), 
         pos = substr(SNP, 7, nchar(SNP))) %>%
  arrange(chrom, -n, p_val_sum) 


#### THIS IS FOR FLAG LEAF AREA, REPRODUCTIVE TILLER NUMBER, HEIGHT, AND YIELD PER SPIKE
## WHERE THERE ARE TONS OF PRIORITY SNPS ON THE SAME CHR.


##### HEIGHT
# chr 4
chr4 <- priority_snps %>% filter(chrom == "04") %>% 
  mutate(dist = abs(as.numeric(pos) - 173993463)) %>% # this is the one most signficant and in the most envs
  filter(dist < bp_threshold) %>%
  arrange(dist)

remove_priority <- chr4$SNP[2:nrow(chr4)]
priority_snps <- priority_snps %>% filter(! SNP %in% remove_priority)

# chr 9
chr9 <- priority_snps %>% filter(chrom == "09") %>% 
  mutate(dist = abs(as.numeric(pos) - 450570468)) %>% # this is the one most signficant and in the most envs
  filter(dist < bp_threshold) %>%
  arrange(dist)

remove_priority <- chr9$SNP[2:nrow(chr9)]
priority_snps <- priority_snps %>% filter(! SNP %in% remove_priority)


##### YIELD PER SPIKE
# chr 5
chr5 <- priority_snps %>% filter(chrom == "05") %>% 
  mutate(dist = abs(as.numeric(pos) - 426524999)) %>% # this is the one most signficant and in the most envs
  filter(dist < bp_threshold)

remove_priority <- chr9$SNP[2:nrow(chr9)]

priority_snps <- priority_snps %>% filter(! SNP %in% remove_priority) %>% select(-pos)
priority_snps

##### FLAG LEAF AREA
# chr 9
chr9 <- priority_snps %>% filter(chrom == "09") %>% 
  mutate(dist = abs(as.numeric(pos) - 450570468)) %>% # this is the one most signficant and in the most envs
  filter(dist < bp_threshold) %>%
  arrange(dist)

remove_priority <- chr9$SNP[2:nrow(chr9)]
priority_snps <- priority_snps %>% filter(! SNP %in% remove_priority)

####


#### Step 4: Run 'mmmer' in Sommer ####
# declare trait here - edit for each 
traits
phenotype <- stem_diameter
gwas <- stem_diameter_gwas

for (i in 1:4) { 

  # identify samples in common
  phen <- phenotype[[i]]$longID
  gen <- colnames(genotype)[-1:-5]
  common <- Reduce(intersect, list(phen, gen))
  
  # filter phenotype and genotype so they match 
  phenotype[[i]] <- phenotype[[i]] %>% filter(longID %in% common)
  genotype_common <- genotype[, colnames(genotype) %in% common]
  print(summary(phenotype[[i]]$longID == colnames(genotype_common))) 
  
  # create the genotype matrix
  genotype_matrix <- t(genotype_common)
  colnames(genotype_matrix) <- id_frame$ID
  genotype_matrix <- data.frame(apply(genotype_matrix, 2, function(x) as.numeric(as.character(x))))
  rownames(genotype_matrix) <- colnames(genotype_common)
  
  # extract markers to test
  snps <- gwas[[i]]$SNP
  
  # create phenotype_genotype dataset
  snps_to_test <- cbind(phenotype[[i]], genotype_matrix[,colnames(genotype_matrix) %in% snps])
  rownames(snps_to_test) <- NULL
  #colnames(snps_to_test)[3] <- c("Chr04_7259337")
  
  # calculate the A matrix
  A <- A.mat(as.matrix(genotype_matrix))
  
  # prepare model terms
  fixed <- as.formula(paste("emmean ~ ",  paste0(snps, collapse = " + "), sep = ""))
  random <- ~vs(longID, Gu=A)
  
  # run the model
  fit <- mmer(fixed = fixed, random = random, data = snps_to_test)
  
  # look at the anova 
  anova <- anova.mmer(fit)
  anova_snps <- anova %>% slice(3:nrow(anova) - 1) %>% # ignore first and last terms
    rownames_to_column("SNP") %>%
    mutate(chrom = substr(SNP, 4, 5)) %>%
    arrange(chrom)
  
  anova_snps_nest <- anova_snps %>% group_by(chrom) %>% nest()
  

  for (j in 1:length(anova_snps_nest$data)) {
    if (nrow(anova_snps_nest$data[[j]]) > 1) { # if there's more than one on a chromosome
      keep <- anova_snps_nest$data[[j]] %>% filter(`Pr(>F)` <= 0.00025) # keep the ones that are significant
      
      # but if none are significant to that point, then keep the most significant one
      if (nrow(keep) == 0) {
        keep <- anova_snps_nest$data[[j]] %>% arrange(`Pr(>F)`) %>% slice(1:1)
      }
      
      remove <- anova_snps_nest$data[[j]] %>% filter(`Pr(>F)` > 0.001) # and remove those that are not
      remove_priority <- remove %>% filter(SNP %in% priority_snps$SNP)
      
      if (length(remove_priority$SNP) > 1) {
        keep <- rbind(keep, remove_priority) # but if a priority SNP was not significant we'll keep it in for ease of reporting
      }
  
      if (nrow(keep) > 1) {
        # go through the keepers and test to see if any are priority SNPs
        p_snps <- keep %>%  
          mutate(chr = substr(SNP, 4, 5), 
                 pos = substr(SNP, 7, nchar(SNP)))
        
        # determine position of priority SNP in remaining snps
        p <- p_snps %>% filter(SNP %in% priority_snps$SNP)

        position <- as.numeric(p$pos)
        
        keep <- p_snps %>% mutate(bp_diff = abs(as.numeric(pos) - position)) %>% filter(! between(bp_diff, 1, bp_threshold))
      }
      anova_snps_nest$data[[j]] <- anova_snps_nest$data[[j]] %>% filter(SNP %in% keep$SNP) # only keep the final SNPs
    }
  }

  
  # unnest the final snps
  anova_snps_final <- anova_snps_nest %>% unnest(cols = c(data))
  
  # re-run sommer for the final variance explained 
  # prepare model terms
  fixed <- as.formula(paste("emmean ~ ",  paste0(anova_snps_final$SNP, collapse = " + "), sep = ""))
  random <- ~vs(longID, Gu=A)
  
  # run the model
  fit <- mmer(fixed = fixed, random = random, data = snps_to_test)
  
  # look at the anova 
  anova <- anova.mmer(fit)
  
  # extract out variance explained
  fixed_beta <- fit$Beta$Estimate
  X <- model.matrix(fixed, snps_to_test)
  var_fixed <- var(matrix(data = fixed_beta, nrow = nrow(X), ncol = ncol(X), byrow = T) * X)
  snp_var <- diag(var_fixed)[-1]
  var <- snp_var / (snp_var + sum(unlist(fit$sigma)))
  var_percent <- as.data.frame(var*100) %>% rownames_to_column("SNP")
  colnames(var_percent)[2] <- "percent_variation_explained"
  
  
  gwas[[i]] <- gwas[[i]] %>% 
    filter(SNP %in% anova_snps_final$SNP) %>% 
    left_join(., var_percent, by = "SNP")

  print(paste(env[[i]], " is complete at ", Sys.time(), sep = ""))
}



## edit trait here
for (i in 1:length(gwas)) {
  gwas[[i]] <- gwas[[i]] %>% mutate(trait = "stem_diameter",
                                    loc = loc[i], 
                                    year = year[i])
}



gwas_final <- do.call("rbind", gwas)
keep <- gwas_final %>% group_by(SNP) %>% tally() %>% arrange(-n) %>% filter(n >= 2) # must occur in more than 1 env

gwas_final <- gwas_final %>% filter(SNP %in% keep$SNP)

gwas_final %>% group_by(SNP) %>% tally()

# change file name here for trait
write.csv(gwas_final, paste(dir, folder, "/output/final_gwas_stem_diameter.csv", sep = ""), row.names = FALSE)

gwas_final %>% arrange(SNP)


bp_threshold
#### Step 5: Create a Table of Results ####
# read in the output

flag_leaf_area_output <- read.csv(paste(dir, folder, "/output/final_gwas_flag_leaf_area.csv", sep = ""), header = T)
floret_site_utilization_output <- read.csv(paste(dir, folder, "/output/final_gwas_floret_site_utilization.csv", sep = ""), header = T)
florets_per_spikelet_output <- read.csv(paste(dir, folder, "/output/final_gwas_florets_per_spikelet.csv", sep = ""), header = T)
height_output <- read.csv(paste(dir, folder, "/output/final_gwas_height.csv", sep = ""), header = T)
reproductive_tiller_ct_output <- read.csv(paste(dir, folder, "/output/final_gwas_reproductive_tiller_ct.csv", sep = ""), header = T)
seeds_per_spike_output <- read.csv(paste(dir, folder, "/output/final_gwas_seeds_per_spike.csv", sep = ""), header = T)
spikelet_density_output <- read.csv(paste(dir, folder, "/output/final_gwas_spikelet_density.csv", sep = ""), header = T)
spikelets_per_spike_output <- read.csv(paste(dir, folder, "/output/final_gwas_spikelets_per_spike.csv", sep = ""), header = T)
stem_diameter_output <- read.csv(paste(dir, folder, "/output/final_gwas_stem_diameter.csv", sep = ""), header = T)
thousand_grain_weight_output <- read.csv(paste(dir, folder, "/output/final_gwas_thousand_grain_weight.csv", sep = ""), header = T)
yield_per_spike_output <- read.csv(paste(dir, folder, "/output/final_gwas_yield_per_spike.csv", sep = ""), header = T)



colnames(flag_leaf_area_output)
colnames(floret_site_utilization_output)
colnames(florets_per_spikelet_output)
colnames(height_output)
colnames(reproductive_tiller_ct_output)
colnames(seeds_per_spike_output)
colnames(spikelet_density_output)
colnames(spikelets_per_spike_output)
colnames(stem_diameter_output)
colnames(thousand_grain_weight_output)
colnames(yield_per_spike_output)

gwas_all <- rbind(flag_leaf_area_output, floret_site_utilization_output, florets_per_spikelet_output, 
                  height_output, reproductive_tiller_ct_output, seeds_per_spike_output, spikelet_density_output, 
                  spikelets_per_spike_output, stem_diameter_output, thousand_grain_weight_output, yield_per_spike_output)
write.csv(gwas_all, "/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Yield Components Mapping/Scripts for Github/GWAS/output/final_gwas_all.csv", row.names = F)


#### summary
gwas_all %>% tally()
gwas_all %>% group_by(SNP) %>% tally()
gwas_all %>% group_by(SNP, trait) %>% tally() #%>% View()
gwas_all %>% group_by(trait) %>% tally() %>% arrange(-n)
gwas_all %>% group_by(Chromosome) %>% tally() %>% arrange(-n)
gwas_all %>% group_by(loc, year) %>% tally() %>% arrange(-n)
gwas_all %>% summarise(min = min(percent_variation_explained), mean = mean(percent_variation_explained), 
                       max = max(percent_variation_explained))
gwas_all %>% filter(percent_variation_explained > 6)
gwas_all %>% arrange(-percent_variation_explained)
gwas_all %>% filter(SNP == "Chr09_450570468")
gwas_all %>% filter(SNP == "Chr05_426524999")


# bring in reference and alternate alleles
colnames(id_frame)[3] <- "SNP"
gwas_all1 <- left_join(gwas_all, id_frame, by = "SNP") %>% 
  mutate(Alleles = paste(REF, ALT, sep = "/"))

# determine segregating famililes
# first create a vector of families
famID <- backbone %>% 
  mutate(famID = substr(famID_plantID3, 1, 5)) %>%
  dplyr::select(famID) %>%
  distinct()
famID <- as.vector(famID[,1])

# calculate segregation types for each family at each significant locus
output <- data.frame(matrix(NA, ncol=3, nrow=length(unique(gwas_all1$SNP))))
colnames(output) <- c("het", "hom_alt", "hom_ref") # create an output dataframe/list
output_list <- replicate(10, output, simplify = FALSE) # one for each family 

for (i in 1:length(famID)) {
  fam_genotype <- cbind(id_frame, genotype[, grepl(famID[i], names(genotype))]) %>% #iterating through families one at a time
    filter(SNP %in% gwas_all1$SNP) %>%
    arrange(SNP) # extract and arrange all the unique significant SNPs
  
  fam_id_frame <- fam_genotype %>% select(CHROM:ALT) # extract the ID to bind back later  
  fam_genotype <- fam_genotype %>% select(-CHROM:-ALT) # and the genotype data to analyze
  
  
  output_list[[i]]$het <- rowSums(fam_genotype == "0") / ncol(fam_genotype)
  output_list[[i]]$hom_alt <- rowSums(fam_genotype == "1") / ncol(fam_genotype)
  output_list[[i]]$hom_ref <- rowSums(fam_genotype == "-1") / ncol(fam_genotype)
}

# now assess which families are segregating 
# create a dataframe for the output
seg_output <- data.frame(fam_id_frame)

for (i in 1:length(output_list)) {
  for (j in 1:nrow(output_list[[i]])) {
    seg_pattern <- output_list[[i]][j,]
    
    if (length(seg_pattern[seg_pattern > 0.05]) > 2) { # if it's not segregating, one of the types would be 1, or 0.95 and 0.05, but if two are over 0.05, then we 
      output_list[[i]]$seg[j] <- substr(famID[i], 4, 5) # know it's segregating 
    }
    else (output_list[[i]]$seg[j] <- NA)
  }
  seg_output[,i] <- output_list[[i]]$seg
}

id_seg <- cbind(fam_id_frame, seg_output)
colnames(id_seg)[-1:-5] <- famID

# combine all the results
id_seg_frame <- id_seg %>% 
  unite("Segregating_Families", WGN07:WGN63, na.rm = TRUE, remove = FALSE, sep = ", ") %>% 
  select(-WGN07:-WGN63)


# there's one QTL on chromsome 15 that's associated with floret socre mean but it technically does not segregate except only within family 39
# het: 0.975409836 hom_alt: 0.01639344 hom_ref: 0.008196721
# every other family has 100% hom_alt.

# join it back with the gwas_all1 and arrange into a nice table
gwas_all2 <- left_join(gwas_all1, id_seg_frame, by = "SNP") %>% 
  select(trait, SNP, Alleles, loc, year, Segregating_Families, maf, P.value, effect, percent_variation_explained) %>%
  mutate(P.value = round(-log10(P.value), digits = 2),
         effect = round(effect, digits = 2),
         percent_variation_explained = as.numeric(round(percent_variation_explained, digits = 2)), 
         maf = round(maf, digits = 2)) %>%
  pivot_wider(names_from = c(loc, year), values_from = c("P.value", "effect", "percent_variation_explained")) %>%
  select(trait, SNP, Alleles, maf, Segregating_Families, P.value_STP_2017, effect_STP_2017, percent_variation_explained_STP_2017,
         P.value_STP_2018, effect_STP_2018, percent_variation_explained_STP_2018,
         P.value_TLI_2017, effect_TLI_2017, percent_variation_explained_TLI_2017,
         P.value_TLI_2018, effect_TLI_2018, percent_variation_explained_TLI_2018) %>%
  separate(SNP, into = c("Chromosome", "Position"), sep = "_") %>%
  arrange(trait, Chromosome, as.numeric(Position)) %>%
  mutate(SNP = paste(Chromosome, Position, sep = "_")) %>%
  select(-Chromosome, -Position) 

 
gwas_all2$trait = factor(gwas_all2$trait, levels = traits)
gwas_all2 <- gwas_all2[order(gwas_all2$trait), ]


View(gwas_all2)

gwas_all3 <- gwas_all2 %>% mutate(chr = substr(SNP, 4, 5)) %>%
  mutate(pos = substr(SNP, 6, nchar(SNP)))


write.csv(gwas_all3, paste(dir, folder, "/output/gwas_final_table.csv", sep = ""), row.names = F)

# remake figs, table, insert into doc, review results/discussion


View(gwas_all1)

# calculate some summary statistics for the results section
gwas_all1 %>% select(SNP) %>% group_by(SNP) %>% tally() # 92 total SNPs detected overall
gwas_all1 %>% group_by(trait) %>% tally() # range between 2 and 19 QTL detected per trait
gwas_all1 %>% select(SNP) %>% group_by(SNP) %>% tally() %>% arrange(-n) # 7 detected in 2 environments, 2 in three or more environment/trait combinations

View(gwas_all1 %>% group_by(trait, loc, year) %>% summarise(mean_pve = mean(percent_variation_explained), 
                                            min_pve = min(percent_variation_explained), 
                                            max_pve = max(percent_variation_explained), 
                                            sum_pve = sum(percent_variation_explained)) %>% arrange(-sum_pve)) # 1.8 and 1.7


gwas_all1 %>% group_by(trait) %>% summarise(mean_effect = (median(abs(effect)))) # 0.324, 0.0647
gwas_all1 %>% select(SNP, Chromosome) %>% distinct() %>% group_by(Chromosome) %>% tally() %>% arrange(-n) # most detected on 5, 2, 6, 21, 9 and 16
gwas_all1 %>% group_by(loc, year) %>% tally() %>% arrange(-n)

# looking at the QTL on Chr 17
gwas_all1 %>% filter(Chromosome == 17)

# iterate through them
for (i in 1:length(snps)) {
  locus <- genotype_data %>% filter(rs == snps[i])
  
  alleles[[i]] <- locus$alleles
  loc[[i]] <- locus$rs
  
  # extract one family at a time
  for (j in 1:length(famID)) {
    fam <- locus[, grepl(famID[j], colnames(locus))]
    output_list[[1]][i, j] <- rowSums(fam == "0") / ncol(fam)
    output_list[[2]][i, j] <- rowSums(fam == "1") / ncol(fam)
    output_list[[3]][i, j] <- rowSums(fam == "-1") / ncol(fam)
    output_list[[4]][i, j] <- rowSums(is.na(fam)) / ncol(fam)
  }
}


##### Step 6: these are some extra bits of code to use for calculating LD if necessary ####   
# how about LD
check <- keep$SNP
ld_test <- genotype_matrix[, names(genotype_matrix) %in% check]
cor_out <- cor(ld_test)^2

# change the diagonals to NA
diag(cor_out) <- 0

#### looking at boxplots of allele states for QTL ####
genotype[1:10, 1:10]

flag_leaf_area

# extract the SNP of interest
# Chr09_450570468
# Chr05_426524999

chr9 <- genotype %>% filter(ID == "Chr09_450570468") %>% 
  select(-CHROM, -POS, -ID, -REF, -ALT) %>%
  t(.) %>%
  as.data.frame(.) %>%
  rownames_to_column("longID") %>%
  rename(Chr09_450570468 = V1)

head(flag_leaf_area[[i]])

for (i in 1:length(reproductive_tiller_ct)) {
  reproductive_tiller_ct[[i]] <- left_join(chr9, reproductive_tiller_ct[[i]], by = "longID") %>% 
    mutate(loc = loc[i], 
           year = year[i])
  
}

flag_leaf_area_dat <- do.call("rbind", reproductive_tiller_ct) %>%
  mutate(famID = substr(longID, 1, 5))

ggplot(flag_leaf_area_dat, aes(x=Chr09_450570468, y=emmean)) + 
  geom_boxplot() + 
  facet_grid(loc ~ year ~ famID)

genotypes_numeric["Chr09_450570468", "WGN59P20_"] 
