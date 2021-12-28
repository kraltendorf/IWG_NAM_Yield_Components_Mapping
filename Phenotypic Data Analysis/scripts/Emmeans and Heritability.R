# Project: IWG_NAM_Yield_Components_Mapping
# Analysis - Phenotypic Data Analysis
# Author: Kayla R. Altendorf
# Date: 4/29/21

# Required Pacakges:
library("dplyr") 
library("tidyr")
library("lme4")
library("lmerTest")
library("emmeans")
library("stringr")
library("reshape")
library("plyr")
library("multcompView")
library("multcomp")


# Declare where you want the output to go 
path <- c("/users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Yield Components Mapping/Scripts for Github/")
folder <- c("Phenotypic Data Analysis")


#### Step 1: Load Phenotypic Data #### 
# this is the same format that is available from the intermediate wheatgrass database
# using this dataset requires a bit of fanagalling, but for posterity it's best to have
# one version of the data in use
dat <- read.table(paste(path, folder, "/data/NAM_Data.txt", sep = ""), header = T, sep = "\t")
unique(dat$trait_id)
# select the traits that will be used for this analysis
traits <- c("SDSPK", "SDAREA")
dat1 <- dat %>% filter(trait_id %in% traits)

# I prefer to work with full length names, so I'll sub them in here
my_names <- c("seeds_per_spike", "seed_area")

trait_names <- data.frame(trait_id = traits, trait_id_full = my_names)
dat2 <- left_join(dat1, trait_names, by = "trait_id")

#### Step 2: Format Phenotypic Data #### 
# get rid of the unnecessary columns
# note: since spikelets per spike was taken three times, the samples need separate names so they can be averaged
dat3 <- dat2 %>% 
  dplyr::rename(year = phenotype_year, # rename cols to match personal preference
                famID = family_name, 
                col = range) %>%
  mutate(loc = substr(experiment_id, 4, 6), # extract location
         trait_id_full = case_when(trait_id_full == "spikelets_per_spike" ~  paste(trait_id_full, sample_number, sep = ""), 
                                   trait_id_full != "spikelets_per_spike" ~ paste(trait_id_full)), # give each subsample of spikelets_per_spike a unique name
         parent = substr(germplasm_id, 6, 6)) %>% # extract parent (C for common, d for donor, p for parent)
  filter(parent != "P") %>% # filter to exclude parents, as this project deals only with progeny
  dplyr::select(famID, germplasm_id, loc, year, rep, trait_id_full, phenotype_value, plant_id) %>% 
  pivot_wider(names_from = trait_id_full, values_from = phenotype_value) %>% 
  dplyr::select(-plant_id) %>% # pivot to wide format
  mutate(loc = str_replace(loc, "SAL", "TLI")) # replace SAL (Salina) with TLI

# create column for merging
dat4 <- dat3 %>% 
  mutate(merge_col = paste(germplasm_id, loc, year, rep, sep = "_")) %>%
  dplyr::select(-famID, -germplasm_id, -loc, -year, -rep)

# the database includes *all* entries, including parents, plants that were identified later as selfs, and so on.
# futhermore, the population itself is not balanced (e.g. unequal numbers of individuals within families, 
# entries per location and year), which causes problems with the ANOVA
# to address this, we will load in a backbone dataset, which is balanced using NA values,
# and we'll format the data to match. 

# read in the backbone csv
backbone <- read.csv(paste(path, folder, "/data/backbone.csv", sep = ""))
# plantID3 is the balanced plantID
# example:
backbone %>% group_by(loc, year, rep, famID) %>% tally() # all have 133 entries

# left join with the data
dat6 <- left_join(backbone, dat4, by = "merge_col") %>% dplyr::select(-merge_col)
dat6 %>% group_by(loc, year, rep, famID) %>% tally() # make sure it's still balanced


# change all selfs to NA
for (i in 1:nrow(dat6)) {
  if (! is.na(dat6$self[i])) { # if self is not NA (e.g. if it is outcross or self)
    dat6[i, 13:14] <- NA # change all phneotype data columns to NA
  }
}

# write out the final dataset
#write.csv(dat6, paste(path, folder, "/data/data_including_parents.csv", sep = "/"), row.names = F)

write.csv(dat6, paste(path, folder, "/data/data.csv", sep = "/"), row.names = F)


# read it back in to convert everything to numerics
dat <- read.csv(paste(path, folder, "/data/data.csv", sep = ""))
dat$year <- as.factor(dat$year)
dat$rep <- as.factor(dat$rep)

#### Step 3: Analysis of Variance for Combined Analysis #### 
traits <- colnames(dat)[13:14] # grab the trait names from column headers

# make data as factor
dat$rep <- as.factor(dat$rep)
dat$year <- as.factor(dat$year)
dat$plantID3 <- as.factor(dat$plantID3)
dat$famID <- as.factor(dat$famID)

# run through each trait
traits

# seed_area
model <- lmer(seed_area ~ famID * loc * year + (1|loc:rep) + (1|loc:famID:plantID3), data = dat)
anova <- anova(model, type="II")
plot(model)
dir.create(paste(path, folder, "/data/seed_area", sep = "/"))
emmeans_loc_year <- emmeans(model, ~  year * loc)
write.table(anova, paste(path, folder, "/data/seed_area/seed_area_anova.txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
write.table(emmeans_loc_year, paste(path, folder, "/data/seed_area/seed_area_emmeans_loc_year.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

# seeds_per_spike
model <- lmer(seeds_per_spike ~ famID * loc * year + (1|loc:rep) + (1|loc:famID:plantID3), data = dat)
anova <- anova(model, type="II")
plot(model)
dir.create(paste(path, folder, "/data/seeds_per_spike", sep = "/"))
emmeans_loc_year <- emmeans(model, ~  year * loc)
write.table(anova, paste(path, folder, "/data/seeds_per_spike/seeds_per_spike_anova.txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
write.table(emmeans_loc_year, paste(path, folder, "/data/seed_area/seeds_per_spike_emmeans_loc_year.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")


#### Step 4: Analysis of Variance within Environments #### 

# filter data by location
stp17 <- filter(dat, loc == "STP" & year == "2017")
stp18 <- filter(dat, loc == "STP" & year == "2018")
tli17 <- filter(dat, loc == "TLI" & year == "2017")
tli18 <- filter(dat, loc == "TLI" & year == "2018")

loc_list <- list(stp17, stp18, tli17, tli18)
loc_names <- c("stp17", "stp18", "tli17", "tli18")

for (j in 1:length(loc_list)) {
  for (i in 1:length(traits)) {
    formula <- paste0(traits[i], " ~ famID/plantID3 + rep", sep = "")
    model <- lm(formula, data = loc_list[[j]])
    an <- as.data.frame(anova(model))
    emmeans_fam <- as.data.frame(emmeans(model, ~ famID), Letters = c(LETTERS))
    emmeans_genet <- as.data.frame(emmeans(model, ~ plantID3|famID))
    write.table(an, paste(path, folder, "/data/", traits[i], "/", traits[i], "_anova_", loc_names[j], ".txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
    write.table(emmeans_fam, paste(path, folder, "/data/", traits[i], "/", traits[i], "_emmeans_fam_", loc_names[j], ".txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
    write.table(emmeans_genet, paste(path, folder, "/data/", traits[i], "/", traits[i], "_emmeans_genet_", loc_names[j], ".txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
    print(paste(Sys.time(), "done with", traits[i], sep = " "))
  }
}

#### calculate heritability for these two traits ####

#### Step 3: Calculate Broad and Narrow Sense Heritability ####

# create a dataframe for the output
herit_df <- data.frame(env = c("STP", "TLI"), broad = NA, narrow = NA, trait = NA)
herit_df_list <- replicate(2, herit_df, simplify = FALSE)

# set vectors for iterating through years and locations
loc <- c("STP", "TLI")

# make important terms as factor
dat$year <- as.factor(dat$year)
dat$rep <- as.factor(dat$rep)

for (i in 1:length(traits)) {
  for (j in 1:length(loc)) {
    
    # extract data 
    dat_loc_year <- dat[dat$loc == loc[j],]
    formula <- paste0(traits[i], " ~ (1|famID:plantID3) + (1|year) + (1|famID:plantID3:year) + (1|rep) + (1|famID:plantID3:rep)", sep = "")
    # calculate broad sense on a genet mean basis
    model <- lmer(formula, data = dat_loc_year)
    
    output <- as.data.frame(VarCorr(model))
    vg <- output$vcov[3] # extracting appropriate variance components
    vgy <- output$vcov[1]
    vgr <- output$vcov[6]
    heritability_broad <- vg / ((vg) + (vgy / 2) + (vgr / 4))
    
    # clean out variables before next iteration
    model <- NA
    output <- NA
    vg <- NA
    vgy <- NA
    vgr <- NA
      
    # calculate narrow sense according to falconer
    # in this case, family and rep are random
    formula <- paste0(traits[i], "~ (1|famID) + (1|year) + (1|famID:year) + (1|rep) + (1|famID:rep)" , sep = "")
    model <- lmer(formula, data = dat_loc_year)
    output <- as.data.frame(VarCorr(model))
    vf <- output$vcov[3] # extracting appropriate variance components
    ve <- output$vcov[6]
    heritability_narrow <- (vf * 4) / ((vf * 4) + (ve))
      
    # clean out variables before next iteration
    model <- NA
    output <- NA
    vf <- NA
    ve <- NA

    # output result into heritability dataframe
    herit_df_list[[i]][j,2] <- heritability_broad
    herit_df_list[[i]][j,3] <- heritability_narrow
    herit_df_list[[i]][j,4] <- traits[i]
  }
}


# format the output into a nice table with means
# and sub in pub_names
herit_df <- do.call("rbind", herit_df_list)
herit_df_wide <- herit_df %>% 
  pivot_wider(values_from = c("broad", "narrow"), names_from = c("env")) %>%
  dplyr::select(trait, broad_STP, narrow_STP, broad_TLI, narrow_TLI) 

# calculate and append averages
herit_df_wide_avg <- herit_df_wide %>% 
  mutate_if(is.numeric, funs(round(., 2)))

# write out result
write.table(herit_df_wide_avg, paste(path, folder, "/output/heritabilities.txt", sep = ""),  quote = F, row.names = F, sep = "\t")


