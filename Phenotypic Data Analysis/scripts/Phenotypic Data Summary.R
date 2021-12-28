# Project: IWG_NAM_Yield_Components_Mapping
# Analysis - Phenotypic Data Summary
# Author: Kayla R. Altendorf
# Date: 4/29/2021

# load required packages
library("dplyr")
library("tidyr")
library("stringr")
library("tibble")
library("lme4")
library("lmerTest")
library("emmeans")
library("gtools")
library("Hmisc")

# location of the github directory
dir <- "/users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Yield Components Mapping/Scripts for Github/"

# folder we're on
folder <- "/Phenotypic Data Analysis"

# declare where you want the output to go
out_path <- paste(dir, folder, "/output", sep = "")
in_path <- paste(dir, folder, "/data", sep = "")


#### Step 1: Read in Emmeans ####
# get the id_frame from data.csv - this is needed to convert the emmeans identity (which was required to be balanced for proper DF calcs) to their NAM identity (not balanced)
dat <- read.csv(paste(out_path, "/data.csv", sep = ""), header = T)
id_frame <- dat %>% 
  dplyr::select(famID, plantID3, longID) %>% 
  mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% 
  distinct() %>%
  filter(! is.na(longID))

# which traits are we dealing with? 
dirs <- list.dirs(path = in_path, full.names = TRUE, recursive = FALSE)

# remove the trait correlations output file if it's already been created 
dirs1 <- dirs[grepl("seed_area", dirs)]
dirs2 <- dirs[grepl("seeds_per_spike", dirs)]

dirs <- c(dirs1, dirs2)

# set locations in order
env <- c("stp17", "stp18", "tli17", "tli18")

# create empty dataframe
emmeans_env <- data.frame(matrix(NA, nrow = 1295, ncol = 2))
emmeans_all <- list()

for (j in 1:length(env)) {
  for (i in 1:length(dirs)) {
    files <- list.files(dirs[i], full.names = TRUE)
    files <- files[grepl(env[j], files)]
    file <- files[grepl("genet", files)]
    
    if (length(file) > 1) { # if there's more than one file that matches the "env" and genet characteristics, go with the untransformed for the sake of reporting the phnenotypic data
      file <- file[! grepl("transformed", file)]
    }
    
    emmeans <- read.table(file, head = T)
    emmeans1 <- emmeans %>% mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% dplyr::select(8, 3)
    emmeans2 <- left_join(id_frame, emmeans1, by = "famID_plantID3")
    emmeans_env[,i] <- emmeans2[,5] 
    colnames(emmeans_env)[i] <- str_split(dirs[i], "/")[[1]][11] 
  }
  emmeans_all[[j]] <- emmeans_env
}


#### Step 2: Phenotypic Data Summary Table ####
# traits we're working with
all_traits <- c("seed_area", "seeds_per_spike")
trait_summary <- list()

for (i in 1:length(emmeans_all)){
  emmeans_all[[i]] <- emmeans_all[[i]] %>% dplyr::select(all_traits)
  trait_min <- as.data.frame(sapply(emmeans_all[[i]], min, na.rm = T))
  trait_max <- as.data.frame(sapply(emmeans_all[[i]], max, na.rm = T))
  trait_mean <- as.data.frame(sapply(emmeans_all[[i]], mean, na.rm = T))
  trait_sd <- as.data.frame(sapply(emmeans_all[[i]], function(x)sd(x, na.rm = T)))
  trait_sum <- cbind(trait_min, trait_max[,1], trait_mean[,1], trait_sd[,1])
  trait_sum1 <- rownames_to_column(trait_sum, "trait")
  trait_sum2 <- trait_sum1 %>% mutate_if(is.numeric, funs(round(., 2)))
  trait_summary[[i]] <- trait_sum2
}

# rbind all environments together and rename column headers
trait_summary <- cbind(trait_summary[[1]][,1], trait_summary[[1]][,-1], trait_summary[[2]][,-1], trait_summary[[3]][,-1], trait_summary[[4]][,-1])
colnames(trait_summary) <- c("trait", "min", "max", "mean", "sd", "min", "max", "mean", "sd", "min", "max", "mean", "sd", "min", "max", "mean", "sd")


# write output
write.table(trait_summary, paste(out_path, "/trait_summary.txt", sep = ""), quote = F, row.names = F, sep = "\t")


