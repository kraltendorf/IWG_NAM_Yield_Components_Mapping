# Project: IWG_NAM_Domestication_Traits
# Analysis - Compiling Within Location ANOVAs
# Author: Kayla R. Altendorf
# Date: 12/23/2020

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
dir <- "/users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Domestication Traits/Scripts for Github/"

# folder we're on
folder <- "/Phenotypic Data Analysis"

# declare where you want the output to go
out_path <- paste(dir, folder, "/output", sep = "")
in_path <- paste(dir, folder, "/data", sep = "")

#### Step 1: Read in ANOVAs ####

# which traits are we dealing with? 
traits <- c("rachis_breaks_mean", "floret_score_mean", "reproductive_tiller_ct", "threshability", "floret_site_utilization", "height", "thousand_grain_weight")
dirs <- list.dirs(path = out_path, full.names = TRUE, recursive = FALSE) # determine their directories

# remove the trait correlations output file if it's already been created 
dirs <- dirs[! grepl("trait_correlations", dirs)]
dirs <- dirs[grepl(paste(traits, collapse = "|"), dirs)]

# set locations in order
env <- c("stp17", "stp18", "tli17", "tli18")

# create empty dataframe
out_file <- paste(dir, folder, "/output/anovas.txt", sep = "")

# iterate through and write to output file
for (i in 1:length(traits)) {
  trait_dir <- dirs[grepl(traits[i], dirs)] # extract the directory associated with the trait
  
  # list the files in that dir
  files <- list.files(trait_dir, full.names = TRUE)
  anovas <- files[grepl("anova", files)] # get the anovas
  anova <- anovas[! grepl("anova.txt", anovas)] # remove the full location anova if present
  
  # write out the trait we're dealing with in the output file
  write.table(str_split(trait_dir, pattern = "/")[[1]][11], out_file, append = T, col.names = F, row.names = F, quote = F, sep = "\t")
  
  for (j in 1:length(anova)) { 
    an <- read.table(anova[j], skip = 1)
    write.table(env[j], out_file, append = T, col.names = F, row.names = F, quote = F, sep = "\t")
    cat("\n", file = out_file, append = T)
    write.table(an, out_file, append = T, col.names = F, row.names = F, quote = F, sep = "\t")
    cat("\n", file = out_file, append = T)
  }
}



