# Project: IWG_NAM_Yield_Components_Mapping
# Analysis - Preparing .qua Files for MapQTL
# Author: Kayla Altendorf 
# Date: 4/30/21

# load required packages
library("dplyr") 

# location of the github directory
dir <- "/users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Yield Components Mapping/Scripts for Github/"

# folder we're on
folder <- "/Preparing .qua Files for MapQTL"

# declare where you want the output to go
out_path <- paste(dir, folder, "/output", sep = "")
in_path <- paste(dir, folder, "/data", sep = "")

# prepare backbone
backbone <- read.csv(paste(dir, "Phenotypic Data Analysis/data/backbone.csv", sep = "/"))
backbone1 <- backbone %>% dplyr::select(famID, parent, plantID, plantID3, longID) %>% 
  mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% 
  distinct() %>%
  dplyr::select(-plantID3, -famID)

# set vectors
famID <- unique(backbone$famID)
traits <- c("flag_leaf_area", "floret_site_utilization", "florets_per_spikelet", "height", "reproductive_tiller_ct", "seed_area", "seeds_per_spike", 
           "spikelet_density", "spikelets_per_spike", "stem_diameter", "thousand_grain_weight", "yield_per_spike")
env_names <- c("stp17", "stp18", "tli17", "tli18")

# create output
output <- list()

# import and prepare files 
for (i in 1:length(traits)) {
  file_names <- list.files(paste(dir, "Phenotypic Data Analysis/data/", traits[i], sep = ""), pattern = "emmeans_genet", full.names = T)
  files <- lapply(file_names, read.table, header = T)
  files <- files[! grepl("transformed", files)]
  
  # create empty list
  files_ed <- list()
  
  for (j in 1:length(files)) {
    file <- files[[j]]
    colnames(file)[3] <- "emmean"
    file <- file %>% mutate(famID_plantID3 = paste(famID, "_", plantID3, sep = "")) %>% 
      dplyr::select(famID_plantID3, emmean)
    file <- inner_join(file, backbone1, by = "famID_plantID3") %>%
      dplyr::select(longID, emmean) %>% filter(! is.na(longID))
    colnames(file)[2] <- paste(substr(traits[i], 1, 12), "_", env_names[j], sep = "") # determines length of trait name included
    files_ed[[j]] <- file
  }
  
  # now combine into one 
  file1 <- left_join(files_ed[[1]], files_ed[[2]], by = "longID")
  file2 <- left_join(file1, files_ed[[3]], by = "longID")
  file3 <- left_join(file2, files_ed[[4]], by = "longID")
  file4 <- arrange(file3, longID)
  all <- file4
  
  output[[i]] <- all 
}


# join together
phenotypes1 <- left_join(output[[1]], output[[2]], by = "longID")
phenotypes2 <- left_join(phenotypes1, output[[3]], by = "longID")
phenotypes3 <- left_join(phenotypes2, output[[4]], by = "longID")
phenotypes4 <- left_join(phenotypes3, output[[5]], by = "longID")
phenotypes5 <- left_join(phenotypes4, output[[6]], by = "longID")
phenotypes6 <- left_join(phenotypes5, output[[7]], by = "longID")
phenotypes7 <- left_join(phenotypes6, output[[8]], by = "longID")
phenotypes8 <- left_join(phenotypes7, output[[9]], by = "longID")
phenotypes9 <- left_join(phenotypes8, output[[10]], by = "longID")
phenotypes10 <- left_join(phenotypes9, output[[11]], by = "longID")
all_phenotypes <- left_join(phenotypes10, output[[12]], by = "longID")


# round to 3 digits
all_phenotypes_round <- all_phenotypes
for (i in 1:nrow(all_phenotypes)) {
  row <- unlist(all_phenotypes[i,-1])
  row_ed <- round(row, digits = 2)
  all_phenotypes_round[i,-1] <- row_ed
}

# create directory for qua files
dir_qua <- c(paste(out_path, "/qua_files/", sep = ""))
dir.create(dir_qua)
setwd(dir_qua)

# retrieve order of individuals from loc files to make sure they match
name_order <- list()

for (j in 1:length(famID)) {
  loc <- read.table(paste(dir, "Preparing .qua Files for MapQTL/data/", famID[j], "_chr1.txt", sep = ""), skip = 7, fill = TRUE)
  line <- which(grepl("individual", loc$V1)) # find where the individual names start
  names <- slice(loc, (line+1):nrow(loc)) %>% dplyr::select(V1) # extract them 
  name_order[[j]] <- names
}

# separate out each file into families, order them, and write them out
for (j in 1:length(famID)) {
  all_phenotypes_fam <- all_phenotypes_round %>% filter(grepl(famID[j], longID)) # just take those from one family at a time
  target <- name_order[[j]] # 
  colnames(target) <- c("longID")
  all_phenotypes_fam_target_order <- left_join(target, all_phenotypes_fam, by = "longID") # reorder using left_join
  print(summary(target$longID == all_phenotypes_fam_target_order$longID)) # make sure they're matching, print to console
  all_phenotypes_fam_target_order[is.na(all_phenotypes_fam_target_order)] <- c("*") # if it's NA, change to *
  write.table(paste("ntrt = ", ncol(all_phenotypes_fam_target_order), sep = "" ), # start writing out the files, with all the required components
              paste(dir_qua, famID[j], ".qua", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
  write.table(paste("nind = ", nrow(all_phenotypes_fam_target_order), sep = "" ),
              paste(dir_qua, famID[j], ".qua", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t", append = T)
  write.table(paste("miss = *", sep = "" ),
              paste(dir_qua, famID[j], ".qua", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t", append = T)
  cat("\n", file= paste(dir_qua, famID[j], ".qua", sep = ""), append=TRUE)
  colnames(all_phenotypes_fam_target_order)[1] <- "nr"
  trait_names <- colnames(all_phenotypes_fam_target_order)
  
  for (i in 1:length(trait_names)) {
    write.table(trait_names[i], paste(dir_qua, famID[j], ".qua", sep = ""), quote = F, row.names = F, col.names = F, sep = "\t", append = T)
  }
  cat("\n", file= paste(dir_qua, famID[j], ".qua", sep = ""), append=TRUE)
  write.table(all_phenotypes_fam_target_order, paste(dir_qua, famID[j], ".qua", sep = ""), quote = F, row.names = F, col.names = F, sep = "\t", append = T)
}

