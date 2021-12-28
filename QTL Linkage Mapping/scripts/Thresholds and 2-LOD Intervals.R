# IWG_Nested_Association_Mapping_Yield Components
# Script: Analysis of Results from Linkage Mapping in MapQTL 
# Kayla Altendorf  
# 8/21/21

### METHOD: 
# determine permutation test significance threshold for combined and within family
# first run interval mapping, select most significant QTL on a LG as co-factors
# output a file each step to sort in R 
# run rMQM, select most significant QTL on a LG as cofactors
# final, output data, identify 2-LOD dropoffs in this R script

### NOTE:
# settings for the permutation test and IM were 10 neighboring markers, 1,000 permutations, 
# algorithm was Mixture Model, mapping step size 1, 20 iterations
# the permutation test had to exclude family 07 (because it gives a singularity error, 
# and the first family cannot yield an error because there's a
# programming issue) and group 18_CP and 18_DP, because they have so much missing data

# load necessary packages
library("dplyr")
library("tidyr")
library("ggplot2")
library("tibble")
library("dplyr")
library("tidyr")
library("stringr")

# set paths
dir <- "/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Yield Components Mapping/Scripts for Github/"
folder <- "QTL Linkage Mapping"


#### Step 1: Summarizing Permutation Tests ####
# read in all permutation test files
files <- list.files(paste(dir, folder, "data", sep = "/"), full.names = T)

# set iteration vectors
famID <- c("WGN07", "WGN15", "WGN26", "WGN36", "WGN38", "WGN39", "WGN45", "WGN46", "WGN55", "WGN63")
trait <- c("flag_leaf_area", "floret_site_utilization", "florets_per_spikelet", "height", "reproductive_tiller_ct", 
           "seed_area", "seeds_per_spikelet", "spikelet_density", "spikelets_per_spike", "stem_diameter", 
           "thousand_grain_weight", "yield_per_spike")
env <- c("stp17", "stp18", "tli17", "tli18")

# make an empty dataframe nested in a list
output <- data.frame(famID = c("combined", famID), stp17 = NA, stp18 = NA, tli17 = NA, tli18 = NA, stringsAsFactors = FALSE)
LOD_thresh_list <- replicate(length(trait), output, simplify = FALSE) 

# make tables of LOD thresholds
for (i in 1:length(trait)) {
  for (e in 1:length(env)) {

    files <- list.files(paste(dir, folder, "data", trait[i], "permutation_test/", sep = "/"), pattern = env[e], full.names = T)
    
    if (length(files == 11)) {
      for (f in 1:length(files)) {
        permutation_test <- read.table(files[f], header = T, fill = T)


        LOD <- permutation_test %>% filter(Group == "GW") %>% 
          mutate(Rel.cum.count = as.numeric(as.character(Rel.cum.count))) %>%
          filter(Rel.cum.count > 0.90 & Rel.cum.count < 0.98) %>%
          dplyr::select(Interval)

          LOD_thresh_list[[i]][f, (e+1)] <- mean(LOD[,1])
      }
    }
  }
}



#### Step 2: Find 2-LOD Peaks for Families ####
# must manually iterate through each trait, then each environment all the way through writing out the results

# set trait
traits <- c("flag_leaf_area", "floret_site_utilization", "florets_per_spikelet", "height", "reproductive_tiller_ct", 
           "seed_area", "seeds_per_spikelet", "spikelet_density", "spikelets_per_spike", "stem_diameter", 
           "thousand_grain_weight", "yield_per_spike")
trait <- "spikelets_per_spike"

# set env
env <- "stp17" # stp17, stp18, tli17, tli18
# set family ID
famID <- c("WGN07", "WGN15", "WGN26", "WGN36", "WGN38", "WGN39", "WGN45", "WGN46", "WGN55", "WGN63")


# load in results
files <- list.files(paste(dir, folder, "/data/", trait, "/rmqm/", sep = ""), full.names = T)
files <- files[grepl(env, files)]
if (length(files) > 1) {file <- files[grepl("Combined", files)]
} else {
  file <- files
}

if (str_detect(file, "csv")) {
  rmqm <- read.csv(file, header = T, na.strings = c(" ", ""))
} else {
  rmqm <- read.delim(file, header = TRUE, sep = "\t", quote = "\"")
}


# prepare the LOD_thresh_list based on what trait is being investigated
# determine which order of the list the trait is
if (trait == "flag_leaf_area") {
  df <- LOD_thresh_list[[1]] }
if (trait == "floret_site_utilization") {
  df <- LOD_thresh_list[[2]] }
if (trait == "florets_per_spikelet") {
  df <- LOD_thresh_list[[3]] }
if (trait == "height") {
  df <- LOD_thresh_list[[4]] }
if (trait == "reproductive_tiller_ct") {
  df <- LOD_thresh_list[[5]] }
if (trait == "seed_area") {
  df <- LOD_thresh_list[[6]] }
if (trait == "seeds_per_spikelet") {
  df <- LOD_thresh_list[[7]] }
if (trait == "spikelet_density") {
  df <- LOD_thresh_list[[8]] }
if (trait == "spikelets_per_spike") {
  df <- LOD_thresh_list[[9]] }
if (trait == "stem_diameter") {
  df <- LOD_thresh_list[[10]] }
if (trait == "thousand_grain_weight") {
  df <- LOD_thresh_list[[11]] }
if (trait == "yield_per_spike") {
  df <- LOD_thresh_list[[12]] }

df_env <- df[,paste(env, sep = "")]

LOD_family <- df_env[-1] # remove combined for this instance
LOD_combined <- df_env[1] # keep only combined 

# select necessary columns
if (str_detect(file, "csv")) {
  rmqm_mapping <- rmqm %>% 
  dplyr::select(Group, Position, Locus, Combined.LOD, WGN07_dh.LOD, WGN15_dh.LOD, WGN26_dh.LOD, WGN36_dh.LOD, 
                WGN39_dh.LOD, WGN46_dh.LOD, WGN55_dh.LOD, WGN45_dh.LOD, WGN38_dh.LOD, WGN63_dh.LOD) %>%
  mutate(Group = as.factor(Group),
         Position = as.numeric(Position)) %>% 
  filter(Position >= 0) 
  } else {
    rmqm_mapping <- rmqm %>% 
      dplyr::select(Group, Position, Locus, Combined..LOD, WGN07..LOD, WGN15..LOD, WGN26..LOD, WGN36..LOD, 
                    WGN39..LOD, WGN46..LOD, WGN55..LOD, WGN45..LOD, WGN38..LOD, WGN63..LOD) %>%
      mutate(Group = as.factor(Group),
             Position = as.numeric(Position)) %>% 
      filter(Position >= 0)
  }

# regardless the file format, rename columns accordingly
colnames(rmqm_mapping) <- c("Group", "Position", "Locus", "Combined.LOD", "WGN07_dh.LOD", "WGN15_dh.LOD", "WGN26_dh.LOD", "WGN36_dh.LOD", 
                            "WGN39_dh.LOD", "WGN46_dh.LOD", "WGN55_dh.LOD", "WGN45_dh.LOD", "WGN38_dh.LOD", "WGN63_dh.LOD")


##########################################
# create output dataframes
groups <- data.frame(lg = rep(1:21, 2), map = c(rep("CP", 21), rep("DP", 21))) %>% mutate(group = paste(lg, map, sep = "_"))
result_list <- list()

# loop
for (j in 1:length(famID)) {
  
  # empty the result frame
  result <- groups %>%
    mutate(left = NA, 
           center = NA, 
           right = NA, 
           max_lod = NA, 
           locus = NA)
  
  # then iterate through groups
  for (i in 1:nrow(groups)) {
    
    group <- rmqm_mapping %>% 
      dplyr::select(Position, Group, Locus, paste(famID[j], "_dh.LOD", sep = "")) %>%
      filter(Group == groups$group[i])
    
    if (nrow(group) > 1) {
      if (max(group[,4], na.rm = T) >= LOD_family[j]) {
        peak <- group[which.max(group[,4]),]
        
        max_lod <- peak[,4]
        locus <- as.character(peak$Locus)
        center <- peak[,1] # record peak position
        LOD_drop <- peak[,4] - 2 # two LOD dropoff
        
        LOD_drop_less <- group %>% filter(.[,4] < LOD_drop & .[,4] > 0) # find positions that are less than the dropoff
        LOD_drop_id <- LOD_drop_less %>% 
          mutate(proximity_abs = abs(Position - center)) %>% 
          mutate(proximity = Position - center) %>% 
          arrange(proximity_abs) # find positions closest to "center"
        
        if (nrow(LOD_drop_id) == 0) { # if there are not any instances where LOD dropoff gets below the threshold
          right <- min(group$Position) # take the min and max positions for that group
          left <- max(group$Position) }
        
        else { # the first row should be one side of the peak, if negative, left, if positive right
          LOD_drop_id_pn <- LOD_drop_id %>% 
          mutate(pos_neg = ifelse(proximity < 0, "neg", ifelse(proximity > 0, "pos", NA))) %>%
          rownames_to_column("row_number")
          
          if (LOD_drop_id_pn$pos_neg[1] == "neg") {
            left <-  LOD_drop_id_pn$Position[1]
            right_id  <- LOD_drop_id_pn %>% filter(pos_neg == "pos") 
            if (length(right_id$row_number) == 0) {right <- center}
            else {right <- right_id[1,2]}
            }
          
          else if (LOD_drop_id_pn$proximity[1] > 0) {
            right <- LOD_drop_id_pn$Position[1]
            left_id  <- LOD_drop_id_pn %>% filter(pos_neg == "neg") 
            if (length(left_id$row_number) == 0) {left <- "0"}
            else {left <- left_id[1,2]}
          }
        }
        
        # export results
        result$left[i] <- left
        result$center[i] <- center
        result$right[i] <- right
        result$max_lod[i] <- max_lod 
        result$locus[i] <- locus
      }
    }
  }
  result_list[[j]] <- result
  result_list[[j]]$famID <- famID[j]
}

results <- do.call("rbind", result_list) %>% filter(! is.na(left))
results


trait


#### Step 3: Find 2-LOD Peaks for Combined ####

# empty out the result frame
result <- groups %>%
  mutate(left = NA, 
         center = NA, 
         right = NA, 
         locus = NA, 
         max_lod = NA)

for (i in 1:nrow(result)) {
  
  # extract groups 
  group  <- rmqm_mapping %>% 
    dplyr::select(Position, Group, Locus, Combined.LOD) %>%
    filter(Group == groups$group[i])

  if (nrow(group) > 1) {
    if (max(group[,4], na.rm = T) >= LOD_combined) {
      peak <- group[which.max(group[,4]),]
    
      locus <- as.character(peak$Locus)
      max_lod <- peak[,4]
      center <- peak[,1] # record peak position
      LOD_drop <- peak[,4] - 2 # two LOD dropoff
    
      LOD_drop_less <- group %>% filter(.[,4] < LOD_drop & .[,4] > 0) # find positions that are less than the dropoff
      LOD_drop_id <- LOD_drop_less %>% 
        mutate(proximity_abs = abs(Position - center)) %>% 
        mutate(proximity = Position - center) %>% 
        arrange(proximity_abs) # find positions closest to "center"
    
      # the first row should be one side of the peak, if negative, left, if positive right
      LOD_drop_id_pn <- LOD_drop_id %>% 
        mutate(pos_neg = ifelse(proximity < 0, "neg", ifelse(proximity > 0, "pos", NA))) %>%
        rownames_to_column("row_number")
    
      if (LOD_drop_id_pn$pos_neg[1] == "neg") {
        left <-  LOD_drop_id_pn$Position[1]
        right_id  <- LOD_drop_id_pn %>% filter(pos_neg == "pos") 
      if (length(right_id$row_number) == 0) {right <- center}
        else {right <- right_id[1,2]}
      }
    
      if (LOD_drop_id_pn$proximity[1] > 0) {
        right <- LOD_drop_id_pn$Position[1]
        left_id  <- LOD_drop_id_pn %>% filter(pos_neg == "neg") 
        if (length(left_id$row_number) == 0) {left <- "0"}
        else {left <- left_id[1,2]}
     }
    
      result$left[i] <- left
      result$center[i] <- center
      result$right[i] <- right
      result$locus[i] <- locus
      result$max_lod[i] <- max_lod
      result$famID <- "combined"
    }
  }
}

result <- result %>% filter(! is.na(left))

final_results <- rbind(result, results) %>% arrange(lg)
final_results






#### Step 4: Extract r2 for the Peak Marker #### 

# add a column for r2
final_results$r2 <- NA

# extract variance explained from the results
id_frame <- rmqm %>% dplyr::select(Group, Position, Locus)
var_explained <- rmqm[, grep("Expl", colnames(rmqm))] # extract the variance explained columns

for (i in 1:nrow(final_results)) {
  if (final_results$famID[i] == "combined") {final_results$r2[i] <- NA } # no r2 for the combined analysis
  else {
    family_name <- final_results$famID[i]
    var <- var_explained[, grep(family_name, colnames(var_explained))] # extract the variance explained columns
    id_var <- cbind(id_frame, var)
    
    r_frame <- id_var %>% 
      filter(Group == final_results$group[i], 
             Position == final_results$center[i])
    final_results$r2[i] <- r_frame$var
  }
}


#### Step 5: Extract the Allele Effects #### 
id_frame <- rmqm %>% dplyr::select(Group, Position, Locus)
mu_A <- rmqm[, grep("mu_A", colnames(rmqm))] # extract mu_A cols
mu_B <- rmqm[, grep("mu_B", colnames(rmqm))] # extract mu_B cols

allele_effects <- cbind(mu_A, mu_B)

# create output cols
final_results$mu <- NA
final_results$alpha <- NA
final_results$gamma <- NA

for (i in 1:nrow(final_results)) {
  if (final_results$famID[i] == "combined") {
    final_results$mu[i] <- NA 
    final_results$gamma[i] <- NA } # there are no allele effects calculated for the combined analysis 

  else {
    family_name <- final_results$famID[i] # extract family name
    al <- allele_effects[, grep(family_name, colnames(allele_effects))] # extract the variance explained columns
    id_al <- cbind(id_frame, al)
    
    al_frame <- id_al %>% 
      filter(Group == final_results$group[i], 
             Position == final_results$center[i])
    
    mu_A <- al_frame[,4]
    mu_B <- al_frame[,5]
    mu_AB <- (mu_A + mu_B) / 2 
    final_results$mu[i] <- mu_AB
    
    if (final_results$map[i] == "CP") { # the calculation varies depending on which parent the QTL was detected in
      final_results$alpha[i] <- (mu_B - mu_A) / 2 # if detected in the common parent, muB - muA and vice versa
      final_results$gamma[i] <- NA }
    
    if (final_results$map[i] == "DP") {
      final_results$gamma[i] <- (mu_A - mu_B) / 2
      final_results$alpha[i] <- NA }
  }
}





# END - format the results with the trait, environment, and year
final_results$trait <- trait
final_results$env <- substr(env, 1, 3)
final_results$year <- paste("20", substr(env, 4, 6), sep = "") 
final_results # view the results

# write out results
write.table(final_results, paste(dir, folder, "/output/", trait, "_", env, "_LOD_intervals.txt", sep = ""), 
            col.names = T, row.names = F, quote = F, sep = "\t")










### REPEAT STEPS 2-5 for each trait and environment combination ####

##############################
##############################
##############################
##############################
##############################





#### Step 6: Create a Table of Family Effects and Summary ####
files <- list.files(paste(dir, folder, "/output/", sep = ""), full.names = TRUE)
all_files <- lapply(files, read.table, header = T, sep="\t")
lod_intervals <- do.call("rbind", all_files)

# arrange the results by trait, year, lg, map and family
lod_intervals_final <- lod_intervals %>% arrange(trait, env, year, as.numeric(lg), map, famID) %>% select(-group)
head(lod_intervals_final)

lod_intervals_final$r2 <- as.numeric(lod_intervals_final$r2)
lod_intervals_final$alpha <- as.numeric(lod_intervals_final$alpha)
lod_intervals_final$gamma <- as.numeric(lod_intervals_final$gamma)



#### Step 7: Invert the Results on the Chromosomes that Require it #### 
# some of the linkage groups were inverted in the map making process - should have caught this earlier, but... here we are
# instead of re-doing the analysis, we can convert the positions here 
# these are the LGs that need inversion
invert <- c(01, 03, 05, 07, 08, 09, 10, 13, 14, 16, 17, 20)

# load in the map and calculate max positions 
map <- read.table("/users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Yield Components Mapping/Scripts for Github/QTL Linkage Mapping/data/no_group_inverted.map", header = F)
colnames(map) <- c("marker", "pos")
max_pos <- map %>% mutate(chrom = as.numeric(substr(marker, 4, 5))) %>% group_by(chrom) %>% summarise(max(pos)) %>% as.data.frame()

head(lod_intervals_final) # watch what happens to 9

for (i in 1:nrow(lod_intervals_final)) {
  if (lod_intervals_final$lg[i] %in% invert) {
    conversion_factor <- max_pos %>% filter(chrom  == lod_intervals_final$lg[i]) %>% as.data.frame()
    conversion_factor <- conversion_factor$`max(pos)`
    lod_intervals_final$left[i] <- abs(conversion_factor - as.numeric(lod_intervals_final$left[i]))
    lod_intervals_final$center[i] <- abs(conversion_factor - as.numeric(lod_intervals_final$center[i]))
    lod_intervals_final$right[i] <- abs(conversion_factor - as.numeric(lod_intervals_final$right[i]))
  }
}

head(lod_intervals_final)

## sometimes the left and rights are out of order now that I reversed the order of the map.
# fix that here 

for (i in 1:nrow(lod_intervals_final)) {
  interval <- c(lod_intervals_final$left[i], lod_intervals_final$center[i], lod_intervals_final$right[i])
  interval_sorted <- sort(interval)
  lod_intervals_final$left[i] <- interval_sorted[1]
  lod_intervals_final$center[i] <- interval_sorted[2]
  lod_intervals_final$right[i] <- interval_sorted[3]
}

#### write out all results
write.csv(lod_intervals_final, paste(dir, folder, "/output", "/lod_intervals_final.csv", sep = ""), row.names = F)


                                                                                                       