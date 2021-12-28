# Project: IWG_NAM_Yield_Components_Mapping
# Visualize Results
# Author: Kayla R. Altendorf
# Date: 01/30/21

# load packages
library("readr")
library("vcfR")
library("LinkageMapView")
library("dplyr")
library("stringr")
library("tibble")
library("stringr")

# location of github directory
dir <- "/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Yield Components Mapping/Scripts for Github/"
# script we're on
folder <- c("/Visualize Results")


#### Step 1: Prepare in SNP Markers ####
# read in vcf file
vcf <- read.vcfR(paste(dir, "GWAS/data/", "NAM_GATK_filtered_maf_selfs.vcf", sep = ""), convertNA = TRUE) # 8003 variants

# extract snp marker information
snp_markers <- as.data.frame(vcf@fix[,1:5]) %>% 
  mutate(ID = paste(CHROM, "_", POS, sep = "")) %>% 
  dplyr::select(ID) %>% 
  rename(marker_name = ID)

# read the inverted map formatted for R
map <- read.table("/users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Yield Components Mapping/Scripts for Github/QTL Linkage Mapping/data/no_group_inverted.map", header = F) %>% 
  mutate(group = substr(V1, 1, 5), 
         bp_pos = as.numeric(substr(V1, 7, nchar(as.character(V1)))), 
         chr = as.numeric(substr(V1, 4, 5)))

colnames(map) <- c("marker_name", "cm_pos", "group", "bp_pos", "chr") # make the column names uniform

# filter snp markers to remove those that are already on the genetic map
snp_markers <- snp_markers %>% filter(! marker_name %in% map$marker_name)

# create a dataframe of snp markers to rbind onto the map file
to_add <- snp_markers %>% 
  mutate(cm_pos = NA) %>%
  mutate(group = substr(marker_name, 1,5)) %>%
  mutate(bp_pos = as.numeric(as.character(substr(marker_name, 7, nchar(as.character(marker_name)))))) %>%
  mutate(chr = as.numeric(substr(group, 4, 5)))

map1 <- rbind(map, to_add) %>% mutate(display = NA) # display is NO here


#### Step 2: Read in and Format Known Candidate Gene Positions ####
# and make it the same format as the map 
View(test)
domestication_genes <- read.csv("/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Domestication Traits/Scripts for Github/Visualize Results/data/candidate_genes.csv", header = T) %>% 
  mutate(group = paste("Chr", str_pad(string = chr, width = 2, side = "left", pad = 0), sep = ""),
         cm_pos = NA,
         display = "gene") %>%
  filter(trait %in% c("Brittle rachis", "Free threshing")) %>%
  dplyr::select(marker_name, cm_pos, group, bp_pos, chr, display)

# rbind these domestication_genes to the map
map2 <- rbind(domestication_genes, map1)


#### Step 3: Read in Significant GWAS Markers ####
gwas_hits <- read.csv("/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Yield Components Mapping/Scripts for Github/GWAS/output/final_gwas_all.csv", header = T) %>%
  dplyr::select(Chromosome, Position, trait, loc, year)%>% 
  mutate(loc_year = paste(loc, substr(year, 3, 4), sep = "")) %>%
  dplyr::select(-loc, -year) %>%
  pivot_wider(names_from = loc_year, values_from = loc_year) %>%
  select(Chromosome, Position, trait, STP17, STP18, TLI17, TLI18) %>%
  unite("envs", STP17:TLI18, sep = ", ", na.rm = T)


# rename traits
trait_names <- data.frame(trait = unique(gwas_hits$trait), trait2 = c("FLA", "FSU", "FPS", "HGT", "RTC", "SESP", "SPKSPK", "STDIA", "TGW", "YPS"))

gwas_hits2 <- left_join(gwas_hits, trait_names, by = "trait") %>%
  mutate(marker_name = paste(trait2, " (", envs, ")", sep = "")) %>%
  dplyr::select(-envs) %>%
  rename(chr = Chromosome, 
         bp_pos = Position) %>%
  mutate(cm_pos = NA, 
         group = paste("Chr", str_pad(string = chr, width = 2, side = "left", pad = "0"), sep = ""),
         bp_pos = as.numeric(bp_pos), 
         display = "gwas_qtl") %>%
  dplyr::select(marker_name, cm_pos, group, bp_pos, chr, display, trait2) %>%
  arrange(trait2, group, bp_pos)

gwas_hits2 %>% group_by(marker_name, chr) %>% tally() %>% arrange(-n)


# write this out, and manually edit it to make sure there are no repeats on a chromosome, 
# if there are, add a "1_" before the marker name to make sure they show up independently on the figure
# rename with "_ed" and read back in 
write.csv(gwas_hits2, "/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Yield Components Mapping/Scripts for Github/Visualize Results/output/gwas_hits_formatted.csv", row.names = F)

gwas_hits_formatted <- read.csv("/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Yield Components Mapping/Scripts for Github/Visualize Results/output/gwas_hits_formatted.csv")

gwas_hits_formatted %>% group_by(marker_name, chr) %>% tally() %>% arrange(-n) #making sure there are no repeat names on the same chrom 

gwas_hits_formatted <- gwas_hits_formatted %>% dplyr::select(-trait2)

# the linkagemapview package is not ideal, and the chromosome maps come out differently
# sized depending on the content of the map. Here I am adding some random markers to Chromosome 1 to make
# it sized more appropriately
#rows_to_add <- gwas_hits_formatted %>% 
#  slice(10:12) %>%
#  mutate(group = "Chr01",
#         chr = "1")

# rbind to current map
map3 <- rbind(map1, gwas_hits_formatted) #rows_to_add) # add the gwas markers and extra markers to the map


#### Step 4: Normalize the Positions ####
# calculate minimums and maximums for each linkage group to normalize everything for the sake of visualization
map_sum <- map3 %>% 
  group_by(group) %>% 
  summarise(bp_min = min(bp_pos, na.rm = T), 
            bp_max = max(bp_pos, na.rm = T),
            cm_min = min(cm_pos, na.rm = T), 
            cm_max = max(cm_pos, na.rm = T))

# add in minimum and maximums to the map
map4 <- left_join(map3, map_sum, by = "group")

# create two copies of this map, one that is physical, the other that is genetic
genetic <- map4 %>% 
  mutate(group = paste("LG", substr(group, 4,5), sep = "")) %>% 
  mutate(map = "genetic") # %>% 
  #filter(! marker_name %in% domestication_genes$marker_name) # remove gene markers from genetic map

physical <- map4 %>% mutate(map = "physical")

# bind them together
all_data <- rbind(genetic, physical)

# make a column for normalized position
all_data$position <- NA

# calculate the normalized positions
for (i in 1:nrow(all_data)) {
  if (! is.na(all_data$map[i])) {
    if (all_data$map[i] == "genetic") {
      all_data$position[i] <- (1000 * (all_data$cm_pos[i] - all_data$cm_min[i]) / (all_data$cm_max[i] - all_data$cm_min[i]))
    }
    else if (all_data$map[i] == "physical") {
      all_data$position[i] <- (1000 * (as.numeric(all_data$bp_pos[i]) - as.numeric(all_data$bp_min[i])) / (as.numeric(all_data$bp_max[i]) - as.numeric(all_data$bp_min[i])))
    }
  }
}



#### Step 5: Create Reference Markers for Displaying cM and Mb the Map #####
# round position to two digits and create a new column for display, which indicates
# whether or not the marker name will be shown on the map
all_data <- all_data %>% mutate(position = round(position, 2))

genetic_beginning <- all_data %>%
  filter(map == "genetic") %>%
  filter(position == 0) %>%
  mutate(marker_name = paste(cm_pos, "cM"))

genetic_end <- all_data %>% 
  filter(map == "genetic") %>% 
  filter(position == 1000) %>%
  mutate(marker_name = paste(round(cm_pos, 0), "cM"))

physical_beginning <- all_data %>% 
  filter(map == "physical") %>% 
  filter(position == 0) %>%
  mutate(marker_name = paste(round(bp_pos/1000000, 0), "Mb"))

physical_end <- all_data %>% 
  filter(map == "physical") %>% 
  filter(position == 1000) %>%
  mutate(marker_name = paste(round(bp_pos/1000000, 0), "Mb"))

# bind them all together and call this series the 'ruler'
display_markers <- rbind(genetic_beginning, genetic_end, physical_end, physical_beginning) %>% 
  mutate(display = "ruler")

# remove any doubles, where the rounding causes several to match at the 1000 position
display_markers_distinct <- display_markers %>% 
  dplyr::select(-bp_pos, -cm_pos) %>% 
  distinct() %>%
  mutate(bp_pos = NA,
         cm_pos = str_replace_all(marker_name, c(" cM" = "", " Mb" = ""))) %>%
  dplyr::select(marker_name, cm_pos, group, bp_pos, chr:position)

# bind them back to all_data
all_data1 <- rbind(all_data, display_markers_distinct) %>% 
  arrange(group, position)

# remove the genetic markers from the map that don't have a cm_pos!
genetic <- all_data1 %>% filter(map == "genetic") %>% filter(! is.na(cm_pos))
physical <- all_data1 %>% filter(map == "physical")
all_data2 <- rbind(genetic, physical)


#### Step 6: Add in a Label Adjuster Dataframe ####
# in some cases the intervals for QTL will collapse on top
# of each other and so we create a blank white marker here to bump them out
# i had to make the map figures several times to identify and resolve and troubleshoot this problem 
# see: for brief tutorial on this: https://github.com/bio-services/LinkageMapView/issues/13
all_data2 %>% filter(group == "LG05") 

combined_qtl
adjuster <- data.frame(marker_name = c("w", "w", "w", "w"), 
                       cm_pos = c(400, 310, 63.489, 47.654), 
                       group = c("LG02", "LG08", "LG05", "LG05"), 
                       bp_pos = rep(NA, 4),
                       chr = c(2, 8, 5, 5), 
                       display = rep("adjuster", 4),
                       bp_min = c(109865, 1162650, 1153145, 1153145), 
                       bp_max = c(430937086, 431292282, 434141578, 434141578), 
                       cm_min = rep(0, 4), 
                       cm_max = c(143.98, 180.245, 146.36, 146.36),
                       map = rep("genetic", 4),
                       position = c(380, 310, 433.79, 325.59))

# add the adjuster markers to the map 
all_data2 <- rbind(all_data2, adjuster)

combined_qtldf





#### Step 7: Load in the Linkage Mapping Results ####
# read in final LOD intervals from "Thresholds and 2-LOD Intervals.R
lod_intervals <- read.csv(paste(dir, "/QTL Linkage Mapping/output/lod_intervals_final.csv", sep = ""))
                          
                          
# select only QTL that were detected in the combined analysis to feature in the figure
lod_intervals_ed <- lod_intervals %>% 
  filter(famID == "combined") %>%
  mutate(env = toupper(env), 
         chr = paste("LG", str_pad(lg, width = 2, side = c("left"), pad = "0"), sep = ""), 
         env_year = paste(env, year, sep = ""), 
         marker_name = NA) %>%
  dplyr::select(lg, map, left, right, famID, trait, env, year, marker_name, chr, env_year) %>%
  arrange(trait, map, lg, left, env_year) 

  
# write out the result and manually combine approximately overlapping intervals and edit marker names
write.csv(lod_intervals_ed, "/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Yield Components Mapping/Scripts for Github/Visualize Results/data/lod_intervals_ed.csv", row.names = F)
# delete markers that do not occur in more than one environment
View(lod_intervals_ed)
# read back in
combined_qtl <- read.csv("/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Yield Components Mapping/Scripts for Github/Visualize Results/data/lod_intervals_for_lmv.csv", header = T)
combined_qtl <- combined_qtl %>% filter(trait != "seed_area")
#combined_qtl <- read.csv("/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Yield Components Mapping/Scripts for Github/Visualize Results/data/lod_intervals_for_lmv_all.csv", header = T)


combined_qtl <- combined_qtl %>%
  filter(., grepl(',', marker_name)) # filter by those with a comma, indicating detected in more than one ev. 
# though I eliminted these manually anyway

# prepare the map summary
map_sum_lg <- map_sum %>%
  mutate(chr = str_replace(string = group, pattern = "Chr", replacement = "LG"))

# left join in the map sum
combined_qtl_map <- left_join(combined_qtl, map_sum_lg, by = "chr")

# calculate the ajusted positions and prepare order of dataframe for qtldf
combined_qtldf <- combined_qtl_map %>% 
  mutate(so = 1000 * (as.numeric(left) - cm_min) / (cm_max - cm_min),
         si = so,
         ei = 1000 * (as.numeric(right) - cm_min) / (cm_max - cm_min),
         eo = ei,
         qtl = marker_name, 
         col = "black") %>%
  dplyr::select(chr, qtl, so, si, eo, ei, col, map)


# change the color based on the environment in which it was detected
for (i in 1:nrow(combined_qtldf)) {
  if (combined_qtldf$map[i] == "DP") {combined_qtldf$col[i] <- "#56B4E9"}
  else if (combined_qtldf$map[i] == "CP") {combined_qtldf$col[i] <- "#154ba1"}
}

# remove the group column
combined_qtldf <- combined_qtldf %>% dplyr::select(-map)

#### Step 8: Make the Figures ####
# select only needed columns
lmv_dat <- all_data2 %>% 
  dplyr::select(group, position, marker_name)

# take all of the physical map, genes, rulers, and GWAS QTL (Linkage Mapping QTL are in a separate dataframe)
# and nest them into vectors and then lists later on as required by the package

# first, all names that will be featured
names <- all_data2 %>% 
  filter(! is.na(display)) %>%
  dplyr::select(marker_name) 
names <- as.vector(as.matrix(names))

# then the genes
gene <- all_data2 %>% 
  filter(display == "gene") %>%
  dplyr::select(marker_name) 
gene <- as.vector(as.matrix(gene))

# the gwas QTL
gwas_qtl <- all_data2 %>% 
  filter(display == "gwas_qtl") %>%
  dplyr::select(marker_name)
gwas_qtl <- as.vector(as.matrix(gwas_qtl))

# the ruler
ruler <- all_data2 %>% 
  filter(display == "ruler") %>%
  dplyr::select(marker_name)
ruler <- as.vector(as.matrix(ruler))

# and the adjuster data
adjuster <- all_data2 %>% 
  filter(display == "adjuster") %>%
 dplyr::select(marker_name)
adjuster <- as.vector(as.matrix(adjuster))


#### Step 9: Create the Format List ####
# create formatting list for marker names
flist <- list()

# gene names will be italic, black font
genes_vector <- gene
cex <- c(1) #size
font <- c(3) 
col <- c("black")
flist[[1]] <- list(locus = genes_vector, font = font, cex = cex, col = col)

# gwas_qtl will be bold and blue
cex <- c(1) #size
font <- c(1) #bold
col <- c("black")
flist[[2]] <- list(locus = gwas_qtl, font = font, cex = cex, col = col)

# ruler
ruler_vector <- ruler
cex <- c(1) #size
font <- c(2) #plain
col <- c("black")
flist[[3]] <- list(locus = ruler_vector, font = font, cex = cex, col = col)

# adjuster values will be white so they're invisible
adjuster_vector <- adjuster
cex <- c(1) #size
font <- c(2) #plain
col <- c("white")
flist[[4]] <- list(locus = adjuster_vector, font = font, cex = cex, col = col)



#### Step 10: Create Each Plot ####
# set the working directory
setwd(paste(dir, folder, "output", sep = "/"))

# plot homeoelogous groups as columns
# ROW 1 - this row gets messed up because it doesn't have as many components and it's
# a different size compared to the others, so we'll add more adjusters

row_to_add_1 <- combined_qtldf %>% 
  slice(1:1) %>%
  mutate(chr = "LG02", 
         col = '#FFFFFF')

row_to_add_2 <- combined_qtldf %>% 
  slice(9:9) %>%
  mutate(chr = "LG08", 
         col = '#FFFFFF')

row_to_add_3 <- combined_qtldf %>% 
  slice(1:3) %>%
  mutate(chr = "LG07", 
         col = '#FFFFFF')


#combined_qtldf_2 <- rbind(combined_qtldf, row_to_add_1, row_to_add_2, row_to_add_3)

setwd("/users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Yield Components Mapping/Scripts for Github/Visualize Results/output/")
# ROW 1
outfile = file.path(file="./row_1.pdf")
mirror <- c(rep(c(FALSE, TRUE), 7))
lg_list <- c('Chr01', 'LG01', 'Chr04', 'LG04', 'Chr07', 'LG07', 'Chr10', 'LG10', 
             'Chr13', 'LG13', 'Chr16', 'LG16', 'Chr19', 'LG19')
lmv.linkage.plot(lmv_dat, outfile, 
                 mapthese=lg_list,
                 showonly=names,ruler=TRUE,posonleft=mirror,labdist=0.1,lgw=0.2,
                 pdf.height=4, pdf.width = 24,
                 pdf.pointsize=10,maxnbrcolsfordups=1,markerformatlist=flist,qtldf = combined_qtldf,
                 lgperrow=14, 
                 par(lwd=0.25),lty.axis=0.5)

# ROW 2
outfile = file.path(file="./row_2.pdf")
mirror <- c(rep(c(FALSE, TRUE), 7))
lg_list <- c('Chr02', 'LG02', 'Chr05', 'LG05', 'Chr08', 'LG08',  'Chr11', 'LG11', 'Chr14', 'LG14', 'Chr17', 'LG17', 'Chr20', 'LG20')
lmv.linkage.plot(lmv_dat, outfile, 
                 mapthese=lg_list,
                 showonly=names,ruler=TRUE,posonleft=mirror,labdist=0.1,lgw=0.2,
                 pdf.height=4, pdf.width = 24,
                 pdf.pointsize=10,maxnbrcolsfordups=1,markerformatlist=flist,qtldf = combined_qtldf,
                 lgperrow=14, 
                 par(lwd=0.25),lty.axis=0.5)

# ROW 3
outfile = file.path(file="./row_3.pdf")
mirror <- c(rep(c(FALSE, TRUE), 7))
lg_list <- c('Chr03', 'LG03', 'Chr06', 'LG06', 'Chr09', 'LG09',  'Chr12', 'LG12', 'Chr15', 'LG15', 'Chr18', 'LG18', 'Chr21', 'LG21')
lmv.linkage.plot(lmv_dat, outfile, 
                 mapthese=lg_list, 
                 showonly=names,ruler=TRUE,posonleft=mirror,labdist=0.1,lgw=0.2,
                 pdf.height=4,pdf.width = 24,
                 pdf.pointsize=10,maxnbrcolsfordups=1,markerformatlist=flist,qtldf = combined_qtldf,
                 lgperrow=14, 
                 par(lwd=0.25),lty.axis=0.5)


### just using the chromosomes of interest
2, 3, 4, 5, 7, 
8, 9, 14, 15, 
16, 18, 20, 21 # 13 chrs


outfile = file.path(file="./chrs_of_interest_1.pdf")
mirror <- c(rep(c(FALSE, TRUE), 5))
lg_list <- c('Chr02', 'LG02', 'Chr03', 'LG03', 'Chr04', 'LG04',  'Chr05', 'LG05', 'Chr07', 'LG07')
lmv.linkage.plot(lmv_dat, outfile, 
                 mapthese=lg_list, 
                 showonly=names,ruler=TRUE,posonleft=mirror,labdist=0.1,lgw=0.2,
                 pdf.height=4,pdf.width = 24,
                 pdf.pointsize=10,maxnbrcolsfordups=1,markerformatlist=flist,qtldf = combined_qtldf,
                 lgperrow=14, 
                 par(lwd=0.25),lty.axis=0.5)


outfile = file.path(file="./chrs_of_interest_2.pdf")
mirror <- c(rep(c(FALSE, TRUE), 4))
lg_list <- c('Chr08', 'LG08', 'Chr09', 'LG09', 'Chr14', 'LG14', 'Chr15', 'LG15')
lmv.linkage.plot(lmv_dat, outfile, 
                 mapthese=lg_list, 
                 showonly=names,ruler=TRUE,posonleft=mirror,labdist=0.1,lgw=0.2,
                 pdf.height=4,pdf.width = 24,
                 pdf.pointsize=10,maxnbrcolsfordups=1,markerformatlist=flist,qtldf = combined_qtldf,
                 lgperrow=14, 
                 par(lwd=0.25),lty.axis=0.5)


outfile = file.path(file="./chrs_of_interest_3.pdf")
mirror <- c(rep(c(FALSE, TRUE), 7))
lg_list <- c('Chr16', 'LG16', 'Chr18', 'LG18', 'Chr20', 'LG20', 'Chr21', 'LG21', 'Chr19', 'LG19', 'Chr17', 'LG17', 'Chr13', 'LG13')
lmv.linkage.plot(lmv_dat, outfile, 
                 mapthese=lg_list, 
                 showonly=names,ruler=TRUE,posonleft=mirror,labdist=0.1,lgw=0.2,
                 pdf.height=4,pdf.width = 24,
                 pdf.pointsize=10,maxnbrcolsfordups=1,markerformatlist=flist,qtldf = combined_qtldf,
                 lgperrow=14, 
                 par(lwd=0.25),lty.axis=0.5)

