# Project: IWG_NAM_Yield Components Mapping
# Visualize Results
# Author: Kayla R. Altendorf
# Date: 09/06/21

# load packages
library('sommer')

# location of github directory
dir <- "/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Yield Components Mapping/Scripts for Github"

# script we're on
folder <- c("/Visualize Results")
traits <- c("flag_leaf_area", "floret_site_utilization", "florets_per_spikelet", "height", "reproductive_tiller_ct", 
            "seed_area", "seeds_per_spike", "spikelet_density", "spikelets_per_spike", "stem_diameter", 
            "thousand_grain_weight", "yield_per_spike")

year <- c("2017", "2018", "2017", "2018")
loc <- c("STP", "STP", "TLI", "TLI")

#### Step 1: Read in GWAS Results ####
files <- list.files(paste(dir, "GWAS/output/GAPIT/flag_leaf_area/", sep = "/"), pattern = "GWAS.Results.csv", full.names = T)
flag_leaf_area <- lapply(files, read.csv)

files <- list.files(paste(dir, "GWAS/output/GAPIT/floret_site_utilization/", sep = "/"), pattern = "GWAS.Results.csv", full.names = T)
floret_site_utilization <- lapply(files, read.csv)

files <- list.files(paste(dir, "GWAS/output/GAPIT/florets_per_spikelet/", sep = "/"), pattern = "GWAS.Results.csv", full.names = T)
florets_per_spikelet <- lapply(files, read.csv)

files <- list.files(paste(dir, "GWAS/output/GAPIT/height/", sep = "/"), pattern = "GWAS.Results.csv", full.names = T)
height <- lapply(files, read.csv)

files <- list.files(paste(dir, "GWAS/output/GAPIT/reproductive_tiller_ct/", sep = "/"), pattern = "GWAS.Results.csv", full.names = T)
reproductive_tiller_ct <- lapply(files, read.csv)

files <- list.files(paste(dir, "GWAS/output/GAPIT/seed_area/", sep = "/"), pattern = "GWAS.Results.csv", full.names = T)
seed_area <- lapply(files, read.csv)

files <- list.files(paste(dir, "GWAS/output/GAPIT/seeds_per_spike/", sep = "/"), pattern = "GWAS.Results.csv", full.names = T)
seeds_per_spike <- lapply(files, read.csv)

files <- list.files(paste(dir, "GWAS/output/GAPIT/spikelet_density/", sep = "/"), pattern = "GWAS.Results.csv", full.names = T)
spikelet_density <- lapply(files, read.csv)

files <- list.files(paste(dir, "GWAS/output/GAPIT/spikelets_per_spike/", sep = "/"), pattern = "GWAS.Results.csv", full.names = T)
spikelets_per_spike <- lapply(files, read.csv)

files <- list.files(paste(dir, "GWAS/output/GAPIT/stem_diameter/", sep = "/"), pattern = "GWAS.Results.csv", full.names = T)
stem_diameter <- lapply(files, read.csv)

files <- list.files(paste(dir, "GWAS/output/GAPIT/thousand_grain_weight/", sep = "/"), pattern = "GWAS.Results.csv", full.names = T)
thousand_grain_weight <- lapply(files, read.csv)

files <- list.files(paste(dir, "GWAS/output/GAPIT/yield_per_spike/", sep = "/"), pattern = "GWAS.Results.csv", full.names = T)
yield_per_spike <- lapply(files, read.csv)

#### Step 2: Make Manhattan Plots ####

# name publication titles
traits
pub_names <- c("Flag Leaf Area", "Floret Site Utilization", "Florets Spikelet-1", "Height", 
               "Reproductive Tiller", "Seed Area", "Seeds Spike-1", "Spikelet Density",
               "Spikelets Spike-1", "Stem Diameter", "Thousand Grain Weight", "Yield Spike-1")

name_frame <- data.frame(traits = traits, pub_names = pub_names)
# edit trait and save each plot separately

# first read in files
for (i in 1:length(traits)) {
  files <- list.files(paste(dir, "GWAS/output/GAPIT/", traits[i], sep = "/"), pattern = "GWAS.Results.csv", full.names = T)
  files <- lapply(files, read.csv)
  
  for (j in 1:length(files)) {
    files[[j]] <- files[[j]] %>% dplyr::select(Chromosome, Position, P.value)
    colnames(files[[j]]) <- c("Chrom", "Position", "p.val")
    
    files[[j]] <- files[[j]] %>% 
      
    # compute chromosome size
    group_by(Chrom) %>%
    summarise(chr_len = max(as.numeric(Position))) %>%
      
    # calculate cumulative position of each chromosome
    mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
    dplyr::select(-chr_len) %>% 
  
    # Add this info to the initial dataset
    left_join(files[[j]], ., by = c("Chrom" = "Chrom")) %>%
  
    # Add a cumulative position of each SNP
    arrange(Chrom, Position) %>%
    mutate(BPcum = Position + tot) %>%
  
    # add loc and year info to each plot
    mutate(loc = loc[j], 
          year = year[j])
  }
  
  dat <- do.call("rbind", files)
  files <- list()
  axisdf = dat %>% group_by(Chrom) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  #Ready to make the plot using ggplot2:
  plot <- ggplot(dat, aes(x=BPcum, y=-log10(p.val))) +
    facet_grid(rows = loc ~ year) +
  
  # show all points
  geom_point( aes(color=as.factor(Chrom)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  geom_hline(yintercept = 3.6, linetype = "dashed", color = "grey", size = 0.5) + 
  
  # custom X axis:
  scale_x_continuous( label = axisdf$Chrom, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  xlab("Chromosome") + 
    
  # custom theme:
  theme_bw() +
    labs(title = name_frame$pub_names[i]) +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
          axis.text.x = element_text(size = 15), 
          axis.text.y = element_text(size = 20),
          axis.title = element_text(size=20, face = "bold"), 
          strip.text.x = element_text(size = 20),
          strip.background.x = element_rect(fill = "white", colour = NA), 
          strip.text.y = element_text(size = 20),
          strip.background.y = element_rect(fill = "white", colour = NA), 
          legend.position = "none", 
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.spacing.y = unit(1, "lines")
          )
  
  ggsave(paste(traits[i], "_gwas.tiff", sep = ""), 
         plot = last_plot(), 
         path = paste(dir, folder, "output", sep = "/"),
         device = "tiff",
         scale = 1, 
         width = 15, 
         height = 10, 
         units = "in")
}

    