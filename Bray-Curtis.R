library(data.table)
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(patchwork)

#--------------------------------------------------
# 1) LOAD DATA
#--------------------------------------------------

taxa <- fread("/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/taxa_lists/23S/23S_combined_taxa_aggregated.tsv")
meta <- fread("/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/Metadata/23S_LEAP22_metadata.tsv")

#--------------------------------------------------
# 2) CLEAN TO GENUS LEVEL
#--------------------------------------------------

# Keep genus + sample columns
sample_cols <- setdiff(names(taxa),
                       c("domain","supergroup","division","subdivision",
                         "class","order","family","species"))

taxa_genus <- taxa[, c(sample_cols), with = FALSE]

# Remove unknown/unassigned
taxa_genus <- taxa_genus[!genus %in% c("Unknown","unknown","Unassigned","unassigned")]

# Collapse duplicate genera
taxa_genus <- taxa_genus %>%
  group_by(genus) %>%
  summarise(across(where(is.numeric), sum), .groups="drop")

# Compute column sums (excluding genus column)
sample_matrix <- taxa_genus[, -1, with = FALSE]
col_totals <- colSums(sample_matrix)

# Identify empty samples
empty_samples <- names(col_totals[col_totals == 0])

# Remove empty samples
if(length(empty_samples) > 0) {
  message("Removing empty samples: ", paste(empty_samples, collapse = ", "))
  taxa_genus <- taxa_genus[, !names(taxa_genus) %in% empty_samples, with = FALSE]
}

# Convert to relative abundance (sum = 100 per sample)
taxa_genus_rel <- taxa_genus
taxa_genus_rel[ , -1] <- sweep(taxa_genus[ , -1], 2,
                               colSums(taxa_genus[ , -1]), "/") * 100

#--------------------------------------------------
# 3) TRANSPOSE TO SAMPLE x GENUS
#--------------------------------------------------

taxa_mat <- as.data.frame(taxa_genus_rel)
rownames(taxa_mat) <- taxa_mat$genus
taxa_mat$genus <- NULL

taxa_mat <- t(taxa_mat)
taxa_mat <- as.data.frame(taxa_mat)
taxa_mat$sample_name <- rownames(taxa_mat)

# Merge with metadata
data_full <- merge(meta, taxa_mat, by="sample_name")

#--------------------------------------------------
# 4) FUNCTION TO RUN BRAY + PCoA PER TIME
#--------------------------------------------------
run_pcoa_time <- function(timepoint) {
  
  df <- data_full %>% filter(time == timepoint)
  
  if(nrow(df) < 3) return(NULL)
  
  genus_cols <- colnames(taxa_mat)[!colnames(taxa_mat) %in% "sample_name"]
  
  comm <- as.matrix(df[, ..genus_cols])
  
  # Replace NA with 0 for Bray-Curtis
  comm[is.na(comm)] <- 0
  
  # Bray-Curtis
  bray <- vegdist(comm, method="bray")
  
  # PCoA
  ord <- cmdscale(bray, eig=TRUE, k=2)
  
  scores <- as.data.frame(ord$points)
  colnames(scores) <- c("PC1","PC2")
  
  df_plot <- cbind(df, scores)
  
  #--------------------------------------------------
  # Nutrient handling
  #--------------------------------------------------
  
  df_plot$nutrient_addition[is.na(df_plot$nutrient_addition)] <- "control"
  
  df_plot$nutrient_addition <- factor(
    as.character(df_plot$nutrient_addition),
    levels = c("control","0","25","50","100","200")
  )
  
  #--------------------------------------------------
  # Shape handling (EXPLICIT replicate usage)
  #--------------------------------------------------
  
  # Start by using replicate column directly
  df_plot$shape_group <- as.character(df_plot$replicate)
  
  # Override blanks
  df_plot$shape_group[df_plot$SampleType == "blank"] <- "blank"
  
  # Override reservoir
  df_plot$shape_group[df_plot$SampleType == "resv"] <- "resv"
  
  #--------------------------------------------------
  # Build connectivity dataframe (UNCHANGED)
  #--------------------------------------------------
  
  connections <- data.frame(
    x = numeric(),
    y = numeric(),
    xend = numeric(),
    yend = numeric(),
    flow = numeric()
  )
  
  for(i in 1:nrow(df_plot)) {
    
    row_i <- df_plot[i,]
    
    if(is.na(row_i$flow_rate) || row_i$flow_rate == 0) next
    
    target <- df_plot %>%
      filter(tank == row_i$connected_tank,
             replicate == row_i$replicate)
    
    if(nrow(target) == 0) next
    
    connections <- rbind(connections,
                         data.frame(
                           x = row_i$PC1,
                           y = row_i$PC2,
                           xend = target$PC1[1],
                           yend = target$PC2[1],
                           flow = row_i$flow_rate
                         ))
  }
  
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  
  p <- ggplot(df_plot, aes(PC1, PC2)) +
    
    geom_segment(data = connections,
                 aes(x=x, y=y, xend=xend, yend=yend, size=flow),
                 arrow = arrow(length=unit(0.15,"cm")),
                 inherit.aes = FALSE) +
    
    geom_point(aes(color = nutrient_addition,
                   shape = shape_group),
               size=3,
               stroke=1) +
    
    scale_color_manual(
      values = c(
        "control" = "#000000",
        "0"  = "#0072B2",
        "25" = "#56B4E9",
        "50" = "#E69F00",
        "100"= "#F0E442",
        "200"= "#009E73"
      ),
      name = "Nutrient Addition"
    ) +
    
    scale_shape_manual(
      values = c(
        "R1" = 16,      # circle
        "R2" = 18,      # diamond
        "blank" = 21,   # hollow circle
        "resv" = 16     # filled circle
      ),
      name = "Replicate"
    ) +
    
    scale_size(range=c(0.5,2), name="Flow Rate") +
    
    labs(x=NULL, y=NULL) +
    
    ggtitle(paste("Time", gsub("T","",timepoint))) +
    
    theme_classic() +
    theme(
      plot.title = element_text(hjust=0.5),
      legend.position = "bottom"
    )
  
  return(p)
}





#--------------------------------------------------
# 5) RUN FOR ALL TIMEPOINTS
#--------------------------------------------------

times <- c("T1","T2","T3","T4","T5")

plots <- lapply(times, run_pcoa_time)

#--------------------------------------------------
# 6) FIX IDENTICAL AXIS LIMITS
#--------------------------------------------------

# Extract all coordinates to calculate global limits
all_scores <- data.frame()

for(t in times){
  df <- data_full %>% filter(time == t)
  if(nrow(df) < 3) next
  
  genus_cols <- colnames(taxa_mat)[!colnames(taxa_mat) %in% "sample_name"]
  
  comm <- as.matrix(df[, ..genus_cols])   # <-- FIXED HERE
  
  bray <- vegdist(comm, method="bray")
  ord <- cmdscale(bray, eig=TRUE, k=2)
  scores <- as.data.frame(ord$points)
  all_scores <- rbind(all_scores, scores)
}

xlim <- range(all_scores$V1)
ylim <- range(all_scores$V2)

plots <- lapply(plots, function(p){
  p + coord_cartesian(xlim=xlim, ylim=ylim)
})

#--------------------------------------------------
# 7) COMBINE INTO ONE ROW
#--------------------------------------------------

final_plot <- wrap_plots(plots, nrow=1, guides="collect") &
  theme(legend.position="bottom")

final_plot

ggsave(
  filename = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/Bray-Curtis_18S.png",
  plot = final_plot,
  width = 14,
  height = 4,
  units = "in",
  dpi = 600
)
