library(tidyverse) 
library(ggbeeswarm)

##### Aivia output was munged using "scripts/munge_AiviaOutput.R" and the master_table was generated)

master_table <- read.table("C:/myPath/2026-07-03-15-22-29 IF_analysis/combined_measurements_2026-07-03-15-22-29 IF_analysis.csv",
                           sep=",", header=TRUE)


##########  Add treatment columns 
master_table$treatment <- NA

master_table$treatment[grepl("ES-", master_table$Image)] <- "ES"
master_table$treatment[grepl("DE-", master_table$Image)] <- "DE"
master_table$treatment[grepl("DE1", master_table$Image)] <- "DE1"
master_table$treatment[grepl("DE2", master_table$Image)] <- "DE2"
master_table$treatment[grepl("DE3", master_table$Image)] <- "DE3"
master_table$treatment[grepl("DE4", master_table$Image)] <- "DE4"
master_table$treatment[grepl("DE5", master_table$Image)] <- "DE5"
master_table$treatment[grepl("DE6", master_table$Image)] <- "DE6"
master_table$treatment[grepl("DE7", master_table$Image)] <- "DE7"
master_table$treatment[grepl("HB", master_table$Image)] <- "HB"
master_table$treatment[grepl("iHep", master_table$Image)] <- "iHep"
master_table$treatment[grepl("mHep", master_table$Image)] <- "mHep"



###### number images from each treatment

image_number <- NA

master_table <- master_table %>%
  mutate(image_number = if_else(
    grepl("ES", treatment),
    as.integer(factor(Image, levels = unique(Image[grepl("ES", treatment)]))),
    image_number  
  ))

master_table <- master_table %>%
  mutate(image_number = if_else(
    grepl("DE", treatment),
    as.integer(factor(Image, levels = unique(Image[grepl("DE", treatment)]))),
    image_number  
  ))

master_table <- master_table %>%
  mutate(image_number = if_else(
    grepl("DE1", treatment),
    as.integer(factor(Image, levels = unique(Image[grepl("DE1", treatment)]))),
    image_number  
  ))

master_table <- master_table %>%
  mutate(image_number = if_else(
    grepl("DE2", treatment),
    as.integer(factor(Image, levels = unique(Image[grepl("DE2", treatment)]))),
    image_number  
  ))

master_table <- master_table %>%
  mutate(image_number = if_else(
    grepl("DE3", treatment),
    as.integer(factor(Image, levels = unique(Image[grepl("DE3", treatment)]))),
    image_number  
  ))

master_table <- master_table %>%
  mutate(image_number = if_else(
    grepl("DE4", treatment),
    as.integer(factor(Image, levels = unique(Image[grepl("DE4", treatment)]))),
    image_number  
  ))

master_table <- master_table %>%
  mutate(image_number = if_else(
    grepl("DE5", treatment),
    as.integer(factor(Image, levels = unique(Image[grepl("DE5", treatment)]))),
    image_number  
  ))

master_table <- master_table %>%
  mutate(image_number = if_else(
    grepl("DE6", treatment),
    as.integer(factor(Image, levels = unique(Image[grepl("DE6", treatment)]))),
    image_number  
  ))

master_table <- master_table %>%
  mutate(image_number = if_else(
    grepl("DE7", treatment),
    as.integer(factor(Image, levels = unique(Image[grepl("DE7", treatment)]))),
    image_number  
  ))

master_table <- master_table %>%
  mutate(image_number = if_else(
    grepl("HB", treatment),
    as.integer(factor(Image, levels = unique(Image[grepl("HB", treatment)]))),
    image_number  
  ))

master_table <- master_table %>%
  mutate(image_number = if_else(
    grepl("iHep", treatment),
    as.integer(factor(Image, levels = unique(Image[grepl("iHep", treatment)]))),
    image_number  
  ))

master_table <- master_table %>%
  mutate(image_number = if_else(
    grepl("mHep", treatment),
    as.integer(factor(Image, levels = unique(Image[grepl("mHep", treatment)]))),
    image_number  
  ))



###### Sort table

master_table <- master_table %>%
  mutate(treatment = factor(treatment, 
                            levels = c("ES", "DE", "DE1","DE2","DE3","DE4","DE5","DE6","DE7","HB","iHep", "mHep")))


master_table <- master_table %>%
  arrange(treatment, image_number)


#### Convert pixel values to um

master_table <- master_table %>%
  mutate(
    vx = if_else(grepl("^DE|HB", treatment), 0.053935186, 0.038449),
    vy = if_else(grepl("^DE|HB", treatment), 0.053935186, 0.038449),
    vz = 0.27,
    Volume..µm.. = Volume..px.. * vx * vy * vz
  ) %>%
  select(-vx, -vy, -vz) %>%
  relocate(Volume..µm.., .after = 3)


long_vol <- master_table[,c('Image','treatment','Volume..µm..')]
long_vol <- long_vol[long_vol$treatment ==  c("ES","DE", "HB", "iHep","mHep"),]
long_vol$treatment <- factor(long_vol$treatment, levels = c("ES","DE","HB","iHep","mHep"))

count_data <- long_vol %>%
  na.omit() %>%
  group_by(treatment) %>%
  summarise(n = n(), y_pos = max(Volume..µm..), .groups = "drop")

## Fig.2b - Speckle volume by treatment (IF)

setwd("C:/myPath/")
pdf("Speckle_vol_by_stage_IF.pdf", height = 2.34, width =5)
ggplot(na.omit(long_vol), aes(treatment, Volume..µm..)) +
  geom_violin(aes(fill = treatment), color = NA) +
  geom_boxplot(outlier.colour = NA, width = 0.15, fill = "white", linewidth = 0.5) +
  geom_text(data = count_data, aes(x = treatment, y = y_pos, label = n),
            size = 6, vjust = -0.5, show.legend = FALSE, color = "black") +
  scale_y_log10(limits = c(-10000, 200)) +
  # scale_y_log10() +
  scale_fill_manual(values = c(
    "ES" = "grey50"
    # add other treatments here, e.g. "TPL/DRB" = "blue"
  )) +
  theme_minimal(base_size = 20) +
  scale_fill_manual(values = c("ES"= "#2C6FB5","DE"="#00A878","HB"="#D4AC0D","iHep"="#EE5A24","mHep"="#5E4327")) +
  labs(x = NULL) +
  ylab(bquote('Volume '(µm^3))) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        axis.ticks = element_line(colour = "black", linewidth = 0.4, lineend = "round"),
        axis.ticks.length = unit(0.1, "cm"),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 30, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")

dev.off()


###########  DE timecourse ############

long_vol <- master_table[,c('Image','treatment','Volume..µm..')]
long_vol <- long_vol[long_vol$treatment ==  c("DE", "DE1","DE2","DE3","DE4","DE5","DE6","DE7","HB"),]

count_data <- long_vol %>%
  na.omit() %>%
  group_by(treatment) %>%
  summarise(n = n(), y_pos = max(Volume..µm..), .groups = "drop")

pal_fun <- colorRampPalette(c("#00A878", "#D4AC0D"))

# e.g. 5 evenly spaced colors between the two
pal <- pal_fun(7)

## Fig.2b - Speckle volume by treatment (IF)

setwd("C:/myPath/")
pdf("Speckle_vol_by_stageDEtoHB_IF.pdf", height = 3.8, width = 9)
ggplot(na.omit(long_vol), aes(treatment, Volume..µm..)) +
  geom_violin(aes(fill = treatment), color = NA) +
  geom_boxplot(outlier.colour = NA, width = 0.15, fill = "white", linewidth = 0.5) +
  geom_text(data = count_data, aes(x = treatment, y = y_pos, label = n),
            size = 6, vjust = -0.5, show.legend = FALSE, color = "black") +
  scale_y_log10(limits = c(-10000, 200)) +
  # scale_y_log10() +
  scale_fill_manual(values = c(
    "ES" = "grey50"
    # add other treatments here, e.g. "TPL/DRB" = "blue"
  )) +
  theme_minimal(base_size = 20) +
  scale_fill_manual(values = c(
    "DE"="#00A878","DE1"="#00A878","DE2"="#23A866","DE3"="#46A954",
    "DE4"="#6AAA42","DE5"="#8DAA30","DE6"="#B0AB1E","DE7"="#D4AC0D",
    "HB"="#D4AC0D")) +
  labs(x = NULL) +
  ylab(bquote('Volume '(µm^3))) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2),
        axis.ticks = element_line(colour = "black", linewidth = 0.2, lineend = "round"),
        axis.ticks.length = unit(0.1, "cm"),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 30, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")
dev.off()



######## Plot number of speckles per cell
w <- as.data.frame(table(master_table$Image))

#w$treatment <- stringr::str_extract(w$Var1, "HB|ES|DE-|DE1|DE2|DE3|DE4|DE5|DE6|DE7|iHep|mHep")
w$treatment <- stringr::str_extract(w$Var1, "HB|ES|DE-|iHep|mHep")
w$treatment <- gsub("-", "", w$treatment)
names(w)[names(w) == "Var1"] <- "cell"

w$treatment <- factor(w$treatment, levels = c("ES","DE","HB","iHep","mHep"))

## Fig.2c - Speckle count by treatment (IF)

setwd("C:/myPath/")
pdf("Speckle_count_by_stage_IF_v10.pdf", height = 2.5, width = 5)
count_data <- na.omit(w) %>%
  group_by(treatment) %>%
  summarise(n = n(), y_pos = max(Freq), .groups = "drop")

ggplot(na.omit(w), aes(treatment, Freq)) +
  geom_violin(aes(fill = treatment), color = NA) +
  geom_boxplot(outlier.colour = NA, width = 0.15, fill = "white", linewidth = 0.5) +
  geom_beeswarm(cex = 2, size = 1, color = "black") +
  
  geom_text(data = count_data, aes(x = treatment, y = y_pos, label = n),
            size = 6, vjust = -0.5, show.legend = FALSE, color = "black") +
  # scale_y_log10(limits = c(-10000, 200)) +
  # scale_y_log10() +
  scale_fill_manual(values = c(
    "ES" = "grey50"
    # add other treatments here, e.g. "TPL/DRB" = "blue"
  )) +
  ylim(0,230) +
  theme_minimal(base_size = 20) +
  scale_fill_manual(values = c("ES"= "#2C6FB5","DE"="#00A878","HB"="#D4AC0D","iHep"="#EE5A24","mHep"="#5E4327")) +
  labs(x = NULL, y = NULL) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        axis.ticks = element_line(colour = "black", linewidth = 0.4, lineend = "round"),
        axis.ticks.length = unit(0.1, "cm"),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 30, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")

dev.off()

w <- as.data.frame(table(master_table$Image))

#w$treatment <- stringr::str_extract(w$Var1, "HB|ES|DE-|DE1|DE2|DE3|DE4|DE5|DE6|DE7|iHep|mHep")
w$treatment <- stringr::str_extract(w$Var1, "DE-|DE1|DE2|DE3|DE4|DE5|DE6|DE7|HB")
w$treatment <- gsub("-", "", w$treatment)
names(w)[names(w) == "Var1"] <- "cell"

w$treatment <- factor(w$treatment, levels = c("DE", "DE1","DE2","DE3","DE4","DE5","DE6","DE7","HB"))

## Fig.2c - Speckle count by treatment (IF)

setwd("C:/myPath/")
pdf("Speckle_count_by_stageDEtoHB_IF_v10.pdf", height = 7.3, width = 18)
count_data <- na.omit(w) %>%
  group_by(treatment) %>%
  summarise(n = n(), y_pos = max(Freq), .groups = "drop")

ggplot(na.omit(w), aes(x = treatment, y = Freq)) +
  geom_violin(aes(fill = treatment), color = NA, trim = FALSE) +
  geom_boxplot(outlier.colour = NA, width = 0.15, fill = "white",
               linewidth = 0.1) +
  geom_beeswarm(cex = 0.5, size = 0.5, color = "black") +
  geom_text(data = count_data, aes(x = treatment, y = y_pos, label = n),
            vjust = -0.5, size = 8, color = "black") +
  scale_fill_manual(values = c("DE"="#00A878","DE1"="#00A878","DE2"="#23A866","DE3"="#46A954","DE4"="#6AAA42","DE5"="#8DAA30","DE6"="#B0AB1E","DE7"="#D4AC0D","HB"="#D4AC0D")) +
  theme_minimal(base_size = 20) +
  
  labs(x = NULL, y = "Speckles per cell") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2),
        axis.ticks = element_line(colour = "black", linewidth = 0.2, lineend = "round"),
        axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(size = 40, colour = "black"),
        axis.title = element_text(size = 30, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")


dev.off()

