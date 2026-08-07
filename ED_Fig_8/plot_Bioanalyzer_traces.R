library(tidyverse)
library(stringr)
library(RColorBrewer)
library(paletteer)

####  Parameters
data_folder <- "C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Fragment Analyzer/Analysis/data/4DN_sample/"

run_dates <- c("20210329","20210607","20210623","20240104","20240607","20250625","20250806","20260714")
samples_ByRun_dates <- c("20210329","20210607","20210623","20240125","20240611","20250324|20250624","20250806","20260616")


###### List file names  
LadderList <- list.files(paste(data_folder, "ladders", sep =""), pattern = ".csv", full.name = TRUE)

SampleList <- list.files(paste(data_folder,"samples", sep=""), pattern = ".txt", full.name = TRUE)



##########################################   Functions  ######################################################
#######  Functions to extract column names from file names
col_name_ladder <- function(x){
  a <- unlist(str_split(x,"/"))  
  b <- tail(a,1)
  ladder_name <- paste(unlist(str_split(b," "))[1], unlist(str_split(b," "))[2],unlist(str_split(b," "))[3], "_ladder", sep="")
}

col_name_sample <- function(x){
  a <- unlist(str_split(x,"/"))  
  b <- tail(a,1)
  sample_name <- unlist(str_split(b,"\\."))[1]
}


####### Function to generate point to point linear models from ladder data
create_piecewise_models <- function(data) {
  # Extract x and y values
  x_vals <- data[[1]]
  y_vals <- data[[2]]
  
  n_points <- length(x_vals)
  
  # Create list to store models
  models <- list()
  
  # Create linear model for each segment
  for (i in 1:(n_points - 1)) {
    # Extract two consecutive points
    x_segment <- x_vals[i:(i+1)]
    y_segment <- y_vals[i:(i+1)]
    
    # Create data frame for this segment
    segment_data <- data.frame(x = x_segment, y = y_segment)
    
    # Fit linear model
    model <- lm(y ~ x, data = segment_data)
    
    # Store model with segment information
    models[[i]] <- list(
      model = model,
      x_range = c(x_segment[1], x_segment[2]),
      segment = i,
      slope = coef(model)[2],
      intercept = coef(model)[1]
    )
    
    # Print segment information
    cat(sprintf("Segment %d: Migration Time = [%.1f, %.1f] sec, Fragment Size = %.2f + %.4f*time\n", 
                i, x_segment[1], x_segment[2], 
                coef(model)[1], coef(model)[2]))
  }
  
  return(models)
}


######## Function to predict fragment size from migration time
predict_fragment_size <- function(migration_times, models) {
  predictions <- numeric(length(migration_times))
  
  for (i in seq_along(migration_times)) {
    time_val <- migration_times[i]
    
    # Find which segment this migration time belongs to
    segment_found <- FALSE
    
    for (j in seq_along(models)) {
      time_range <- models[[j]]$x_range
      
      # Check if time_val is within this segment's range
      if (time_val >= time_range[1] && time_val <= time_range[2]) {
        # Predict using this segment's model
        predictions[i] <- predict(models[[j]]$model, newdata = data.frame(x = time_val))
        segment_found <- TRUE
        break
      }
    }
    
    # Handle extrapolation cases
    if (!segment_found) {
      if (time_val < models[[1]]$x_range[1]) {
        # Extrapolate using first segment
        predictions[i] <- predict(models[[1]]$model, newdata = data.frame(x = time_val))
        warning(paste("Extrapolating below calibration range for migration time =", time_val))
      } else if (time_val > models[[length(models)]]$x_range[2]) {
        # Extrapolate using last segment
        predictions[i] <- predict(models[[length(models)]]$model, newdata = data.frame(x = time_val))
        warning(paste("Extrapolating above calibration range for migration time =", time_val))
      }
    }
  }
  
  return(predictions)
}


predict_migration_time <- function(fragment_sizes, models) {
  n_predictions <- length(fragment_sizes)
  n_segments <- length(models)
  
  # Pre-calculate all segment parameters to avoid repeated calculations
  segment_params <- matrix(nrow = n_segments, ncol = 6)
  colnames(segment_params) <- c("slope", "intercept", "x_min", "x_max", "y_min", "y_max")
  
  for (j in seq_len(n_segments)) {
    slope <- models[[j]]$slope
    intercept <- models[[j]]$intercept
    x_range <- models[[j]]$x_range
    
    # Calculate y-range once
    y_start <- intercept + slope * x_range[1]
    y_end <- intercept + slope * x_range[2]
    
    segment_params[j, ] <- c(slope, intercept, x_range[1], x_range[2], 
                             min(y_start, y_end), max(y_start, y_end))
  }
  
  # Vectorized approach for finding segments
  predictions <- numeric(n_predictions)
  
  # Create matrices for vectorized comparisons
  size_matrix <- matrix(fragment_sizes, nrow = n_predictions, ncol = n_segments)
  y_min_matrix <- matrix(segment_params[, "y_min"], nrow = n_predictions, ncol = n_segments, byrow = TRUE)
  y_max_matrix <- matrix(segment_params[, "y_max"], nrow = n_predictions, ncol = n_segments, byrow = TRUE)
  
  # Find which segments contain each fragment size
  in_range <- (size_matrix >= y_min_matrix) & (size_matrix <= y_max_matrix)
  
  # Find the first matching segment for each fragment size
  segment_indices <- apply(in_range, 1, function(x) {
    idx <- which(x)
    if (length(idx) > 0) idx[1] else NA
  })
  
  # Process predictions where segments were found
  valid_idx <- !is.na(segment_indices)
  valid_segments <- segment_indices[valid_idx]
  valid_sizes <- fragment_sizes[valid_idx]
  
  # Vectorized reverse prediction for valid segments
  slopes <- segment_params[valid_segments, "slope"]
  intercepts <- segment_params[valid_segments, "intercept"]
  
  # Check for near-zero slopes
  near_zero <- abs(slopes) < 1e-10
  
  # Calculate predictions
  predictions[valid_idx] <- ifelse(near_zero, 
                                   (segment_params[valid_segments, "x_min"] + segment_params[valid_segments, "x_max"]) / 2,
                                   (valid_sizes - intercepts) / slopes)
  
  # Warn about near-zero slopes
  if (any(near_zero)) {
    zero_indices <- which(valid_idx)[near_zero]
    for (idx in zero_indices) {
      warning(paste("Near-zero slope in segment", valid_segments[sum(valid_idx[1:idx])], 
                    "for fragment size =", fragment_sizes[idx]))
    }
  }
  
  # Handle extrapolation cases
  extrap_idx <- is.na(segment_indices)
  if (any(extrap_idx)) {
    extrap_sizes <- fragment_sizes[extrap_idx]
    
    # Determine extrapolation direction
    first_y_min <- segment_params[1, "y_min"]
    first_y_max <- segment_params[1, "y_max"]
    last_y_min <- segment_params[n_segments, "y_min"]
    last_y_max <- segment_params[n_segments, "y_max"]
    
    # Find overall y-range
    overall_y_min <- min(first_y_min, last_y_min)
    overall_y_max <- max(first_y_max, last_y_max)
    
    below_range <- extrap_sizes < overall_y_min
    above_range <- extrap_sizes > overall_y_max
    
    # Extrapolate below range using first segment
    if (any(below_range)) {
      below_indices <- which(extrap_idx)[below_range]
      slope_1 <- segment_params[1, "slope"]
      intercept_1 <- segment_params[1, "intercept"]
      
      if (abs(slope_1) < 1e-10) {
        predictions[below_indices] <- segment_params[1, "x_min"]
        for (idx in below_indices) {
          warning(paste("Near-zero slope in first segment for fragment size =", fragment_sizes[idx]))
        }
      } else {
        predictions[below_indices] <- (extrap_sizes[below_range] - intercept_1) / slope_1
      }
      
      for (idx in below_indices) {
        warning(paste("Extrapolating below calibration range for fragment size =", fragment_sizes[idx]))
      }
    }
    
    # Extrapolate above range using last segment
    if (any(above_range)) {
      above_indices <- which(extrap_idx)[above_range]
      slope_last <- segment_params[n_segments, "slope"]
      intercept_last <- segment_params[n_segments, "intercept"]
      
      if (abs(slope_last) < 1e-10) {
        predictions[above_indices] <- segment_params[n_segments, "x_max"]
        for (idx in above_indices) {
          warning(paste("Near-zero slope in last segment for fragment size =", fragment_sizes[idx]))
        }
      } else {
        predictions[above_indices] <- (extrap_sizes[above_range] - intercept_last) / slope_last
      }
      
      for (idx in above_indices) {
        warning(paste("Extrapolating above calibration range for fragment size =", fragment_sizes[idx]))
      }
    }
  }
  
  return(predictions)
}





#################################### Load and munge ladder tables  ##########################

#####  Load ladder tables  #####   NOTE: these tables were obtained by clicking File > Export Data and selecting Size Calibration data (for each run)
ladder_names <- lapply(LadderList, col_name_ladder)
ladder_names <- unlist(ladder_names) 
ladder_Data <- lapply(1:length(LadderList), function(x) { data.frame(read_csv(LadderList[[x]], col_names = TRUE))})


for (i in c(1:length(ladder_Data))) {
  names(ladder_Data[[i]]) <- c("ladder", paste("time_",ladder_names[[i]],sep=""))
}
names(ladder_Data) <- c(ladder_names)


ladder_Data <- lapply(ladder_Data, function(df) {
  df[, c(2, 1)] # Reorder columns by index
})


piecewise_models <- lapply(ladder_Data, create_piecewise_models)

###################################### Load and munge sample tables   #########################

#####  Load sample tables  #####

sample_names <- lapply(SampleList, col_name_sample)
sample_names <- unlist(sample_names) 
sample_Data <- lapply(1:length(SampleList), function(x) { read.delim(SampleList[[x]], header = TRUE)})

sample_Data <- lapply(sample_Data, function(df) {
  df[,1:2]
})


for (i in c(1:length(sample_Data))) {
  names(sample_Data[[i]]) <- c("time","rfu")
}
names(sample_Data) <- c(sample_names)


##### Subset tables and group by electrophoresis run (same FA run)
for (i in c(1:length(run_dates))) {
  assign(paste("samples_",run_dates[i],sep=""), sample_Data[grep(samples_ByRun_dates[i], names(sample_Data), value = TRUE)])
}

samples_ByRun_names <- c()
for (i in c(1:length(run_dates))) {
  samples_ByRun_names <- append(samples_ByRun_names, paste("samples_", run_dates[i],sep=""))
}

samples_ByRun <- mget(samples_ByRun_names)


################################ Run predicted fragment sizes  ########################

####### use piece-wise linear models from calibration curves to predict fragment size

for (i in c(1:length(samples_ByRun))) {
  for (j in c(1:length(samples_ByRun[[i]]))) {
    a <- predict_fragment_size(samples_ByRun[[i]][[j]][,1], piecewise_models[[i]])  
    samples_ByRun[[i]][[j]][["Fragment_size"]] <- a
  }
}




##### Make merged calibration data in order to later calculate average migration times from calibration curves used in comparison  ######

##### Merge ladder files by fragment size  ####
ladder_merge <- ladder_Data[[1]][2]

for (i in c(1:length(ladder_Data))) {
  ladder_merge <- merge(ladder_merge, ladder_Data[[i]], by=c("ladder"), all.x= TRUE)
}










############################# Normalize data from different runs and axis scaling ####################################

run_dates <- "20210329|20210607|20210623|20240104|20240607|20250625|20250806|20260714"


ladder_names <- names(ladder_Data)[grep((run_dates), names(ladder_Data))]
s_dat_names <- names(samples_ByRun)[grep((run_dates), names(samples_ByRun))]
s_dat_nameslong <- samples_ByRun[grep((run_dates), names(samples_ByRun))]
s_dat_nameslong <- names(unlist(s_dat_nameslong, recursive = FALSE))
s_dat_nameslong <- substr(s_dat_nameslong, 18, nchar(s_dat_nameslong))

s_dat <- samples_ByRun[grep((run_dates), names(samples_ByRun))]
s_dat <- unlist(s_dat, recursive = FALSE)
names(s_dat) <- c(s_dat_nameslong) 

#### Average time measurements from relevant calibration curves and render new table called dat
dat <- ladder_Data[[1]][2]
ladder_time_names <- unlist(lapply(ladder_names, function(x) paste0("time_", x)))
dat$mean_time <- rowMeans(ladder_merge[, c(ladder_time_names)])
dat <- dat[, c(2, 1)]

#####  Scale time for ladder file  #####
dat[["scaled_time"]] <- (dat[,1] - min(dat[,1])) / (head(tail(dat[,1], n=1),n=1) - min(dat[,1]))

#### Calculate piece-wise linear models from Averaged calibration curve
plot_models <- create_piecewise_models(dat)

#### Convert predicted FA to predicted time using model (averaged calibration curves)
#### Calculate normalized migration time

for (i in c(1:length(s_dat))) {
  a <- predict_migration_time(s_dat[[i]][,3], plot_models) 
  s_dat[[i]][["timePlot"]] <- a
}

##### Scale RFU values by main DpnII peak  ######
s_dat_names <- names(s_dat) 

s_dat <- lapply(names(s_dat), function(name) {
  df <- s_dat[[name]]
  new_col_name <- "rfu_scaled"
  df[[new_col_name]] <- (df[,2] - min(df[,2])) / (max(df[1200:2500,2]) - min(df[,2]))
  return(df)
})
names(s_dat) <- c(s_dat_names) 

##### Scale RFU values by large band  ######

s_dat <- lapply(names(s_dat), function(name) {
  df <- s_dat[[name]]
  new_col_name <- "rfu_upper_scaled"
  df[[new_col_name]] <- (df[,2] - min(df[,2])) / (max(df[2500:3000,2]) - min(df[,2]))
  return(df)
})
names(s_dat) <- c(s_dat_names) 


##### Min/max scaling of time using upper and lower bands (TimePLot)

s_dat <- lapply(names(s_dat), function(name) {
  df <- s_dat[[name]]
  subset_df <- df[df[,1] >= 900 & df[,1] <= 1300, ]
  small_marker <- subset_df[which.max(subset_df[,2]),1] +1
  subset_df <- df[df[,1] >= 2600 & df[,1] <= 3200, ]
  big_marker <- subset_df[which.max(subset_df[,2]),1] +1
  new_col_name <- "timePlot_scaled"
  df[[new_col_name]] <- (df[,4] - df[,4][as.numeric(small_marker)]) / (df[,4][as.numeric(big_marker)] - df[,4][as.numeric(small_marker)])
  return(df)
})
names(s_dat) <- c(s_dat_names) 


##### Min/max scaling using upper and lower bands (raw_time)

s_dat <- lapply(names(s_dat), function(name) {
  df <- s_dat[[name]]
  subset_df <- df[df[,1] >= 900 & df[,1] <= 1300, ]
  small_marker <- subset_df[which.max(subset_df[,2]),1] +1
  subset_df <- df[df[,1] >= 2600 & df[,1] <= 3200, ]
  big_marker <- subset_df[which.max(subset_df[,2]),1] +1
  new_col_name <- "time_scaled"
  df[[new_col_name]] <- (df[,1] - df[,1][as.numeric(small_marker)]) / (df[,1][as.numeric(big_marker)] - df[,1][as.numeric(small_marker)])
  return(df)
})
names(s_dat) <- c(s_dat_names) 




############################################ Plot ############################################################## 



##############  U2OStimecourse Grey ###########

setwd("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2025_chromatin_hardening/Results/Fig_5/Fragment_analyzer")
pdf("FA_traces_U2OSno3h_timecourse_LMnorm_upperMarkerScale_TimePlotScaling_grey.pdf", height = 6, width = 8)

colors <- c("U2OS_1h" = "#000000", "U2OS_2h"= "#888888", "U2OS_4h" = "#E0E0E0")

ggplot() + 
  geom_line(data=s_dat$'20210329_U2OS1h-D', aes(x=timePlot_scaled, y = rfu_scaled, color = "U2OS_1h")) +
  geom_line(data=s_dat$'20210329_U2OS2h-D', aes(x=timePlot_scaled, y = rfu_scaled, color = "U2OS_2h")) +
  geom_line(data=s_dat$'20210329_U2OS4h-D', aes(x=timePlot_scaled, y = rfu_scaled, color = "U2OS_4h")) +
  labs(x = "Fragent Size (bp)",
       y = "scaled RFU",
       color = "Legend") +
  ylim(0,1)+
  scale_color_manual(values = colors) +
  theme_classic() + 
  theme(axis.text = element_text(size=22),
        text = element_text(size=22),
        legend.position = "none",
        axis.title = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        axis.line = element_line(linewidth=0.2),
        axis.ticks = element_line(linewidth=0.2),
        axis.text.x=element_text(angle = 45, hjust = 0.6, vjust = 0.9,size=22)) +
  scale_x_continuous(breaks=dat$scaled_time, labels= dat$ladder, limits = c(-0.01, 1.05)) +
  theme(legend.position = "none")


dev.off()



##############  HepG2 timecourse Grey ###########

setwd("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2025_chromatin_hardening/Results/Fig_5/Fragment_analyzer")
pdf("FA_traces_HepG2_timecourse_LMnorm_upperMarkerScale_TimePlotScaling_grey.pdf", height = 6, width = 8)

colors <- c("HepG2_1h_R2"= "#000000", "HepG2_2h_R1"= "#888888", "HepG2_4h_R2"= "#E0E0E0")

ggplot() + 
  geom_line(data=s_dat$'20250624_HepG2_1h_D', aes(x=timePlot_scaled, y = rfu_scaled, color = "HepG2_1h_R2")) +
  geom_line(data=s_dat$'20250624_HepG2_2h_D', aes(x=timePlot_scaled, y = rfu_scaled, color = "HepG2_2h_R1")) +
  geom_line(data=s_dat$'20250624_HepG2_4h_D', aes(x=timePlot_scaled, y = rfu_scaled, color = "HepG2_4h_R2")) +
  labs(x = "Fragent Size (bp)",
       y = "scaled RFU",
       color = "Legend") +
  ylim(0,1)+
  scale_color_manual(values = colors) +
  theme_classic() + 
  theme(axis.text = element_text(size=22),
        text = element_text(size=22),
        legend.position = "none",
        axis.title = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        axis.line = element_line(linewidth=0.2),
        axis.ticks = element_line(linewidth=0.2),
        axis.text.x=element_text(angle = 45, hjust = 0.6, vjust = 0.9,size=22)) +
  scale_x_continuous(breaks=dat$scaled_time, labels= dat$ladder, limits = c(-0.01, 1.05)) +
  theme(legend.position = "none")


dev.off()




##############  H1ESC timecourse R2 ###########

setwd("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2025_chromatin_hardening/Results/Fig_5/Fragment_analyzer")
pdf("FA_traces_H1R2_timecourse_LMnorm_upperMarkerScale_TimePlotScaling.pdf", height = 6, width = 8)

colors <- c("H1_p25_1h"= "#08519c", "H1_p25_2h"= "#3182bd", "H1_p25_4h"= "#bdd7e7")

ggplot() + 
  geom_line(data=s_dat$'20250806_H11hp25-D', aes(x=timePlot_scaled, y = rfu_scaled, color = "H1_p25_1h")) +
  geom_line(data=s_dat$'20250806_H12h-D', aes(x=timePlot_scaled, y = rfu_scaled, color = "H1_p25_2h")) +
  geom_line(data=s_dat$'20250806_H14h-D', aes(x=timePlot_scaled, y = rfu_scaled, color = "H1_p25_4h")) +
  labs(x = "Fragent Size (bp)",
       y = "scaled RFU",
       color = "Legend") +
  ylim(0,1)+
  scale_color_manual(values = colors) +
  theme_classic() + 
  theme(axis.text = element_text(size=22),
        text = element_text(size=22),
        legend.position = "none",
        axis.title = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        axis.line = element_line(linewidth=0.2),
        axis.ticks = element_line(linewidth=0.2),
        axis.text.x=element_text(angle = 45, hjust = 0.6, vjust = 0.9,size=22)) +
  scale_x_continuous(breaks=dat$scaled_time, labels= dat$ladder, limits = c(-0.01, 1.05)) +
  theme(legend.position = "none")


dev.off()



#######
###### Stages only###

setwd("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2025_chromatin_hardening/Results/Fig_5/Fragment_analyzer")
pdf("FA_traces_Stages_only_2h_LMnorm_upperMarkerScale_TimePlotScaling.pdf", height = 6, width = 8)

colors <- c("H1_p25_2h"= "#3182bd", "HB_2h_R2" = "#D4AC0D", "HB_2h_R3" = "#D4AC0D", "mHep_2h_R2" = "#5E4327", "DE_2h_R2" = "#00A878")


ggplot() + 
  geom_line(data=s_dat$'20250806_H12h-D', aes(x=timePlot_scaled, y = rfu_scaled, color = "H1_p25_2h")) +
  geom_line(data=s_dat$'20250324_HB_2h_D', aes(x=timePlot_scaled, y = rfu_scaled, color = "HB_2h_R2")) +
  geom_line(data=s_dat$'20250324_mHep_2h_D', aes(x=timePlot_scaled, y = rfu_scaled, color = "mHep_2h_R2")) +
  geom_line(data=s_dat$'20260616_DER2-D', aes(x=timePlot_scaled, y = rfu_scaled, color = "DE_2h_R2")) +
  
  labs(x = "Fragent Size (bp)",
       y = "scaled RFU",
       color = "Legend") +
  ylim(0,1)+
  scale_color_manual(values = colors) +
  theme_classic() + 
  theme(axis.text = element_text(size=22),
        text = element_text(size=22),
        legend.position = "none",
        axis.title = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        axis.line = element_line(linewidth=0.2),
        axis.ticks = element_line(linewidth=0.2),
        axis.text.x=element_text(angle = 45, hjust = 0.6, vjust = 0.9,size=22)) +
  scale_x_continuous(breaks=dat$scaled_time, labels= dat$ladder, limits = c(-0.01, 1.05)) +
  theme(legend.position = "none")

dev.off()


#######


##############  H1ESC titration ###########

setwd("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2025_chromatin_hardening/Results/Fig_5/Fragment_analyzer")
pdf("FA_traces_H1_titration_LMnorm_upperMarkerScale_TimePlotScaling.pdf", height = 6, width = 8)

colors <- c("H1_p05_1h"= "#08306b", "H1_p1_1h" = "#08519c", "H1_p15_1h"= "#2171b5", "H1_p25_1h"= "#6baed6")

ggplot() + 
  geom_line(data=s_dat$'20250806_H11hp05-D', aes(x=timePlot_scaled, y = rfu_scaled, color = "H1_p05_1h")) +
  geom_line(data=s_dat$'20250806_H11hp1-D', aes(x=timePlot_scaled, y = rfu_scaled, color = "H1_p1_1h")) +
  geom_line(data=s_dat$'20250806_H11hp15-D', aes(x=timePlot_scaled, y = rfu_scaled, color = "H1_p15_1h")) +
  geom_line(data=s_dat$'20250806_H11hp25-D', aes(x=timePlot_scaled, y = rfu_scaled, color = "H1_p25_1h")) +
  labs(x = "Fragent Size (bp)",
       y = "scaled RFU",
       color = "Legend") +
  ylim(0,1)+
  scale_color_manual(values = colors) +
  theme_classic() + 
  theme(axis.text = element_text(size=22),
        text = element_text(size=22),
        legend.position = "none",
        axis.title = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        axis.line = element_line(linewidth=0.2),
        axis.ticks = element_line(linewidth=0.2),
        axis.text.x=element_text(angle = 45, hjust = 0.6, vjust = 0.9,size=22)) +
  scale_x_continuous(breaks=dat$scaled_time, labels= dat$ladder, limits = c(-0.01, 1.05)) +
  theme(legend.position = "none")

dev.off()



###### DE vs comparable H1 ###

setwd("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2025_chromatin_hardening/Results/Fig_5/Fragment_analyzer")
pdf("FA_traces_H1_comparable_vs_DE_LMnorm_upperMarkerScale_TimePlotScaling.pdf", height = 6, width = 8)

colors <- c("H1_p25_1h"= "#2C6FB5", "DE_2h_R2" = "#00A878")


ggplot() + 
  geom_line(data=s_dat$'20250806_H11hp25-D', aes(x=timePlot_scaled, y = rfu_scaled, color = "H1_p25_1h")) +
  geom_line(data=s_dat$'20260616_DER2-D', aes(x=timePlot_scaled, y = rfu_scaled, color = "DE_2h_R2")) +
  
  labs(x = "Fragent Size (bp)",
       y = "scaled RFU",
       color = "Legend") +
  ylim(0,1)+
  scale_color_manual(values = colors) +
  theme_classic() + 
  theme(axis.text = element_text(size=22),
        text = element_text(size=22),
        legend.position = "none",
        axis.title = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        axis.line = element_line(linewidth=0.2),
        axis.ticks = element_line(linewidth=0.2),
        axis.text.x=element_text(angle = 45, hjust = 0.6, vjust = 0.9,size=22)) +
  scale_x_continuous(breaks=dat$scaled_time, labels= dat$ladder, limits = c(-0.01, 1.05)) +
  theme(legend.position = "none")

dev.off()




###### HB vs comparable H1 ###

setwd("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2025_chromatin_hardening/Results/Fig_5/Fragment_analyzer")
pdf("FA_traces_H1_comparable_vs_HB_LMnorm_upperMarkerScale_TimePlotScaling.pdf", height = 6, width = 8)

colors <- c("H1_p1_1h" = "#2C6FB5", "HB_2h_R2" = "#D4AC0D")


ggplot() + 
  geom_line(data=s_dat$'20250806_H11hp1-D', aes(x=timePlot_scaled, y = rfu_scaled, color = "H1_p1_1h")) +
  geom_line(data=s_dat$'20250324_HB_2h_D', aes(x=timePlot_scaled, y = rfu_scaled, color = "HB_2h_R2")) +
  labs(x = "Fragent Size (bp)",
       y = "scaled RFU",
       color = "Legend") +
  ylim(0,1)+
  scale_color_manual(values = colors) +
  theme_classic() + 
  theme(axis.text = element_text(size=22),
        text = element_text(size=22),
        legend.position = "none",
        axis.title = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        axis.line = element_line(linewidth=0.2),
        axis.ticks = element_line(linewidth=0.2),
        axis.text.x=element_text(angle = 45, hjust = 0.6, vjust = 0.9,size=22)) +
  scale_x_continuous(breaks=dat$scaled_time, labels= dat$ladder, limits = c(-0.01, 1.05)) +
  theme(legend.position = "none")

dev.off()

#######

###### HLC vs comparable H1 ###

setwd("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2025_chromatin_hardening/Results/Fig_5/Fragment_analyzer")
pdf("FA_traces_H1_comparable_vs_HLC_LMnorm_upperMarkerScale_TimePlotScaling.pdf", height = 6, width = 8)

colors <- c("H1_p05_1h"= "#2C6FB5", "mHep_2h_R2" = "#5E4327")


ggplot() + 
  geom_line(data=s_dat$'20250806_H11hp05-D', aes(x=timePlot_scaled, y = rfu_scaled, color = "H1_p05_1h")) +
  geom_line(data=s_dat$'20250324_mHep_2h_D', aes(x=timePlot_scaled, y = rfu_scaled, color = "mHep_2h_R2")) +
  labs(x = "Fragent Size (bp)",
       y = "scaled RFU",
       color = "Legend") +
  ylim(0,1)+
  scale_color_manual(values = colors) +
  theme_classic() + 
  theme(axis.text = element_text(size=22),
        text = element_text(size=22),
        legend.position = "none",
        axis.title = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        axis.line = element_line(linewidth=0.2),
        axis.ticks = element_line(linewidth=0.2),
        axis.text.x=element_text(angle = 45, hjust = 0.6, vjust = 0.9,size=22)) +
  scale_x_continuous(breaks=dat$scaled_time, labels= dat$ladder, limits = c(-0.01, 1.05)) +
  theme(legend.position = "none")

dev.off()

#######



