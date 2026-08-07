library(tidyverse)
library(gridExtra)



df_residuals <- read.table("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2025_chromatin_hardening/analysis/dataframes/df_residuals_20260806_50kb.bed",
                sep="\t", header=TRUE)


######  Convert low DpnII-seq values to NA
df_residuals <- df_residuals %>%
  mutate(across(contains("signal"), 
                ~replace(., . %in% 0:20, NA)))

###### Remove SPIN and IPG columns
df_residuals[c("ipg", "SPIN")] <- list(NULL)



###############  Replicate scatters LOS

panels <- data.frame(
  xcol = c("LOS_range2Mb_H1ESC_60mPDp15_DpnII_R3_20250820", "LOS_range2Mb_DE_120mPD_DpnII_R2_20260616", "LOS_range2Mb_HB_2hPD_DpnII_R2_20240326","LOS_range2Mb_mHep_2hPD_DpnII_R2_20240326"),
  ycol = c("LOS_range2Mb_H1ESC_nPI120mPD_DpnII_R1_20240125", "LOS_range2Mb_DE_120mPD_DpnII_R1_20260616", "LOS_range2Mb_HB_2hPD_DpnII_R3_20250624","LOS_range2Mb_mHep_2hPD_DpnII_R2_20250624"),
  stringsAsFactors = FALSE
)

plots <- Map(function(xcol, ycol) {
  ggplot(df_residuals, aes(.data[[xcol]], .data[[ycol]])) +
    geom_point(size = 0.25) +
    labs(x = xcol, y = ycol) +
    theme_classic()
}, panels$xcol, panels$ycol)

setwd("C:/MyPath")
pdf("Replicate_scatters.pdf", height = 4, width = 20)
grid.arrange(grobs = plots, ncol = 5)
dev.off()

