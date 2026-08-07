library(tidyverse)


df_residuals <- read.table("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2025_chromatin_hardening/analysis/dataframes/df_residuals_20260806_50kb.bed",
                sep="\t", header=TRUE)


######  Convert low DpnII-seq values to NA
df_residuals <- df_residuals %>%
  mutate(across(contains("signal"), 
                ~replace(., . %in% 0:20, NA)))


df_residuals_150kb <- read.table("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2025_chromatin_hardening/analysis/dataframes/df_residuals_20260806_150kb.bed",
                           sep="\t", header=TRUE)

ipt_150kb <- read.table("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2025_chromatin_hardening/analysis/data/IPG_IPT/IPT_150kb_rebinned.bed",
                        sep="\t", header=TRUE)


ipt_levels <- c("A1", "A2", "A3", "B1", "Quies", "Inactive", "B4")

df_residuals_150kb <- df_residuals_150kb %>%
  mutate(ipt = factor(ipt, levels = ipt_levels))



##### Adding segment columns for tracks #####

###  IPTs
df_residuals$seg_1 <- -1
df_residuals <- cbind(df_residuals, replicate(7, df_residuals$seg_1))
n <- ncol(df_residuals)
names(df_residuals)[(n-6):n] <- c("A1_ipt_seg", "A2_ipt_seg", "A3_ipt_seg", "B1_ipt_seg", "Quies_ipt_seg","Inactive_ipt_seg","B4_ipt_seg")


df_residuals$A1_ipt_seg <- ifelse(df_residuals$ipt == "A1",df_residuals$A1_ipt_seg, 10)
df_residuals$A2_ipt_seg <- ifelse(df_residuals$ipt == "A2",df_residuals$A2_ipt_seg, 10)
df_residuals$A3_ipt_seg <- ifelse(df_residuals$ipt == "A3",df_residuals$A3_ipt_seg, 10)
df_residuals$B1_ipt_seg <- ifelse(df_residuals$ipt == "B1",df_residuals$B1_ipt_seg, 10)
df_residuals$Quies_ipt_seg <- ifelse(df_residuals$ipt == "Quies",df_residuals$Quies_ipt_seg, 10)
df_residuals$Inactive_ipt_seg <- ifelse(df_residuals$ipt == "Inactive",df_residuals$Inactive_ipt_seg, 10)
df_residuals$B4_ipt_seg <- ifelse(df_residuals$ipt == "B4",df_residuals$B4_ipt_seg, 10)



######### Add compartment calls for each stage  ###########
df_residuals$H1_comp <- ifelse(df_residuals$PC1_H130 > 0,"A","B") 
df_residuals$DE_comp <- ifelse(df_residuals$PC1_DE30 > 0,"A","B")
df_residuals$HB_comp <- ifelse(df_residuals$PC1_HB30 > 0,"A","B") 
df_residuals$mHep_comp <- ifelse(df_residuals$PC1_mHep30 > 0,"A","B") 

###  H1 compartment segs
df_residuals$seg_3 <- -3
df_residuals <- cbind(df_residuals, replicate(2, df_residuals$seg_3))
n <- ncol(df_residuals)
names(df_residuals)[(n-1):n] <- c("A_H1_seg", "B_H1_seg")

df_residuals$A_H1_seg <- ifelse(df_residuals$H1_comp == "A",df_residuals$A_H1_seg, 10)
df_residuals$B_H1_seg <- ifelse(df_residuals$H1_comp == "B",df_residuals$B_H1_seg, 10)

###  DE compartment segs
df_residuals$seg_3 <- -3
df_residuals <- cbind(df_residuals, replicate(2, df_residuals$seg_3))
n <- ncol(df_residuals)
names(df_residuals)[(n-1):n] <- c("A_DE_seg", "B_DE_seg")

df_residuals$A_DE_seg <- ifelse(df_residuals$DE_comp == "A",df_residuals$A_DE_seg, 10)
df_residuals$B_DE_seg <- ifelse(df_residuals$DE_comp == "B",df_residuals$B_DE_seg, 10)

###  HB compartment segs
df_residuals$seg_3 <- -3
df_residuals <- cbind(df_residuals, replicate(2, df_residuals$seg_3))
n <- ncol(df_residuals)
names(df_residuals)[(n-1):n] <- c("A_HB_seg", "B_HB_seg")

df_residuals$A_HB_seg <- ifelse(df_residuals$HB_comp == "A",df_residuals$A_HB_seg, 10)
df_residuals$B_HB_seg <- ifelse(df_residuals$HB_comp == "B",df_residuals$B_HB_seg, 10)

###  mHep compartment segs
df_residuals$seg_3 <- -3
df_residuals <- cbind(df_residuals, replicate(2, df_residuals$seg_3))
n <- ncol(df_residuals)
names(df_residuals)[(n-1):n] <- c("A_mHep_seg", "B_mHep_seg")

df_residuals$A_mHep_seg <- ifelse(df_residuals$mHep_comp == "A",df_residuals$A_mHep_seg, 10)
df_residuals$B_mHep_seg <- ifelse(df_residuals$mHep_comp == "B",df_residuals$B_mHep_seg, 10)



##################################### PLOT #############################################

#### Fig. 5 IPT and H1comp bars

starts <- 0
ends <- 120000000
seg_size <- 9



setwd("C:/MyPath/")


x <-df_residuals[df_residuals$chrom == "chr1", ]
x <- x[x$start>starts & x$end<ends,]

pdf("H1_NICEeigenBars_IPTsegs_chr1p.pdf", height = 1.5, width = 5)
par(mfrow=c(1,1), oma=c(1,0,1,0), mar=c(0.5, 4, 0.5, 0) + 0.1, tcl=-0.25)

plot(x$start/1000000, x$seg_1, type="n",
     xlab = "", ylab= "", axes=FALSE, ylim=c(-4,1))

A1 <- x$A1_ipt_seg
segments(x0 = x$end/1000000,                   
         y0 = A1,
         x1 = x$start/1000000,
         y1 = A1,
         lend=1,lwd=seg_size, col="#e23838")


A2 <- x$A2_ipt_seg
segments(x0 = x$end/1000000,                   
         y0 = A2,
         x1 = x$start/1000000,
         y1 = A2,
         lend=1,lwd=seg_size, col="#f78200")

A3 <- x$A3_ipt_seg
segments(x0 = x$end/1000000,                   
         y0 = A3,
         x1 = x$start/1000000,
         y1 = A3,
         lend=1,lwd=seg_size, col="#ffb900")

B1 <- x$B1_ipt_seg
segments(x0 = x$end/1000000,                   
         y0 = B1,
         x1 = x$start/1000000,
         y1 = B1,
         lend=1,lwd=seg_size, col="#009cdf")

Quies <- x$Quies_ipt_seg
segments(x0 = x$end/1000000,                   
         y0 = Quies,
         x1 = x$start/1000000,
         y1 = Quies,
         lend=1,lwd=seg_size, col="#B0C4DE")

Inactive <- x$Inactive_ipt_seg
segments(x0 = x$end/1000000,                   
         y0 = Inactive,
         x1 = x$start/1000000,
         y1 = Inactive,
         lend=1,lwd=seg_size, col="#A9A9A9")

B4 <- x$B4_ipt_seg
segments(x0 = x$end/1000000,                   
         y0 = B4,
         x1 = x$start/1000000,
         y1 = B4,
         lend=1,lwd=seg_size, col="#973999")


A <- x$A_H1_seg
segments(x0 = x$end/1000000,                   
         y0 = A,
         x1 = x$start/1000000,
         y1 = A,
         lend=1,lwd=seg_size, col="red")


B <- x$B_H1_seg
segments(x0 = x$end/1000000,                   
         y0 = B,
         x1 = x$start/1000000,
         y1 = B,
         lend=1,lwd=seg_size, col="blue")

axis(1, lwd=2, cex.axis=1, labels=FALSE) 

dev.off()


##### DE compartment bars ####

x <-df_residuals[df_residuals$chrom == "chr1", ]
x <- x[x$start>starts & x$end<ends,]

pdf("DE_NICEeigenBars_chr1p.pdf", height = 1.5, width = 5)
par(mfrow=c(1,1), oma=c(1,0,1,0), mar=c(0.5, 4, 0.5, 0) + 0.1, tcl=-0.25)

plot(x$start/1000000, x$seg_1, type="n",
     xlab = "", ylab= "", axes=FALSE, ylim=c(-4,1))

A <- x$A_DE_seg
segments(x0 = x$end/1000000,                   
         y0 = A,
         x1 = x$start/1000000,
         y1 = A,
         lend=1,lwd=seg_size, col="red")


B <- x$B_DE_seg
segments(x0 = x$end/1000000,                   
         y0 = B,
         x1 = x$start/1000000,
         y1 = B,
         lend=1,lwd=seg_size, col="blue")

axis(1, lwd=2, cex.axis=1, labels=FALSE) 

dev.off()


##### HB compartment bars ####

x <-df_residuals[df_residuals$chrom == "chr1", ]
x <- x[x$start>starts & x$end<ends,]

pdf("HB_NICEeigenBars_chr1p.pdf", height = 1.5, width = 5)
par(mfrow=c(1,1), oma=c(1,0,1,0), mar=c(0.5, 4, 0.5, 0) + 0.1, tcl=-0.25)

plot(x$start/1000000, x$seg_1, type="n",
     xlab = "", ylab= "", axes=FALSE, ylim=c(-4,1))

A <- x$A_HB_seg
segments(x0 = x$end/1000000,                   
         y0 = A,
         x1 = x$start/1000000,
         y1 = A,
         lend=1,lwd=seg_size, col="red")


B <- x$B_HB_seg
segments(x0 = x$end/1000000,                   
         y0 = B,
         x1 = x$start/1000000,
         y1 = B,
         lend=1,lwd=seg_size, col="blue")

axis(1, lwd=2, cex.axis=1, labels=FALSE) 

dev.off()


##### mHep compartment bars ####

x <-df_residuals[df_residuals$chrom == "chr1", ]
x <- x[x$start>starts & x$end<ends,]

pdf("mHep_NICEeigenBars_chr1p.pdf", height = 1.5, width = 5)
par(mfrow=c(1,1), oma=c(1,0,1,0), mar=c(0.5, 4, 0.5, 0) + 0.1, tcl=-0.25)

plot(x$start/1000000, x$seg_1, type="n",
     xlab = "", ylab= "", axes=FALSE, ylim=c(-4,1))

A <- x$A_mHep_seg
segments(x0 = x$end/1000000,                   
         y0 = A,
         x1 = x$start/1000000,
         y1 = A,
         lend=1,lwd=seg_size, col="red")


B <- x$B_mHep_seg
segments(x0 = x$end/1000000,                   
         y0 = B,
         x1 = x$start/1000000,
         y1 = B,
         lend=1,lwd=seg_size, col="blue")

axis(1, lwd=2, cex.axis=1, labels=FALSE) 

dev.off()
###########


##### LOS tracks ESC vs DE

x <-df_residuals[df_residuals$chrom == "chr1", ]
x <- x[x$start>starts & x$end<ends,]
pdf("H11h_vs_DER2_LOS_chr1_zoom_multiLine.pdf", height = 1.5, width = 5)
par(mfrow=c(1,1), oma=c(1,0,1,0), mar=c(0.5, 4, 0.5, 0) + 0.1, tcl=-0.25)
plot(x$start/1000000, x$LOS_range2Mb_DE_120mPD_DpnII_R2_20260616, type="l", col="#00A878", 
     ylim=c(0.73, 1), xlab = "", ylab= "", axes=FALSE, lwd=0.5)
lines(x$start/1000000, x$LOS_range2Mb_H1ESC_60mPDp25_DpnII_R3_20250820,
      col = "#2C6FB5",  
      lwd = 0.5)
axis(1, lwd=0.7, cex.axis=0.9, labels=TRUE)
axis(2, lwd=0.7, cex.axis=0.9)

dev.off()



##### LOS tracks ESC vs HB
x <-df_residuals[df_residuals$chrom == "chr1", ]
x <- x[x$start>starts & x$end<ends,]
pdf("H1p1_vs_HBR2_LOS_chr1_zoom_multiLine.pdf", height = 1.5, width = 5)
par(mfrow=c(1,1), oma=c(1,0,1,0), mar=c(0.5, 4, 0.5, 0) + 0.1, tcl=-0.25)
plot(x$start/1000000, x$LOS_range2Mb_H1ESC_60mPDp10_DpnII_R3_20250820, type="l", col="#2C6FB5", 
     ylim=c(0.73, 1), xlab = "", ylab= "LOS", axes=FALSE, lwd=0.5)
lines(x$start/1000000, x$LOS_range2Mb_HB_2hPD_DpnII_R2_20240326,
      col = "#D4AC0D",         
      lwd = 0.5)
axis(1, lwd=0.7, cex.axis=0.9, labels=TRUE)
axis(2, lwd=0.7, cex.axis=0.9)

dev.off()




##### LOS tracks ESC vs HLC
x <-df_residuals[df_residuals$chrom == "chr1", ]
x <- x[x$start>starts & x$end<ends,]
pdf("H1p05_vs_HLCR1_LOS_chr1_zoom_multiLine.pdf", height = 1.5, width = 5)
par(mfrow=c(1,1), oma=c(1,0,1,0), mar=c(0.5, 4, 0.5, 0) + 0.1, tcl=-0.25)
plot(x$start/1000000, x$LOS_range2Mb_H1ESC_60mPDp05_DpnII_R3_20250820, type="l", col="#2C6FB5", 
     ylim=c(0.73, 1), xlab = "", ylab= "LOS", axes=FALSE, lwd=0.5)
lines(x$start/1000000, x$LOS_range2Mb_mHep_2hPD_DpnII_R2_20240326,
      col = "#5E4327",        
      lwd = 0.5)
axis(1, lwd=0.7, cex.axis=0.9, labels=TRUE)
axis(2, lwd=0.7, cex.axis=0.9)

dev.off()



######### Dodged LOS violins by stage by IPT - H1 vs DE ##########

setwd("C:/MyPath/")

LOS_cols <- c('LOS_range2Mb_H1ESC_60mPDp25_DpnII_R3_20250820',
              'LOS_range2Mb_DE_120mPD_DpnII_R2_20260616')

NT <- df_residuals_150kb[, c('chrom','start','end','ipt', LOS_cols)]
NT <- NT[!is.na(NT$ipt), ]


NT_long <- NT %>%
  pivot_longer(cols = all_of(LOS_cols),
               names_to = "sample", values_to = "LOS") %>%
  filter(!is.na(LOS))


pdf("LOS_by_ipt_violin_combinedH1vsDE.pdf", width=15, height=6)
ggplot(NT_long, aes(x = ipt, y = LOS, fill = sample)) +
  geom_violin(trim=FALSE, position=position_dodge(0.65), linewidth=0.15) +
  geom_boxplot(aes(group = interaction(ipt, sample)),
               width=0.1, fill="white", linewidth=0.15,
               position=position_dodge(0.65)) +
  scale_fill_manual(values = c("#2C6FB5","#00A878"),
                    labels = c("H1ESC 60mPD", "DE 2hPD")) +
  ylim(0.5, 0.97) +
  geom_hline(yintercept = 0, color="black", linetype=2, linewidth=0.3) +
  theme_classic() +
  theme(axis.text = element_text(size=38),
        text = element_text(size=38),
        legend.position = "none",
        axis.title = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        axis.line = element_line(linewidth=0.2),
        axis.ticks = element_line(linewidth=0.2),
        axis.text.x = element_blank())
dev.off()



######### Dodged LOS violins by stage by IPT - H1 vs HB ##########

LOS_cols <- c('LOS_range2Mb_H1ESC_60mPDp10_DpnII_R3_20250820',
              'LOS_range2Mb_HB_2hPD_DpnII_R2_20240326')

NT <- df_residuals_150kb[, c('chrom','start','end','ipt', LOS_cols)]
NT <- NT[!is.na(NT$ipt), ]


NT_long <- NT %>%
  pivot_longer(cols = all_of(LOS_cols),
               names_to = "sample", values_to = "LOS") %>%
  filter(!is.na(LOS))

pdf("LOS_by_ipt_violin_combinedH1vsHB.pdf", width=15, height=6)
ggplot(NT_long, aes(x = ipt, y = LOS, fill = sample)) +
  geom_violin(trim=FALSE, position=position_dodge(0.65), linewidth=0.15) +
  geom_boxplot(aes(group = interaction(ipt, sample)),
               width=0.1, fill="white", linewidth=0.15,
               position=position_dodge(0.65)) +
  scale_fill_manual(values = c("#2C6FB5", "#D4AC0D"),
                    labels = c("H1ESC 60mPD","HB 2hPD")) +
  ylim(0.5, 0.97) +
  geom_hline(yintercept = 0, color="black", linetype=2, linewidth=0.3) +
  theme_classic() +
  theme(axis.text = element_text(size=38),
        text = element_text(size=38),
        legend.position = "none",
        axis.title = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        axis.line = element_line(linewidth=0.2),
        axis.ticks = element_line(linewidth=0.2),
        axis.text.x = element_blank())
dev.off()




######### Dodged LOS violins by stage by IPT - H1 vs HLC ##########

LOS_cols <- c('LOS_range2Mb_H1ESC_60mPDp10_DpnII_R3_20250820',
              'LOS_range2Mb_mHep_2hPD_DpnII_R2_20240326')

NT <- df_residuals_150kb[, c('chrom','start','end','ipt', LOS_cols)]
NT <- NT[!is.na(NT$ipt), ]


NT_long <- NT %>%
  pivot_longer(cols = all_of(LOS_cols),
               names_to = "sample", values_to = "LOS") %>%
  filter(!is.na(LOS))

pdf("LOS_by_ipt_violin_combinedH1vsHLC.pdf", width=15, height=6)
ggplot(NT_long, aes(x = ipt, y = LOS, fill = sample)) +
  geom_violin(trim=FALSE, position=position_dodge(0.65), linewidth=0.15) +
  geom_boxplot(aes(group = interaction(ipt, sample)),
               width=0.1, fill="white", linewidth=0.15,
               position=position_dodge(0.65)) +
  scale_fill_manual(values = c("#2C6FB5", "#5E4327"),
                    labels = c("H1ESC 60mPD","HB 2hPD")) +
  ylim(0.5, 0.97) +
  geom_hline(yintercept = 0, color="black", linetype=2, linewidth=0.3) +
  theme_classic() +
  theme(axis.text = element_text(size=38),
        text = element_text(size=38),
        legend.position = "none",
        axis.title = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        axis.line = element_line(linewidth=0.2),
        axis.ticks = element_line(linewidth=0.2),
        axis.text.x = element_blank())
dev.off()




######### Dodged LOSres violins by stage by IPT - H1 vs DE ##########


LOS_cols <- c('LOS_residuals_range2Mb_H1ESC_60mPDp25_DpnII_R3_20250820',
              'LOS_residuals_range2Mb_DE_120mPD_DpnII_R2_20260616')

NT <- df_residuals_150kb[, c('chrom','start','end','ipt', LOS_cols)]
NT <- NT[!is.na(NT$ipt), ]


NT_long <- NT %>%
  pivot_longer(cols = all_of(LOS_cols),
               names_to = "sample", values_to = "LOS") %>%
  filter(!is.na(LOS))

pdf("LOSres_by_ipt_violin_combinedH1vsDE.pdf", width=15, height=6)
ggplot(NT_long, aes(x = ipt, y = LOS, fill = sample)) +
  geom_violin(trim=FALSE, position=position_dodge(0.65), linewidth=0.15) +
  geom_boxplot(aes(group = interaction(ipt, sample)),
               width=0.1, fill="white", linewidth=0.15,
               position=position_dodge(0.65)) +
  scale_fill_manual(values = c("#2C6FB5","#00A878"),
                    labels = c("H1ESC 60mPD", "DE 2hPD")) +
  ylim(-0.3, 0.1) +
  geom_hline(yintercept = 0, color="black", linetype=2, linewidth=0.3) +
  theme_classic() +
  theme(axis.text = element_text(size=38),
        text = element_text(size=38),
        legend.position = "none",
        axis.title = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        axis.line = element_line(linewidth=0.2),
        axis.ticks = element_line(linewidth=0.2),
        axis.text.x = element_blank())
dev.off()



######### Dodged LOSres violins by stage by IPT - H1 vs HB ##########

LOS_cols <- c('LOS_residuals_range2Mb_H1ESC_60mPDp10_DpnII_R3_20250820',
              'LOS_residuals_range2Mb_HB_2hPD_DpnII_R2_20240326')

NT <- df_residuals_150kb[, c('chrom','start','end','ipt', LOS_cols)]
NT <- NT[!is.na(NT$ipt), ]


NT_long <- NT %>%
  pivot_longer(cols = all_of(LOS_cols),
               names_to = "sample", values_to = "LOS") %>%
  filter(!is.na(LOS))

pdf("LOSres_by_ipt_violin_combinedH1vsHB.pdf", width=15, height=6)
ggplot(NT_long, aes(x = ipt, y = LOS, fill = sample)) +
  geom_violin(trim=FALSE, position=position_dodge(0.65), linewidth=0.15) +
  geom_boxplot(aes(group = interaction(ipt, sample)),
               width=0.1, fill="white", linewidth=0.15,
               position=position_dodge(0.65)) +
  scale_fill_manual(values = c("#2C6FB5", "#D4AC0D"),
                    labels = c("H1ESC 60mPD","HB 2hPD")) +
  ylim(-0.3, 0.1) +
  geom_hline(yintercept = 0, color="black", linetype=2, linewidth=0.3) +
  theme_classic() +
  theme(axis.text = element_text(size=38),
        text = element_text(size=38),
        legend.position = "none",
        axis.title = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        axis.line = element_line(linewidth=0.2),
        axis.ticks = element_line(linewidth=0.2),
        axis.text.x = element_blank())
dev.off()




######### Dodged LOSres violins by stage by IPT - H1 vs HLC ##########

LOS_cols <- c('LOS_residuals_range2Mb_H1ESC_60mPDp10_DpnII_R3_20250820',
              'LOS_residuals_range2Mb_mHep_2hPD_DpnII_R2_20240326')

NT <- df_residuals_150kb[, c('chrom','start','end','ipt', LOS_cols)]
NT <- NT[!is.na(NT$ipt), ]


NT_long <- NT %>%
  pivot_longer(cols = all_of(LOS_cols),
               names_to = "sample", values_to = "LOS") %>%
  filter(!is.na(LOS))

pdf("LOSres_by_ipt_violin_combinedH1vsHLC.pdf", width=15, height=6)
ggplot(NT_long, aes(x = ipt, y = LOS, fill = sample)) +
  geom_violin(trim=FALSE, position=position_dodge(0.65), linewidth=0.15) +
  geom_boxplot(aes(group = interaction(ipt, sample)),
               width=0.1, fill="white", linewidth=0.15,
               position=position_dodge(0.65)) +
  scale_fill_manual(values = c("#2C6FB5", "#5E4327"),
                    labels = c("H1ESC 60mPD","HB 2hPD")) +
  ylim(-0.3, 0.1) +
  geom_hline(yintercept = 0, color="black", linetype=2, linewidth=0.3) +
  theme_classic() +
  theme(axis.text = element_text(size=38),
        text = element_text(size=38),
        legend.position = "none",
        axis.title = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        axis.line = element_line(linewidth=0.2),
        axis.ticks = element_line(linewidth=0.2),
        axis.text.x = element_blank())
dev.off()
#####



