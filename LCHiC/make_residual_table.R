library(stringr)
library(yaml)




##### Load and pre-process tables
resolution <- "1000kb"

data_folder <- "C:/MyDataFolderPath/"


bins <- read.table(paste0("C:/MyPath/hg38_", resolution,".bed"),
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end"))
bins <- bins[ - grep("chrY", bins$chrom),]
bins <- bins[ - grep("chrM", bins$chrom),]

bins$chrom <- factor(bins$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))
bins <- bins[order(bins$chrom, bins$start,bins$end),]
rownames(bins) <- NULL

IPT <- read.table(paste("C:/MyPath/",resolution,"/IPT/IPT_", resolution,"_rebinned.bed", sep=""),
                  sep="\t", header=TRUE)# col.names=c("chrom", "start", "end", "IPG"))
IPT$ipt   <- replace(IPT$ipt, IPT$ipt == '<NA>', NA)


##### Functions to extract column name from filename (would have to adapt theses functions for altered file name)

col_name_LOS <- function(x){
  a <- unlist(str_split(x,"/"))  
  b <- tail(a,1)
  LOS_name <- paste("LOS_",unlist(str_split(b,"_"))[4], "_", unlist(str_split(b,"-"))[4],"_",unlist(str_split(b,"-"))[5],"_",unlist(str_split(b,"-"))[6],"_", unlist(str_split(b,"-"))[7],"_",unlist(str_split(b,"-"))[2], sep="")
}

col_name_DpnIIseq <- function(x){
  a <- unlist(str_split(x,"/"))  
  b <- tail(a,1)
  DpnII_name <- paste("signal_", unlist(str_split(b,"-"))[4],"_",unlist(str_split(b,"-"))[5],"_",unlist(str_split(b,"-"))[6],"_",unlist(str_split(b,"-"))[2],sep = "")
}

col_name_RNAseq <- function(x){
  a <- unlist(str_split(x,"/"))  
  b <- tail(a,1)
  DpnII_name <- paste("CPM_", unlist(str_split(b,"_"))[1],"_",unlist(str_split(b,"_"))[2],sep = "")
}

##### Core Functions

get_moving_average <- function(x, y, n, s, dec) {
  ## Args:
  #  x: unordered numeric vector to slide windows across
  #  y: unordered numeric values to mean over x window
  #  n: window size
  #  s: step size
  #  dec: number of decimal places to round for merge
  # Returns:
  # df: dataframe with columns:
  #   w: window centers
  #   mu: moving average
  if ((sum(is.na(x)) > 0) | (sum(is.na(y)) > 0)) {
    stop("Remove NAs")
  }
  if (length(x)!=length(y)) {
    stop("Unequal vector lengths")
  }
  xmin <- round(min(x, na.rm=TRUE), dec)
  xmax <- round(max(x, na.rm=TRUE), dec)
  w <- seq(xmin, xmax, by=s)
  mu <- c()
  for (i in w) {
    m <- mean(y[x > (i-(n/2)) & x < (i + (n/2))])
    mu <- c(mu, m)
  }
  df <- data.frame(w, mu)
  return(df)
}

get_ma_residuals <- function(o, df, dec) {
  # Get moving average residuals
  ## Args:
  #  o: original dataframe with cols:
  #   column1: chrom
  #   column2: start
  #   column3: end
  #   column4: numeric vector original x-axis
  #   column5:  numeric vector original y-axis
  #  df: dataframe output by get_moving_average
  #      function
  #  dec: number of decimal places to round for merge
  # Returns:
  #  r: dataframe of residuals and locations
  colnames(o) <- c("chrom", "start", "end", "w", "y")
  o["w"] <- round(o["w"], dec)
  df["w"] <- round(df["w"], dec)
  m <- merge(o, df, by="w")
  m_resid <- m$y-m$mu
  r <- cbind(m[c("chrom", "start", "end")], m_resid)
  return(r)
}

##### list data from each folder (cis_coverage, cis_percent, LOS, smooth_LOS)


LOSList <- list.files(paste(data_folder, "/LOS", sep =""), pattern = "LOS", full.name = TRUE)

DpnIIList <- list.files(paste(data_folder,"DpnII_seq", sep=""), pattern = "copy_correct_coverage_", full.name = TRUE)

RNAseqList <- list.files(paste(data_folder,"rnaseq", sep=""), pattern = ".cpm", full.name = TRUE)



##### Load and munge LOS tables

LOS_names <- lapply(LOSList, col_name_LOS)
LOS_names <- unlist(LOS_names) 
LOS_Data <- lapply(1:length(LOSList), function(x) { read.delim(LOSList[[x]], header = TRUE)})


for (i in c(1:length(LOS_Data))) {
  names(LOS_Data[[i]]) <- c("chrom","start","end",LOS_names[[i]])
}
names(LOS_Data) <- c(LOS_names)

for (i in c(1:length(LOS_Data))) {
  LOS_Data[[i]][2] <- LOS_Data[[i]][2] + 1
}

for (i in c(1:length(LOS_Data))) {
  LOS_Data[[i]][!grepl("chrY", LOS_Data[[i]]$chrom),]
}


##### Load and munge DpnII-seq tables
signal_names <- lapply(DpnIIList, col_name_DpnIIseq)
signal_names <- unlist(signal_names)
signal_Data <- lapply(1:length(DpnIIList), function(x) { read.delim(DpnIIList[[x]], header = FALSE,col.names= c("chrom", "start", "end", "signal"))})
for (i in c(1:length(signal_Data))) {
  names(signal_Data[[i]]) <- c("chrom","start","end",signal_names[[i]])
}
names(signal_Data) <- c(signal_names)



##### Merge LOS and DpnII-seq data 

d <- data.frame(bins[,1:3], IPT["ipt"])


for (i in c(1:length(signal_Data))) {
  d <- merge(d, signal_Data[[i]], by=c("chrom", "start", "end"), all.x= TRUE)
}

for (i in c(1:length(LOS_Data))) {
  d <- merge(d, LOS_Data[[i]], by=c("chrom", "start", "end"), all.x= TRUE)
}



d$chrom <- factor(d$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))
d <- d[order(d$chrom, d$start,d$end),]
rownames(d) <- NULL

##### Calculate LOS residuals, merge with LOS table and create individual LOS residual bed files
df_residuals <- d

# Read the YAML configuration file
config <- yaml.load_file("C:/myPath/LOS_signal_pairs_1000kb.yaml")
setwd("C:/MyPath/")
# Loop through the key-value pairs in the YAML file
for (pair in names(config)) {
  LOS_name <- config[[pair]]$LOS
  signal_name <- config[[pair]]$signal

  d_clean <-  na.omit(df_residuals[, c("chrom", "start", "end", LOS_name, signal_name)])
  
  los_sample <- d_clean[,LOS_name]
  signal_sample <- d_clean[,signal_name]

  n <- mean(na.omit(d_clean[,signal_name]))/4
  
  ma <- get_moving_average(signal_sample, los_sample, n, 1, 0)
  los_r <- get_ma_residuals(d_clean[c("chrom", "start", "end", signal_name, LOS_name)], ma, 0)
  names(los_r)[4] <- paste("LOS_residuals_", unlist(str_split(LOS_name,"_"))[2], "_", unlist(str_split(LOS_name,"_"))[3], "_",unlist(str_split(LOS_name,"_"))[4], "_",unlist(str_split(LOS_name,"_"))[5], "_",unlist(str_split(LOS_name,"_"))[6], "_",tail(unlist(str_split(LOS_name,"_"))[-1], n=1), sep="")
  
  df_residuals <- merge(df_residuals,los_r, by=c("chrom", "start", "end"), all.x= TRUE)
  los_r_clean <- merge(bins, los_r, by=c("chrom", "start", "end"), all.x= TRUE)
  write.table(los_r_clean, paste("LOS_residuals_", unlist(str_split(LOS_name,"_"))[2], "_", unlist(str_split(LOS_name,"_"))[3], "_",unlist(str_split(LOS_name,"_"))[4], "_",unlist(str_split(LOS_name,"_"))[5], "_",unlist(str_split(LOS_name,"_"))[6],"_",tail(unlist(str_split(LOS_name,"_"))[-1], n=1),".bedGraph", sep=""), col.names = TRUE, row.names = FALSE, sep = '\t',quote = FALSE)
}
  
setwd("C:/myPath/")
write.table(df_residuals, paste0("df_residuals_20260731_",resolution,".bed"), col.names = TRUE, row.names = FALSE, sep = '\t',quote = FALSE)


#######



