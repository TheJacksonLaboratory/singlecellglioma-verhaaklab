##############################################
# Determine the clonality of copy number events
# Updated: 2021.05.17
# Author: Kevin J.
##################################################

# Working directory for this analysis in the SCGP-analysis project. 
mybasedir = "/Users/johnsk/mnt/verhaak-lab/scgp"
setwd(mybasedir)

########################################
# Necessary packages:
library(tidyverse)
library(openxlsx)
library(RColorBrewer)
library(parallel)
library(data.table)
library(grid)
library(gridExtra)
library(gtable)
library(egg)
########################################
## Generate plot theme.
plot_theme    <- theme_bw(base_size = 12) + theme(axis.title = element_text(size = 12),
                                                  axis.text = element_text(size = 12),
                                                  panel.background = element_rect(fill = "transparent"),
                                                  axis.line = element_blank(),
                                                  strip.background = element_blank(),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(), 
                                                  panel.border = element_blank(),
                                                  axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                                  axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))

#### Load the aneuploidy values from bulk WGS.
SM001 <- read.table("/Users/johnsk/mnt/verhaak-lab/SCGP-analysis/results/cnv/callsegments/SCGP-SM-001-R1-01D-WGS-FC04O9.called.seg", sep="\t", header=T, skip=88, stringsAsFactors = F)
SM001$SAMPLE <- rep("SCGP-SM-001", dim(SM001)[1])
SM001$range <- SM001$END-SM001$START
SM001_filt = SM001 %>% 
  filter(CALL!=0)

SM002 <- read.table("/Users/johnsk/mnt/verhaak-lab/SCGP-analysis/results/cnv/callsegments/SCGP-SM-002-R1-01D-WGS-8CSURD.called.seg", sep="\t", header=T, skip=88, stringsAsFactors = F)
SM002$SAMPLE <- rep("SCGP-SM-002", dim(SM002)[1])
SM002$range <- SM002$END-SM002$START
SM002_filt = SM002 %>% 
  filter(CALL!=0)

SM004 <- read.table("/Users/johnsk/mnt/verhaak-lab/SCGP-analysis/results/cnv/callsegments/SCGP-SM-004-TP-01D-WGS-4WK1HC.called.seg", sep="\t", header=T, skip=88, stringsAsFactors = F)
SM004$SAMPLE <- rep("SCGP-SM-004", dim(SM004)[1])
SM004$range <- SM004$END-SM004$START
SM004_filt = SM004 %>% 
  filter(CALL!=0)

SM006 <- read.table("/Users/johnsk/mnt/verhaak-lab/SCGP-analysis/results/cnv/callsegments/SCGP-SM-006-TP-01D-WGS-YSLCR9.called.seg", sep="\t", header=T, skip=88, stringsAsFactors = F)
SM006$SAMPLE <- rep("SCGP-SM-006", dim(SM006)[1])
SM006$range <- SM006$END-SM006$START
SM006_filt = SM006 %>% 
  filter(CALL!=0)

SM008 <- read.table("/Users/johnsk/mnt/verhaak-lab/SCGP-analysis/results/cnv/callsegments/SCGP-SM-008-R1-01D-WGS-JX9XB5.called.seg", sep="\t", header=T, skip=88, stringsAsFactors = F)
SM008$SAMPLE <- rep("SCGP-SM-008", dim(SM008)[1])
SM008$range <- SM008$END-SM008$START
SM008_filt = SM008 %>% 
  filter(CALL!=0)

SM011 <- read.table("/Users/johnsk/mnt/verhaak-lab/SCGP-analysis/results/cnv/callsegments/SCGP-SM-011-R1-01D-WGS-HSFBM4.called.seg", sep="\t", header=T, skip=88, stringsAsFactors = F)
SM011$SAMPLE <- rep("SCGP-SM-011", dim(SM011)[1])
SM011$range <- SM011$END-SM011$START
SM011_filt = SM011 %>% 
  filter(CALL!=0)

SM012 <- read.table("/Users/johnsk/mnt/verhaak-lab/SCGP-analysis/results/cnv/callsegments/SCGP-SM-012-TP-01D-WGS-MDRUO5.called.seg", sep="\t", header=T, skip=88, stringsAsFactors = F)
SM012$SAMPLE <- rep("SCGP-SM-012", dim(SM012)[1])
SM012$range <- SM012$END-SM012$START
SM012_filt = SM012 %>% 
  filter(CALL!=0)

SM015 <- read.table("/Users/johnsk/mnt/verhaak-lab/SCGP-analysis/results/cnv/callsegments/SCGP-SM-015-TP-01D-WGS-KXLLQO.called.seg", sep="\t", header=T, skip=88, stringsAsFactors = F)
SM015$SAMPLE <- rep("SCGP-SM-015", dim(SM015)[1])
SM015$range <- SM015$END-SM015$START
SM015_filt = SM015 %>% 
  filter(CALL!=0)

SM017 <- read.table("/Users/johnsk/mnt/verhaak-lab/SCGP-analysis/results/cnv/callsegments/SCGP-SM-017-TP-01D-WGS-ZA1161.called.seg", sep="\t", header=T, skip=88, stringsAsFactors = F)
SM017$SAMPLE <- rep("SCGP-SM-017", dim(SM017)[1])
SM017$range <- SM017$END-SM017$START
SM017_filt = SM017 %>% 
  filter(CALL!=0)

SM018 <- read.table("/Users/johnsk/mnt/verhaak-lab/SCGP-analysis/results/cnv/callsegments/SCGP-SM-018-TP-01D-WGS-XUKWEY.called.seg", sep="\t", header=T, skip=88, stringsAsFactors = F)
SM018$SAMPLE <- rep("SCGP-SM-018", dim(SM018)[1])
SM018$range <- SM018$END-SM018$START
SM018_filt = SM018 %>% 
  filter(CALL!=0)

UC917 <- read.table("/Users/johnsk/mnt/verhaak-lab/SCGP-analysis/results/cnv/callsegments/SCGP-UC-917-TP-01D-WGS-GEHXUI.called.seg", sep="\t", header=T, skip=88, stringsAsFactors = F)
UC917$SAMPLE <- rep("SCGP-UC-917", dim(UC917)[1])
UC917$range <- UC917$END-UC917$START
UC917_filt = UC917 %>% 
  filter(CALL!=0)

##################################
## TITAN calls w/ manual review.
##################################
titan_SM001 <- read.table("/Users/johnsk/mnt/verhaak-lab/SCGP-analysis/results/cnv/titan/SCGP-SM-001-R1-01-NB-01D-WGS/ploidy2/SCGP-SM-001-R1-01-NB-01D-WGS_cluster02.segs.txt", sep="\t", header=T, stringsAsFactors = F)
titan_SM002 <- read.table("/Users/johnsk/mnt/verhaak-lab/SCGP-analysis/results/cnv/titan/SCGP-SM-002-R1-01-NB-01D-WGS/ploidy2/SCGP-SM-002-R1-01-NB-01D-WGS_cluster01.segs.txt", sep="\t", header=T, stringsAsFactors = F)
titan_SM004 <- read.table("/Users/johnsk/mnt/verhaak-lab/SCGP-analysis/results/cnv/titan/SCGP-SM-004-TP-01-NB-01D-WGS/ploidy2/SCGP-SM-004-TP-01-NB-01D-WGS_cluster01.segs.txt", sep="\t", header=T, stringsAsFactors = F)
titan_SM006 <- read.table("/Users/johnsk/mnt/verhaak-lab/SCGP-analysis/results/cnv/titan/SCGP-SM-006-TP-01-NB-01D-WGS/ploidy2/SCGP-SM-006-TP-01-NB-01D-WGS_cluster03.segs.txt", sep="\t", header=T, stringsAsFactors = F)
titan_SM008 <- read.table("/Users/johnsk/mnt/verhaak-lab/SCGP-analysis/results/cnv/titan/SCGP-SM-008-R1-01-NB-01D-WGS/ploidy2/SCGP-SM-008-R1-01-NB-01D-WGS_cluster04.segs.txt", sep="\t", header=T, stringsAsFactors = F)
titan_SM011 <- read.table("/Users/johnsk/mnt/verhaak-lab/SCGP-analysis/results/cnv/titan/SCGP-SM-011-R1-01-NB-01D-WGS/ploidy2/SCGP-SM-011-R1-01-NB-01D-WGS_cluster02.segs.txt", sep="\t", header=T, stringsAsFactors = F)
titan_SM012 <- read.table("/Users/johnsk/mnt/verhaak-lab/SCGP-analysis/results/cnv/titan/SCGP-SM-012-TP-01-NB-01D-WGS/ploidy2/SCGP-SM-012-TP-01-NB-01D-WGS_cluster03.segs.txt", sep="\t", header=T, stringsAsFactors = F)
titan_SM015 <- read.table("/Users/johnsk/mnt/verhaak-lab/SCGP-analysis/results/cnv/titan/SCGP-SM-015-TP-01-NB-01D-WGS/ploidy2/SCGP-SM-015-TP-01-NB-01D-WGS_cluster01.segs.txt", sep="\t", header=T, stringsAsFactors = F)
titan_SM017 <- read.table("/Users/johnsk/mnt/verhaak-lab/SCGP-analysis/results/cnv/titan/SCGP-SM-017-TP-01-NB-01D-WGS/ploidy2/SCGP-SM-017-TP-01-NB-01D-WGS_cluster02.segs.txt", sep="\t", header=T, stringsAsFactors = F)
titan_SM018 <- read.table("/Users/johnsk/mnt/verhaak-lab/SCGP-analysis/results/cnv/titan/SCGP-SM-018-TP-01-NB-01D-WGS/ploidy2/SCGP-SM-018-TP-01-NB-01D-WGS_cluster01.segs.txt", sep="\t", header=T, stringsAsFactors = F)
titan_UC917 <- read.table("/Users/johnsk/mnt/verhaak-lab/SCGP-analysis/results/cnv/titan/SCGP-UC-917-TP-01-NB-01D-WGS/ploidy2/SCGP-UC-917-TP-01-NB-01D-WGS_cluster02.segs.txt", sep="\t", header=T, stringsAsFactors = F)



##################################################
## Intersect data between WGS and TITAN
##################################################
## Objective: Intersect the WGS data w/ TITAN calls.
### SM001 ###
copynum_calls = data.table(
  chr = as.character(SM001_filt$CONTIG),
  start = as.numeric(SM001_filt$START), 
  end =  as.numeric(SM001_filt$END),
  sample_id = SM001_filt$SAMPLE,
  copy_ratio = SM001_filt$MEAN_LOG2_COPY_RATIO,
  CALL = SM001_filt$CALL)
setkey(copynum_calls, chr, start, end)

# The TITAN coordinates.
titan_coord = data.table(
  chr = as.character(titan_SM001$Chromosome),
  start = as.numeric(titan_SM001$`Start_Position.bp.`), 
  end =  as.numeric(titan_SM001$`End_Position.bp.`),
  cellular_prevalence = titan_SM001$Cellular_Prevalence,
  clonal_cluster = titan_SM001$Clonal_Cluster)
setkey(titan_coord, chr, start, end)

# Find overlaps for between the TITAN coordinates and WGS CNV data.
titan_wgs_overlap <- foverlaps(copynum_calls, titan_coord, type = "any", nomatch = 0)
titan_wgs_overlap_df <- as.data.frame(titan_wgs_overlap)

titan_filt = titan_wgs_overlap_df %>% 
  select(clonal_cluster, chr, i.start, i.end) %>% 
  distinct() %>% 
  mutate(range = i.end-i.start) %>% 
  group_by(clonal_cluster) %>% 
  summarise(clone_seg_length = sum(range))
## SM001 clonality: 0.7892959
547888719/(547888719+146259962)

#### SM002 ####
## Try to do this the data.table way.
copynum_calls = data.table(
  chr = as.character(SM002_filt$CONTIG),
  start = as.numeric(SM002_filt$START), 
  end =  as.numeric(SM002_filt$END),
  sample_id = SM002_filt$SAMPLE,
  copy_ratio = SM002_filt$MEAN_LOG2_COPY_RATIO,
  CALL = SM002_filt$CALL)
setkey(copynum_calls, chr, start, end)

# The TITAN coordinates.
titan_coord = data.table(
  chr = as.character(titan_SM002$Chromosome),
  start = as.numeric(titan_SM002$`Start_Position.bp.`), 
  end =  as.numeric(titan_SM002$`End_Position.bp.`),
  cellular_prevalence = titan_SM002$Cellular_Prevalence,
  clonal_cluster = titan_SM002$Clonal_Cluster)
setkey(titan_coord, chr, start, end)

# Find overlaps for between the TITAN coordinates and WGS CNV data.
titan_wgs_overlap <- foverlaps(copynum_calls, titan_coord, type = "any", nomatch = 0)
titan_wgs_overlap_df <- as.data.frame(titan_wgs_overlap)

titan_filt = titan_wgs_overlap_df %>% 
  select(clonal_cluster, chr, i.start, i.end) %>% 
  distinct() %>% 
  mutate(range = i.end-i.start) %>% 
  group_by(clonal_cluster) %>% 
  summarise(clone_seg_length = sum(range))
## SM002 clonality: 1.0
126249970/(126249970)

#### SM004 ####
## Try to do this the data.table way.
copynum_calls = data.table(
  chr = as.character(SM004_filt$CONTIG),
  start = as.numeric(SM004_filt$START), 
  end =  as.numeric(SM004_filt$END),
  sample_id = SM004_filt$SAMPLE,
  copy_ratio = SM004_filt$MEAN_LOG2_COPY_RATIO,
  CALL = SM004_filt$CALL)
setkey(copynum_calls, chr, start, end)

# The TITAN coordinates.
titan_coord = data.table(
  chr = as.character(titan_SM004$Chromosome),
  start = as.numeric(titan_SM004$`Start_Position.bp.`), 
  end =  as.numeric(titan_SM004$`End_Position.bp.`),
  cellular_prevalence = titan_SM004$Cellular_Prevalence,
  clonal_cluster = titan_SM004$Clonal_Cluster)
setkey(titan_coord, chr, start, end)

# Find overlaps for between the TITAN coordinates and WGS CNV data.
titan_wgs_overlap <- foverlaps(copynum_calls, titan_coord, type = "any", nomatch = 0)
titan_wgs_overlap_df <- as.data.frame(titan_wgs_overlap)

titan_filt = titan_wgs_overlap_df %>% 
  select(clonal_cluster, chr, i.start, i.end) %>% 
  distinct() %>% 
  mutate(range = i.end-i.start) %>% 
  group_by(clonal_cluster) %>% 
  summarise(clone_seg_length = sum(range))
## SM004 clonality: 1.0
153029927/(153029927)


#### SM006 ####
## Try to do this the data.table way.
copynum_calls = data.table(
  chr = as.character(SM006_filt$CONTIG),
  start = as.numeric(SM006_filt$START), 
  end =  as.numeric(SM006_filt$END),
  sample_id = SM006_filt$SAMPLE,
  copy_ratio = SM006_filt$MEAN_LOG2_COPY_RATIO,
  CALL = SM006_filt$CALL)
setkey(copynum_calls, chr, start, end)

# The TITAN coordinates.
titan_coord = data.table(
  chr = as.character(titan_SM006$Chromosome),
  start = as.numeric(titan_SM006$`Start_Position.bp.`), 
  end =  as.numeric(titan_SM006$`End_Position.bp.`),
  cellular_prevalence = titan_SM006$Cellular_Prevalence,
  clonal_cluster = titan_SM006$Clonal_Cluster)
setkey(titan_coord, chr, start, end)

# Find overlaps for between the TITAN coordinates and WGS CNV data.
titan_wgs_overlap <- foverlaps(copynum_calls, titan_coord, type = "any", nomatch = 0)
titan_wgs_overlap_df <- as.data.frame(titan_wgs_overlap)

titan_trim = titan_wgs_overlap_df %>% 
  select(clonal_cluster, chr, i.start, i.end) %>% 
  distinct()
to_drop_2 = which(titan_trim$clonal_cluster==2 & titan_trim$chr%in%c(7,10,17,19))
to_drop_3 = which(titan_trim$clonal_cluster==3 & titan_trim$chr%in%c(7,10,17,19))
titan_filt = titan_trim[-c(to_drop_2, to_drop_3), ]
titan_filt_res = titan_filt %>% 
  mutate(range = i.end-i.start) %>% 
  group_by(clonal_cluster) %>% 
  summarise(clone_seg_length = sum(range))
## Clonality: 0.5827958
c(126189966+149419887)/(197299951+126189966+149419887)

#### SM008 ####
## Try to do this the data.table way.
copynum_calls = data.table(
  chr = as.character(SM008_filt$CONTIG),
  start = as.numeric(SM008_filt$START), 
  end =  as.numeric(SM008_filt$END),
  sample_id = SM008_filt$SAMPLE,
  copy_ratio = SM008_filt$MEAN_LOG2_COPY_RATIO,
  CALL = SM008_filt$CALL)
setkey(copynum_calls, chr, start, end)

# The TITAN coordinates.
titan_coord = data.table(
  chr = as.character(titan_SM008$Chromosome),
  start = as.numeric(titan_SM008$`Start_Position.bp.`), 
  end =  as.numeric(titan_SM008$`End_Position.bp.`),
  cellular_prevalence = titan_SM008$Cellular_Prevalence,
  clonal_cluster = titan_SM008$Clonal_Cluster)
setkey(titan_coord, chr, start, end)

# Find overlaps for between the TITAN coordinates and WGS CNV data.
titan_wgs_overlap <- foverlaps(copynum_calls, titan_coord, type = "any", nomatch = 0)
titan_wgs_overlap_df <- as.data.frame(titan_wgs_overlap)

titan_filt = titan_wgs_overlap_df %>% 
  select(clonal_cluster, chr, i.start, i.end) %>% 
  distinct()
#to_drop_2 = which(titan_trim$clonal_cluster==2 & titan_trim$chr%in%c(7,10,17,19))
#to_drop_3 = which(titan_trim$clonal_cluster==3 & titan_trim$chr%in%c(7,10,17,19))
#titan_filt = titan_trim[-c(to_drop_2, to_drop_3), ]
titan_filt_res = titan_filt %>% 
  mutate(range = i.end-i.start) %>% 
  group_by(clonal_cluster) %>% 
  summarise(clone_seg_length = sum(range))
## Combine clones 1+2: 0.78 
c(47249820+69559871)/(47249820+69559871+13269975+18359975)

#### SM011 ####
## Try to do this the data.table way.
copynum_calls = data.table(
  chr = as.character(SM011_filt$CONTIG),
  start = as.numeric(SM011_filt$START), 
  end =  as.numeric(SM011_filt$END),
  sample_id = SM011_filt$SAMPLE,
  copy_ratio = SM011_filt$MEAN_LOG2_COPY_RATIO,
  CALL = SM011_filt$CALL)
setkey(copynum_calls, chr, start, end)

# The TITAN coordinates.
titan_coord = data.table(
  chr = as.character(titan_SM011$Chromosome),
  start = as.numeric(titan_SM011$`Start_Position.bp.`), 
  end =  as.numeric(titan_SM011$`End_Position.bp.`),
  cellular_prevalence = titan_SM011$Cellular_Prevalence,
  clonal_cluster = titan_SM011$Clonal_Cluster)
setkey(titan_coord, chr, start, end)

# Find overlaps for between the TITAN coordinates and WGS CNV data.
titan_wgs_overlap <- foverlaps(copynum_calls, titan_coord, type = "any", nomatch = 0)
titan_wgs_overlap_df <- as.data.frame(titan_wgs_overlap)

titan_filt = titan_wgs_overlap_df %>% 
  select(clonal_cluster, chr, i.start, i.end) %>% 
  distinct() %>% 
  mutate(range = i.end-i.start) %>% 
  group_by(clonal_cluster) %>% 
  summarise(clone_seg_length = sum(range))
## 0.3522107 clonal.
50809865/(50809865+93449976)

#### SM012 ####
## Try to do this the data.table way.
copynum_calls = data.table(
  chr = as.character(SM012_filt$CONTIG),
  start = as.numeric(SM012_filt$START), 
  end =  as.numeric(SM012_filt$END),
  sample_id = SM012_filt$SAMPLE,
  copy_ratio = SM012_filt$MEAN_LOG2_COPY_RATIO,
  CALL = SM012_filt$CALL)
setkey(copynum_calls, chr, start, end)

# The TITAN coordinates.
titan_coord = data.table(
  chr = as.character(titan_SM012$Chromosome),
  start = as.numeric(titan_SM012$`Start_Position.bp.`), 
  end =  as.numeric(titan_SM012$`End_Position.bp.`),
  cellular_prevalence = titan_SM012$Cellular_Prevalence,
  clonal_cluster = titan_SM012$Clonal_Cluster)
setkey(titan_coord, chr, start, end)

# Find overlaps for between the TITAN coordinates and WGS CNV data.
titan_wgs_overlap <- foverlaps(copynum_calls, titan_coord, type = "any", nomatch = 0)
titan_wgs_overlap_df <- as.data.frame(titan_wgs_overlap)

titan_trim = titan_wgs_overlap_df %>% 
  select(clonal_cluster, chr, i.start, i.end) %>% 
  distinct()
to_drop_2 = which(titan_trim$clonal_cluster==2 & titan_trim$chr%in%c(3,4,7,10,13,14,19,21,22))
to_drop_3 = which(titan_trim$clonal_cluster==3 & titan_trim$chr%in%c(3,4,7,10,13,14,19,21,22))
titan_filt = titan_trim[-c(to_drop_2, to_drop_3), ]
titan_filt_res = titan_filt %>% 
  mutate(range = i.end-i.start) %>% 
  group_by(clonal_cluster) %>% 
  summarise(clone_seg_length = sum(range))
## 0.6278271 clonal
472979296/(472979296+81749968+198629885)

#### SM015 ####
## Try to do this the data.table way.
copynum_calls = data.table(
  chr = as.character(SM015_filt$CONTIG),
  start = as.numeric(SM015_filt$START), 
  end =  as.numeric(SM015_filt$END),
  sample_id = SM015_filt$SAMPLE,
  copy_ratio = SM015_filt$MEAN_LOG2_COPY_RATIO,
  CALL = SM015_filt$CALL)
setkey(copynum_calls, chr, start, end)

# The TITAN coordinates.
titan_coord = data.table(
  chr = as.character(titan_SM015$Chromosome),
  start = as.numeric(titan_SM015$`Start_Position.bp.`), 
  end =  as.numeric(titan_SM015$`End_Position.bp.`),
  cellular_prevalence = titan_SM015$Cellular_Prevalence,
  clonal_cluster = titan_SM015$Clonal_Cluster)
setkey(titan_coord, chr, start, end)

# Find overlaps for between the TITAN coordinates and WGS CNV data.
titan_wgs_overlap <- foverlaps(copynum_calls, titan_coord, type = "any", nomatch = 0)
titan_wgs_overlap_df <- as.data.frame(titan_wgs_overlap)

titan_filt = titan_wgs_overlap_df %>% 
  select(clonal_cluster, chr, i.start, i.end) %>% 
  distinct() %>% 
  mutate(range = i.end-i.start) %>% 
  group_by(clonal_cluster) %>% 
  summarise(clone_seg_length = sum(range))
## SM015 clonality:
285879908/(285879908)


#### SM017 ####
## Try to do this the data.table way.
copynum_calls = data.table(
  chr = as.character(SM017_filt$CONTIG),
  start = as.numeric(SM017_filt$START), 
  end =  as.numeric(SM017_filt$END),
  sample_id = SM017_filt$SAMPLE,
  copy_ratio = SM017_filt$MEAN_LOG2_COPY_RATIO,
  CALL = SM017_filt$CALL)
setkey(copynum_calls, chr, start, end)

# The TITAN coordinates.
titan_coord = data.table(
  chr = as.character(titan_SM017$Chromosome),
  start = as.numeric(titan_SM017$`Start_Position.bp.`), 
  end =  as.numeric(titan_SM017$`End_Position.bp.`),
  cellular_prevalence = titan_SM017$Cellular_Prevalence,
  clonal_cluster = titan_SM017$Clonal_Cluster)
setkey(titan_coord, chr, start, end)

# Find overlaps for between the TITAN coordinates and WGS CNV data.
titan_wgs_overlap <- foverlaps(copynum_calls, titan_coord, type = "any", nomatch = 0)
titan_wgs_overlap_df <- as.data.frame(titan_wgs_overlap)

titan_filt = titan_wgs_overlap_df %>% 
  select(clonal_cluster, chr, i.start, i.end) %>% 
  distinct() %>% 
  mutate(range = i.end-i.start) %>% 
  group_by(clonal_cluster) %>% 
  summarise(clone_seg_length = sum(range))
## 0.83 clonal
409369707/(409369707+81109984)

#### SM018 ####
## Try to do this the data.table way.
copynum_calls = data.table(
  chr = as.character(SM018_filt$CONTIG),
  start = as.numeric(SM018_filt$START), 
  end =  as.numeric(SM018_filt$END),
  sample_id = SM018_filt$SAMPLE,
  copy_ratio = SM018_filt$MEAN_LOG2_COPY_RATIO,
  CALL = SM018_filt$CALL)
setkey(copynum_calls, chr, start, end)

# The TITAN coordinates.
titan_coord = data.table(
  chr = as.character(titan_SM018$Chromosome),
  start = as.numeric(titan_SM018$`Start_Position.bp.`), 
  end =  as.numeric(titan_SM018$`End_Position.bp.`),
  cellular_prevalence = titan_SM018$Cellular_Prevalence,
  clonal_cluster = titan_SM018$Clonal_Cluster)
setkey(titan_coord, chr, start, end)

# Find overlaps for between the TITAN coordinates and WGS CNV data.
titan_wgs_overlap <- foverlaps(copynum_calls, titan_coord, type = "any", nomatch = 0)
titan_wgs_overlap_df <- as.data.frame(titan_wgs_overlap)

titan_filt = titan_wgs_overlap_df %>% 
  select(clonal_cluster, chr, i.start, i.end) %>% 
  distinct() %>% 
  mutate(range = i.end-i.start) %>% 
  group_by(clonal_cluster) %>% 
  summarise(clone_seg_length = sum(range))
## 1.0 clonal
363069956/(363069956)


#### UC917 ####
## Try to do this the data.table way.
copynum_calls = data.table(
  chr = as.character(UC917_filt$CONTIG),
  start = as.numeric(UC917_filt$START), 
  end =  as.numeric(UC917_filt$END),
  sample_id = UC917_filt$SAMPLE,
  copy_ratio = UC917_filt$MEAN_LOG2_COPY_RATIO,
  CALL = UC917_filt$CALL)
setkey(copynum_calls, chr, start, end)

# The TITAN coordinates.
titan_coord = data.table(
  chr = as.character(titan_UC917$Chromosome),
  start = as.numeric(titan_UC917$`Start_Position.bp.`), 
  end =  as.numeric(titan_UC917$`End_Position.bp.`),
  cellular_prevalence = titan_UC917$Cellular_Prevalence,
  clonal_cluster = titan_UC917$Clonal_Cluster)
setkey(titan_coord, chr, start, end)

# Find overlaps for between the TITAN coordinates and WGS CNV data.
titan_wgs_overlap <- foverlaps(copynum_calls, titan_coord, type = "any", nomatch = 0)
titan_wgs_overlap_df <- as.data.frame(titan_wgs_overlap)

titan_filt = titan_wgs_overlap_df %>% 
  select(clonal_cluster, chr, i.start, i.end) %>% 
  distinct() %>% 
  mutate(range = i.end-i.start) %>% 
  group_by(clonal_cluster) %>% 
  summarise(clone_seg_length = sum(range))
## 0.29 clonal
143609720/(143609720+351329957)

## Create stacked barplot based on SCNA burden and CNV clonality.
# Combine all aneuploidy values.
sum(SM001$range[SM001$CALL!=0])/sum(SM001$range)
sum(SM002$range[SM002$CALL!=0])/sum(SM002$range)
sum(SM004$range[SM004$CALL!=0])/sum(SM004$range)
sum(SM006$range[SM006$CALL!=0])/sum(SM006$range)
sum(SM008$range[SM008$CALL!=0])/sum(SM008$range)
sum(SM011$range[SM011$CALL!=0])/sum(SM011$range)
sum(SM012$range[SM012$CALL!=0])/sum(SM012$range)
sum(SM015$range[SM015$CALL!=0])/sum(SM015$range)
sum(SM017$range[SM017$CALL!=0])/sum(SM017$range)
sum(SM018$range[SM018$CALL!=0])/sum(SM018$range)
sum(UC917$range[UC917$CALL!=0])/sum(UC917$range)
## Clonality values:
## SM001 clonality: 0.7892959
547888719/(547888719+146259962)
## SM002 clonality: 1.0
126249970/(126249970)
## SM004 clonality: 1.0
153029927/(153029927)
## SM006 clonality: 0.5827958
c(126189966+149419887)/(197299951+126189966+149419887)
## SM008 clonality (clones 1+2): 0.78 
c(47249820+69559871)/(47249820+69559871+13269975+18359975)
## SM011 clonality: 0.3522107 
50809865/(50809865+93449976)
## SM012 clonality: 0.6278271
472979296/(472979296+81749968+198629885)
## SM015 clonality: 1.0
285879908/(285879908)
## SM017 clonality: 0.83 
409369707/(409369707+81109984)
## SM018 clonality: 1.0
363069956/(363069956)
## UC917 clonality: 0.29
143609720/(143609720+351329957)

case_barcode <- c("SM001", "SM002", "SM004","SM006", "SM008", "SM011", "SM012", "SM015", "SM017", "SM018", "UC917")
aneuploidy_values = c(0.196, 0.0420, 0.0515, 0.155, 0.122, 0.062, 0.312,
                      0.101, 0.144, 0.128, 0.143)
clonality_values <- c(0.7892959, 1, 1, 0.5827958, 0.78, 1, 0.6278271, 1.0, 0.83, 0.98, 0.29)
scgp_aneuploidy <- as.data.frame(cbind(case_barcode, aneuploidy_values, clonality_values))
scgp_aneuploidy$aneuploidy_values <- as.numeric(as.character(scgp_aneuploidy$aneuploidy_values))
scgp_aneuploidy$clonality_values <- as.numeric(as.character(scgp_aneuploidy$clonality_values))
scgp_aneuploidy$clonal <- scgp_aneuploidy$aneuploidy_values*scgp_aneuploidy$clonality_values
scgp_aneuploidy$subclonal <- scgp_aneuploidy$aneuploidy_values-scgp_aneuploidy$clonal

scgp_aneuploidy_grp = scgp_aneuploidy %>% 
  gather(key = "clonality", value = "aneuploidy_prop", clonal, subclonal)

clonality_order <- c("subclonal", "clonal")
scgp_aneuploidy_grp <- scgp_aneuploidy_grp %>% mutate(clonality = factor(clonality, levels = clonality_order))

case_order <- c("SM004", "SM001", "SM015", "UC917", "SM002", "SM008", "SM006", "SM012", "SM017", "SM018", "SM011")
scgp_aneuploidy_grp <- scgp_aneuploidy_grp %>% mutate(case_barcode = factor(case_barcode, levels = case_order))

scgp_aneuploidy_grp$idh_status <- ifelse(scgp_aneuploidy_grp$case_barcode%in%c("SM004", "SM001", "SM015", "UC917", "SM002", "SM008"), "IDHmut", "IDHwt")
scgp_aneuploidy_grp$evo_mode <- ifelse(scgp_aneuploidy_grp$case_barcode%in%c("SM004", "SM001", "SM017", "SM018", "UC917", "SM011"), "linear", "branched")

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig6/SuppFig9-scna-clonality.pdf", height = 4, width = 6, bg = "transparent", useDingbats = FALSE)
ggplot(scgp_aneuploidy_grp, aes(x=case_barcode, fill = clonality, y=aneuploidy_prop)) + geom_bar(position="stack", stat="identity") + 
  plot_theme +
  theme(panel.spacing.x = unit(1.5, "lines"),
        axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="", y="SCNA burden", fill="SCNA clonality") +
  scale_fill_manual(values = c("clonal" = "#2078b4", 
                                "subclonal" = "#B47846")) +
  ylim(0,0.35) +
  facet_grid(. ~ idh_status, scales = "free_x", space = "free")
dev.off()

### END ###

