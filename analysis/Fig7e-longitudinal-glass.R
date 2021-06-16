##############################################
# Generate a tiled heatmap for GLASS methylation data
# Updated: 2020.06.07
# Author: Kevin J.
##################################################

## Working directory for this analysis.
mybasedir = "/Users/johnsk/github/"
setwd(mybasedir)

########################################
# Necessary packages:
library(tidyverse)
library(grid)
library(gridExtra)
library(gtable)
library(egg)
library(RColorBrewer)
library(openxlsx)
library(DBI)
library(viridis)
########################################
## Establish connection with GLASS database (publication version).
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB2")


######### CLINICAL ########
# Load additional tables from the database.
subtypes = dbReadTable(con,  Id(schema="clinical", table="subtypes"))  
cases = dbReadTable(con,  Id(schema="clinical", table="cases"))  
subclonalselection = dbReadTable(con,  Id(schema="analysis", table="subclonalselection"))
seqz_params = dbReadTable(con,  Id(schema="variants", table="seqz_params"))
pairs = dbReadTable(con,  Id(schema="analysis", table="pairs"))
silver_set = dbReadTable(con,  Id(schema="analysis", table="silver_set"))  
clinical_tumor_pairs_query = read_file("misc/clinical-tumor-pairs-db2.sql")
clinical_tumor_pairs <- dbGetQuery(con, clinical_tumor_pairs_query)

## Create an IDHmut vs. IDHwt subtype-level.
subtypes = subtypes %>% 
  mutate(idh_subtype = ifelse(idh_codel_subtype=="IDHwt", "IDHwt", "IDHmut"))

######### Meth Changes ########
summary_freq = readRDS("data/glass-beta-discord-freq.rds")
summary_freq = summary_freq %>% 
  mutate(idh_subtype = ifelse(idh_codel_subtype=="IDHwt", "IDHwt", "IDHmut"))

####### Purity   ########
## Load purity based on Sequenza (WXS|WGS) estimates:
seqz_params_filt = seqz_params %>% 
  inner_join(pairs, by="pair_barcode") %>% 
  mutate(sample_barcode = substr(tumor_barcode, 1, 15)) %>% 
  select(sample_barcode, cellularity) %>% 
  group_by(sample_barcode) %>% 
  summarise(avg_purity = mean(cellularity)) %>% 
  ungroup()

## Subset purity data.
puritydata = clinical_tumor_pairs %>% 
  filter(tumor_pair_barcode%in%silver_set$tumor_pair_barcode) %>% 
  mutate(sample_barcode_a = substr(tumor_barcode_a, 1, 15),
         sample_barcode_b = substr(tumor_barcode_b, 1, 15)) %>% 
  left_join(seqz_params_filt, by=c("sample_barcode_a"= "sample_barcode")) %>% 
  left_join(seqz_params_filt, by=c("sample_barcode_b"= "sample_barcode")) %>% 
  select(tumor_pair_barcode:tumor_barcode_b, purity_a = avg_purity.x, purity_b = avg_purity.y) %>% 
  filter(case_barcode%in%unique(summary_freq$case_barcode)) %>% 
  left_join(subtypes, by="case_barcode")


## GENETIC-BASED DATA ##
###### Mutations ########
mutfdata <- dbGetQuery(con, read_file("misc/heatmap_mf.sql"))
mut_freq_prop_case = mutfdata %>% 
  mutate(P = (count_a-intersection_ab)/union_ab,
         R =  (count_b-intersection_ab)/union_ab,
         S = intersection_ab/union_ab, 
         idh_subtype = ifelse(idh_codel_subtype=="IDHwt", "IDHwt", "IDHmut")) %>% 
  select(P, R, S, case_barcode, union_ab, idh_codel_subtype, idh_subtype) %>% 
  gather(mutation_type, mutation_percent, c(P, R, S), -case_barcode, -union_ab, -idh_codel_subtype, -idh_subtype) %>%
  mutate(mutation_type = factor(mutation_type, levels = c("R", "P", "S"))) %>% 
  filter(case_barcode%in%unique(summary_freq$case_barcode))

####### Clinical ########
clin_data = clinical_tumor_pairs %>% 
  filter(case_barcode%in%unique(summary_freq$case_barcode)) %>% 
  left_join(subtypes, by="case_barcode") %>% 
  filter(tumor_pair_barcode%in%unique(mutfdata$tumor_pair_barcode)) 
  
####### CNVs  ########
anpldata <- dbGetQuery(con, read_file("misc/heatmap_aneuploidy.sql"))
anpldata = anpldata %>% 
  mutate(idh_subtype = ifelse(idh_codel_subtype=="IDHwt", "IDHwt", "IDHmut")) %>% 
  filter(case_barcode%in%unique(summary_freq$case_barcode))


## METHYLATION-BASED DATA ##
####### MNP subtype switch  #######
il450k_g_filt_meta_long = readRDS(file ="data/glass-Mvals-annot.rds")
## Retrieve a list of paired IDAT files
il450k_g_filt_meta_wide = il450k_g_filt_meta_long %>% 
  select(case_barcode, timepoint, idat) %>% 
  mutate(timepoint = recode(timepoint, `0` = "initial", `1` = "recurrence")) %>% 
  spread(timepoint, idat) 
mnp_class  <- read.table("data/glass-mnp-classification-20190408.txt", sep="\t", header=T, stringsAsFactors = F)
mnpdata = il450k_g_filt_meta_wide %>% 
  left_join(mnp_class, by=c("initial"= "idat")) %>% 
  left_join(mnp_class, by=c("recurrence"= "idat")) %>% 
  select(case_barcode, mnp_a = mnp_classification.x, mnp_b = mnp_classification.y) %>%
  mutate(subtype_stability = ifelse(mnp_a!=mnp_b, "discordant", "concordant"))  %>% 
  left_join(subtypes, by="case_barcode")
  
#### RNA-based subtypes ######
glass_rna  <- read.table("data/glass-analysis_transcriptional_subtypes.tsv", sep="\t", header=T, stringsAsFactors = F)

## Create a joinable sample_barcode.
glass_rna_filt = glass_rna %>% 
  mutate(sample_barcode = substr(aliquot_barcode, 1, 15)) %>% 
  select(aliquot_barcode, transcriptional_subtype = subtype, sample_barcode)
  
## Use silver set and filter to samples for which we have DNA methylation data.
rna_data = clinical_tumor_pairs %>% 
  filter(tumor_pair_barcode%in%silver_set$tumor_pair_barcode) %>% 
  mutate(sample_barcode_a = substr(tumor_barcode_a, 1, 15),
         sample_barcode_b = substr(tumor_barcode_b, 1, 15)) %>% 
  left_join(glass_rna_filt, by=c("sample_barcode_a"="sample_barcode")) %>% 
  left_join(glass_rna_filt, by=c("sample_barcode_b"="sample_barcode")) %>% 
  select(tumor_pair_barcode:tumor_barcode_b, subtype_a = transcriptional_subtype.x, subtype_b = transcriptional_subtype.y) %>% 
  filter(case_barcode%in%unique(summary_freq$case_barcode)) %>% 
  left_join(subtypes, by="case_barcode")

####### mDNAsi    #######
mDNAsi  <- read.table("data/glass-mDNAsi-values.csv", sep=",", header=T, stringsAsFactors = F)
mDNAsidata = il450k_g_filt_meta_wide %>% 
  left_join(mDNAsi, by=c("initial"= "idat")) %>% 
  left_join(mDNAsi, by=c("recurrence"= "idat")) %>% 
  select(case_barcode, mDNAsi_a = mDNAsi.x, mDNAsi_b = mDNAsi.y)  %>% 
  left_join(subtypes, by="case_barcode")

####### epiTOC    #######
epi_TOC = read.table("data/glass-epiTOC-values.csv", sep = ",", row.names = 1, header = T)
epi_TOCdata = il450k_g_filt_meta_wide %>% 
  left_join(epi_TOC, by=c("initial"= "idat")) %>% 
  left_join(epi_TOC, by=c("recurrence"= "idat")) %>% 
  select(case_barcode, epiTOC_a = epiTOC.x, epiTOC_b = epiTOC.y)  %>% 
  left_join(subtypes, by="case_barcode")

####### Horvath Clock    #######
horvath = read.table("data/full_glass_horvath_betas_unnormalized.output.csv", sep = ",", row.names = 1, header = T)
horvath$idat <- gsub("^X", "", rownames(horvath))
horvathdata = il450k_g_filt_meta_wide %>% 
  left_join(horvath, by=c("initial"= "idat")) %>% 
  left_join(horvath, by=c("recurrence"= "idat")) %>% 
  select(case_barcode, DNAmAge_a = DNAmAge.x, DNAmAge_b = DNAmAge.y)  %>% 
  left_join(subtypes, by="case_barcode")

####### eITH    #######
eITH = read.table("data/glass-eITH-values.csv", sep = ",", row.names = 1, header = T)
eITHdata = il450k_g_filt_meta_wide %>% 
  left_join(eITH, by=c("initial"= "idat")) %>% 
  left_join(eITH, by=c("recurrence"= "idat")) %>% 
  select(case_barcode, eITH_a = eITH.x, eITH_b = eITH.y) %>% 
  left_join(subtypes, by="case_barcode")

##########################
## Order all of the data:
sort_df <- summary_freq %>%
  filter(concordance == "no_change") %>% 
  arrange(desc(freq)) 
case_order <- unique(sort_df$case_barcode)
## Now order all data by the discordant methylation.
summary_freq <- summary_freq %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
puritydata <- puritydata %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
anpldata  <- anpldata %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
mut_freq_prop_case <- mut_freq_prop_case %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
mnpdata <- mnpdata %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
mDNAsidata <- mDNAsidata %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
epi_TOCdata <- epi_TOCdata %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
horvathdata <- horvathdata %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
eITHdata <- eITHdata %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
rna_data <- rna_data %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
clin_data <- clin_data %>% mutate(case_barcode = factor(case_barcode, levels = case_order))

########################################################
## Tests for correlations between proportion methylation differences and abs(change) in other variables
df_red = summary_freq %>% 
  filter(concordance%in%c("meth_gain", "meth_loss")) %>% 
  group_by(case_barcode) %>% 
  summarise(delta_meth_prop = sum(freq)) %>% 
  ungroup() %>% 
  inner_join(mDNAsidata, by="case_barcode") %>% 
  inner_join(epi_TOCdata, by="case_barcode") %>%
  inner_join(horvathdata, by="case_barcode") %>%
  inner_join(eITHdata, by="case_barcode") %>% 
  inner_join(puritydata, by="case_barcode") %>% 
  inner_join(anpldata, by="case_barcode") %>% 
  inner_join(rna_data, by="case_barcode") %>% 
  inner_join(clin_data, by="case_barcode")
## Delta purity was NOT associated with delta methylation prop.
cor.test(df_red$delta_meth_prop, abs(df_red$purity_b-df_red$purity_a), method = "spearman")

## ***Delta copy number was associated with change in DNA methylation proportion.
cor.test(df_red$delta_meth_prop, abs(df_red$aneuploidy_b-df_red$aneuploidy_a), method = "spearman")
cor.test(df_red$delta_meth_prop, df_red$aneuploidy_b-df_red$aneuploidy_a, method = "spearman")

## The abs(mDNAsi) score was weakly associated with delta methylation change.
cor.test(df_red$delta_meth_prop, abs(df_red$mDNAsi_b-df_red$mDNAsi_a), method = "spearman") # 0.01
cor.test(df_red$delta_meth_prop, abs(df_red$epiTOC_b-df_red$epiTOC_a), method = "spearman")
cor.test(df_red$delta_meth_prop, abs(df_red$DNAmAge_b-df_red$DNAmAge_a), method = "spearman")
## The epigenetic heterogeneity / DNA methylation instability was associated with the eITH score.
cor.test(df_red$delta_meth_prop, abs(df_red$eITH_b-df_red$eITH_a), method = "spearman")

## The epigenetic heterogeneity / DNA methylation instability was associated with the treatment?
wilcox.test(df_red$delta_meth_prop~df_red$received_tmz)
wilcox.test(df_red$delta_meth_prop~df_red$received_rt)
wilcox.test(df_red$delta_meth_prop~df_red$received_treatment)
wilcox.test(df_red$eITH_b-df_red$eITH_a~df_red$received_tmz)
wilcox.test(df_red$eITH_b-df_red$eITH_a~df_red$received_rt)
wilcox.test(df_red$eITH_b-df_red$eITH_a~df_red$received_treatment)
## There appears to be a larger difference between tumor_a and tumor_b upon treatment (non-directional).
ggplot(df_red, aes(x=received_treatment, y=abs(df_red$eITH_b-df_red$eITH_a))) + geom_boxplot()

##################################################
## TEST PLOT FUNCTION
##################################################

testPlot <- function(gg, grid = TRUE) {
  if(grid)
    gg + plot_theme + plot_grid + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  else
    gg + plot_theme + theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

######################## 
## Common plotting elements
########################
## Extract legend
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


plot_grid     <- facet_grid(. ~ idh_subtype, scales = "free_x", space = "free")
plot_theme    <- theme_bw(base_size = 12) + theme(axis.title = element_text(size = 12),
                                                  axis.text = element_text(size = 12),
                                                  panel.background = element_rect(fill = "transparent"),
                                                  axis.line = element_blank(),
                                                  strip.background = element_blank(),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(), 
                                                  panel.border = element_blank(),
                                                  panel.spacing.x = unit(1.5, "lines"),
                                                  axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                                  axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))

null_legend   <- theme(legend.position = 'none')
null_x        <- theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())
null_y        <- theme(axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())
bottom_x      <- theme(axis.text.x=element_blank())
null_facet    <- theme(strip.background = element_blank(),
                       strip.text.x = element_blank())
top_margin    <- theme(plot.margin= unit(c(1, 1, 0.1, 1), "lines")) ## Top, Right, Bottom, Left
middle_margin <- theme(plot.margin= unit(c(0, 1, 0.1, 1), "lines"))
bottom_margin <- theme(plot.margin= unit(c(0, 1, 1, 1), "lines"))

########################
## View test plot
########################

testPlot <- function(gg, grid = TRUE) {
  if(grid)
    gg + plot_theme + plot_grid + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  else
    gg + plot_theme + theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

gg_cbind <- function(..., widths = NULL) {
  if(length(match.call()) - 2 != length(widths))
    message("Number of widths does not match number of columns")
  gg <- gtable_cbind(...)
  panels = gg$layout$t[grep("panel", gg$layout$name)]
  gg$widths[panels] <- unit(widths, "null")
  return(gg)
}


gg_rbind <- function(..., heights = NULL, ncol = 2) {
  if(length(match.call()) - 3 != length(heights))
    message("Number of heights does not match number of rows")
  gg <- gtable_rbind(...)
  panels = gg$layout$t[grep("panel", gg$layout$name)]
  gg$heights[panels] <- unit(rep(heights,each = ncol), "null")
  return(gg)
}

## Extract legend
gg_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

########################
## Blank plot
########################

gg_blank <-
  ggplot(data.frame()) +
  geom_blank()

#########################

########################
## Plot discordant CpGs
########################
summary_freq$concordance = gsub("cnv_driven_change", "CNV-driven meth change", summary_freq$concordance)
summary_freq$concordance = gsub("meth_gain", "Meth. gain only", summary_freq$concordance)
summary_freq$concordance = gsub("meth_loss", "Meth. loss only", summary_freq$concordance)
summary_freq$concordance = gsub("no_change", "No meth. change", summary_freq$concordance)

gg_stacked <- ggplot(summary_freq, aes(x=case_barcode, fill = factor(concordance), y=freq)) +
  geom_bar(stat = "identity") +
  labs(y = "Fraction of longitudinal \n methylation differences", fill="Methylation change") +
  scale_fill_manual(values = c("#d7191c", "#fdae61", "#2b83ba", "gray70"))

testPlot(gg_stacked) 

########################
## Plot proportions shared/private per patient
########################
gg_mut_freq_prop_case <-
  ggplot(mut_freq_prop_case, aes(x = case_barcode, y = mutation_percent, fill=mutation_type)) +
  geom_bar(stat="identity") +
  labs(y = "% mutation") +
  scale_fill_manual(values=c("R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F")) +
  plot_grid

testPlot(gg_mut_freq_prop_case)


########################
### Purity
#######################
gg_puritydata <-
  ggplot(puritydata, aes(x=case_barcode)) +
  geom_bar(aes(y =  purity_b - purity_a), fill = "#e34a33", alpha = 1, stat = "identity") +
  coord_cartesian(ylim=c(-0.5,0.5)) +
  labs(y = "Tumor puity\ndifference")

testPlot(gg_puritydata)


########################
## Copy Number Variants
########################
gg_anpl_data <-
  ggplot(anpldata, aes(x=case_barcode)) +
  geom_bar(aes(y = aneuploidy_b - aneuploidy_a), fill = "#beaed4", alpha = 1, stat = "identity") +
  coord_cartesian(ylim=c(-0.35,0.35)) +
  labs(y = "SCNA burden\ndifference")

testPlot(gg_anpl_data)


#######################
### RNA
#######################
subtype_order <- c("Proneural", "Classical", "Mesenchymal")

## Need to profive additional ordering for MNP classes
rnadata_ordered = rna_data %>%
  gather(key = "type", value = "value", subtype_a, subtype_b) %>%
  mutate(type = factor(type,
                       levels = c("subtype_b", "subtype_a"),
                       labels = c("Recurrence","Primary"))) 
rnadata_ordered <- rnadata_ordered %>% mutate(value = factor(value, levels = subtype_order))

gg_rna <-
  rnadata_ordered %>%
  ggplot(aes(x=case_barcode)) +
  geom_tile(aes(fill = factor(value), y = type)) +
  scale_fill_manual(values = c("Proneural" = "#00458A",
                               "Classical" = "#008A22",
                               "Mesenchymal" = "#8A0000")) +
  labs(y = "RNA subtype")
testPlot(gg_rna)


#######################
### MNP
#######################
class_order <- c("CONTR, HEMI", "CONTR, INFLAM",
                 "O IDH", "A IDH", "A IDH, HG",
                 "GBM, RTK I", "GBM, RTK II", "GBM, MID", "GBM, MES",
                 "HGNET, BCOR")

## Need to profive additional ordering for MNP classes
  mnpdata_ordered = mnpdata %>%
  select(-subtype_stability) %>% 
  gather(key = "type", value = "value", mnp_a, mnp_b) %>%
  mutate(type = factor(type,
                       levels = c("mnp_b", "mnp_a"),
                       labels = c("Recurrence","Primary"))) 
  mnpdata_ordered <- mnpdata_ordered %>% mutate(value = factor(value, levels = class_order))
  
  gg_mnp <-
    mnpdata_ordered %>%
    ggplot(aes(x=case_barcode)) +
    geom_tile(aes(fill = factor(value), y = type)) +
    scale_fill_manual(values = c("O IDH" = "#9ecae1", "A IDH" = "#4292c6", "A IDH, HG" = "#08306b",
                               "CONTR, HEMI" = "#d9d9d9", "CONTR, INFLAM" = "#737373",
                               "GBM, MES" = "#a50026", "GBM, MID" = "#d73027", "GBM, RTK I" = "#fee090", "GBM, RTK II"="#fdae61",
                               "HGNET, BCOR" = "#a6d96a")) +
    labs(y = "Methylation\nclassification")
testPlot(gg_mnp)


########################
## Plot clinical
########################

gg_clinical <-
  clin_data %>% 
  mutate(grade_change = ifelse(grade_change=="Grade up", 1, ifelse(grade_change=="Grade stable", 0, "NA"))) %>% 
  gather(key = "type", value = "value", grade_change, received_treatment) %>%
  mutate(type = factor(type,
                       levels = c("grade_change", "received_treatment"),
                       labels = c("Grade increase", "Received Treatment"))) %>%
  ggplot(aes(x=case_barcode)) +
  geom_tile(aes(fill = factor(value), y = type)) +
  scale_fill_manual(values=c("white", "#4d4d4d"), na.value = "gray90") +
  labs(y="", fill = "Event")

testPlot(gg_clinical)


########################
## mDNAsi
########################
gg_mDNAsidata <-
  ggplot(mDNAsidata, aes(x=case_barcode)) +
  geom_bar(aes(y = mDNAsi_b - mDNAsi_a), fill = "#2b8cbe", alpha = 1, stat = "identity") +
  coord_cartesian(ylim=c(-0.5,0.5)) +
  labs(y = "Stemness\ndifference")

testPlot(gg_mDNAsidata)


########################
## epiTOC
########################
gg_epi_TOCdata <-
  ggplot(epi_TOCdata, aes(x=case_barcode)) +
  geom_bar(aes(y = epiTOC_b - epiTOC_a), fill = "#2b8cbe", alpha = 1, stat = "identity") +
  coord_cartesian(ylim=c(-0.25,0.25)) +
  labs(y = "epiTOC\ndifference")

testPlot(gg_epi_TOCdata)

########################
## Horvath clock
########################
gg_horvathdata <-
  ggplot(horvathdata, aes(x=case_barcode)) +
  geom_bar(aes(y = DNAmAge_b - DNAmAge_a), fill = "#2b8cbe", alpha = 1, stat = "identity") +
  coord_cartesian(ylim=c(-50,50)) +
  labs(y = "DNAmAge\ndifference")

testPlot(gg_horvathdata)

########################
## eITH
########################

gg_eITHdata <-
  eITHdata %>%
  gather(key = "type", value = "value", eITH_a, eITH_b) %>%
  mutate(type = factor(type,
                       levels = c("eITH_b", "eITH_a"),
                       labels = c("Recurrence", "Initial"))) %>%
  ggplot(aes(x=case_barcode)) +
  geom_tile(aes(fill = value, y = type)) +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(y = "DNAme disorder")

testPlot(gg_eITHdata)



##################################################
## Final combined plots
##################################################

## Legends
gleg1 = g_legend(gg_stacked)
gleg2 = g_legend(gg_eITHdata)
gleg3 = g_legend(gg_anpl_data) 


## Plot v1:
g1 = ggplotGrob(gg_stacked + plot_grid + plot_theme + null_legend + null_x + top_margin)  %>% gtable_frame()
g2 = ggplotGrob(gg_eITHdata + plot_grid + plot_theme + null_legend + null_x + null_facet + middle_margin)  %>% gtable_frame()
g3 = ggplotGrob(gg_anpl_data + plot_grid + plot_theme + null_legend + null_x + null_facet + bottom_margin)  %>% gtable_frame()


g = gtable_rbind(g1, g2, g3)
gleg = gtable_rbind(gleg1, gleg2)

## Adjust relative height of panels
panels = g$layout$t[grep("panel", g$layout$name)]
g$heights[panels] <- unit(c(0.6, 0.2, 0.2), "null")

plot(g)
plot(gleg)


pdf(file = "results/Fig7/Fig7e-longitudinal-DNAme.pdf", width = 6, height = 6, useDingbats = FALSE, bg="transparent")
grid.newpage()
grid.draw(g)
dev.off()

### END ###