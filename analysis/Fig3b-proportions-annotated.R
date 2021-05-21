##############################################
# Generate a heatmap with bulk and single-cell RNA data
# Updated: 2020.05.31
# Author: Kevin J.
##################################################

########################################
# Necessary packages:
library(tidyverse)
library(grid)
library(gridExtra)
library(gtable)
library(egg)
library(tidyverse)
library(RColorBrewer)
library(openxlsx)
library(DBI)
library(vegan)
########################################


######### CLINICAL ##########
## Load in clinical data (subtype, age, treatment, hypermutation_status).
# Supply metadata so that 10X filenames and samples can be linked together.
metadata = readWorkbook("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
metadata$`10X_id_short` <- as.character(metadata$`10X_id_short`)
clin_data = metadata %>% 
  mutate(case_barcode = gsub("-", "", substr(subject_id, 6, 11))) %>% 
  mutate(idh_status = ifelse(subtype == "IDHwt", "IDHwt", "IDHmut")) %>% 
  select(case_barcode,idh_status, idh_codel_subtype = subtype, grade = who_grade, timepoint = initial_recurrence,  is_hypermutator = hypermutation)


########### SINGLE CELL #############
## Load in proportions of cell types as defined by all captured AND by Suva classifications.
neftel_class = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/10X/cell-state-labels-IDHwt.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
neftel_prop = neftel_class %>% 
  select(case_barcode = sample_id, class) %>% 
  left_join(clin_data, by = "case_barcode") %>% 
  select(case_barcode, idh_status, class) %>% 
  group_by(case_barcode, idh_status, class) %>% 
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n))

venteicher_class = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/10X/cell-state-labels-IDHmut.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
venteicher_prop = venteicher_class %>% 
  select(case_barcode, class) %>% 
  left_join(clin_data, by = "case_barcode") %>% 
  select(case_barcode, idh_status, class) %>% 
  group_by(case_barcode, idh_status, class) %>% 
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n))

## All subtypes together:
cell_states <- rbind(neftel_prop, venteicher_prop) %>% ungroup()

## scRNAseq ##
# Load the 10X data for all samples.
load("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/10X-aggregation-20190905/RV18001-2-3-RV19001-2-4-5-6-7-8-9_20190827.Rds")

# Strip the rownames of the sample-specific identifier. It's numeric and assumed to be in the same order as the 10X samples.
tsne.data$sample_id = sapply(strsplit(rownames(tsne.data), "-"), "[[", 3)

# Create a new sample identifier based on the marker gene expression.
clust_annot = tsne.data %>% 
  mutate(cell_type = recode(dbCluster, `1` = "differentiated_tumor",  `2` = "myeloid", `3` = "stemcell_tumor",
                            `4` = "oligodendrocyte", `5` = "prolif_stemcell_tumor", `6` = "granulocyte", `7` = "endothelial",
                            `8` = "t_cell", `9` = "pericyte", `10` = "fibroblast", `11` = "b_cell", `12` = "dendritic_cell")) 
cells_to_keep = which(clust_annot$cell_type%in%c("differentiated_tumor", "stemcell_tumor", "prolif_stemcell_tumor"))
clust_annot <- clust_annot[cells_to_keep, ]


######### Shannon diversity for JAX states ###########
clust_annot = clust_annot %>%  
  group_by(sample_id, cell_type) %>%
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>% 
  left_join(metadata, by = c("sample_id"="10X_id_short")) %>% 
  mutate(case_barcode = gsub("-", "", substr(subject_id, 6, 11))) %>% 
  mutate(idh_status = ifelse(subtype == "IDHwt", "IDHwt", "IDHmut")) %>% 
  ungroup() %>% 
  select(case_barcode, idh_status, cell_type, freq)

cell_type_summary <- clust_annot %>% 
  spread(cell_type, freq) 
sample_names <- cell_type_summary$case_barcode

## Tabulate the Shannon diversity index for cell type composition.
shannon_div <- diversity(cell_type_summary[,3:5], index = "shannon", MARGIN = 1)
names(shannon_div) <- sample_names
shannon_div_df = cell_type_summary[, 1]
shannon_div_df$entropy = as.numeric(shannon_div)


######## Epimutation rate ###########
epimut_cpg <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/methylation/Samples_pdr_score_table-all_single_cells.txt", sep="\t", header=T, stringsAsFactors = F)

## Load in the quality control data for these samples.
rrbs_qc <- read.table(file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/scgp_cnv_status.txt", sep="\t", header=T, stringsAsFactors = F)

## Restrict to the samples that pass the general QC guidelines of CpG count, bisulfite conversion, and cell number.
rrbs_qc_pass <- rrbs_qc %>% 
  filter(cell_num == 1, cpg_unique > 40000, conversion_rate > 95, tumor_cnv == 1) %>% 
  # remove the single-end libraries as these may have founding effects.
  filter(!library_id %in% c("SCGP-UC-917-1-1", "SCGP-UC-917-1-2", "SCGP-UC-917-1-3", "SCGP-UC-917-2-2")) %>% 
  filter(!grepl("SCGP-HF-", sample_barcode)) %>%
  mutate(case_barcode = gsub("-", "", substr(sample, 6, 11))) %>% 
  left_join(clin_data, by="case_barcode") 

## Calculate the mean epimutation for each tumor.
epimut_mean <- epimut_cpg %>% 
  filter(Sample%in%rrbs_qc_pass$sample_barcode) %>% 
  group_by(Patient) %>% 
  summarise(epimut_mean = mean(PDR)) %>% 
  ungroup() %>% 
  select(case_barcode = Patient, epimutation = epimut_mean)


######## Mutation Freq. ###########
epi_mut_hmapdata = read.table("/Users/johnsk/mnt/verhaak-lab/scgp/results/mutation_vs_epimutation/Samples_WGS_epimutation_and_mutation_rates.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
epi_mut_hmapdata = epi_mut_hmapdata %>% 
  select(case_barcode = Patient, idh_codel_subtype = subtype, timepoint = initial_recurrence, mut_freq = global)

######### Aneuploidy CNVs #########
aneuploidy_hmapdata = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/wgs/scgp_aneuploidy_hmp_data.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)



########### Combine annotation meta data ############
meta_comb = clin_data %>% 
  left_join(epi_mut_hmapdata, by=c("case_barcode", "idh_codel_subtype")) %>% 
  left_join(aneuploidy_hmapdata, by=c("case_barcode", "idh_codel_subtype")) %>% 
  inner_join(shannon_div_df, by="case_barcode") %>% 
  inner_join(epimut_mean, by="case_barcode") %>% 
  mutate(idh_status = ifelse(idh_codel_subtype == "IDHwt", "IDHwt", "IDHmut")) %>% 
  select(case_barcode:grade,idh_status, timepoint = timepoint.x,  entropy, epimutation, mut_freq, aneuoploidy_value)



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


plot_grid     <- facet_grid(. ~ idh_status, scales = "free_x", space = "free")
## Generate plot theme.
plot_theme    <- theme_bw(base_size = 12) + theme(axis.title = element_text(size = 12),
                                                  axis.text = element_text(size = 12),
                                                  panel.background = element_rect(fill = "transparent"),
                                                  axis.text.x = element_text(angle=45, hjust=1),
                                                  axis.line = element_blank(),
                                                  strip.background = element_blank(),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(), 
                                                  panel.border = element_blank(),
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

######################
## Clinical data  ####
######################
case_order <- c("SM004", "SM001", "SM015", "UC917", "SM002", "SM008", "SM006", "SM012", "SM017", "SM018", "SM011")
clin_data <- clin_data %>% mutate(case_barcode = factor(case_barcode, levels = case_order))

########################
## Plot clinical
########################

gg_clinical <-
  clin_data %>% 
  gather(key = "type", value = "value", grade, timepoint, is_hypermutator) %>%
  mutate(type = factor(type,
                       levels = c("grade","timepoint",  "is_hypermutator"),
                       labels = c("WHO Grade", "Timepoint", "Hypermutator"))) %>%
  ggplot(aes(x=case_barcode)) +
  geom_tile(aes(fill = factor(value), y = type)) +
  scale_fill_manual(values = c("0" = "white", "1" = "#377eb8", "initial" = "#CA2F66", "recurrence" = "#2FB3CA",
                               "II" = "#fee0d2", "III" = "#fc9272", "IV" = "#de2d26")) +
  labs(y="", fill = "Clinical")

## Relabel some of the variables.
gg_clinical$data$value <- factor(gg_clinical$data$value, levels = c("0", "1", "initial", "recurrence", "II",
                                                                    "III", "IV"))


testPlot(gg_clinical)

#########################
## Plot subtype
#########################
epimut_mean <- epimut_mean %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
epimut_mean = epimut_mean %>% 
  inner_join(clin_data, by="case_barcode")
#### Plot Simplicity scores.
gg_epimut = ggplot() +
  geom_tile(data = epimut_mean, aes(x = case_barcode, y = 1, fill = epimutation), color = "black") +
  labs(y = "Epimutation", fill = "Epimutation") +
  scale_fill_distiller(palette = "Blues", direction = 1, na.value = "gray70")

testPlot(gg_epimut)


#########################
## Plot cell states (Suva).
#########################
cell_states_data <- cell_states %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
cell_states_data$class <- factor(cell_states_data$class, levels = c("Oligo-like", "Astro-like", "Undifferentiated",
                                                                        "AC", "MES", "NPC", "OPC"))

gg_cell_states <-
  ggplot(cell_states_data, aes(x = case_barcode, y = freq, fill=class)) +
  geom_bar(stat="identity") +
  labs(y = "% Tumor Cell states") +
  scale_fill_manual(values=c("Oligo-like" = "#a6dba0", "Astro-like" ="#008837", "Undifferentiated"="#7b3294",
                             "MES" = "#c5168a", "AC" = "#fcc5c0", "NPC" = "#fdc086", "OPC" = "#018571")) +
  plot_grid

testPlot(gg_cell_states)

#########################
## Plot cell states (SCGP).
#########################
clust_annot_data <- clust_annot %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
clust_annot_data$cell_type <- factor(clust_annot_data$cell_type, levels = c("differentiated_tumor", "stemcell_tumor", "prolif_stemcell_tumor"))

gg_clust_annot <-
  ggplot(clust_annot_data, aes(x = case_barcode, y = freq, fill=cell_type)) +
  geom_bar(stat="identity") +
  labs(y = "Proportion tumor cells") +
  scale_fill_manual(values=c("stemcell_tumor" = "#fb6a4a", "differentiated_tumor" = "#fcbba1", "prolif_stemcell_tumor" = "#a50f15")) +
  plot_grid

testPlot(gg_clust_annot)

## Calculate the proportion of stem-like and diff-like per subtype:
clust_annot_data %>% 
  group_by(idh_status, cell_type) %>% 
  summarise(med_freq = median(freq))

## Test for statistical significance.
wilcox.test(clust_annot_data$freq[clust_annot_data$cell_type=="prolif_stemcell_tumor"]~clust_annot_data$idh_status[clust_annot_data$cell_type=="prolif_stemcell_tumor"])

###########################
### Mutation frequency
###########################
meta_comb <- meta_comb %>% mutate(case_barcode = factor(case_barcode, levels = case_order))

meta_comb$mut_freq[meta_comb$case_barcode=="SM011"] <- 15
gg_mf = ggplot() +
  geom_tile(data = meta_comb, aes(x = case_barcode, y = 1, fill = mut_freq), color = "black") +
  labs(y = "", x="", fill = "Mutations/Mb") +
  scale_fill_distiller(palette = "Blues", direction = 1)
testPlot(gg_mf)

###########################
### Chromosomal instability
###########################
gg_cin = ggplot() +
  geom_tile(data = meta_comb, aes(x = case_barcode, y = 1, fill = aneuoploidy_value), color = "black") +
  labs(y = "", fill = "SCNA burden") +
  scale_fill_distiller(palette = "Reds", direction = 1)
testPlot(gg_cin)

###########################
### Shannon entropy value
###########################
gg_shannon = ggplot() +
  geom_tile(data = meta_comb, aes(x = case_barcode, y = 1, fill = entropy), color = "black") +
  labs(y = "", fill = "Cell type diversity") +
  scale_fill_distiller(palette = "Oranges", direction = 1)
testPlot(gg_shannon)

###########################
### Epimutation rate
###########################
gg_epimut = ggplot() +
  geom_tile(data = meta_comb, aes(x = case_barcode, y = 1, fill = epimutation), color = "black") +
  labs(y = "", fill = "Cell type diversity") +
  scale_fill_distiller(palette = "Greens", direction = 1)
testPlot(gg_epimut)

##################################################
## Final combined plots
##################################################

# ## Legends
gleg1 = g_legend(gg_cell_states) 
gleg2 = g_legend(gg_clust_annot)
gleg3 = g_legend(gg_shannon) 
gleg4 = g_legend(gg_epimut) 
gleg5 = g_legend(gg_mf) 
gleg6 = g_legend(gg_cin) 

## Plots
g1 = ggplotGrob(gg_clust_annot + plot_grid + plot_theme + null_legend + null_x + top_margin)  %>% gtable_frame()
#g2 = ggplotGrob(gg_cell_states + plot_grid + plot_theme + null_legend + null_x + null_facet + middle_margin) %>% gtable_frame()
g3 = ggplotGrob(gg_shannon + plot_grid + plot_theme + null_legend + null_x + null_y + null_facet + middle_margin) %>% gtable_frame()
g4 = ggplotGrob(gg_epimut + plot_grid + plot_theme + null_legend + null_x + null_y + null_facet + middle_margin) %>% gtable_frame()
g5 = ggplotGrob(gg_cin + plot_grid + plot_theme + null_legend + null_x + null_y + null_facet + middle_margin) %>% gtable_frame()
g6 = ggplotGrob(gg_mf + plot_grid + plot_theme + null_legend + null_y + null_facet + bottom_margin)  %>% gtable_frame()

g = gtable_rbind(g1, g3, g4, g5, g6)
gleg = gtable_rbind(gleg1, gleg2, gleg3, gleg4, gleg5, gleg6)

## Adjust relative height of panels
panels = g$layout$t[grep("panel", g$layout$name)]
g$heights[panels] <- unit(c(0.5, 0.05, 0.05, 0.05, 0.05), "null")

plot(g)
plot(gleg)

pdf(file = "github/results/Fig3/Fig3b-cell-states-annotated.pdf", width = 7, height = 5, useDingbats = FALSE, bg="transparent")
grid.newpage()
grid.draw(g)
dev.off()

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/manuscript/figure-drafts/Fig2/cell-states-annotated-legend.pdf", width = 4, height = 11)
grid.newpage()
grid.draw(gleg)
dev.off()

################################
# Calculate the correlations
################################
# Split into IDHmut vs. IDHwt
meta_comb$idh_status <- ifelse(meta_comb$idh_codel_subtype=="IDHwt", "IDHwt", "IDHmut")
meta_comb_idh = meta_comb %>% filter(idh_status=="IDHmut")
meta_comb_wt = meta_comb %>% filter(idh_status=="IDHwt")

## Shannon entropy vs. other summary metrics (all).
cor.test(meta_comb$entropy, meta_comb$epimutation, method="s")
cor.test(meta_comb$entropy, meta_comb$aneuoploidy_value, method="s")
cor.test(meta_comb$entropy, meta_comb$mut_freq, method="s")

## Shannon entropy vs. other summary metrics (IDHmut).
cor.test(meta_comb_idh$entropy, meta_comb_idh$epimutation, method="s")
cor.test(meta_comb_idh$entropy, meta_comb_idh$aneuoploidy_value, method="s")
cor.test(meta_comb_idh$entropy, meta_comb_idh$mut_freq, method="s")

## Shannon entropy vs. other summary metrics (IDHwt).
cor.test(meta_comb_wt$entropy, meta_comb_wt$epimutation, method="s")
cor.test(meta_comb_wt$entropy, meta_comb_wt$aneuoploidy_value, method="s")
cor.test(meta_comb_wt$entropy, meta_comb_wt$mut_freq, method="s")

### END ###