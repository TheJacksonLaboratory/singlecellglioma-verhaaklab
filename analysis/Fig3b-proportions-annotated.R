##############################################
# Generate a heatmap with bulk and single-cell RNA data
# Updated: 2020.05.31
# Author: Kevin J.
##################################################

## Project working directory
mybasedir = "/Users/johnsk/github/"
setwd(mybasedir)

########################################
# Necessary packages:
library(tidyverse)
library(grid)
library(gridExtra)
library(gtable)
library(egg)
library(tidyverse)
library(RColorBrewer)
library(DBI)
library(vegan)
########################################


######### CLINICAL ##########
## Load in clinical data (subtype, age, treatment, hypermutation_status).
meta_data = read.csv("data/clinical_metadata.csv", sep = ",", header = TRUE)
clin_data = meta_data %>% 
  mutate(idh_status = ifelse(idh_codel_subtype == "IDHwt", "IDHwt", "IDHmut"),
         is_hypermutator = ifelse(case_barcode == "SM011", 1, 0)) %>% 
  select(case_barcode, idh_status, idh_codel_subtype, grade = who_grade, timepoint = time_point,  is_hypermutator)


########### SINGLE CELL #############
## Load in proportions of cell types as defined by all captured AND by Neftel classifications.
neftel_class = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/10X/cell-state-labels-IDHwt.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
neftel_prop = neftel_class %>% 
  select(case_barcode = sample_id, class) %>% 
  left_join(clin_data, by = "case_barcode") %>% 
  select(case_barcode, idh_status, class) %>% 
  group_by(case_barcode, idh_status, class) %>% 
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n))

venteicher_class = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/10X/cell-state-labels-IDHmut.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
venteicher_class$case_barcode <- gsub("UC917", "SM019", venteicher_class$case_barcode)
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
## 2D UMAP coordinates.
umap_coords_2d <- read.csv("data/analysis_scRNAseq_tumor_metadata.csv", sep = ",", header = TRUE)
clust_annot <- umap_coords_2d
cells_to_keep = which(clust_annot$cell_state%in%c("Diff.-like", "Stem-like", "Prolif. stem-like"))
clust_annot <- clust_annot[cells_to_keep, ]


######### Shannon diversity for JAX states ###########
clust_annot = clust_annot %>%  
  group_by(case_barcode, cell_state) %>%
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>% 
  left_join(clin_data, by = "case_barcode") %>% 
  ungroup() %>% 
  select(case_barcode, idh_status, cell_state, freq)

cell_state_summary <- clust_annot %>% 
  spread(cell_state, freq) 
sample_names <- cell_state_summary$case_barcode

## Tabulate the Shannon diversity index for cell type composition.
shannon_div <- diversity(cell_state_summary[,3:5], index = "shannon", MARGIN = 1)
names(shannon_div) <- sample_names
shannon_div_df = cell_state_summary[, 1]
shannon_div_df$entropy = as.numeric(shannon_div)


######## DNAme disorder rate ###########
disorder_cpg <- read.table(file="data/analysis_scRRBS_context_specific_DNAme_disorder.csv", sep = ",", header = TRUE)

## Additional information about single-cells passing QC.
rrbs_qc <- read.table(file="data/analysis_scRRBS_sequencing_qc.csv", sep = ",", header = TRUE)

## Restrict to the samples that pass the general QC guidelines of CpG count, bisulfite conversion, and tumor status (inferred CNVs).
rrbs_qc_pass <- rrbs_qc %>% 
  filter(cpg_unique > 40000, bisulfite_conversion_rate > 95, tumor_status == 1) %>% 
  left_join(clin_data, by="case_barcode") 

## Calculate the mean disorder for each tumor.
disorder_mean <- disorder_cpg %>% 
  filter(cell_barcode%in%rrbs_qc_pass$cell_barcode) %>% 
  group_by(case_barcode) %>% 
  summarise(disorder_mean = mean(global_PDR)) %>% 
  ungroup() %>% 
  select(case_barcode, disorder = disorder_mean)


######## Mutation Freq. ###########
me_disorder_hmapdata = read.table("/Users/johnsk/mnt/verhaak-lab/scgp/results/mutation_vs_epimutation/Samples_WGS_epimutation_and_mutation_rates.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
me_disorder_hmapdata = me_disorder_hmapdata %>% 
  select(case_barcode = Patient, idh_codel_subtype = subtype, timepoint = initial_recurrence, mut_freq = global)
#me_disorder_hmapdata$case_barcode <- gsub("UC917", "SM019", me_disorder_hmapdata$case_barcode)

######### Aneuploidy CNVs #########
aneuploidy_hmapdata = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/wgs/scgp_aneuploidy_hmp_data.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
aneuploidy_hmapdata$case_barcode <- gsub("UC917", "SM019", aneuploidy_hmapdata$case_barcode)


########### Combine annotation meta data ############
meta_comb = clin_data %>% 
  left_join(me_disorder_hmapdata, by=c("case_barcode", "idh_codel_subtype")) %>% 
  left_join(aneuploidy_hmapdata, by=c("case_barcode", "idh_codel_subtype")) %>% 
  inner_join(shannon_div_df, by="case_barcode") %>% 
  inner_join(disorder_mean, by="case_barcode") %>% 
  dplyr::select(case_barcode,idh_status, timepoint = timepoint.x,  entropy, disorder, mut_freq, aneuoploidy_value)



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
case_order <- c("SM004", "SM001", "SM015", "SM019", "SM002", "SM008", "SM006", "SM012", "SM017", "SM018", "SM011")
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
disorder_mean <- disorder_mean %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
disorder_mean = disorder_mean %>% 
  inner_join(clin_data, by="case_barcode")
#### Plot Simplicity scores.
gg_me_disorder = ggplot() +
  geom_tile(data = disorder_mean, aes(x = case_barcode, y = 1, fill = disorder), color = "black") +
  labs(y = "disorder", fill = "disorder") +
  scale_fill_distiller(palette = "Blues", direction = 1, na.value = "gray70")

testPlot(gg_me_disorder)


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
clust_annot_data$cell_state <- factor(clust_annot_data$cell_state, levels = c("Diff.-like", "Stem-like", "Prolif. stem-like"))

gg_clust_annot <-
  ggplot(clust_annot_data, aes(x = case_barcode, y = freq, fill=cell_state)) +
  geom_bar(stat="identity") +
  labs(y = "Proportion tumor cells") +
  scale_fill_manual(values=c("Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
  plot_grid

testPlot(gg_clust_annot)

## Calculate the proportion of stem-like and diff-like per subtype:
clust_annot_data %>% 
  group_by(idh_status, cell_state) %>% 
  summarise(med_freq = median(freq))

## Test for statistical significance.
wilcox.test(clust_annot_data$freq[clust_annot_data$cell_state=="Prolif. stem-like"]~clust_annot_data$idh_status[clust_annot_data$cell_state=="Prolif. stem-like"])

###########################
### Mutation frequency
###########################
meta_comb <- meta_comb %>% mutate(case_barcode = factor(case_barcode, levels = case_order))

## Set an upper limit so that hypermutator doesn't make other comparisons impossible.
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
### disorder rate
###########################
gg_me_disorder = ggplot() +
  geom_tile(data = meta_comb, aes(x = case_barcode, y = 1, fill = disorder), color = "black") +
  labs(y = "", fill = "Cell type diversity") +
  scale_fill_distiller(palette = "Greens", direction = 1)
testPlot(gg_me_disorder)

##################################################
## Final combined plots
##################################################

# ## Legends
gleg1 = g_legend(gg_cell_states) 
gleg2 = g_legend(gg_clust_annot)
gleg3 = g_legend(gg_shannon) 
gleg4 = g_legend(gg_me_disorder) 
gleg5 = g_legend(gg_mf) 
gleg6 = g_legend(gg_cin) 

## Plots
g1 = ggplotGrob(gg_clust_annot + plot_grid + plot_theme + null_legend + null_x + top_margin)  %>% gtable_frame()
#g2 = ggplotGrob(gg_cell_states + plot_grid + plot_theme + null_legend + null_x + null_facet + middle_margin) %>% gtable_frame()
g3 = ggplotGrob(gg_shannon + plot_grid + plot_theme + null_legend + null_x + null_y + null_facet + middle_margin) %>% gtable_frame()
g4 = ggplotGrob(gg_me_disorder + plot_grid + plot_theme + null_legend + null_x + null_y + null_facet + middle_margin) %>% gtable_frame()
g5 = ggplotGrob(gg_cin + plot_grid + plot_theme + null_legend + null_x + null_y + null_facet + middle_margin) %>% gtable_frame()
g6 = ggplotGrob(gg_mf + plot_grid + plot_theme + null_legend + null_y + null_facet + bottom_margin)  %>% gtable_frame()

g = gtable_rbind(g1, g3, g4, g5, g6)
gleg = gtable_rbind(gleg1, gleg2, gleg3, gleg4, gleg5, gleg6)

## Adjust relative height of panels
panels = g$layout$t[grep("panel", g$layout$name)]
g$heights[panels] <- unit(c(0.5, 0.05, 0.05, 0.05, 0.05), "null")

plot(g)
plot(gleg)

pdf(file = "results/Fig3/Fig3b-cell-states-annotated.pdf", width = 7, height = 5, useDingbats = FALSE, bg="transparent")
grid.newpage()
grid.draw(g)
dev.off()

pdf(file = "results/Fig3/cell-states-annotated-legend.pdf", width = 4, height = 11)
grid.newpage()
grid.draw(gleg)
dev.off()

################################
# Calculate the correlations
################################
# Split into IDHmut vs. IDHwt
meta_comb_idh = meta_comb %>% filter(idh_status=="IDHmut")
meta_comb_wt = meta_comb %>% filter(idh_status=="IDHwt")

## Shannon entropy vs. other summary metrics (all).
cor.test(meta_comb$entropy, meta_comb$disorder, method="s")
cor.test(meta_comb$entropy, meta_comb$aneuoploidy_value, method="s")
cor.test(meta_comb$entropy, meta_comb$mut_freq, method="s")

## Shannon entropy vs. other summary metrics (IDHmut).
cor.test(meta_comb_idh$entropy, meta_comb_idh$disorder, method="s")
cor.test(meta_comb_idh$entropy, meta_comb_idh$aneuoploidy_value, method="s")
cor.test(meta_comb_idh$entropy, meta_comb_idh$mut_freq, method="s")

## Shannon entropy vs. other summary metrics (IDHwt).
cor.test(meta_comb_wt$entropy, meta_comb_wt$disorder, method="s")
cor.test(meta_comb_wt$entropy, meta_comb_wt$aneuoploidy_value, method="s")
cor.test(meta_comb_wt$entropy, meta_comb_wt$mut_freq, method="s")

### END ###