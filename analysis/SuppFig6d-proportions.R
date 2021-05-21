##############################################
# Generate stacked barplots for supplementary materials
# Updated: 2021.02.06
# Author: Kevin J.
##################################################

# Working directory for this analysis in the SCGP-analysis project. 
mybasedir = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/"
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
library(openxlsx)
########################################
## Generate plot theme.
plot_theme    <- theme_bw(base_size = 12) + theme(axis.title = element_text(size = 12),
                                                  axis.text = element_text(size = 12),
                                                  axis.text.x = element_text(angle=45, hjust=1),
                                                  panel.background = element_rect(fill = "transparent"),
                                                  axis.line = element_blank(),
                                                  strip.background = element_blank(),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(), 
                                                  panel.border = element_blank(),
                                                  axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                                  axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))


######### CLINICAL ##########
## Load in clinical data (subtype, age, treatment, hypermutation_status).
# Supply metadata so that 10X filenames and samples can be linked together.
metadata = readWorkbook("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
metadata$`10X_id_short` <- as.character(metadata$`10X_id_short`)
clin_data = metadata %>% 
  mutate(case_barcode = gsub("-", "", substr(subject_id, 6, 11)),
         idh_status = ifelse(subtype=="IDHwt", "IDHwt", "IDHmut")) %>% 
  select(case_barcode, idh_status, grade = who_grade, timepoint = initial_recurrence,  is_hypermutator = hypermutation)


########### SINGLE CELL #############
## Load in proportions of cell types as defined by all captured AND by Suva classifications.
neftel_class = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/10X/cell-state-labels-IDHwt.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
neftel_class = neftel_class %>% 
  mutate(cell_type = recode(cell_type, "Prolif.-Stem-like" = "Prolif. stem-like"))

neftel_prop = neftel_class %>% 
  mutate(class_neftel = paste(class, "like", sep="-")) %>% 
  select(case_barcode = sample_id, class_neftel) %>% 
  left_join(clin_data, by = "case_barcode") %>% 
  select(case_barcode, idh_status, class_neftel) %>% 
  group_by(case_barcode, idh_status, class_neftel) %>% 
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n))

cell_state_order <- c("OPC-like", "NPC-like", "MES-like", "AC-like")
neftel_prop <- neftel_prop %>% mutate(class_neftel = factor(class_neftel, levels = cell_state_order))

ggplot(neftel_prop, aes(x = case_barcode, y = freq, fill=class_neftel)) +
  geom_bar(stat="identity") +
  labs(x="", y = "Proportion tumor\ncell state", fill="Neftel\ncell state") +
  scale_fill_manual(values=c("MES-like" = "#d7191c", "AC-like" = "#fdae61",
                             "NPC-like" = "#abd9e9", "OPC-like" = "#2c7bb6")) +
  plot_theme 


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

neftel_class_merge <- neftel_class %>% 
  mutate(suva_class = paste(class, "like", sep="-")) %>% 
  select(cell_name, sample_id, suva_class, pan_glioma = cell_type)
venteicher_class_merge <- venteicher_class %>% 
  select(cell_name, sample_id, suva_class = class, cell_type) %>% 
  mutate(pan_glioma = recode(cell_type, "differentiated_tumor" = "Diff.-like", "prolif_stemcell_tumor" = "Prolif. stem-like", "stemcell_tumor" = "Stem-like"))
all_states <- bind_rows(neftel_class_merge, venteicher_class_merge)

pan_glioma_class_prop <- all_states %>% 
  group_by(sample_id, pan_glioma) %>%
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>% 
  left_join(clin_data, by = c("sample_id"="case_barcode")) %>% 
  ungroup() %>% 
  select(case_barcode = sample_id, idh_status, pan_glioma, n , freq)


suva_class_prop <- all_states %>% 
  group_by(sample_id, suva_class) %>%
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>% 
  left_join(clin_data, by = c("sample_id"="case_barcode")) %>% 
  ungroup() %>% 
  select(case_barcode = sample_id, idh_status, suva_class, n , freq)


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
## Plot cell states (Suva).
#########################
suva_class_prop <- suva_class_prop %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
suva_class_prop$suva_class <- factor(suva_class_prop$suva_class, levels = c("Astro-like", "Oligo-like", "Undifferentiated",
                                                                    "MES-like", "AC-like", "OPC-like", "NPC-like"))

gg_cell_states <-
  ggplot(suva_class_prop, aes(x = case_barcode, y = freq, fill=suva_class)) +
  geom_bar(stat="identity") +
  labs(y = "Proportion tumor state") +
  scale_fill_manual(values=c("Oligo-like" = "#a6dba0", "Astro-like" ="#008837", "Undifferentiated"="#7b3294",
                             "MES-like" = "#d7191c", "AC-like" = "#fdae61",
                             "NPC-like" = "#abd9e9", "OPC-like" = "#2c7bb6")) +
  plot_grid

testPlot(gg_cell_states)

#########################
## Plot cell states (SCGP).
#########################
pan_glioma_class_prop <- pan_glioma_class_prop %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
state_levels <- c("Diff.-like", "Stem-like", "Prolif. stem-like")
pan_glioma_class_prop$pan_glioma <- factor(pan_glioma_class_prop$pan_glioma, levels = state_levels)

gg_clust_annot <-
  ggplot(pan_glioma_class_prop, aes(x = case_barcode, y = freq, fill=pan_glioma)) +
  geom_bar(stat="identity") +
  labs(y = "Proportion tumor state") +
  scale_fill_manual(values=c("Stem-like" = "#fb6a4a", 
                             "Diff.-like" = "#fcbba1", 
                             "Prolif. stem-like" = "#a50f15")) +
  plot_grid

testPlot(gg_clust_annot)

### Create extra plot enumerating proportions.
pdf("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/10X/tumor-scgp-cellstates.pdf", width = 6, height = 4)
 testPlot(gg_clust_annot) +
  geom_text(aes(label = round(freq, 2)), size = 2, hjust = 0.5, vjust = 3, position ="stack") 
dev.off()



##################################################
## Final combined plots
##################################################

## Legends
gleg1 = g_legend(gg_clust_annot) 
gleg2 = g_legend(gg_cell_states) 


## Plots
g1 = ggplotGrob(gg_clust_annot + plot_grid + plot_theme + null_legend + null_x + top_margin)  %>% gtable_frame()
g2 = ggplotGrob(gg_cell_states + plot_grid + plot_theme + null_legend + null_facet + bottom_margin)  %>% gtable_frame()

g = gtable_rbind(g1, g2)
gleg = gtable_rbind(gleg1, gleg2)

## Adjust relative height of panels
panels = g$layout$t[grep("panel", g$layout$name)]
g$heights[panels] <- unit(c(0.5,0.5), "null")

plot(g)
plot(gleg)

pdf(file = "github/results/Fig3/SuppFig6d-cell-states-compare.pdf", width = 6, height = 5, useDingbats = FALSE, bg="transparent")
grid.newpage()
grid.draw(g)
dev.off()

pdf(file = "github/results/Fig3/SuppFig6d-cell-states-compare-legends.pdf", width = 6, height = 5, useDingbats = FALSE, bg="transparent")
grid.newpage()
grid.draw(gleg)
dev.off()

### END ###