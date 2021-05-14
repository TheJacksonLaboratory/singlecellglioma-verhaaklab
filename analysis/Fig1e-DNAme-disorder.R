##################################
# Generate boxplots for sample-specific global DNAme disorder (annot. with bulk features).
# Updated: 2021.05.12
# Author: Kevin J.
###################################

# Working directory for this analysis in the SCGP project. 
mybasedir = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/"
setwd(mybasedir)

###################################
## Load the essential packages.
library(tidyverse)
library(grid)
library(gridExtra)
library(gtable)
library(egg)
library(openxlsx)
library(EnvStats)
###################################

## Load in the data with the total number of epialleles.
## Load the SCGP subject-level metadata.
full_meta_data = readWorkbook("/Users/johnsk/mnt/verhaak-lab/scgp/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

## Need to extract subtype and IDHmut status.
meta_data = full_meta_data %>% 
  select(case_barcode = subject_id, subtype) %>% 
  mutate(idh_status = ifelse(subtype == "IDHwt", "IDHwt", "IDHmut")) %>% 
  mutate(case_barcode_short = gsub("-", "", substr(case_barcode, 6, 11)))


######### CLINICAL ##########
## Load in clinical data (subtype, age, treatment, hypermutation_status).
# Supply metadata so that 10X filenames and samples can be linked together.
metadata = readWorkbook("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
metadata$`10X_id_short` <- as.character(metadata$`10X_id_short`)
clin_data = metadata %>% 
  mutate(case_barcode = gsub("-", "", substr(subject_id, 6, 11))) %>% 
  mutate(idh_status = ifelse(subtype == "IDHwt", "IDHwt", "IDHmut")) %>% 
  select(case_barcode, idh_codel_subtype = subtype, idh_status, grade = who_grade, timepoint = initial_recurrence,  is_hypermutator = hypermutation)

######## Mutation Freq. ###########
epi_mut_hmapdata = read.table("/Users/johnsk/mnt/verhaak-lab/scgp/results/mutation_vs_epimutation/Samples_WGS_epimutation_and_mutation_rates.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
epi_mut_hmapdata = epi_mut_hmapdata %>% 
  mutate(idh_status = ifelse(subtype == "IDHwt", "IDHwt", "IDHmut")) %>% 
  select(case_barcode = Patient, idh_codel_subtype = subtype, idh_status, timepoint = initial_recurrence, mut_freq = global)

######### Aneuploidy CNVs #########
aneuploidy_hmapdata = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/wgs/scgp_aneuploidy_hmp_data.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
aneuploidy_hmapdata = aneuploidy_hmapdata %>% 
  mutate(idh_status = ifelse(idh_codel_subtype == "IDHwt", "IDHwt", "IDHmut")) 

########### Combine annotation meta data ############
meta_comb = clin_data %>% 
  left_join(epi_mut_hmapdata, by=c("case_barcode", "idh_codel_subtype")) %>% 
  left_join(aneuploidy_hmapdata, by=c("case_barcode", "idh_codel_subtype")) %>% 
  mutate(`1p_19q_codel` = ifelse(idh_codel_subtype=="IDHmut_codel", 1, 0)) %>% 
  select(case_barcode:grade, `1p_19q_codel`,  idh_status, timepoint = timepoint.x,  mut_freq, aneuoploidy_value)
                          

##### Epimutation ##############
epiallele_info <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/Samples-passQC_single_cells_global_and_context-specific_pdr_score_table.txt", sep="\t", header=T, stringsAsFactors = F)

## Create a new variable to indicate "non_tumor" or shortened case barcode.
gg_epiallele = epiallele_info %>% 
  #mutate(case_normal_barcode = ifelse(tumor_cnv==0, "Non-tumor", gsub("-","", substr(sample, 6, 11)))) %>% 
  mutate(case_barcode = gsub("-","", substr(sample, 6, 11))) %>% 
  left_join(meta_data, by=c("sample"="case_barcode")) %>% 
  #  mutate(idh_codel_subtype = ifelse(case_normal_barcode=="Non-tumor", "Non-tumor", subtype)) %>% 
  select(sample_barcode, case_barcode, PDR, idh_status, idh_codel_subtype = subtype)



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
plot_theme    <- theme_bw(base_size = 12) + theme(axis.title = element_text(size = 12),
                                                 axis.text = element_text(size=12),
                                                 axis.line = element_blank(),
                                                 panel.background = element_rect(fill = "transparent"),
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

########################
## Plot clinical
########################
subtype_order <- c("IDHmut", "IDHwt")
case_order <- c("SM004", "SM001", "SM015", "UC917", "SM002", "SM008", "SM006", "SM012", "SM017", "SM018", "SM011")
meta_comb <- meta_comb %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
meta_comb <- meta_comb %>% mutate(idh_status = factor(idh_status, levels = subtype_order))


gg_timepoint <-
  meta_comb %>% 
  gather(key = "type", value = "value", timepoint) %>%
  mutate(type = factor(type,
                       levels = c("timepoint"),
                       labels = c("Timepoint"))) %>%
  ggplot(aes(x=case_barcode)) +
  geom_tile(aes(fill = factor(value), y = type)) +
  scale_fill_manual(values = c("initial" = "#CA2F66", "recurrence" = "#2FB3CA")) +
  labs(y="", fill = "Clinical")

## Relabel some of the variables.
gg_timepoint$data$value <- factor(gg_timepoint$data$value, levels = c("initial", "recurrence"))
testPlot(gg_timepoint)

gg_codel <-
  meta_comb %>% 
  gather(key = "type", value = "value", `1p_19q_codel`) %>%
  mutate(type = factor(type,
                       levels = c("1p_19q_codel"),
                       labels = c("Chr1p/19q codel"))) %>%
  ggplot(aes(x=case_barcode)) +
  geom_tile(aes(fill = factor(value), y = type)) +
  scale_fill_manual(values = c( "0" = "gray90", "1" = "#377eb8")) +
  labs(y="", fill = "Chr1p/19q Codeletion")

## Relabel some of the variables.
gg_codel$data$value <- factor(gg_codel$data$value, levels = c("1", "0"))

testPlot(gg_codel)


###########################
### Mutation frequency
###########################
meta_comb$mut_freq[meta_comb$case_barcode=="SM011"] <- 15
gg_mf = ggplot() +
  geom_tile(data = meta_comb, aes(x = case_barcode, y = 1, fill = mut_freq), color = "black") +
  labs(y = "", fill = "Mutations/Mb") +
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

############################
### Epimutation boxplots
############################
gg_epiallele <- gg_epiallele %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
gg_epiallele <- gg_epiallele %>% mutate(idh_codel_subtype = factor(idh_codel_subtype, levels = subtype_order))


## Color by idh_codel_subtype.
gg_epimut = ggplot(gg_epiallele, aes(x=case_barcode, y = PDR, fill=idh_status)) +
  geom_boxplot(outlier.shape = NA) + 
  plot_theme +
  scale_fill_manual(name='Glioma subtype', values=c("IDHmut"= "#AF8DC3", "IDHwt"="#7FBF7B")) +
  labs(y = "DNAme disorder (PDR)", x="") + 
  guides(alpha = FALSE, color = FALSE) +
  facet_grid(~idh_codel_subtype, scales = "free", space="free_x") +
  ylim(0.25, 0.45) 

testPlot(gg_epimut)


############################
### Epimutation variation per sample
############################
gg_epimut_mad <- gg_epiallele %>% 
  group_by(case_barcode) %>% 
  summarise(epimut_mad = mad(PDR),
            epimut_median = median(PDR))

gg_epimut_mad = gg_epimut_mad %>% 
  inner_join(meta_comb, by="case_barcode")

gg_epimut_var = ggplot() +
  geom_tile(data = gg_epimut_mad, aes(x = case_barcode, y = 1, fill = epimut_mad), color = "black") +
  labs(y = "", fill = "Epimutation variability") +
  scale_fill_distiller(palette = "Greens", direction = 1)
testPlot(gg_epimut_var)

##################################################
## Final combined plots
##################################################
# ## Legends
gleg1 = g_legend(gg_epimut) 
gleg2 = g_legend(gg_codel)
gleg3 = g_legend(gg_timepoint)
gleg4 = g_legend(gg_mf) 
gleg5 = g_legend(gg_cin)


## Plots
g1 = ggplotGrob(gg_epimut + plot_grid + plot_theme + null_legend + null_x + top_margin)  %>% gtable_frame()
g2 = ggplotGrob(gg_codel + plot_grid + plot_theme + null_legend + null_x + null_y  + null_facet + middle_margin) %>% gtable_frame()
g3 = ggplotGrob(gg_timepoint + plot_grid + plot_theme + null_legend + null_x + null_y + null_facet + middle_margin) %>% gtable_frame()
g4 = ggplotGrob(gg_mf + plot_grid + plot_theme + null_legend + null_x + null_y + null_facet + middle_margin) %>% gtable_frame()
g5 = ggplotGrob(gg_cin + plot_grid + plot_theme + null_legend + null_x + null_y + null_facet + bottom_margin)  %>% gtable_frame()

g = gtable_rbind(g1, g2, g3, g4, g5)
gleg = gtable_rbind(gleg1, gleg2, gleg3)

## Adjust relative height of panels
panels = g$layout$t[grep("panel", g$layout$name)]
g$heights[panels] <- unit(c(4, 0.25, 0.25, 0.25, 0.25), "null")

plot(g)
plot(gleg)

pdf(file = "github/results/Fig1/Fig1e-DNAme-disorder.pdf", height = 4, width = 6, useDingbats = FALSE, bg="transparent")
grid.newpage()
grid.draw(g)
dev.off()

pdf(file = "github/results/Fig1/Fig1e-DNAme-disorder-legend.pdf", width = 4, height = 6, useDingbats = FALSE)
grid.newpage()
grid.draw(gleg)
dev.off()


### Calculate some statistics for the epimut
meta_comb_stats = meta_comb %>% 
  inner_join(gg_epimut_mad, by="case_barcode") 

## Test association between median epimutation value and "chromosomal instability"
cor.test(meta_comb_stats$epimut_median, meta_comb_stats$aneuoploidy_value, method="s") 
## Test association between median epimutation value and "Mutations/Mb".
meta_comb_stats$mut_freq[meta_comb_stats$case_barcode=="SM011"] <- 269.8537783
cor.test(meta_comb_stats$epimut_median, log10(meta_comb_stats$mut_freq), method="s") 

## Test association between median epimutation value and MAD epimutation value.
cor.test(meta_comb_stats$epimut_median, log10(meta_comb_stats$epimut_mad), method="s") 

### END ###