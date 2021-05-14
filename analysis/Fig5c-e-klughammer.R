##################################
# Investigate relationships between promoter DNAme and other variables (Klughammer et al. Nature Medicine).
# Updated: 2021.05.12
# Author: Kevin J.
###################################

# Working directory for this analysis.
mybasedir = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/"
setwd(mybasedir)

###################################
# Necessary packages:
library(tidyverse)
library(DBI)
library(ggpubr)
library(openxlsx)
library(survminer)
library(survival)
library(EnvStats)
###################################
## ggplot theme
plot_theme    <- theme_bw(base_size = 12) + theme(axis.title = element_text(size = 12),
                                                  axis.text = element_text(size=12),
                                                  panel.background = element_rect(fill = "transparent"),
                                                  axis.line = element_blank(),
                                                  strip.background =element_rect(fill="white"),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank())

## Retrieve the chromosome arm-data for aneuploidy value estimates from database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB3")
chr_arms = dbReadTable(con,  Id(schema="ref", table="chr_arms"))

## Load the metadata available from Klughammer (http://www.medical-epigenomics.org/papers/GBMatch/).
klughammer_dat_all <- read.delim("/Users/johnsk/Documents/Single-Cell-DNAmethylation/public-rrbs/GBMatch-master/GBMatch_sampleAnnotation.tsv", sep="\t", header=T, stringsAsFactors = F)
klughammer_dict <- read.delim("/Users/johnsk/Documents/Single-Cell-DNAmethylation/public-rrbs/GBMatch-master/GBMatch_columnAnnotation.tsv", sep="\t", header=T, stringsAsFactors = F)

# Restrict to tumor tissue by removing non-tumor tissue from epilepsy patients. Add some additional filtering criteria.
klughammer_dat = klughammer_dat_all %>% 
  # A few tumors do not have IDH status and a single sample appears to have an epimutation of 0.
  filter(tissue != "White matter",  mean_pdr > 0, !is.na(IDH))

# "IDH" status separates based on IDHmut status.
klughammer_idh_mut = klughammer_dat %>% 
  filter(IDH == "mut")
klughammer_idh_wt = klughammer_dat %>% 
  filter(IDH == "wt")

# Create a new variable to account for the different data types as they might relate with the idh_codel_subtype project.
klughammer_dat$idh_codel_subtype <- c()
klughammer_dat$idh_codel_subtype[klughammer_dat$IDH=="mut"] <- "IDHmut_noncodel"
klughammer_dat$idh_codel_subtype[klughammer_dat$WHO2016_classification%in%c("Anaplastic Oligodendroglioma", "Oligo II")] <- "IDHmut_codel"
klughammer_dat$idh_codel_subtype[klughammer_dat$IDH=="wt"] <- "IDHwt"

## There are clear differences across the subtypes with IDHmut tumors having a higher promoter level, but likely also greater tumor purity.
kruskal.test(klughammer_dat$mean_pdr~klughammer_dat$idh_codel_subtype)

# Overall, what are some factors that are related with DNAme disorder? Promoter DNAme disorder is nearly normally distributed.
epi_mut_model <- glm(mean_pdr ~ Age + idh_codel_subtype + material + cohort, family = gaussian, data = klughammer_dat)
summary(epi_mut_model) # There's a slight association between PDR and 1) IDHmut, 2) FFPE, and 3) cohort.


##################################################
### Somatic Copy Number Alteration (SCNA) burden
##################################################
## Investigate the relationship between SCNAs and DNAme disorder,
## First, examine general SCNA burden.
klughammer_cnv = klughammer_dat %>% 
  select(patID, surgery, id, Age, MIB, idh_codel_subtype, mean_pdr, ends_with("_amplification"), ends_with("_deletion")) %>% 
  # Remove sex chromosomes.
    select(-c(X_p_amplification,X_q_amplification, Y_p_amplification, Y_q_amplification,
            X_p_deletion, X_q_deletion, Y_p_deletion, Y_q_deletion))

## Calculate the total bases affected by copy number alteration.
klughammer_fga = klughammer_cnv %>% 
  ## Retrieve the amplification estimates for 22 autosomes.
  mutate(amp_bp = apply(klughammer_cnv[, 8:49], 1, sum, na.rm=TRUE),
         ## Retrieve the deletion estimates for 22 autosomes.
         del_bp = apply(klughammer_cnv[, 50:91], 1, sum, na.rm=TRUE)) %>% 
  select(patID, surgery, id, Age, MIB, idh_codel_subtype, mean_pdr, amp_bp, del_bp) 

# What's the max number of bases altered for each chromosome. We can consider this the theoretical maximum MEASURED.
chrom_amp_max <- apply(klughammer_cnv[, 8:49], 2, max, na.rm=TRUE)
chrom_del_max <- apply(klughammer_cnv[, 50:91], 2, max, na.rm=TRUE)

## There were 46 samples that had no amplifications or deletions. 
## This would be VERY unusual in glioma so filtering out as these are likely low purity or quality.
sum(klughammer_fga$del_bp==0 & klughammer_fga$amp_bp==0)

## Remove these samples from analysis, as they likely represent contaminating normal or low quality segmentations.
klughammer_fga_filtered = klughammer_fga %>% 
  # Calculate the fraction of the genome with amplifications and deletions.
  mutate(amp_fga = amp_bp/sum(chrom_amp_max),
         del_fga = del_bp/sum(chrom_del_max),
         # To derive total fraction of the genome measured with alterations add both features together.
         total_fga = amp_fga+del_fga) %>% 
  filter(total_fga!=0, surgery%in%c(1, 2)) %>% 
  mutate(idh_status = ifelse(idh_codel_subtype=="IDHwt", "IDHwt", "IDHmut"),
         timepoint = ifelse(surgery=="1", "Initial", "Recurrence"))
  
## Restrict to IDHwt samples that have the highest sample size needed for bulk cohort validation analyses.
klughammer_fga_filtered_wt = klughammer_fga_filtered %>% 
  filter(idh_status=="IDHwt")

# Perform a statistical test to determine the association between genetic (SCNA) and epigenetic (DNAme disorder) instability.
# Add back some additional covariates. SCNA seems generally positively associated with promoter DNAme disorder.
disorder_fga_model <- glm(mean_pdr~ Age  + total_fga, family = gaussian, data = klughammer_fga_filtered_wt)
summary(disorder_fga_model) 
disorder_amp_model <- glm(mean_pdr~ Age + as.factor(surgery) + MIB  + amp_fga, family = gaussian, data = klughammer_fga_filtered_wt)
summary(disorder_amp_model)
disorder_del_model <- glm(mean_pdr~ Age + as.factor(surgery) + MIB  + del_fga, family = gaussian, data = klughammer_fga_filtered_wt)
summary(disorder_del_model)
disorder_fga_surgery_model <- glm(mean_pdr~ Age + as.factor(surgery) + MIB  + total_fga, family = gaussian, data = klughammer_fga_filtered_wt)
summary(disorder_fga_surgery_model)

## How many samples? (384 across all conditions)
dim(klughammer_fga_filtered_wt)

## Create a bubble chart for each of these variables. Color = direction of relationship, Height -log10(p-value).
variables <- c("Age", "Timepoint", "Proliferation\n(MIB staining)", "SCNA")
## P-values from disorder_fga_surgery_model results 
pvalues <-  c(summary(disorder_fga_surgery_model)$coefficients[2,4], 
              summary(disorder_fga_surgery_model)$coefficients[3,4], 
              summary(disorder_fga_surgery_model)$coefficients[4,4], 
              summary(disorder_fga_surgery_model)$coefficients[5,4])
sign <- c("-", "+", "-", "+")
disorder_assoc <- data.frame(pvalues, sign, variables)
disorder_assoc$variables <- factor(disorder_assoc$variables, levels = c("Age", "Timepoint", "Proliferation\n(MIB staining)", "SCNA"))

## Can demonstrate that copy number alterations are associated with DNAme disorder.
pdf(file = "results/methylation/public/klughammer-disorder-vars.pdf", height = 5, width = 8, bg = "transparent", useDingbats = FALSE)
ggplot(disorder_assoc, aes(x=variables,  y=-log10(pvalues), size = -log10(pvalues), color = sign)) + 
  geom_point() + 
  geom_hline(yintercept = -log10(0.05), alpha=0.8, linetype=2) +
  labs(x="", y = "-log10(P-value)", size="-log10(P-value)", color = "Direction of assoc.") +
  plot_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


## It appears that the relationship is stronger for the initial IDHwt samples. This could be a bias of sample size.
pdf(file = "github/results/Fig5c-disorder-SCNA.pdf", height = 5, width = 8, bg = "transparent", useDingbats = FALSE)
ggplot(klughammer_fga_filtered_wt, aes(x=total_fga,  y=mean_pdr, color=idh_status)) + 
  geom_point() + 
  ylab("Promoter DNAme disorder") + xlab("SCNA burden") +
  plot_theme + 
  geom_smooth(method = "lm") +   
  xlim(0, 1) +
  scale_color_manual(name='Glioma subtype', values=c('IDHmut'='#AF8DC3', "IDHwt"="#7FBF7B")) +
  facet_grid(. ~ timepoint, scales = "free_x", space = "free") +
  stat_cor(method = "spearman", label.x = 0.4) +
  guides(color=FALSE)
dev.off()

##############################################################
### Analyze longitudinal changes in SCNA and DNAme disorder. 
### How do they correlated?
##############################################################
## Create a longitudinal data frame for fraction of the genome with copy number alterations.
## Determine the longitudinal differences in FGA.
long_fga = klughammer_fga_filtered %>% 
  group_by(patID, surgery) %>%
  summarise(avg_total_fga = mean(total_fga),
            avg_mean_pdr = mean(mean_pdr)) %>% 
  select(-avg_mean_pdr, -surgery) %>% 
  mutate(grouped_id = row_number()) %>% 
  spread(grouped_id, avg_total_fga) %>% 
  filter(!is.na(`2`)) %>% 
  mutate(aneuploidy_diff = `2`-`1`) %>% 
  select(patID, fga_1 = `1`, fga_2 = `2`, aneuploidy_diff)

# Determine the longitudinal differences in epimutation rate.
long_epi_mut = klughammer_fga_filtered %>% 
  group_by(patID, surgery) %>%
  summarise(avg_total_fga = mean(total_fga),
            avg_mean_pdr = mean(mean_pdr)) %>% 
  select(-avg_total_fga, -surgery) %>% 
  mutate(grouped_id = row_number()) %>% 
  spread(grouped_id, avg_mean_pdr) %>% 
  filter(!is.na(`2`)) %>% 
  mutate(disorder_diff = `2`-`1`) %>% 
  select(patID, epi_mut_1 = `1`, epi_mut_2 = `2`, disorder_diff)

## Be able to stitch together the metadata.
# meta_data for all filtered samples.
long_meta_data = klughammer_fga_filtered %>% 
  select(patID, idh_codel_subtype) %>% 
  distinct()

## There should be 110 distinct ind. for whose genomic data passed QC.
long_instability = long_fga %>%
  inner_join(long_epi_mut, by="patID") %>% 
  inner_join(long_meta_data, by="patID")
n_distinct(long_instability$patID)

# Separate by subtype. IDHmut does not have much information.
long_instability_idh_wt = long_instability %>% 
  filter(idh_codel_subtype =="IDHwt")


## Restrict only to the IDHwt tumors to aid in clarity.
pdf(file = "results/methylation/public/klughammer-longitudinal-SCNA-disorder.pdf", height = 5, width = 8, bg = "transparent", useDingbats = FALSE)
ggplot(long_instability_idh_wt, aes(x=aneuploidy_diff,  y=disorder_diff, color=idh_codel_subtype)) + 
  geom_point() + ylab("Delta DNAme disorder\n (Recurrence - Initial)") + xlab("Delta SCNA burden\n(Recurrence - Initial)") +
  theme_bw() + geom_smooth(method = "lm") +   
  ylim(-0.1, 0.1) + xlim(-0.5, 0.5) +
  guides(color=FALSE) +
  scale_color_manual(name='Glioma subtype', values=c("IDHwt"="#7FBF7B"))  + 
  plot_theme +
  stat_cor(method = "spearman", label.x = -0.5, label.y = 0.08)
dev.off()

######################
### Survival
######################
## Retrieve survival information for these samples.
klughammer_clin = klughammer_dat %>% 
  select(patID, Age, IDH, timeToFirstProg, timeToSecSurg, VitalStatus, DateOfDeath_LastFollow.up, 
         Follow.up_years, treatment, Sex) %>% 
  mutate(patient_vital = ifelse(VitalStatus=="alive", 0, 1),
         patient_vital_surgical_int = 1) 
## Only retain the first instance for a particular patient.
klughammer_clin_filt = klughammer_clin[match(unique(klughammer_clin$patID), klughammer_clin$patID), ]

## Create a set of data.frames on which to perform survival analyses.
## It seems that `pat123` and `pat124` were post-treatment at initial disease.
long_instability_clin_wt = long_instability %>% 
  left_join(klughammer_clin_filt, by="patID") %>% 
  ungroup() %>% 
  ## Restrict to IDHwt only.
  filter(IDH=="wt") %>% 
  mutate(disorder_diff_group = ifelse(disorder_diff>median(disorder_diff), "1", "0"),
         ## Convert days to months.
         timeToSecSurg.mo = timeToSecSurg/30.4167,
         timeToFirstProg.mo = timeToFirstProg/30.4167,
         primary_disorder_group = ifelse(epi_mut_1>median(epi_mut_1), "1", "0"),
         recurrence_disorder_group = ifelse(epi_mut_2>median(epi_mut_2), "1", "0"))

## Perform Cox-proportional hazards models with epimutation rate + Age.
primary_pdr_cox <- coxph(Surv(Follow.up_years, patient_vital) ~ Age + Sex + primary_disorder_group, data = long_instability_clin_wt)
summary(primary_pdr_cox)
recurrent_pdr_cox <- coxph(Surv(Follow.up_years, patient_vital) ~ Age + Sex + recurrence_disorder_group, data = long_instability_clin_wt)
summary(recurrent_pdr_cox)

## Delta DNAme disorder with overall survival.
## Continuous variable.
diff_pdr_cox <- coxph(Surv(Follow.up_years, patient_vital) ~ Age + Sex  + disorder_diff, data = long_instability_clin_wt)
summary(diff_pdr_cox)
## Treated as a discrete variable.
diff_pdr_cox <- coxph(Surv(Follow.up_years, patient_vital) ~ Age + Sex + disorder_diff_group, data = long_instability_clin_wt)
summary(diff_pdr_cox)

## Add in the delta SCNA value to account for potential shifts in tumor purity (i.e., loss of CNVs may represent decrease in purity).
diff_pdr_cox <- coxph(Surv(Follow.up_years, patient_vital) ~ Age + Sex + disorder_diff_group + aneuploidy_diff, data = long_instability_clin_wt)
summary(diff_pdr_cox)

## Generate visualizations based on groupings:
fit_wt_os <- survfit(Surv(Follow.up_years, patient_vital) ~ disorder_diff_group,
                      data = long_instability_clin_wt)
pdf(file = "results/methylation/public/klughammer-disorder-diff-overall-survival.pdf", height = 5, width = 8, bg = "transparent", useDingbats = FALSE)
ggsurvplot(fit_wt_os, data = long_instability_clin_wt, risk.table = FALSE, pval= TRUE, pval.coord = c(6, 0.75),
           palette = c("royalblue4", "tomato3"),
           ylab = "Overall survival\n probability", xlab = "Time (years)")
dev.off()

## Examining only the primary tumor samples.
fit_wt_os_primary <- survfit(Surv(Follow.up_years, patient_vital) ~ primary_disorder_group,
                     data = long_instability_clin_wt)
pdf(file = "results/methylation/public/klughammer-disorder-group-primary-overall-survival.pdf", height = 5, width = 8, bg = "transparent", useDingbats = FALSE)
ggsurvplot(fit_wt_os_primary, data = long_instability_clin_wt, risk.table = FALSE, pval= TRUE, pval.coord = c(6, 0.75),
           palette = c("royalblue4", "tomato3"),
           ylab = "Overall survival\n probability", xlab = "Time (years)")
dev.off()

## Examining only the recurrent tumor samples.
fit_wt_os_recurence <- survfit(Surv(Follow.up_years, patient_vital) ~ recurrence_disorder_group,
                             data = long_instability_clin_wt)
pdf(file = "results/methylation/public/klughammer-disorder-group-recurrence-overall-survival.pdf", height = 5, width = 8, bg = "transparent", useDingbats = FALSE)
ggsurvplot(fit_wt_os_recurence, data = long_instability_clin_wt, risk.table = FALSE, pval= TRUE, pval.coord = c(6, 0.75),
           palette = c("royalblue4", "tomato3"),
           ylab = "Overall survival\n probability", xlab = "Time (years)")
dev.off()


##############################################
### Assess the "time to second surgery"    ###
##############################################
## Make sure to use the separate patient vital variable for surgical interval.
table(long_instability_clin_wt$patient_vital_surgical_int)
pdr_surg_cox_con <- coxph(Surv(timeToSecSurg, patient_vital_surgical_int) ~ Age + Sex + disorder_diff , data = long_instability_clin_wt)
summary(pdr_surg_cox_con)

## Second surgery analysis treated as a group. Doesn't change much.
pdr_surg_cox_bi <- coxph(Surv(timeToSecSurg, patient_vital_surgical_int) ~ Age + Sex + disorder_diff_group, data = long_instability_clin_wt)
summary(pdr_surg_cox_bi)

## The relationship gets stronger when accounting for SCNA differences.
pdr_surg_cox_bi <- coxph(Surv(timeToSecSurg, patient_vital_surgical_int) ~ Age + Sex + disorder_diff_group + aneuploidy_diff, data = long_instability_clin_wt)
summary(pdr_surg_cox_bi)

## What percentage of these two tissues are the same for the two different progression variables.
sum(long_instability_clin_wt$timeToSecSurg==long_instability_clin_wt$timeToFirstProg)/length(long_instability_clin_wt$timeToFirstProg)

## Assess the time to first progression (assuming MRI). More subjective measure.
pdr_prog_cox <- coxph(Surv(timeToFirstProg, patient_vital_surgical_int) ~ Age + Sex + disorder_diff, data = long_instability_clin_wt)
summary(pdr_prog_cox)
pdr_prog_cox <- coxph(Surv(timeToFirstProg, patient_vital_surgical_int) ~ Age + Sex + disorder_diff_group, data = long_instability_clin_wt)
summary(pdr_prog_cox)

## Again, the survival relationship gets stronger when including a purity proxy.
pdr_prog_cox <- coxph(Surv(timeToFirstProg, patient_vital_surgical_int) ~ Age + Sex + disorder_diff + aneuploidy_diff, data = long_instability_clin_wt)
summary(pdr_prog_cox)
pdr_prog_cox <- coxph(Surv(timeToFirstProg, patient_vital_surgical_int) ~ Age + Sex + disorder_diff_group + aneuploidy_diff, data = long_instability_clin_wt)
summary(pdr_prog_cox)

## Generate visualizations based on groupings:
fit_wt_pdr <- survfit(Surv(timeToSecSurg.mo, patient_vital_surgical_int) ~ disorder_diff_group,
                             data = long_instability_clin_wt)
pdf(file = "results/methylation/public/klughammer-disorder-diff-secondsurg-revised.pdf", height = 5, width = 8, bg = "transparent", useDingbats = FALSE)
ggsurvplot(fit_wt_pdr, data = long_instability_clin_wt, risk.table = FALSE, pval= TRUE, pval.coord = c(30, 0.75),
           palette = c("royalblue4", "tomato3"),
           ylab = "Time to Second Surgery\n probability", xlab = "Time (months)")
dev.off()

## Visualize the time to first progression (MRI) in months. 
fit_wt_pdr <- survfit(Surv(timeToFirstProg.mo, patient_vital_surgical_int) ~ disorder_diff_group,
                      data = long_instability_clin_wt)
pdf(file = "results/methylation/public/klughammer-epimut-diff-pfs.pdf", height = 6, width = 8, bg = "transparent", useDingbats = FALSE)
ggsurvplot(fit_wt_pdr, data = long_instability_clin_wt, risk.table = FALSE, pval= TRUE, pval.coord = c(30, 0.75),
           palette = c("tomato3", "royalblue4"),
           ylab = "Progression Free Survival\n probability", xlab = "Time (months)")
dev.off()

#### END #####

