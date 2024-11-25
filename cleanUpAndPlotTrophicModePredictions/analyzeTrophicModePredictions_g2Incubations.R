library(tidyverse)
library(readxl)

setwd("~/Dropbox/grad/research/g2Incubation/")

trop <- read_csv("g2Inc_trophicModePredictions_updatedMarferret_marmicroDb2023_noOutliers_fall2023")
meta <- read_csv("g2Inc_tpm_updatedMarferret_marmicroDb2023_sampleTaxa_noOutliers_fixedTPM_fall2023.csv")
prob <- read_csv("g2Inc_trophicModePredictionsProbabilities_updatedMarferret_marmicroDb2023_noOutliers_fall2023")

merged <- cbind(meta, trop, prob)

#Mix, Het, Phot
#1, 0, 2,
merged <- merged %>% mutate(xg_pred = ifelse(xg_pred == 1, "Mix", ifelse(xg_pred == 0, "Het", "Phot")))

colnames(merged)[2]
colnames(merged)[2] <- "taxa"

merged <- merged %>% mutate(xg_pred = str_c(xg_pred, "otrophic"))
merged <- merged %>% mutate(xg_pred = str_replace(xg_pred, "Hetotrophic", "Heterotrophic"))
merged %>% distinct(xg_pred)

colnames(merged)[1]
colnames(merged)[1] <- "sample"

merged <- merged %>% mutate(`Sample name` = str_extract(sample, "MS[0-9]{1,}"))

sample <- read_csv("G2.RR_exp.metadata.csv")

merged %>% anti_join(sample, by = c("Sample name"))

nrow(merged)
merged <- merged %>% left_join(sample, by = c("Sample name"))
nrow(merged)

#exclude nonspecies and nonprotists
merged <- merged %>% filter(!taxa %in% c("cellular organisms", "Micromonas <green algae>", 
                 "Calanus finmarchicus", "Eucyclops serrulatus", 
                 "Hydra vulgaris", "Oikopleura dioica", "Paracyclopina nana", 
                 "unclassified Acanthoeca", "unclassified Micromonas", 
                 "unclassified Phaeocystis", "unclassified Chrysochromulina") & 
                   str_detect(taxa, " "))

totalSample <- merged %>% group_by(taxa, Timepoint, Treatment, Expt) %>% summarize(total = n()) 

exclude <- merged %>% group_by(taxa, Timepoint, Treatment, Expt, xg_pred) %>% summarize(numSamples = n()) %>% 
  left_join(totalSample, by = c("taxa", "Timepoint", "Treatment", "Expt")) %>% mutate(prop = numSamples/total*100) %>% 
  filter(xg_pred %in% c("Heterotrophic", "Phototrophic")) %>% ungroup() %>% 
  filter(prop > 25) %>%
  group_by(taxa, Timepoint, Treatment, Expt) %>% distinct(xg_pred) %>% summarize(n = n()) %>% 
  filter(n > 1) %>% ungroup() %>% distinct(taxa, Timepoint, Treatment, Expt)

merged %>% semi_join(exclude, by = c("taxa", "Timepoint", "Treatment", "Expt")) %>% nrow()

merged %>% semi_join(exclude, by = c("taxa", "Timepoint", "Treatment", "Expt")) %>% 
  write_csv("excluded.csv")

merged <- merged %>% anti_join(exclude, by = c("taxa", "Timepoint", "Treatment", "Expt"))

#write cleaned up trophic mode predictions to csv file
merged %>% write_csv("g2Inc_trophicPredictions_cleanedUp_updatedMarferret_marmicroDb2023_noOutliers.csv")

mixed <- read_csv("../g1G2G3/mixedSpeciesStrict_dielIncubations.csv")

#At the northernmost REXP1 station at 41.42°N, in situ dissolved
#inorganic nitrogen (DIN) was 2 μM and iron (Fe) was 0.3 nM resulting in a molar N:Fe of 6.6 ×
#103 and log10N:Fe of 3.8. At the transition zone REXP2 station at 37.00°N, in situ nutrient
#concentrations were 0.06 μM DIN and 0.51 nM Fe resulting in a log10N:Fe of 2.1 At the
#southernmost REXP3 station at 32.93°N, in situ concentrations were 0.01 μM DIN and 0.22 nM
#Fe, resulting in a log10N:Fe of 1.7 (Fig 3.4 a, b).
merged <- merged %>% mutate(latitude = ifelse(Expt == "REXP1", 41.42, NA))
merged <- merged %>% mutate(latitude = ifelse(Expt == "REXP2", 37.00, latitude))
merged <- merged %>% mutate(latitude = ifelse(Expt == "REXP3", 32.93, latitude))

order <- c("32.93_Ctrl_96", "32.93_LoNP_96", "32.93_HiNP_96", "32.93_NPFe_96", "37_Ctrl_96", "37_Fe_96", "37_NP_96", "37_NPFe_96", "41.42_Ctrl_96", "41.42_LFe_96", "41.42_HFe_96", "41.42_NPFe_96")

merged %>% semi_join(mixed, by = c("taxa")) %>%
  mutate(latTreat = str_c(latitude, "_", Treatment, "_", Timepoint)) %>%
  mutate(latitude = str_c(latitude, " °N")) %>%
  mutate(xg_pred = str_replace(xg_pred, "ic", "y")) %>%
  group_by(taxa, latitude, Timepoint, Treatment, xg_pred, latTreat) %>%
  summarize(n = n(), avgProb = mean(probability)) %>%
  mutate(latTreat = factor(latTreat, levels = order)) %>%
  mutate(taxa = factor(taxa, levels = c("Azadinium spinosum", "Karlodinium veneficum", "Pelagodinium beii", "Prorocentrum minimum"))) %>%
  mutate(xg_pred = factor(xg_pred, levels = rev(c("Heterotrophy", "Mixotrophy", "Phototrophy")))) %>%
  filter(taxa %in% c("Karlodinium veneficum", "Pelagodinium beii")) %>%
  filter(!str_detect(latTreat, "0|in situ")) %>%
  ggplot(aes(x = latTreat, y = n, fill = xg_pred)) + geom_bar(stat = 'identity') +  
  facet_grid(rows = vars(taxa), cols = vars(latitude), scales = "free_x") +
  theme_bw() + 
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.y = element_text(size = 18, color = 'black'), strip.text.x = element_text(size = 20, color = 'black'), 
        axis.text.x = element_text(size = 10, color = 'black'), axis.text.y = element_text(size = 18, color = 'black'), 
        axis.title = element_text(size = 26, color = 'black')) +
  theme(legend.title=element_text(size=16), legend.text=element_text(size=16)) +
  labs(x = "", y = "Number of predictions", fill = "Trophic mode prediction") + 
  scale_fill_manual(values = c("Phototrophy" = "deepskyblue2", "Heterotrophy" = "red", "Mixotrophy" = "black"))+
  guides(color = "none") + 
  scale_y_continuous(limits = c(0,6), breaks = c(0,2,4,6)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("g2IncTrophicPredictions.png", dpi = 600, height = 7.5, width = 17)


