library(tidyverse)
library(readxl)
library(lubridate)

setwd("~/Dropbox/grad/research/g3Diel")

trop_diel <- read_csv("G3_diel_trophicModePredictions_updatedMarferret_marmicroDb2023_noOutliers_fall2023")
meta_diel <- read_csv("g3Diel_surface_tpm_updatedMarferret_marmicroDb2023_sampleTaxa_noOutliers_fixedTPM_fall2023.csv")
prob_diel <- read_csv("G3_diel_trophicModePredictionsProbabilities_updatedMarferret_marmicroDb2023_noOutliers_fall2023")

merged <- cbind(meta_diel, trop_diel, prob_diel)

#Mix, Het, Phot
#1, 0, 2,
merged <- merged %>% mutate(xg_pred = ifelse(xg_pred == 1, "Mix", ifelse(xg_pred == 0, "Het", "Phot")))

colnames(merged)[2]
colnames(merged)[2] <- "taxa"

merged <- merged %>% mutate(xg_pred = str_c(xg_pred, "otrophic"))
merged <- merged %>% mutate(xg_pred = str_replace(xg_pred, "Hetotrophic", "Heterotrophic"))

#downloaded from https://docs.google.com/spreadsheets/d/1lKwvlNjUV8dXHKO5HkMNa2m70QQYydrTf-Fb6PRQ9oE/edit#gid=0
sample <- read_csv("../metaT sample log - metaT.csv")

sample <- sample %>% filter(Cruise == "Gradients3" | Cruise == "Gradients 3")

sample <- sample %>% filter(Type == "polyA")

colnames(merged)[1] <- "sample"

merged <- merged %>% mutate(Sample.ID = str_extract(sample, "UW[0-9]{1,}"))
merged <- merged %>% mutate(Sample.ID = ifelse(is.na(Sample.ID), sample, Sample.ID))
merged <- merged %>% mutate(Sample.ID = str_replace(Sample.ID, "G3PA.diel.", ""))
merged <- merged %>% select(-sample)
merged %>% filter(is.na(Sample.ID))

sample <- sample %>% mutate(Sample.ID = ifelse(Exp == "Diel", str_c("S",Station,"C",as.character(Cast), ".", Replicate), Sample.ID))

sample <- sample %>% select(Cruise:Exp)

merged %>% anti_join(sample, by = c("Sample.ID"))

nrow(merged)
merged <- merged %>% left_join(sample %>% select(-c(Date.Time, Replicate)) %>% distinct(), by = c("Sample.ID"))
nrow(merged)

#exclude nonspecies and nonprotists
merged <- merged %>% filter(!taxa %in% c("cellular organisms", "Micromonas <green algae>", "Hydra vulgaris", 
                 "Paracyclopina nana", "Oikopleura dioica", "Calanus finmarchicus", 
                 "Eucyclops serrulatus", "unclassified Acanthoeca", 
                 "unclassified Micromonas", "unclassified Phaeocystis", 
                 "unclassified Chrysochromulina", "Homarus americanus", 
                 "Neocalanus flemingeri") & 
                   str_detect(taxa, " "))

totalSample <- merged %>% group_by(taxa, Date) %>% summarize(total = n()) 

exclude <- merged %>% group_by(taxa, xg_pred, Date) %>% summarize(numSamples = n()) %>% 
  left_join(totalSample, by = c("taxa", "Date")) %>% mutate(prop = numSamples/total*100) %>% 
  filter(xg_pred %in% c("Heterotrophic", "Phototrophic")) %>% ungroup() %>% 
  filter(prop > 25) %>%
  group_by(taxa, Date) %>% distinct(xg_pred) %>% summarize(n = n()) %>% 
  filter(n > 1) %>% ungroup() %>% distinct(taxa, Date)

nrow(exclude)

merged <- merged %>% anti_join(exclude, by = c("taxa", "Date"))

merged$hourMin <- format(as.POSIXct(merged$Time, format="%R", tz="UTC"), format = "%H:%M:%S")
merged$hour <- round_date(as.POSIXct(as.character(merged$hourMin), format="%R", tz="UTC"), unit="hours")

merged <- merged %>% mutate(justHour = str_replace(hour, ".* ", ""))
merged <- merged %>% mutate(justHour = str_replace_all(justHour, ":00:00", ""))
merged <- merged %>% mutate(justHour = as.numeric(justHour))

merged %>% distinct(Date)
merged %>% filter(Date == "4/17/19") %>% head()

merged <- merged %>% mutate(justHour = ifelse(Date == "4/17/19", justHour + 24, justHour))
merged <- merged %>% mutate(justHour = ifelse(Date == "4/18/19", justHour + 24 + 24, justHour))


merged %>% 
  group_by(justHour, taxa, xg_pred) %>%
  summarize(n = n()) %>% ungroup() %>% arrange(desc(n))

#for shading purposes the night of G3-diel goes from 7PM-6AM, for Diel1-2015 the night is 6PM-6AM
merged %>% 
  group_by(justHour, taxa, xg_pred) %>%
  summarize(n = n()) %>%
  mutate(taxa = factor(taxa, levels = c("Bathycoccus prasinos", "Triparma pacifica",
                                        "Oxytricha trifallax", 
                                        "Karlodinium veneficum"))) %>%
  mutate(xg_pred = str_replace(xg_pred, "ic", "y")) %>%
  ggplot(aes(x = justHour, y = xg_pred)) + 
  geom_point(aes(size = n, color = xg_pred)) +  
  geom_point(aes(size = n),pch=21, color = "black") + 
  annotate("rect", xmin = 2.9, xmax = 6, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4) +
  annotate("rect", xmin = 19, xmax = 30, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4) +
  annotate("rect", xmin = 43, xmax = 54, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4) +
  facet_wrap(~taxa) +
  theme_bw() + 
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.y = element_text(size = 14, color = 'black'), strip.text.x = element_text(size = 20, color = 'black'), 
        axis.text.y = element_text(size = 18, color = 'black'), axis.title = element_text(size = 26, color = 'black'), 
        axis.text.x = element_text(size = 18, color = 'black')) +
  theme(strip.background =element_rect(fill="white")) + 
  theme(legend.title=element_text(size=22), legend.text=element_text(size=20)) +
  labs(x = "Time", y = "", color = "", size = "Number of predictions") + 
  scale_size_continuous(limits = c(1,3), breaks = c(1,2,3)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  guides(size = guide_legend(order = 1))+
  scale_x_continuous(limits = c(2.9, 60), breaks = seq(3, 60, by = 3)) +
  scale_color_manual(values = c("Phototrophy" = "deepskyblue2", "Heterotrophy" = "red", "Mixotrophy" = "black")) +
  guides(color = "none")

ggsave("g3_allOrganismsTrophicPredictions_dotPlot_updatedMarferret_marmicroDb2023_notGroupedBySize_diel_noOutliers.png", height = 5.5, width = 18)

#write cleaned up trophic mode predictions to csv file
merged %>% write_csv("G3_trophicPredictions_cleanedUp_updatedMarferret_marmicroDb2023_diel_noOutliers.csv")

