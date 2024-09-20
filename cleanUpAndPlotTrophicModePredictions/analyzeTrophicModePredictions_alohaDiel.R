library(tidyverse)
library(readxl)

setwd("~/Dropbox/grad/research/alohaDiel")

trop_diel <- read_csv("G1_diel_trophicModePredictions_updatedMarferret_marmicroDb2023_noOutliers_fall2023")
meta_diel <- read_csv("alohaDiel_surface_tpm_updatedMarferret_marmicroDb2023_sampleTaxa_noOutliers_fixedTPM_fall2023.csv")
prob_diel <- read_csv("G1_diel_trophicModePredictionsProbabilities_updatedMarferret_marmicroDb2023_noOutliers_fall2023")

merged <- cbind(meta_diel, trop_diel, prob_diel)

#Mix, Het, Phot
#1, 0, 2,
merged <- merged %>% mutate(xg_pred = ifelse(xg_pred == 1, "Mix", ifelse(xg_pred == 0, "Het", "Phot")))

colnames(merged)[2]
colnames(merged)[2] <- "taxa"

merged <- merged %>% mutate(xg_pred = str_c(xg_pred, "otrophic"))
merged <- merged %>% mutate(xg_pred = str_replace(xg_pred, "Hetotrophic", "Heterotrophic"))
merged %>% distinct(xg_pred)

merged <- merged %>% filter(str_detect(taxa, " "))

merged <- merged %>% filter(!(taxa %in% c("cellular organisms", "Micromonas <green algae>")))
merged %>% distinct(taxa)

colnames(merged)[1] <- "sample"
merged <- merged %>% mutate(time = str_extract(sample, "[0-9]{1,}$"))

merged <- merged %>% mutate(time = as.numeric(time))

merged %>% filter(is.na(time))

meta <- read_excel("SCOPE_cruise_sample_log.xlsx")

colnames(meta) <- meta[5,]
meta <- meta[-c(1:5),]

colnames(meta)[20:22] <- c("date", "cast", "time")
meta <- meta %>% select(20:22)

meta <- meta %>% filter(!is.na(date))

merged <- merged %>% mutate(stationCast = str_extract(sample, "[A-Z,0-9]*"))
merged %>% filter(is.na(stationCast))

merged <- merged %>% mutate(stationCast = ifelse(stationCast == "S06C1", "S6C1", stationCast))
merged <- merged %>% mutate(stationCast = ifelse(stationCast == "S07C1", "S7C1", stationCast))
merged <- merged %>% mutate(stationCast = ifelse(stationCast == "S08C1", "S8C1", stationCast))

meta <- meta %>% mutate(stationCast = str_replace_all(cast, " ", "_"))
meta <- meta %>% mutate(stationCast = str_replace(stationCast, ":", ""))
meta <- meta %>% mutate(stationCast = str_replace(stationCast, "S6", "S06"))
meta <- meta %>% mutate(stationCast = str_replace(stationCast, "S7", "S07"))
meta <- meta %>% mutate(stationCast = str_replace(stationCast, "S8", "S08"))

merged %>% anti_join(meta, by = c("sample" = "stationCast")) %>% distinct(sample)

nrow(merged)
merged <- merged %>% left_join(meta %>% select(-cast), by = c("sample" = "stationCast"))
nrow(merged)

merged %>% filter(is.na(date))

merged %>% distinct(date)

merged %>% filter(date == as.Date("2015-07-27"))

#fix times to be continuous
merged <- merged %>% mutate(time.x = ifelse(date == as.Date("2015-07-27"), time.x+2400, time.x))
merged <- merged %>% mutate(time.x = ifelse(date == as.Date("2015-07-28"), time.x+2400+2400, time.x))
merged <- merged %>% mutate(time.x = ifelse(date == as.Date("2015-07-29"), time.x+2400+2400+2400, time.x))
merged <- merged %>% mutate(time.x = ifelse(date == as.Date("2015-07-30"), time.x+2400+2400+2400+2400, time.x))

#exclude non species and non protists
merged <- merged %>% filter(!taxa %in% c("Calanus finmarchicus", "Eucyclops serrulatus", "Hydra vulgaris", 
               "Oikopleura dioica", "Paracyclopina nana", 
               "unclassified Acanthoeca", "unclassified Micromonas", 
               "unclassified Phaeocystis", "unclassified Chrysochromulina", 
               "Homarus americanus", "Neocalanus flemingeri"))

totalSample <- merged %>% group_by(taxa, date) %>% summarize(total = n()) 

exclude <- merged %>% group_by(taxa, xg_pred, date) %>% summarize(numSamples = n()) %>% 
  left_join(totalSample, by = c("taxa", "date")) %>% mutate(prop = numSamples/total*100) %>% 
  filter(xg_pred %in% c("Heterotrophic", "Phototrophic")) %>% ungroup() %>% 
  filter(prop > 25) %>%
  group_by(taxa, date) %>% distinct(xg_pred) %>% summarize(n = n()) %>% 
  filter(n > 1) %>% ungroup() %>% distinct(taxa, date)

merged %>% semi_join(exclude, by = c("taxa", "date")) %>% nrow()

merged %>% semi_join(exclude, by = c("taxa", "date")) %>% 
  write_csv("excluded.csv")

merged <- merged %>% anti_join(exclude, by = c("taxa", "date"))

#write cleaned up file to csv
merged %>% write_csv("G1_trophicPredictions_cleanedUp_updatedMarferret_marmicroDb2023_diel_noOutliers.csv")

#fix time format
merged <- merged %>% mutate(fixTime = as.character(time.x))
merged %>% distinct(fixTime)
merged <- merged %>% mutate(fixTime = str_replace(fixTime, "00$", ":00"))
merged %>% distinct(fixTime)

order <- merged %>% distinct(fixTime)
order <- as.list(order)

merged$fixTime <- factor(merged$fixTime, levels = order$fixTime)

merged %>% 
  group_by(fixTime, taxa, xg_pred) %>%
  summarize(n = n(), avgProb = mean(probability)) %>% ungroup() %>% 
  arrange(n)

merged <- merged %>% mutate(time.x = time.x/100)

merged$xg_pred <- factor(merged$xg_pred, levels = c("Heterotrophic", "Mixotrophic", "Phototrophic"))

#for shading purposes the night of G3-diel goes from 7PM-6AM, for Diel1-2015 the night is 6PM-6AM
merged %>% 
  group_by(time.x, taxa, xg_pred) %>%
  summarize(n = n()) %>%
  arrange(time.x) %>%
  mutate(taxa = str_replace(taxa, "Chrysochromulina sp. ", "Chrysochromulina sp.\n")) %>%
  mutate(taxa = str_replace(taxa, "Gymnodinium catenatum ", "Gymnodinium catenatum\n")) %>%
  ggplot(aes(x = time.x, y = xg_pred)) + 
  annotate("rect", xmin = 5.9, xmax = 6, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4) +
  annotate("rect", xmin = 18, xmax = 30, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4) +
  annotate("rect", xmin = 42, xmax = 54, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4) +
  annotate("rect", xmin = 66, xmax = 78, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4) +
  annotate("rect", xmin = 90, xmax = 98, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4) +
  geom_point(aes(size = n, color = xg_pred)) +  
  geom_point(aes(size = n),pch=21, color = "black") + 
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
  guides(size = guide_legend(order = 1)) +
  scale_x_continuous(limits = c(5.9, 98), breaks = seq(6, 98, by = 4)) +
  scale_color_manual(values = c("Phototrophic" = "deepskyblue2", "Heterotrophic" = "red", "Mixotrophic" = "black")) +
  guides(color = "none")

ggsave("g1_allOrganismsTrophicPredictions_dotPlot_updatedMarferret_marmicroDb2023_notGroupedBySize_diel_noOutliers.png", height = 5, width = 22)
