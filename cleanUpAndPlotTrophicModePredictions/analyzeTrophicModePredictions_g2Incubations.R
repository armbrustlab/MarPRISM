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


###plot figures of g2 incubation trophic predictions

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


#For REXP1, the in situ log10N:Fe of 3.8 was
#shifted to 2.9 (HiFe), 3.5 (LoFe) or 3.7 (NPFe) through amendments with different ratios of N, P,
#and Fe. In REXP2, the in situ log10N:Fe of 2.1 was shifted to 1.6 (Fe), 3.5 (NPFe) or 4.0 (NP). In
#REXP3, the in situ log10N:Fe of 1.7 was shifted to 3.4 (NP), 3.8 (NPFe), or 4.4 (NP)
merged <- merged %>% mutate(NtoFe = ifelse(Expt == "REXP1" & Treatment == "Ctrl", 3.759058, NA))
merged <- merged %>% mutate(NtoFe = ifelse(Expt == "REXP1" & Treatment == "HFe", 2.886808, NtoFe))
merged <- merged %>% mutate(NtoFe = ifelse(Expt == "REXP1" & Treatment == "LFe", 3.46509, NtoFe))
merged <- merged %>% mutate(NtoFe = ifelse(Expt == "REXP1" & Treatment == "NPFe", 3.707533, NtoFe))

merged <- merged %>% mutate(NtoFe = ifelse(Expt == "REXP2" & Treatment == "Ctrl", 1.144683, NtoFe))
merged <- merged %>% mutate(NtoFe = ifelse(Expt == "REXP2" & Treatment == "Fe", 0.6228152, NtoFe))
merged <- merged %>% mutate(NtoFe = ifelse(Expt == "REXP2" & Treatment == "NPFe", 3.544155, NtoFe))
merged <- merged %>% mutate(NtoFe = ifelse(Expt == "REXP2" & Treatment == "NP", 4.066022, NtoFe))

merged <- merged %>% mutate(NtoFe = ifelse(Expt == "REXP3" & Treatment == "Ctrl", 0.9586073, NtoFe))
merged <- merged %>% mutate(NtoFe = ifelse(Expt == "REXP3" & Treatment == "HiNP", 4.356547, NtoFe))
merged <- merged %>% mutate(NtoFe = ifelse(Expt == "REXP3" & Treatment == "NPFe", 3.841638, NtoFe))
merged <- merged %>% mutate(NtoFe = ifelse(Expt == "REXP3" & Treatment == "LoNP", 3.356548, NtoFe))

merged %>% semi_join(mixed, by = c("taxa")) %>%
  filter(Timepoint != 0) %>%
  mutate(latitude = round(NtoFe)) %>%
  mutate(latitude = str_c(as.character(latitude), " °N")) %>%
  mutate(xg_pred = str_replace(xg_pred, "ic", "y")) %>%
  mutate(taxa = str_replace(taxa, " ", "\n")) %>%
  group_by(taxa, latitude, Timepoint, NtoFe, xg_pred) %>%
  summarize(n = n(), avgProb = mean(probability)) %>%
  ggplot(aes(x = NtoFe, y = n, fill = xg_pred)) + 
  geom_bar(stat = 'identity', width = .1) + 
  facet_grid(rows = vars(taxa), cols = vars(latitude)) +
  theme_bw() + 
  theme(strip.background =element_rect(fill="white"))  +
  theme(strip.text.y = element_text(size = 16, color = 'black'), strip.text.x = element_text(size = 20, color = 'black'), 
        axis.text.y = element_text(size = 18, color = 'black'), axis.title = element_text(size = 26, color = 'black'), 
        axis.text.x = element_text(size = 18, color = 'black')) + 
  theme(legend.title=element_text(size=22), legend.text=element_text(size=20)) +
  labs(x = expression(log[10](paste("NO"[3], "_", "NO"[2], "/Fe")))) +
  labs(y = "", size = "Number of predictions") + 
  scale_fill_manual(values = c("Phototrophy" = "deepskyblue2", "Heterotrophy" = "red", "Mixotrophy" = "black"))+
  guides(color = "none") 

#REXP1
#nitrate + nitrite (N+N) was 1.78 μM,
#phosphate was and 0.42 μM
#iron was 0.31 nM
log10(1.78/(0.31*.001))
3.759058

#REXP2 in situ N+N and phosphate
#nitrate 0.006 μM 
#phosphate 0.13 μM
#iron 0.43 nM
log10(0.006/(0.43*.001))
1.144683

#REXP3
#nitrate 0.002 μM
#and phosphate 0.07 μM
log10(0.002/(0.22*.001))
0.9586073


#Two
#iron amendments were given: a low-iron pulse (LoFe; 0.3 nM Fe), a high-iron pulse (HiFe; 2 nM
#Fe), and a combined N+P+Fe pulse (NPFe; 2 nM Fe + 10 μM NO3
#- + 1 μM PO4
#3- added). 


#For REXP1, the in situ log10N:Fe of 3.8 was
#shifted to 2.9 (HiFe), 3.5 (LoFe) or 3.7 (NPFe) through amendments with different ratios of N, P,
#and Fe.

#rexp1 
#hife
log10(1.78/(0.31*.001+2*.001))
2.886808

#lofe
log10(1.78/(0.31*.001+0.3*.001))
3.46509

#NPFe
log10((1.78+10)/(0.31*.001+2*.001))
3.707533

#In REXP2, the in situ log10N:Fe of 2.1 was shifted to 1.6 (Fe), 3.5 (NPFe) or 4.0 (NP).

#Resource ratio experiment 2 (REXP2) used water collected from Station 11 (latitude 37.00°N,
#ambient water temperature 15.2°C) and were given one iron amendment (Fe; 1 nM Fe), an N+P
#amendment (NP; 0.5 μM P + 5 μM NO3) and an N+P+Fe amendment (NPFe; 1 nM Fe + 5 μM
#NO3- + 0.5 μM PO43- added). 

log10(0.006/(0.43*.001))

#Fe
log10(0.006/(0.43*.001+1*.001))
0.6228152

#NP
log10((0.006+5)/(0.43*.001))
4.066022

#NPFe
log10((0.006+5)/(0.43*.001+1*.001))
3.544155


#REXP3, the in situ log10N:Fe of 1.7 was shifted to 3.4 (NP), 3.8 (NPFe), or 4.4 (NP) (Fig 3.4b).

log10(0.002/(0.22*.001))

#low N+P pulse (LoNP, 0.5 μM NO3- + 0.05 μM PO43-), a high N+P pulse (HiNP, 5 μM NO3- + 0.5 μM PO43-
#) and a combined N+P+Fe pulse (NPFe, 5 μM NO3- + 0.5 μM PO43- 0.5 nM + Fe added). 

#LoNP
log10(0.002+0.5/(0.22*.001))
3.356548

#HiNP
log10(0.002+5/(0.22*.001))
4.356547

#NPFe
log10(0.002+5/(0.22*.001+0.5*.001))
3.841638

merged <- merged %>% mutate(nit = ifelse(Expt == "REXP1" & Treatment == "Ctrl", 0.25042, NA))
merged <- merged %>% mutate(nit = ifelse(Expt == "REXP1" & Treatment == "HFe", 0.25042, nit))
merged <- merged %>% mutate(nit = ifelse(Expt == "REXP1" & Treatment == "LFe", 0.25042, nit))
merged <- merged %>% mutate(nit = ifelse(Expt == "REXP1" & Treatment == "NPFe", 1.071145, nit))

merged <- merged %>% mutate(nit = ifelse(Expt == "REXP2" & Treatment == "Ctrl", -2.221849, nit))
merged <- merged %>% mutate(nit = ifelse(Expt == "REXP2" & Treatment == "Fe", -2.221849, nit))
merged <- merged %>% mutate(nit = ifelse(Expt == "REXP2" & Treatment == "NPFe", 0.6994908, nit))
merged <- merged %>% mutate(nit = ifelse(Expt == "REXP2" & Treatment == "NP", 0.6994908, nit))

merged <- merged %>% mutate(nit = ifelse(Expt == "REXP3" & Treatment == "Ctrl", -2.69897, nit))
merged <- merged %>% mutate(nit = ifelse(Expt == "REXP3" & Treatment == "HiNP", 0.6991437, nit))
merged <- merged %>% mutate(nit = ifelse(Expt == "REXP3" & Treatment == "NPFe", 0.6991437, nit))
merged <- merged %>% mutate(nit = ifelse(Expt == "REXP3" & Treatment == "LoNP", -0.2992963, nit))


#REXP1
#nitrate + nitrite (N+N) was 1.78 μM,
#phosphate was and 0.42 μM
#iron was 0.31 nM
log10(1.78)
0.25042

#REXP2 in situ N+N and phosphate
#nitrate 0.006 μM 
#phosphate 0.13 μM
#iron 0.43 nM
log10(0.006)
-2.221849

#REXP3
#nitrate 0.002 μM
#and phosphate 0.07 μM
log10(0.002)
-2.69897


#For REXP1, the in situ log10N:Fe of 3.8 was
#shifted to 2.9 (HiFe), 3.5 (LoFe) or 3.7 (NPFe) through amendments with different ratios of N, P,
#and Fe.

#rexp1 
#hife
log10(1.78)
0.25042

#lofe
log10(1.78)
0.25042

#NPFe
log10((1.78+10))
1.071145

#In REXP2, the in situ log10N:Fe of 2.1 was shifted to 1.6 (Fe), 3.5 (NPFe) or 4.0 (NP).

#Resource ratio experiment 2 (REXP2) used water collected from Station 11 (latitude 37.00°N,
#ambient water temperature 15.2°C) and were given one iron amendment (Fe; 1 nM Fe), an N+P
#amendment (NP; 0.5 μM P + 5 μM NO3) and an N+P+Fe amendment (NPFe; 1 nM Fe + 5 μM
#NO3- + 0.5 μM PO43- added). 

#Fe
log10(0.006)
-2.221849

#NP
log10((0.006+5))
0.6994908

#NPFe
log10((0.006+5))
0.6994908


#REXP3, the in situ log10N:Fe of 1.7 was shifted to 3.4 (NP), 3.8 (NPFe), or 4.4 (NP) (Fig 3.4b).
log10(0.002/(0.22*.001))

#low N+P pulse (LoNP, 0.5 μM NO3- + 0.05 μM PO43-), a high N+P pulse (HiNP, 5 μM NO3- + 0.5 μM PO43-
#) and a combined N+P+Fe pulse (NPFe, 5 μM NO3- + 0.5 μM PO43- 0.5 nM + Fe added). 

#LoNP
log10(0.002+0.5)
-0.2992963

#HiNP
log10(0.002+5)
0.6991437

#NPFe
log10(0.002+5)
0.6991437

merged %>% semi_join(mixed, by = c("taxa")) %>%
  mutate(latitude = round(latitude)) %>%
  group_by(taxa, latitude, Timepoint, nit, xg_pred) %>%
  summarize(n = n(), avgProb = mean(probability)) %>%
  ggplot(aes(x = nit, size = n, y = xg_pred, color = xg_pred)) + geom_point() +  
  facet_grid(rows = vars(taxa), cols = vars(latitude)) +
  theme_bw() + 
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.y = element_text(size = 16), strip.text.x = element_text(size = 10), 
        axis.text = element_text(size = 16), axis.title = element_text(size = 16)) + 
  theme(legend.title=element_text(size=16), legend.text=element_text(size=14)) +
  labs(x = "log10(NO3_NO2)", y = "Prediction", size = "Number of predictions") + 
  scale_color_manual(values = c("Phototrophic" = "deepskyblue2", "Heterotrophic" = "red", "Mixotrophic" = "black"))+
  guides(color = "none")

#from plotAllTrophicModePredictions.R
insitu <- read_csv("g2SurfaceInsituPredictions.csv")

head(insitu)
insitu %>% filter(is.na(probability))
insitu <- insitu %>% group_by(taxa, latitude, xg_pred) %>% summarize(n = n())

merged %>% distinct(latitude) 
insitu %>% distinct(latitude)

insitu <- insitu %>% mutate(Treatment = "in situ")
insitu <- insitu %>% mutate(latitude = round(latitude))

insitu <- insitu %>% mutate(latTreat = str_c(latitude, "_", Treatment))

insitu <- insitu %>% mutate(latitude = str_c(as.character(latitude), " °N"))

insitu <- insitu %>% mutate(xg_pred = str_replace(xg_pred, "ic", "y"))

insitu <- insitu %>% mutate(Timepoint = 1)

order <- c("33_in situ", "33_Ctrl_0", "33_Ctrl_96", "33_LoNP_96", "33_HiNP_96", "33_NPFe_96", "37_in situ", "37_Ctrl_0", "37_Ctrl_96", "37_Fe_96", "37_NP_96", "37_NPFe_96", "41_in situ", "41_Ctrl_0", "41_Ctrl_96", "41_LFe_96", "41_HFe_96", "41_NPFe_96")
order_noTime <- str_replace(order, "_[0-9]{1,}", "")
order_noTime <- unique(order_noTime)

order <- c("33_in situ", "33_Ctrl_0", "33_Ctrl_96", "33_LoNP_96", "33_HiNP_96", "33_NPFe_96", "37_in situ", "37_Ctrl_0", "37_Ctrl_96", "37_Fe_96", "37_NP_96", "37_NPFe_96", "41_in situ", "41_Ctrl_0", "41_Ctrl_96", "41_LFe_96", "41_HFe_96", "41_NPFe_96")

merged %>% semi_join(mixed, by = c("taxa")) %>%
  mutate(latitude = round(latitude)) %>%
  filter(Timepoint != 0) %>%
  mutate(latTreat = str_c(latitude, "_", Treatment, "_", Timepoint)) %>%
  mutate(latitude = str_c(latitude, " °N")) %>%
  mutate(xg_pred = str_replace(xg_pred, "ic", "y")) %>%
  group_by(taxa, latitude, Timepoint, Treatment, xg_pred, latTreat) %>%
  summarize(n = n(), avgProb = mean(probability)) %>% ungroup() %>% arrange(desc(n))

merged %>% semi_join(mixed, by = c("taxa")) %>%
  mutate(latitude = round(latitude)) %>%
  mutate(latTreat = str_c(latitude, "_", Treatment, "_", Timepoint)) %>%
  mutate(latitude = str_c(latitude, " °N")) %>%
  mutate(xg_pred = str_replace(xg_pred, "ic", "y")) %>%
  group_by(taxa, latitude, Timepoint, Treatment, xg_pred, latTreat) %>%
  summarize(n = n(), avgProb = mean(probability)) %>%
  bind_rows(insitu %>% semi_join(mixed, by = c("taxa"))) %>%
  mutate(latTreat = factor(latTreat, levels = order)) %>%
  mutate(taxa = factor(taxa, levels = c("Azadinium spinosum", "Karlodinium veneficum", "Pelagodinium beii", "Prorocentrum minimum"))) %>%
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
  scale_y_continuous(limits = c(0,6), breaks = c(0,2,4,6))

ggsave("g2IncTrophicPredictions.png", dpi = 600, height = 12, width = 17)


merged %>% anti_join(mixed, by = c("taxa")) %>%
  mutate(latitude = round(latitude)) %>%
  mutate(latTreat = str_c(latitude, "_", Treatment)) %>%
  mutate(latitude = str_c(latitude, " °N")) %>%
  mutate(xg_pred = str_replace(xg_pred, "ic", "y")) %>%
  group_by(taxa, latitude, Timepoint, Treatment, xg_pred, latTreat) %>%
  summarize(n = n(), avgProb = mean(probability)) %>%
  bind_rows(insitu %>% anti_join(mixed, by = c("taxa")) %>% semi_join(merged, by = c("taxa"))) %>%
  mutate(latTreat = factor(latTreat, levels = order_noTime)) %>%
  mutate(taxa = str_replace(taxa, " ", "\n")) %>%
  mutate(taxa = factor(taxa, levels = c("Dictyocha\nspeculum", "Aureococcus\nanophagefferens", "Triparma\npacifica", 
                                        "Pelagomonas\ncalceolata", "Bathycoccus\nprasinos", 
                                        "Prasinoderma\nsingulare", "Karenia\nbrevis", "Emiliania\nhuxleyi", "Oxytricha\ntrifallax"))) %>%
  ggplot(aes(x = latTreat, y = n, fill = xg_pred)) + geom_bar(stat = 'identity') +  
  facet_grid(rows = vars(taxa), cols = vars(latitude), scales = "free_x") +
  theme_bw() + 
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.y = element_text(size = 14, color = 'black'), strip.text.x = element_text(size = 20, color = 'black'), 
        axis.text = element_text(size = 18, color = 'black'), axis.title = element_text(size = 26, color = 'black')) +
  theme(legend.title=element_text(size=16), legend.text=element_text(size=16)) +
  labs(x = "", y = "Number of predictions", fill = "Trophic mode prediction") + 
  scale_fill_manual(values = c("Phototrophy" = "deepskyblue2", "Heterotrophy" = "red", "Mixotrophy" = "black"))+
  guides(color = "none")

ggsave("g2IncTrophicPredictions_nonMixed.png", dpi = 600, height = 15, width = 14)


#REXP1
#nitrate + nitrite (N+N) was 1.78 μM,
#phosphate was and 0.42 μM
#iron was 0.31 nM
log10(0.42)
-0.3767507

#REXP2 in situ N+N and phosphate
#nitrate 0.006 μM 
#phosphate 0.13 μM
#iron 0.43 nM
log10(0.13)
-0.8860566

#REXP3
#nitrate 0.002 μM
#and phosphate 0.07 μM
log10(0.07)
-1.154902

merged <- merged %>% mutate(p = ifelse(Expt == "REXP1" & Treatment == "Ctrl", -0.3767507, NA))
merged <- merged %>% mutate(p = ifelse(Expt == "REXP1" & Treatment == "HFe", -0.3767507, p))
merged <- merged %>% mutate(p = ifelse(Expt == "REXP1" & Treatment == "LFe", -0.3767507, p))
merged <- merged %>% mutate(p = ifelse(Expt == "REXP1" & Treatment == "NPFe", 0.1522883, p))

merged <- merged %>% mutate(p = ifelse(Expt == "REXP2" & Treatment == "Ctrl", -0.8860566, p))
merged <- merged %>% mutate(p = ifelse(Expt == "REXP2" & Treatment == "Fe", -0.8860566, p))
merged <- merged %>% mutate(p = ifelse(Expt == "REXP2" & Treatment == "NPFe", -0.2006595, p))
merged <- merged %>% mutate(p = ifelse(Expt == "REXP2" & Treatment == "NP", -0.2006595, p))

merged <- merged %>% mutate(p = ifelse(Expt == "REXP3" & Treatment == "Ctrl", -1.154902, p))
merged <- merged %>% mutate(p = ifelse(Expt == "REXP3" & Treatment == "HiNP", -0.2441251, p))
merged <- merged %>% mutate(p = ifelse(Expt == "REXP3" & Treatment == "NPFe", -0.2441251, p))
merged <- merged %>% mutate(p = ifelse(Expt == "REXP3" & Treatment == "LoNP", -0.9208188, p))

#Two iron amendments were given: a low-iron pulse (LoFe; 0.3 nM Fe), a high-iron pulse (HiFe; 2 nM
#Fe), and a combined N+P+Fe pulse (NPFe; 2 nM Fe + 10 μM NO3
#- + 1 μM PO4
#3- added).

#For REXP1, the in situ log10N:Fe of 3.8 was
#shifted to 2.9 (HiFe), 3.5 (LoFe) or 3.7 (NPFe) through amendments with different ratios of N, P,
#and Fe.

#rexp1 
#hife
log10(0.42)
-0.3767507

#lofe
log10(0.42)
-0.3767507

#NPFe
log10((0.42+1))
0.1522883

#In REXP2, the in situ log10N:Fe of 2.1 was shifted to 1.6 (Fe), 3.5 (NPFe) or 4.0 (NP).

#Resource ratio experiment 2 (REXP2) used water collected from Station 11 (latitude 37.00°N,
#ambient water temperature 15.2°C) and were given one iron amendment (Fe; 1 nM Fe), an N+P
#amendment (NP; 0.5 μM P + 5 μM NO3) and an N+P+Fe amendment (NPFe; 1 nM Fe + 5 μM
#NO3- + 0.5 μM PO43- added). 

#Fe
log10(0.13)
-0.8860566

#NP
log10((0.13+.5))
-0.2006595

#NPFe
log10((0.13+.5))
-0.2006595

#REXP3, the in situ log10N:Fe of 1.7 was shifted to 3.4 (NP), 3.8 (NPFe), or 4.4 (NP) (Fig 3.4b).

log10(0.002/(0.22*.001))

#low N+P pulse (LoNP, 0.5 μM NO3- + 0.05 μM PO43-), a high N+P pulse (HiNP, 5 μM NO3- + 0.5 μM PO43-
#) and a combined N+P+Fe pulse (NPFe, 5 μM NO3- + 0.5 μM PO43- 0.5 nM + Fe added). 

#LoNP
log10(0.07+0.05)
-0.9208188

#HiNP
log10(0.07+0.5)
-0.2441251

#NPFe
log10(0.07+ 0.5)
-0.2441251

merged %>% semi_join(mixed, by = c("taxa")) %>%
  mutate(latitude = round(latitude)) %>%
  group_by(taxa, latitude, Timepoint, p, xg_pred) %>%
  summarize(n = n(), avgProb = mean(probability)) %>%
  ggplot(aes(x = p, size = n, y = xg_pred, color = xg_pred)) + geom_point() +  
  facet_grid(rows = vars(taxa), cols = vars(latitude)) +
  theme_bw() + 
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.y = element_text(size = 16), strip.text.x = element_text(size = 10), 
        axis.text = element_text(size = 16), axis.title = element_text(size = 16)) + 
  theme(legend.title=element_text(size=16), legend.text=element_text(size=14)) +
  labs(x = "log10(NO3_NO2)", y = "Prediction", size = "Number of predictions") + 
  scale_color_manual(values = c("Phototrophic" = "deepskyblue2", "Heterotrophic" = "red", "Mixotrophic" = "black"))+
  guides(color = "none")


