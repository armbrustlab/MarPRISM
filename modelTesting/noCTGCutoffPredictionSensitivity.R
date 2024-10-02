library(tidyverse)

setwd("~/Dropbox/grad/research/g1Surface/")

trop <- read_csv("G1_trophicModePredictions_surface_updatedMarferret_marmicroDb2023_noOutliers_fall2023_noCTGCutoff")
meta <- read_csv("G1_surface_allSamples_processed_updatedMarferret_marmicroDb2023_noOutliers_tpm_noNAPfam_fall2023_noNATaxa.csv")

colnames(meta)[1:5]

merged <- cbind(meta[1:2], trop)

#Mix, Het, Phot
#1, 0, 2,
merged <- merged %>% mutate(xg_pred = ifelse(xg_pred == 1, "Mix", ifelse(xg_pred == 0, "Het", "Phot")))

colnames(merged)
colnames(merged)[3] <- "xg_pred_mine"

trop_ben <- read_csv("G1_trophicModePredictions_surface_updatedMarferret_marmicroDb2023_noOutliers_fall2023_lambert_noCTGCutoff")

trop_ben <- trop_ben %>% mutate(xg_pred = ifelse(xg_pred == 1, "Mix", ifelse(xg_pred == 0, "Het", "Phot")))

colnames(trop_ben)
colnames(trop_ben) <- "xg_pred_ben"

merged <- cbind(merged, trop_ben)

colnames(merged)[1]
colnames(merged)[1] <- "sample"

merged <- merged %>% mutate(Size = ifelse(str_detect(sample, "0.2um"), "0.2um", NA))
merged <- merged %>% mutate(Size = ifelse(str_detect(sample, "3um"), "3um", Size))
merged %>% filter(is.na(Size))

sample <- read_csv("ctd_g1_g2_15m_source_Stephen.csv")
sample <- sample %>% select(SAMPLE_ID, DEPTH, LATITUDE, DATETIME)
sample <- sample %>% as.data.frame()
sample <- sample %>% mutate(SAMPLE_ID = str_replace(SAMPLE_ID, "0_2um", "0.2um"))

sample <- sample %>% select(-DATETIME)

merged %>% anti_join(sample, by = c("sample" = "SAMPLE_ID"))

nrow(merged)
merged <- merged %>% left_join(sample, by = c("sample" = "SAMPLE_ID"))
nrow(merged)

head(merged)

#exclude non species and non protists
merged <- merged %>% filter(!tax_name %in% c("Calanus finmarchicus", "Hydra vulgaris", "Neocalanus flemingeri", 
                                             "Oikopleura dioica", "Paracyclopina nana", 
                                             "Micromonas <green algae>", "unclassified Acanthoeca", 
                                             "unclassified Chrysochromulina", "unclassified Micromonas", 
                                             "unclassified Phaeocystis", "cellular organisms", 
                                             "Pleuromamma xiphias", "uncultured eukaryote") & 
                              str_detect(tax_name, " "))

exclude <- c("Apostichopus", "Caulerpa", "Acartia", "Adineta", "Amphimedon", "Anisakis", "Aplysia", "Bugula", "Calanus", "Capitella", "Chara", "Chlorokybus", "Chondrus", 
             "Cladosiphon", "Compsopogon", "Crassostrea", "Daphnia", "Debaryomyces", "Dibothriocephalus", "Ectocarpus", "Eucyclops", 
             "Eurytemora", "Gracilariopsis", "Helobdella", "Hippoglossus", "Homarus", "Klebsormidium", "Labidocera", "Lepeophtheirus", 
             "Lingula", "Nemacystus", "Nematostella", "Neopyropia", "Octopus", "Opisthorchis", "Porphyra", "Saccharina", "Strongylocentrotus", 
             "Tigriopus", "Trichuris", "Tursiops", "Ulva", "Undaria", "Vaucheria")

merged <- merged %>% mutate(genus = str_extract(tax_name, "[A-Z][a-z]{1,}"))

merged <- merged %>% filter(!(genus %in% exclude))

allProtists <- read_csv("../g1G2G3/absoluteCountsAllTaxaGroups.csv")

merged %>% anti_join(allProtists, by = c("tax_name")) %>% distinct(tax_name)

merged %>% filter(str_detect(tax_name, "unclassified")) %>% distinct(tax_name)

merged <- merged %>% filter(!str_detect(tax_name, "unclassified"))

merged %>% filter(str_detect(tax_name, "<")) %>% distinct(tax_name)

merged <- merged %>% filter(!str_detect(tax_name, "<"))

merged %>% filter(str_detect(tax_name, "group")) %>% distinct(tax_name)

merged <- merged %>% filter(!str_detect(tax_name, "group"))

merged %>% filter(str_detect(tax_name, "clade")) %>% distinct(tax_name)

merged <- merged %>% filter(!str_detect(tax_name, "Chlorella clade"))
merged <- merged %>% filter(!str_detect(tax_name, "CS clade"))
merged <- merged %>% filter(!str_detect(tax_name, "PX clade"))
merged <- merged %>% filter(!str_detect(tax_name, "Elliptochloris clade"))

merged <- merged %>% filter(!str_detect(tax_name, "core chlorophytes"))

merged %>% distinct(tax_name) %>% arrange(tax_name) %>% View()

ncol(meta)

colnames(meta)[1:3]
colnames(meta)[7939]
meta <- meta %>% gather(PF00001:PF20154, key = "pfam", value = "tpm")

head(meta)
meta <- meta %>% filter(tpm > 0)

core <- read_csv("../MarFERReT.v1.core_genes.csv")

head(meta)
head(core)

core <- core %>% mutate(pfam = str_replace(pfam_id, "\\.[0-9]{1,}", ""))

core <- core %>% filter(lineage == "Eukaryota")

meta <- meta %>% semi_join(core, by = c("pfam"))

head(meta)
meta <- meta %>% group_by(var, tax_name) %>% distinct(pfam) %>% summarize(numCore = n())
meta <- meta %>% ungroup()

meta <- meta %>% mutate(propCore = numCore/605)

head(merged)
head(meta)

nrow(merged)
merged <- merged %>% left_join(meta, by = c("sample" = "var", "tax_name"))
nrow(merged)

head(merged)
merged %>% group_by(tax_name, LATITUDE) %>% summarize(n = n()) %>% ungroup() %>% arrange(desc(n))


head(merged)
merged <- merged %>% gather(xg_pred_mine:xg_pred_ben, key = "model", value = "xg_pred")

head(merged)

total <- merged %>% 
  mutate(propCore = round(propCore, 2)) %>%
  group_by(model, propCore) %>% summarize(total = n())

merged %>% 
  mutate(propCore = round(propCore, 2)) %>%
  group_by(model, propCore, xg_pred) %>% summarize(n = n()) %>% 
  left_join(total, by = c("model", "propCore")) %>% 
  mutate(prop = n/total) %>% ungroup() %>% 
  group_by(model, propCore) %>% 
  summarize(prop = sum(prop)) %>% ungroup() %>% 
  distinct(prop)

merged %>% 
  mutate(xg_pred = ifelse(xg_pred == "Het", "Heterotrophic", xg_pred)) %>% 
  mutate(xg_pred = ifelse(xg_pred == "Phot", "Phototrophic", xg_pred)) %>% 
  mutate(xg_pred = ifelse(xg_pred == "Mix", "Mixotrophic", xg_pred)) %>% 
  mutate(propCore = round(propCore, 2)) %>%
  group_by(model, propCore, xg_pred) %>% summarize(n = n()) %>% 
  left_join(total, by = c("model", "propCore")) %>% 
  mutate(model = ifelse(model == "xg_pred_ben", "Previous model\nLambert et al. 2022", "MarPRISM")) %>%
  mutate(model = factor(model, levels = c("MarPRISM", "Previous model\nLambert et al. 2022"))) %>%
  mutate(prop = n/total) %>%
  ggplot(aes(x = propCore, y = prop, color = xg_pred)) + geom_point() + facet_wrap(~model) + 
  geom_smooth() + 
  scale_color_manual(values = c("Phototrophic" = "deepskyblue2", "Heterotrophic" = "red", "Mixotrophic" = "black")) +
  scale_y_continuous(limits = c(-0.1,1.1), breaks = c(0,.25,.5,.75,1)) + 
  xlim(0,1) + 
  labs(x = "Proportion of CTGs in species bin", y = "Proportion of predictions", color = "Trophic mode\nprediction") + 
  theme_classic() +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(axis.text.x = element_text(size = 22, color = 'black'))  + 
  theme(axis.text.y = element_text(size = 22, color = 'black')) + 
  theme(axis.title.y = element_text(size = 23, color = 'black')) + 
  theme(legend.text = element_text(size = 20, color = 'black')) +
  theme(legend.title = element_text(size = 22, color = 'black')) +
  theme(axis.title.x = element_text(size = 23, color = 'black')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(strip.text.x = element_text(size = 23, color = 'black')) +
  geom_vline(xintercept = .7, linetype = 'dashed')
  
ggsave("propPredictionsForEachTrophicModeVsPropCTGs.png", dpi = 600, height = 7, width = 12)

total <- merged %>% group_by(tax_name, LATITUDE, model) %>% summarize(total = n())

total %>% ungroup() %>% arrange(desc(total))

dat <- merged %>% group_by(tax_name, LATITUDE, model, xg_pred) %>% summarize(n = n()) %>% 
  left_join(total, by = c("tax_name", "LATITUDE", "model")) %>% 
  filter(total != 1) %>%
  mutate(prop = n/total) %>% 
  ungroup() %>%
  group_by(tax_name, LATITUDE, model) %>% arrange(desc(prop)) %>% 
  slice(1)

head(dat)

head(merged)
merged %>% group_by(tax_name, LATITUDE) %>% summarize(n = n()) %>% arrange(desc(n))
merged <- merged %>% filter(model == "xg_pred_mine") %>% group_by(tax_name, LATITUDE) %>% summarize(propCore = mean(propCore))
merged <- merged %>% ungroup()


nrow(dat)
dat <- dat %>% left_join(merged, by = c("tax_name", "LATITUDE"))
nrow(dat)

dat %>% ungroup() %>% 
  mutate(model = ifelse(model == "xg_pred_ben", "Previous model\nLambert et al. 2022", "MarPRISM")) %>%
  mutate(model = factor(model, levels = c("MarPRISM", "Previous model\nLambert et al. 2022"))) %>%
  ggplot(aes(x = propCore, y = prop, color = model)) + geom_jitter(alpha = .4) + geom_smooth() + 
  geom_vline(xintercept = .7, linetype = 'dashed') + 
  scale_y_continuous(limits = c(0,1.1), breaks = c(0,.25,.5,.75,1)) + 
  xlim(0,1) + 
  labs(x = "Proportion of CTGs in species bin", y = "Proportion of predictions in agreement\nwith top trophic mode prediction", color = "Model") + 
  theme_classic() +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(axis.text.x = element_text(size = 22, color = 'black'))  + 
  theme(axis.text.y = element_text(size = 22, color = 'black')) + 
  theme(axis.title.y = element_text(size = 23, color = 'black')) + 
  theme(legend.text = element_text(size = 20, color = 'black')) +
  theme(legend.title = element_text(size = 22, color = 'black')) +
  theme(axis.title.x = element_text(size = 23, color = 'black')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("propPredictionsInAgreementVsPropCTGs.png", dpi = 600, height = 7, width = 10)

head(dat)
total <- dat %>% ungroup() %>% 
  mutate(propCore = round(propCore, 2)) %>%
  group_by(propCore, model) %>% 
  summarize(total = n())

dat %>% ungroup() %>% 
  mutate(propCore = round(propCore, 2)) %>%
  filter(prop == 1) %>%
  group_by(propCore, model) %>% 
  summarize(n = n()) %>% 
  full_join(total, by = c("propCore", "model")) %>% 
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  mutate(prop = n/total) %>% 
  mutate(model = ifelse(model == "xg_pred_ben", "Previous model\nLambert et al. 2022", "MarPRISM")) %>%
  mutate(model = factor(model, levels = c("MarPRISM", "Previous model\nLambert et al. 2022"))) %>%
  ggplot(aes(x = propCore, y = prop, color = model)) + 
  geom_point(alpha = .4) + geom_smooth() + 
  geom_vline(xintercept = .7, linetype = 'dashed') + 
  scale_y_continuous(limits = c(0,1.27), breaks = c(0,.25,.5,.75,1)) + 
  scale_x_continuous(limits = c(0,1), breaks = c(0,.25,.5,.75,1)) + 
  labs(x = "Proportion of CTGs in species bin", y = "Proportion of species bins in\ncomplete agreement among replicates", color = "Model") + 
  theme_classic() +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(axis.text.x = element_text(size = 22, color = 'black'))  + 
  theme(axis.text.y = element_text(size = 22, color = 'black')) + 
  theme(axis.title.y = element_text(size = 23, color = 'black')) + 
  theme(legend.text = element_text(size = 20, color = 'black')) +
  theme(legend.title = element_text(size = 22, color = 'black')) +
  theme(axis.title.x = element_text(size = 23, color = 'black')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_manual(values = c("Purple", "Orange"))
  
ggsave("propPredictionsInTotalAgreementVsPropCTGs.png", dpi = 600, height = 7, width = 13)


dat %>% ungroup() %>% 
  mutate(propCore = round(propCore, 2)) %>%
  filter(prop == 1) %>%
  group_by(propCore, model) %>% 
  summarize(n = n()) %>% 
  full_join(total, by = c("propCore", "model")) %>% 
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  mutate(prop = n/total) %>% 
  filter(propCore >= .7) %>%
  mutate(model = ifelse(model == "xg_pred_ben", "Previous model\nLambert et al. 2022", "MarPRISM")) %>%
  mutate(model = factor(model, levels = c("MarPRISM", "Previous model\nLambert et al. 2022"))) %>%
  ggplot(aes(x = model, y = prop, color = model)) + 
  geom_point(position = position_jitter(width = .2, height = 0), alpha = .4) +
  geom_boxplot(fill = NA) + 
  labs(y = "Proportion of species bins in\ncomplete agreement among replicates", x = "") + 
  theme_classic() +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(axis.text.x = element_text(size = 22, color = 'black'))  + 
  theme(axis.text.y = element_text(size = 22, color = 'black')) + 
  theme(axis.title.y = element_text(size = 23, color = 'black')) + 
  theme(legend.text = element_text(size = 20, color = 'black')) +
  theme(legend.title = element_text(size = 22, color = 'black')) +
  theme(axis.title.x = element_text(size = 23, color = 'black')) + 
  ylim(0,1) + 
  theme(legend.position = "") + 
  scale_color_manual(values = c("Purple", "Orange"))

ggsave("propPredictionsInTotalAgreement_atLeast70PercCTGs.png", dpi = 600, height = 6.5, width = 7)



