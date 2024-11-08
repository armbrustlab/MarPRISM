library(tidyverse)
library(readxl)
library(ggrepel)

setwd("~/Dropbox/grad/research/g2Incubation/")

dat <- read_csv("g2Inc_ofInt_28.csv")

mixed <- read_csv("../g1G2G3/mixedSpeciesStrict_dielIncubations.csv") %>% mutate(type = "mixed")
phot <- read_csv("../g1G2G3/phototrophicSpeciesStrict_dielIncubations.csv") %>% mutate(type = "phot")
het <- read_csv("../g1G2G3/heterotrophicSpeciesStrict_dielIncubations.csv") %>% mutate(type = "het")

type <- bind_rows(mixed, phot, het)

type <- type %>% select(taxa, type)

dat %>% anti_join(type, by = c("tax_name" = "taxa"))

nrow(dat)
dat <- dat %>% left_join(type, by = c("tax_name" = "taxa"))
nrow(dat)

sample <- read_csv("G2.RR_exp.metadata.csv")

dat <- dat %>% mutate(`Sample name` = str_extract(sample, "MS[0-9]{1,}"))

dat %>% anti_join(sample, by = c("Sample name"))

nrow(dat)
dat <- dat %>% left_join(sample, by = c("Sample name"))
nrow(dat)

#downloaded from https://github.com/armbrustlab/armbrust-metat/blob/main/Count_standards_workflow/Normalization_factors/G2PA_REXP_norm_factors.csv
stan <- read_csv("G2PA_REXP_norm_factors.csv")

stan <- stan %>% mutate(`Sample name` = str_replace(sample_name, ".sam", ""))

dat %>% anti_join(stan, by = c("Sample name"))

nrow(dat)
dat <- dat %>% left_join(stan %>% select(`Sample name`, NORM_FACTOR), by = c("Sample name"))
nrow(dat)

str(dat)

dat <- dat %>% mutate(absoluteCounts = est_counts*NORM_FACTOR)

colnames(dat)

dat <- dat %>% select(tax_name, type, `Size Fraction`, Expt, Timepoint, Treatment, Rep, absoluteCounts)

head(dat)
dat <- dat %>% spread(key = `Size Fraction`, value = absoluteCounts, fill = 0)

head(dat)
dat <- dat %>% mutate(absoluteCounts = `0.2um`+`3um`)

head(dat)

dat <- dat %>% mutate(NtoFe = ifelse(Expt == "REXP1" & Treatment == "Ctrl", 3.759058, NA))
dat <- dat %>% mutate(NtoFe = ifelse(Expt == "REXP1" & Treatment == "HFe", 2.886808, NtoFe))
dat <- dat %>% mutate(NtoFe = ifelse(Expt == "REXP1" & Treatment == "LFe", 3.46509, NtoFe))
dat <- dat %>% mutate(NtoFe = ifelse(Expt == "REXP1" & Treatment == "NPFe", 3.707533, NtoFe))

dat <- dat %>% mutate(NtoFe = ifelse(Expt == "REXP2" & Treatment == "Ctrl", 1.144683, NtoFe))
dat <- dat %>% mutate(NtoFe = ifelse(Expt == "REXP2" & Treatment == "Fe", 0.6228152, NtoFe))
dat <- dat %>% mutate(NtoFe = ifelse(Expt == "REXP2" & Treatment == "NPFe", 3.544155, NtoFe))
dat <- dat %>% mutate(NtoFe = ifelse(Expt == "REXP2" & Treatment == "NP", 4.066022, NtoFe))

dat <- dat %>% mutate(NtoFe = ifelse(Expt == "REXP3" & Treatment == "Ctrl", 0.9586073, NtoFe))
dat <- dat %>% mutate(NtoFe = ifelse(Expt == "REXP3" & Treatment == "HiNP", 4.356547, NtoFe))
dat <- dat %>% mutate(NtoFe = ifelse(Expt == "REXP3" & Treatment == "NPFe", 3.841638, NtoFe))
dat <- dat %>% mutate(NtoFe = ifelse(Expt == "REXP3" & Treatment == "LoNP", 3.356548, NtoFe))

dat <- dat %>% mutate(latitude = ifelse(Expt == "REXP1", 41.42, NA))
dat <- dat %>% mutate(latitude = ifelse(Expt == "REXP2", 37.00, latitude))
dat <- dat %>% mutate(latitude = ifelse(Expt == "REXP3", 32.93, latitude))


order <- c("32.93_Ctrl_0", "32.93_Ctrl_96", "32.93_LoNP_96", "32.93_HiNP_96", "32.93_NPFe_96", "37_Ctrl_0", "37_Ctrl_96", "37_Fe_96", "37_NP_96", "37_NPFe_96", "41.42_Ctrl_0", "41.42_Ctrl_96", "41.42_LFe_96", "41.42_HFe_96", "41.42_NPFe_96")
order_noTime <- str_replace(order, "_[0-9]{1,}", "")
order_noTime <- unique(order_noTime)


dat_summary <- dat %>% ungroup() %>% 
  filter(Timepoint != 0) %>%
  mutate(latTreat = str_c(latitude, "_", Treatment)) %>%
  group_by(latitude, Timepoint, Treatment, latTreat, type, tax_name) %>% 
  summarize(sd = sd(absoluteCounts), absoluteCounts = mean(absoluteCounts), numSamples = n())

dat_summary <- dat_summary %>% ungroup()

dat_summary <- dat_summary %>% mutate(se = sd/sqrt(numSamples))

dat_summary <- dat_summary %>% mutate(cv = se/absoluteCounts)

dat_summary <- dat_summary %>% mutate(cv_sq = cv^2)

dat_summary <- dat_summary %>% select(-sd, -numSamples, -se, -cv)

dat_summary <- dat_summary %>% ungroup() %>% select(-Timepoint, -Treatment)


dat_summary_mean <- dat_summary %>% select(-cv_sq) %>% spread(key = latTreat, value = absoluteCounts, fill = 0)

dat_summary_cv_sq <- dat_summary %>% select(-absoluteCounts) %>% spread(key = latTreat, value = cv_sq, fill = NA)

dat_summary <- dat_summary_mean


dat_summary <- dat_summary %>% mutate(diff_32.93_LoNP = `32.93_LoNP`/`32.93_Ctrl`)

dat_summary <- dat_summary %>% mutate(diff_32.93_HiNP = `32.93_HiNP`/`32.93_Ctrl`)

dat_summary <- dat_summary %>% mutate(diff_32.93_NPFe = `32.93_NPFe`/`32.93_Ctrl`)


dat_summary <- dat_summary %>% mutate(diff_37_Fe = `37_Fe`/`37_Ctrl`)

dat_summary <- dat_summary %>% mutate(diff_37_NP = `37_NP`/`37_Ctrl`)

dat_summary <- dat_summary %>% mutate(diff_37_NPFe = `37_NPFe`/`37_Ctrl`)


dat_summary <- dat_summary %>% mutate(diff_41.42_LFe = `41.42_LFe`/`41.42_Ctrl`)

dat_summary <- dat_summary %>% mutate(diff_41.42_HFe = `41.42_HFe`/`41.42_Ctrl`)

dat_summary <- dat_summary %>% mutate(diff_41.42_NPFe = `41.42_NPFe`/`41.42_Ctrl`)


###calculate error bars

dat_summary_cv_sq <- dat_summary_cv_sq %>% mutate(diff_32.93_LoNP = sqrt(`32.93_LoNP`+`32.93_Ctrl`))

dat_summary_cv_sq <- dat_summary_cv_sq %>% mutate(diff_32.93_HiNP = sqrt(`32.93_HiNP`+`32.93_Ctrl`))

dat_summary_cv_sq <- dat_summary_cv_sq %>% mutate(diff_32.93_NPFe = sqrt(`32.93_NPFe`+`32.93_Ctrl`))


dat_summary_cv_sq <- dat_summary_cv_sq %>% mutate(diff_37_Fe = sqrt(`37_Fe`+`37_Ctrl`))

dat_summary_cv_sq <- dat_summary_cv_sq %>% mutate(diff_37_NP = sqrt(`37_NP`+`37_Ctrl`))

dat_summary_cv_sq <- dat_summary_cv_sq %>% mutate(diff_37_NPFe = sqrt(`37_NPFe`+`37_Ctrl`))


dat_summary_cv_sq <- dat_summary_cv_sq %>% mutate(diff_41.42_LFe = sqrt(`41.42_LFe`+`41.42_Ctrl`))

dat_summary_cv_sq <- dat_summary_cv_sq %>% mutate(diff_41.42_HFe = sqrt(`41.42_HFe`+`41.42_Ctrl`))

dat_summary_cv_sq <- dat_summary_cv_sq %>% mutate(diff_41.42_NPFe = sqrt(`41.42_NPFe`+`41.42_Ctrl`))



dat_summary <- dat_summary %>% select(latitude, type, tax_name, 16:24)

dat_summary_cv_sq <- dat_summary_cv_sq %>% select(latitude, type, tax_name, 16:24)


dat_summary <- dat_summary %>% gather(diff_32.93_LoNP:diff_41.42_NPFe, key = "treatment", value = "diff")

dat_summary_cv_sq <- dat_summary_cv_sq %>% gather(diff_32.93_LoNP:diff_41.42_NPFe, key = "treatment", value = "mult")



nrow(dat_summary)
dat_summary <- dat_summary %>% left_join(dat_summary_cv_sq, by = c("tax_name", "treatment", "latitude", "type"))
nrow(dat_summary)

dat_summary <- dat_summary %>% mutate(se = diff*mult)


dat_summary %>% filter(!is.na(diff)) %>% 
  distinct(type, tax_name) %>% 
  group_by(type) %>% summarize(n = n())

dat_summary %>% group_by(type) %>% distinct(tax_name) %>% summarize(n = n())

dat_summary %>% filter(!is.na(diff)) %>%
  mutate(latitude = str_c(as.character(latitude), " °N")) %>%
  ggplot(aes(x = treatment, y = diff, color = type)) + 
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 123), alpha = .5) + 
  geom_linerange(aes(ymin=diff-se, ymax=diff+se),
                position = position_jitter(width = 0.15, height = 0, seed = 123), alpha = .5) +
  geom_hline(yintercept = 1, linetype = "dashed") + 
  geom_text_repel(data = dat_summary %>% filter(diff > 5 | (diff > 4 & latitude > 41)) %>%  mutate(latitude = str_c(as.character(latitude), " °N")), 
                  aes(label = tax_name), max.overlaps = 20) + 
  facet_wrap(~latitude, scales = "free_x") + 
  scale_color_manual(values = c("orange", "purple", "forestgreen")) +
  theme_classic() +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size = 26, color = 'black')) + 
  theme(strip.text.y = element_text(size = 26, color = 'black')) + 
  theme(axis.text.x = element_text(size = 22, color = 'black'))  + 
  theme(axis.text.y = element_text(size = 22, color = 'black')) + 
  theme(axis.title.y = element_text(size = 26, color = 'black')) + 
  theme(legend.text = element_text(size = 20, color = 'black')) +
  theme(legend.title = element_text(size = 22, color = 'black')) +
  theme(axis.title.x = element_text(size = 26, color = 'black')) +  
  scale_y_continuous(breaks = c(0,1,5,10,15,20,25), limits = c(0,27.5)) + 
  labs(x = "Treatment", y = "Mean ratio of transcripts in treatment:control", shape = "", color = "")

ggsave("responseBySpecies.png", dpi = 300, height = 7, width = 12)

            