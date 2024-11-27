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

dat <- dat %>% select(tax_name, type, `Size Fraction`, Expt, Timepoint, Treatment, Rep, absoluteCounts)

dat <- dat %>% spread(key = `Size Fraction`, value = absoluteCounts, fill = 0)

dat <- dat %>% mutate(absoluteCounts = `0.2um`+`3um`)

order <- c("32.93_Ctrl_0", "32.93_Ctrl_96", "32.93_LoNP_96", "32.93_HiNP_96", "32.93_NPFe_96", "37_Ctrl_0", "37_Ctrl_96", "37_Fe_96", "37_NP_96", "37_NPFe_96", "41.42_Ctrl_0", "41.42_Ctrl_96", "41.42_LFe_96", "41.42_HFe_96", "41.42_NPFe_96")
order_noTime <- str_replace(order, "_[0-9]{1,}", "")
order_noTime <- unique(order_noTime)

#At the northernmost REXP1 station at 41.42°N, in situ dissolved
#inorganic nitrogen (DIN) was 2 μM and iron (Fe) was 0.3 nM resulting in a molar N:Fe of 6.6 ×
#103 and log10N:Fe of 3.8. At the transition zone REXP2 station at 37.00°N, in situ nutrient
#concentrations were 0.06 μM DIN and 0.51 nM Fe resulting in a log10N:Fe of 2.1 At the
#southernmost REXP3 station at 32.93°N, in situ concentrations were 0.01 μM DIN and 0.22 nM
#Fe, resulting in a log10N:Fe of 1.7 (Fig 3.4 a, b).
dat <- dat %>% mutate(latitude = ifelse(Expt == "REXP1", 41.42, NA))
dat <- dat %>% mutate(latitude = ifelse(Expt == "REXP2", 37.00, latitude))
dat <- dat %>% mutate(latitude = ifelse(Expt == "REXP3", 32.93, latitude))

dat_summary <- dat %>% ungroup() %>% 
  filter(Timepoint != 0) %>%
  mutate(latTreat = str_c(latitude, "_", Treatment)) %>%
  group_by(latitude, Timepoint, Treatment, latTreat, type, tax_name) %>% 
  summarize(sd = sd(absoluteCounts), absoluteCounts = mean(absoluteCounts), numSamples = n())

dat_summary <- dat_summary %>% ungroup()

dat_summary %>% distinct(numSamples)

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

#calculate se of ratio (treatment:controll expression) to plot as error bars

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

dat_summary <- dat_summary %>% filter(!is.na(diff)) 

#calculate p-value for whether ratios (treatment:control expression) are more than 1
dat_summary <- dat_summary %>% mutate(n = 3)

dat_summary <- dat_summary %>%
  rowwise() %>%
  mutate(
    t_stat = (diff - 1) / se,
    df = n - 1,
    p_value = pt(t_stat, df = df, lower.tail = FALSE)
  ) %>%
  ungroup()

dat_summary %>% filter(p_value < .05) %>% arrange(treatment, type)

dat_summary %>% 
  mutate(treatment = factor(treatment, levels = c("diff_32.93_LoNP", "diff_32.93_HiNP", "diff_32.93_NPFe", 
                                                  "diff_37_Fe", "diff_37_NP", "diff_37_NPFe", 
                                                  "diff_41.42_LFe", "diff_41.42_HFe", "diff_41.42_NPFe"))) %>%
  mutate(latitude = str_c(as.character(latitude), " °N")) %>%
  ggplot(aes(x = treatment, y = diff, color = type)) + 
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 123), alpha = .5) + 
  geom_linerange(aes(ymin=diff-se, ymax=diff+se),
                 position = position_jitter(width = 0.15, height = 0, seed = 123), alpha = .5) +
  geom_hline(yintercept = 1, linetype = "dashed") + 
  geom_text_repel(data = dat_summary %>% filter(p_value < .05 | diff > 15) %>% mutate(latitude = str_c(as.character(latitude), " °N")), 
                  aes(label = tax_name), max.overlaps = 20, size = 3) + 
  facet_wrap(~latitude, scales = "free_x") + 
  scale_color_manual(values = c("orange", "purple", "forestgreen")) +
  theme_classic() +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size = 26, color = 'black')) + 
  theme(strip.text.y = element_text(size = 26, color = 'black')) + 
  theme(axis.text.x = element_text(size = 10, color = 'black'))  + 
  theme(axis.text.y = element_text(size = 22, color = 'black')) + 
  theme(axis.title.y = element_text(size = 26, color = 'black')) + 
  theme(legend.text = element_text(size = 20, color = 'black')) +
  theme(legend.title = element_text(size = 22, color = 'black')) +
  theme(axis.title.x = element_text(size = 26, color = 'black')) +  
  scale_y_continuous(breaks = c(0,1,5,10,15,20,25), limits = c(0,27.5)) + 
  labs(x = "Treatment", y = "Mean ratio of transcripts in treatment:control", shape = "", color = "")

ggsave("responseBySpecies.png", dpi = 300, height = 9, width = 14)
ggsave("responseBySpecies.svg", dpi = 300, height = 9, width = 14)
            
dat_summary <- dat %>% ungroup() %>% 
  mutate(latTreat = str_c(latitude, "_", Treatment)) %>%
  group_by(latitude, Timepoint, Treatment, type, tax_name) %>% 
  summarize(sd = sd(absoluteCounts), absoluteCounts = mean(absoluteCounts), numSamples = n())

dat_summary %>% ggplot(aes(x = Treatment, y = absoluteCounts, fill = type)) + geom_bar(stat = 'identity') + 
  facet_wrap(~latitude, scales = 'free')

#chlorophyll data from Nick Hawco 
chlor <- read_excel("Calculating Chl Trilogy Gradients.xlsx", sheet = "RR experiments")

chlor <- chlor %>% select(Station, Depth, `...17`, `...18`)
chlor <- chlor %>% filter(!is.na(`...17`))
chlor <- chlor %>% mutate(se = `...18`/sqrt(3))

#fill in missing station IDs
chlor <- chlor %>% mutate(Station = ifelse(is.na(Station), lag(Station), Station))
chlor <- chlor %>% mutate(Station = ifelse(is.na(Station), lag(Station), Station))
chlor <- chlor %>% mutate(Station = ifelse(is.na(Station), lag(Station), Station))
chlor <- chlor %>% mutate(Station = ifelse(is.na(Station), lag(Station), Station))

#fix station name
chlor <- chlor %>% mutate(Station = ifelse(Station == "TZ", "RR tZ", Station))

#make names for variables more clear
chlor <- chlor %>% select(-`...18`)
colnames(chlor) <- c("station", "treatment", "mean", "se")

#get rid of time zero hours
chlor <- chlor %>% filter(!str_detect(treatment, "T0"))

#make dataframe for control treatments and noncontrol treatments
control <- chlor %>% filter(str_detect(treatment, "control|Control"))
chlor <- chlor %>% filter(!str_detect(treatment, "control|Control"))

control <- control %>% select(-treatment)
colnames(control) <- c("station", "control_mean", "control_se")

#add control treatments as more columns
chlor <- chlor %>% left_join(control, by = c("station"))

#calculate  ratio of treatment:control mean chlorophyll
chlor <- chlor %>% mutate(diff = mean/control_mean)

#calculate se of ratio of treatment:control mean chlorophyll
chlor <- chlor %>% mutate(cv = se/mean)

chlor <- chlor %>% mutate(cv_sq = cv^2)

chlor <- chlor %>% mutate(control_cv = control_se/control_mean)

chlor <- chlor %>% mutate(control_cv_sq = control_cv^2)

chlor <- chlor %>% mutate(se_mult = sqrt(cv_sq+control_cv_sq))

chlor <- chlor %>% mutate(diff_se = diff*se_mult)

#fix station names 
chlor <- chlor %>% mutate(station = ifelse(station == "Gyre T0", "32.93 °N", station))
chlor <- chlor %>% mutate(station = ifelse(station == "RR tZ", "37 °N", station))
chlor <- chlor %>% mutate(station = ifelse(station == "HNLC", "41.42 °N", station))

chlor %>% arrange(station, treatment)
chlor <- chlor %>% mutate(treatment = ifelse(treatment == "FeNP A", "NPFe", treatment))
chlor <- chlor %>% mutate(treatment = ifelse(treatment == "HNP A", "HiNP", treatment))
chlor <- chlor %>% mutate(treatment = ifelse(treatment == "LNP A", "LoNP", treatment))

chlor <- chlor %>% mutate(treatment = ifelse(treatment == "Fe high", "HiFe", treatment))
chlor <- chlor %>% mutate(treatment = ifelse(treatment == "Fe low", "LoFe", treatment))

#put treatments in order
chlor$treatment <- factor(chlor$treatment, levels = c("LoNP", "HiNP", "Fe", "NP", "LoFe", "HiFe", "NPFe"))

chlor %>% ggplot(aes(x = treatment, y = diff)) + geom_bar(stat = 'identity', fill = 'green') + 
  geom_errorbar(aes(ymin=diff-diff_se, ymax=diff+diff_se), width = .5) +
  geom_hline(yintercept = 1, linetype = "dashed") + 
  theme_classic() +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size = 26, color = 'black')) + 
  theme(strip.text.y = element_text(size = 26, color = 'black')) + 
  theme(axis.text.x = element_text(size = 10, color = 'black'))  + 
  theme(axis.text.y = element_text(size = 22, color = 'black')) + 
  theme(axis.title.y = element_text(size = 26, color = 'black')) + 
  theme(legend.text = element_text(size = 20, color = 'black')) +
  theme(legend.title = element_text(size = 22, color = 'black')) +
  theme(axis.title.x = element_text(size = 26, color = 'black')) +  
  facet_wrap(~station, scales = "free_x") + 
  labs(x = "", y = "Mean ratio of chlorophyll\nin treatment:control") + 
  scale_y_continuous(breaks = c(0,1,2,4,6)) 

ggsave("responseInChlorophyll.png", dpi = 300, height = 4.2, width = 14)







