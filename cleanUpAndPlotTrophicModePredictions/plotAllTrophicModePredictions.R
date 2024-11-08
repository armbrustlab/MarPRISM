library(tidyverse)
library(readxl)
library(ggpubr)
library(ggpmisc)

setwd("~/Dropbox/grad/research/")

###clean up and merge datasets

g1 <- read_csv("g1Surface/G1_surface_trophicPredictions_marFerret2023_cleanedUp_noOutliers.csv")
g2 <- read_csv("g2Surface/G2_surface_trophicPredictions_cleanedUp_updatedMarferret_marmicroDb2023_noOutliers.csv")
g3 <- read_csv("g3/G3_trophicPredictions_cleanedUp_updatedMarferret_marmicroDb2023_noOutliers.csv")

colnames(g1) <- c("sample", "taxa", "xg_pred", "probability", "size", "depth", "latitude", "station") 

g2 <- g2 %>% select(-Station, -`Time Start (HST)`)

colnames(g2) <- c("sample", "taxa", "xg_pred", "probability", "size", "station", "cast", "depth", "latitude") 

colnames(g3) <- c("sample", "taxa", "xg_pred", "probability", "Sample.ID", "Cruise", "Cruise.ID", "Type", "station", "cast", "size", "depth", "latitude", "longitude", "exp") 

g3 <- g3 %>% select(sample, station, taxa, xg_pred, probability, size, depth, latitude)

g1$size <- parse_number(g1$size)
g2$size <- parse_number(g2$size)

g3$station <- as.character(g3$station)

g <- bind_rows(g1 %>% mutate(cruise = "g1"), g2 %>% mutate(cruise = "g2"), g3 %>% mutate(cruise = "g3"))

alohaDiel <- read_csv("alohaDiel/G1_trophicPredictions_cleanedUp_updatedMarferret_marmicroDb2023_diel_noOutliers.csv") %>% 
  select(-time.y)

colnames(alohaDiel) <- c("sample", "taxa", "xg_pred", "probability", "time", "station", "date")

g3Diel <- read_csv("g3Diel/G3_trophicPredictions_cleanedUp_updatedMarferret_marmicroDb2023_diel_noOutliers.csv") %>% 
  select(-Cruise, -Type, -Station, -Cast, -Longitude, -Exp, -Cruise.ID)

colnames(g3Diel) <- c("taxa", "xg_pred", "probability", "sample", "size", "depth", "date", "time", "latitude", "hourMin", 
                      "hour", "justHour")

g$depth <- as.character(g$depth)
g3Diel$depth <- as.character(g3Diel$depth)

g3Diel$time <- as.character(g3Diel$time)

g3Diel <- g3Diel %>% mutate(time = str_replace(time, "\\:00", ""))
g3Diel <- g3Diel %>% mutate(time = str_replace(time, "\\:", ""))
g3Diel %>% distinct(time)

alohaDiel <- alohaDiel %>% mutate(date = str_extract(date, "[0-9\\-]{1,}"))
alohaDiel %>% distinct(date)

g3Diel$time <- as.numeric(g3Diel$time)
g3Diel %>% distinct(time)

g <- bind_rows(g, alohaDiel %>% mutate(cruise = "alohaDiel") %>% mutate(xg_pred = str_replace(xg_pred, "y$", "ic")), 
               g3Diel %>% mutate(cruise = "g3Diel"))

g %>% distinct(xg_pred)

g3Depth <- read_csv("g3Depth/g3Depth_trophicPredictions_cleanedUp_updatedMarferret_marmicroDb2023_noOutliers.csv")

g3Depth$sample <- as.character(g3Depth$sample)
g3Depth$depth <- as.character(g3Depth$depth)

g <- bind_rows(g, g3Depth %>% mutate(cruise = "g3Depth"))

type1 <- read_excel("g1Surface/typeOfTaxaInPredictions_noOutliers.xlsx")
type1 <- type1 %>% mutate(species = str_extract(tax_name, "[A-Z]{1,}[ a-z\\.0-9\\-A-Z']{1,}$"))

type2 <- read_excel("g2Surface/typeOfTaxaInPredictions_noOutliers.xlsx")
type2 <- type2 %>% mutate(species = str_extract(taxa, "[A-Z]{1,}[ a-z\\.0-9\\-A-Z]{1,}$"))

type3 <- read_excel("g3/typeOfTaxaInPredictions_noOutliers.xlsx")
type3 <- type3 %>% mutate(species = str_extract(taxa, "[A-Z]{1,}[ a-z\\.0-9\\-A-Z]{1,}$"))

colnames(type2)[1] <- "tax_name"
colnames(type3)[1] <- "tax_name"

type <- bind_rows(type1, type2, type3)

g %>% anti_join(type, by = c("taxa" = "species")) %>% distinct(taxa) %>% arrange(taxa)

nrow(g)
g <- g %>% left_join(type %>% distinct(species, group), by = c("taxa" = "species"))
nrow(g)

g <- g %>% mutate(group = ifelse(taxa == "Calcidiscus leptoporus", "Haptophyte", group))
g %>% filter(is.na(group))

g <- g %>% mutate(typeCruise = str_c(group, "_", cruise))
g <- g %>% mutate(taxaCruise = str_c(taxa, "_", cruise))

g %>% distinct(sample, taxa, xg_pred, cruise) %>% write_csv("g1G2G3/allSpeciesTrophicModePredictions_depthDiel.csv")

g2Inc <- read_csv("g2Incubation/g2Inc_trophicPredictions_cleanedUp_updatedMarferret_marmicroDb2023_noOutliers.csv")

g2Inc <- g2Inc %>% mutate(station = str_c(Expt, Treatment, Timepoint))

g2Inc <- g2Inc %>% select(taxa, xg_pred, probability, `Sample name`, station)

colnames(g2Inc)[4] <- "sample"

g2Inc$cruise <- "G2 incubations"

nrow(g2Inc)
g2Inc <- g2Inc %>% left_join(type %>% distinct(species, group), by = c("taxa" = "species"))
nrow(g2Inc)

g2Inc %>% filter(is.na(group))

g <- bind_rows(g, g2Inc)

g %>% write_csv("g1G2G3/allSpeciesTrophicModePredictions_depthDielIncubations.csv")

g %>% distinct(taxa) %>% write_csv("g1G2G3/28SpeciesWithTrophicPredictions.csv")

g %>% distinct(sample) %>% nrow()

###plot 

g %>% filter(is.na(group))

g$groupNum <- factor(g$group)
g$groupNum <- as.numeric(g$groupNum)

g %>% distinct(group, groupNum)

g$groupNum <- factor(g$groupNum, levels = c(1,4,8,2,5,6,3,7))

taxaList <- g %>% ungroup() %>% arrange(groupNum, taxa) %>% distinct(taxa, group) %>% distinct(taxa)
taxaList <- taxaList$taxa
taxaList <- list(taxaList)
taxaList <- taxaList[[1]]

g %>% group_by(taxa) %>% summarize(n = n()) %>% arrange(desc(n))

g %>% 
  mutate(xg_pred = factor(xg_pred, levels = c("Phototrophic", "Mixotrophic", "Heterotrophic"))) %>%
  group_by(taxa, groupNum, xg_pred) %>% summarize(numPredictions = n()) %>% 
  arrange(groupNum) %>%
  mutate(taxa = factor(taxa, levels = rev(taxaList)))  %>%
  ggplot(aes(y = taxa, x = numPredictions, fill = xg_pred)) + 
  geom_bar(stat = 'identity') +
  theme_classic() + labs(y = "", x = "Number of predictions", fill = "Trophic mode prediction") + 
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1, size = 22, color = "black"))  + 
  scale_fill_manual(values = c("Phototrophic" = "deepskyblue2", "Heterotrophic" = "red", "Mixotrophic" = "black"))  + 
  theme(axis.text.x = element_text(size = 22, color = 'black')) + 
  theme(axis.title.x = element_text(size = 26, color = 'black')) + 
  theme(legend.text = element_text(size = 22, color = 'black')) +
  theme(legend.title = element_text(size = 26, color = 'black')) +
  theme(axis.text.y=element_text(angle=0)) + scale_x_continuous(limits = c(0,300))

ggsave("g1G2G3/trophicModePreditionsBySpeciesOverG1AndG2AndG3_groupedAcrossCruises_noOutliers_withDepth_tall_dielIncubations.png", dpi = 300, width = 16, height = 16)
ggsave("g1G2G3/trophicModePreditionsBySpeciesOverG1AndG2AndG3_groupedAcrossCruises_noOutliers_withDepth_tall_dielIncubations.svg", dpi = 300, width = 16, height = 16)

mixed_strict <- g 

mixed_strict <- mixed_strict %>% group_by(taxa, xg_pred) %>% summarize(n = n())

mixed_strict_total <- mixed_strict %>% group_by(taxa) %>% summarize(total = sum(n))

nrow(mixed_strict)
mixed_strict <- mixed_strict %>% left_join(mixed_strict_total, by = c("taxa"))
nrow(mixed_strict)

mixed_strict <- mixed_strict %>% mutate(prop = n/total)

mixedPreds <- mixed_strict %>% 
  ungroup() %>%
  group_by(taxa) %>% arrange(desc(prop)) %>% slice(1) %>% filter(prop < .90)

mixedPreds %>% 
  ungroup() %>%
  distinct(taxa)

g %>% filter(!is.na(time)) %>% distinct(cruise)

reps <- g

reps_total <- reps %>% group_by(taxa, latitude, depth, station, cruise, time, date) %>% summarize(numSamples = n())

reps <- reps %>% group_by(taxa, latitude, depth, station, xg_pred, cruise, time, date) %>% summarize(numPredictions = n())

nrow(reps)
reps <- reps %>% left_join(reps_total, by = c("taxa", "latitude", "depth", "station", "cruise", "time", "date"))
nrow(reps)

reps %>% arrange(desc(numSamples))
reps %>% arrange(numSamples)

reps <- reps %>% ungroup() 

reps %>% filter(numPredictions != numSamples)
reps %>% filter(numPredictions == numSamples)

reps %>% distinct(taxa)
reps <- reps %>% filter(numPredictions == numSamples)

mixedPreds %>% 
  ungroup() %>%
  distinct(taxa) 

mixedPreds <- mixedPreds %>% 
  ungroup() %>%
  distinct(taxa) %>%
  semi_join(reps %>% group_by(taxa) %>% distinct(xg_pred) %>% summarize(n = n()) %>% 
              filter(n > 1), by = c("taxa"))

mixedPreds %>%
  write_csv("~/Dropbox/grad/research/g1G2G3/mixedSpeciesStrict_dielIncubations.csv")

photPreds <- mixed_strict %>% 
  ungroup() %>%
  anti_join(mixedPreds, by = c("taxa")) %>%
  arrange(desc(prop)) %>% 
  group_by(taxa) %>% 
  slice(1) %>% 
  filter(xg_pred == "Phototrophic") %>% 
  distinct(taxa)

photPreds %>%
  write_csv("~/Dropbox/grad/research/g1G2G3/phototrophicSpeciesStrict_dielIncubations.csv")

hetPreds <- 
  mixed_strict %>% 
  ungroup() %>%
  anti_join(mixedPreds, by = c("taxa")) %>%
  arrange(desc(prop)) %>% group_by(taxa) %>% 
  slice(1) %>% filter(xg_pred == "Heterotrophic") %>% 
  distinct(taxa) 

hetPreds %>%
  write_csv("~/Dropbox/grad/research/g1G2G3/heterotrophicSpeciesStrict_dielIncubations.csv")

nrow(mixedPreds)
nrow(hetPreds)
nrow(photPreds)

g %>% distinct(taxa) %>% nrow()

g <- g %>% mutate(xg_pred = str_replace(xg_pred, "ic", "y"))

dat <- g %>% 
  filter(cruise %in% c("g1", "g2", "g3", "g3Depth")) %>%
  semi_join(mixedPreds, by = c("taxa")) %>%
  mutate(cruise = ifelse(cruise == "g3Depth", "g3", cruise)) %>%
  group_by(xg_pred, taxa, cruise, latitude) %>% summarize(numPredictions = n()) %>% 
  ungroup()  %>%
  mutate(taxa = str_replace(taxa, " ", "\n")) %>%
  mutate(cruise = str_replace(cruise, "g", "G")) %>% 
  mutate(cruise = ifelse(cruise == "G1", "Gradients1: 2016", cruise)) %>%
  mutate(cruise = ifelse(cruise == "G2", "Gradients2: 2017", cruise)) %>%
  mutate(cruise = ifelse(cruise == "G3", "Gradients3: 2019", cruise))

dashed <- dat %>% select(xg_pred, taxa, cruise, numPredictions, latitude) %>% distinct() %>% 
  spread(key = xg_pred, value = numPredictions, fill = 0) %>% 
  gather(Heterotrophy:Phototrophy, key = xg_pred, value = numPredictions)

dashed <- dashed %>% filter(cruise == "Gradients1: 2016")

dat$taxa <- factor(dat$taxa, levels = 
                     c("Karlodinium\nveneficum",
                       "Prorocentrum\nminimum", 
                       "Bolidomonas\nsp. 1657", 
                       "Chrysochromulina\nsp. KB-HA01", 
                       "Pelagodinium\nbeii",   
                     "Azadinium\nspinosum",
                     "Scrippsiella\ntrochoidea", 
                     "Tripos\nfusus", 
                     "Prymnesium\npolylepis"))

dashed$taxa <- factor(dashed$taxa, levels = 
                        c("Karlodinium\nveneficum",
                          "Prorocentrum\nminimum", 
                          "Bolidomonas\nsp. 1657", 
                          "Chrysochromulina\nsp. KB-HA01", 
                          "Pelagodinium\nbeii",   
                          "Azadinium\nspinosum",
                          "Scrippsiella\ntrochoidea", 
                          "Tripos\nfusus", 
                          "Prymnesium\npolylepis"))

dat %>% arrange(desc(numPredictions)) %>% head()

dat %>%
  arrange(taxa) %>%
  mutate(cruise = ifelse(cruise == "G1", "Gradients1: 2016", cruise)) %>%
  mutate(cruise = ifelse(cruise == "G2", "Gradients2: 2017", cruise)) %>%
  mutate(cruise = ifelse(cruise == "G3", "Gradients3: 2019", cruise)) %>%
  ggplot(aes(x = latitude, y = xg_pred)) + 
  geom_point(aes(size = numPredictions, color = xg_pred)) + 
  geom_point(aes(size = numPredictions),pch=21, color = "black") + 
  facet_grid(rows = vars(taxa), cols = vars(cruise)) + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white")) +
  theme(strip.text.y = element_text(size = 12, color = 'black'), strip.text.x = element_text(size = 20, color = 'black'), 
        axis.text = element_text(size = 18, color = 'black'), axis.title = element_text(size = 26, color = 'black')) +
  theme(strip.background =element_rect(fill="white")) + 
  theme(legend.title=element_text(size=22), legend.text=element_text(size=20)) +
  labs(x = "Latitude (°N)", y = "", size = "Number of predictions") +
  geom_vline(data=filter(dashed, cruise=="Gradients1: 2016"), aes(xintercept=32.15), colour="red", linetype="dashed") + 
  geom_vline(data=filter(dashed, cruise=="Gradients1: 2016"), aes(xintercept=33), colour="black", linetype="dashed") + 
  
  geom_vline(data=filter(dashed %>% mutate(cruise = "Gradients2: 2017"), cruise=="Gradients2: 2017"), aes(xintercept=32.5), colour="red", linetype="dashed") + 
  geom_vline(data=filter(dashed %>% mutate(cruise = "Gradients2: 2017"), cruise=="Gradients2: 2017"), aes(xintercept=36.2), colour="black", linetype="dashed") + 
  
  geom_vline(data=filter(dashed %>% mutate(cruise = "Gradients3: 2019"), cruise=="Gradients3: 2019"), aes(xintercept=32.45), colour="red", linetype="dashed") + 
  geom_vline(data=filter(dashed %>% mutate(cruise = "Gradients3: 2019"), cruise=="Gradients3: 2019"), aes(xintercept=35), colour="black", linetype="dashed") + 
  scale_color_gradient(low = "white", high = "red", limits = c(.3, 1), breaks = c(.3, .65, 1), 
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))  +
  scale_size_continuous(breaks = c(2,4,6,8), limits = c(1,8)) + 
  scale_color_manual(values = c("Phototrophy" = "deepskyblue2", "Heterotrophy" = "red", "Mixotrophy" = "black")) + 
  guides(color="none") +
  scale_x_continuous(breaks = c(25,30,35,40))

ggsave("g1G2G3/trophicModePredictionsBySpeciesTypeOverG1AndG2AndG3_mixedPredictions_noOutliers_withG3Depth.png", dpi = 600, width = 15, height = 14)

#write g2 preditions to a csv file for plotting with g2 incubations
insitu <- g %>%
  filter(cruise == "g2") %>% 
  filter(latitude %in% c(32.93, 37.0, 41.42))

insitu %>% write_csv("g2Incubation/g2SurfaceInsituPredictions.csv")

##fix g3 sample names 

g3sample <- read_csv("metaT sample log - metaT.csv")
g3sample <- g3sample %>% filter(Cruise == "Gradients3" | Cruise == "Gradients 3")
g3sample <- g3sample %>% filter(Type == "polyA")
g3sample <- g3sample %>% filter(str_detect(Sample.ID, "UW"))
g3sample <- g3sample %>% select(Sample.ID, Filter, Replicate, Latitude)
g3sample <- g3sample %>% mutate(Sample.ID = str_c("G3PA.", Sample.ID))
g3sample <- g3sample %>% mutate(newSample = str_c(Latitude, "_", Filter, ".", Replicate))

nrow(g)
g <- g %>% left_join(g3sample %>% select(-Replicate, -Filter), by = c("sample" = "Sample.ID"))
nrow(g)

g3_sample <- g %>% filter(cruise == "g3") %>% distinct(sample, newSample)

g <- g %>% mutate(sample = ifelse(cruise == "g3", newSample, sample))

g <- g %>% select(-newSample)

g %>% filter(cruise == "g3") %>% select(sample) %>% head()

g <- g %>% mutate(station = ifelse(cruise == "g3", as.character(latitude), station))
g <- g %>% mutate(station = ifelse(cruise == "g3Diel", str_extract(sample, "S[0-9]{1,}C[0-9]{1,}"), station))
g <- g %>% mutate(station = ifelse(cruise == "g3Depth", str_c(station, "_", Depth), station))

##fix cruise names 
g <- g %>% mutate(cruise = ifelse(cruise == "alohaDiel", "Aloha diel", cruise)) %>% 
  mutate(cruise = ifelse(cruise == "g1", "G1", cruise)) %>% 
  mutate(cruise = ifelse(cruise == "g2", "G2", cruise)) %>% 
  mutate(cruise = ifelse(cruise == "g3", "G3", cruise)) %>% 
  mutate(cruise = ifelse(cruise == "g3Diel", "G3 diel", cruise)) %>%
  mutate(cruise = ifelse(cruise == "g3Depth", "G3 depth", cruise)) 

aloha <- read_csv("alohaDiel/pfamSummary.csv")
g1 <- read_csv("g1Surface/corePfamSummary.csv")
colnames(g1)[1] <- "sample"

g2 <- read_csv("g2Surface/pfamSummary.csv")
colnames(g2)[1] <- "sample"

g3 <- read_csv("g3/pfamSummary.csv")
g3_diel <- read_csv("g3Diel/pfamSummary.csv")
g3_depth <- read_csv("g3Depth/corePfamSummary.csv")
colnames(g3_depth)[2] <- "tax_id"

g2_inc <- read_csv("g2Incubation/pfamSummary.csv")

corePfam <- bind_rows(aloha %>% mutate(cruise = "Aloha diel"), 
                      g1 %>% mutate(cruise = "G1"), 
                      g2 %>% mutate(cruise = "G2"), 
                      g3 %>% mutate(cruise = "G3"), 
                      g3_diel %>% mutate(cruise = "G3 diel"), 
                      g3_depth %>% mutate(cruise = "G3 depth"), 
                      g2_inc %>% mutate(cruise = "G2 incubations"))

corePfam <- corePfam %>% mutate(sample = ifelse(cruise == "G2 incubations", str_extract(sample, "MS[0-9]{1,}"), sample))
corePfam <- corePfam %>% mutate(sample = ifelse(cruise == "G3 diel", str_replace(sample, "G3PA.diel.", ""), sample))

nrow(corePfam)
corePfam <- corePfam %>% left_join(g3_sample, by = c("sample"))
nrow(corePfam)

corePfam <- corePfam %>% mutate(sample = ifelse(cruise == "G3", newSample, sample)) %>% select(-newSample)

g %>% anti_join(corePfam, by = c("sample", "taxa" = "tax_name", "cruise")) %>% distinct(cruise)

nrow(g)
g <- g %>% left_join(corePfam %>% select(-tax_id, -numCorePfams), by = c("sample", "taxa" = "tax_name", "cruise"))
nrow(g)

total <- g %>% group_by(taxa, station, cruise) %>% summarize(total = n())

dat <- g %>% group_by(taxa, station, cruise, xg_pred) %>% summarize(n = n()) %>% 
  left_join(total, by = c("taxa", "station", "cruise")) %>% 
  filter(total != 1) %>%
  mutate(prop = n/total) %>% 
  ungroup() %>%
  group_by(taxa, station, cruise) %>% arrange(desc(prop)) %>% 
  slice(1)

dat <- dat %>% ungroup() 

dat %>% group_by(xg_pred) %>% summarize(prop = mean(prop))

dat %>% 
  ggplot(aes(fill = xg_pred, x = prop)) + 
  geom_bar(position = 'dodge') +
  labs(x = "Proportion of predictions in agreement\nwith top trophic mode prediction", 
       fill = "Top trophic\nmode prediction", y = "Count") +
  theme_classic() +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(axis.text.x = element_text(size = 22, color = 'black'))  + 
  theme(axis.text.y = element_text(size = 22, color = 'black')) + 
  theme(axis.title.y = element_text(size = 23, color = 'black')) + 
  theme(legend.text = element_text(size = 20, color = 'black')) +
  theme(legend.title = element_text(size = 22, color = 'black')) +
  theme(axis.title.x = element_text(size = 23, color = 'black')) +
  scale_fill_manual(values = c("Phototrophy" = "deepskyblue2", "Heterotrophy" = "red", "Mixotrophy" = "black")) + 
  xlim(0,1.02)

ggsave("g1G2G3/numberOfTrophicModePredictionsByPropInAgreement.png", dpi = 600, height = 6, width = 10)

dat %>% group_by(taxa) %>% summarize(prop = mean(prop)) %>% arrange(prop)

nrow(dat)
dat <- dat %>% left_join(g %>% group_by(taxa, station, cruise) %>% summarize(meanPropCore = mean(propCorePfams)), by = c("taxa", "station", "cruise"))
nrow(dat)

dat %>% ggplot(aes(x = meanPropCore, y = prop)) + geom_jitter()+
  labs(x = "Mean proportion of CTGs in species bin", y = "Proportion of predictions in agreement\nwith top trophic mode prediction") + 
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
  geom_vline(xintercept = .7, linetype = 'dashed') + 
  ylim(0,1.02)

ggsave("g1G2G3/propInAgreementByPropCorePfams.png", dpi = 600, height = 7, width = 10)


g %>% 
  ggplot(aes(x = propCorePfams, fill = xg_pred)) + geom_bar(width = .0016, alpha = .8) +
  scale_fill_manual(values = c("Phototrophy" = "deepskyblue2", "Heterotrophy" = "red", "Mixotrophy" = "black")) +
  labs(x = "Proportion of CTGs in species bin", y = "Number of predictions", fill = "Trophic mode\nprediction") + 
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
  geom_vline(xintercept = .7, linetype = 'dashed') + 
  ylim(0,25)

ggsave("g1G2G3/numberOfTrophicModePredictionsByPropCTGs.png", dpi = 600, height = 6, width = 10)

table <- g %>% ungroup() %>% 
  group_by(cruise, taxa, group) %>% summarize(n = n()) %>% 
  ungroup()

colnames(table) <- c("Cruise", "Species", "Taxonomic group", "Number of trophic predictions")

table %>% distinct(Cruise)

table %>% 
  mutate(Cruise = str_replace(Cruise, "G", "Gradients")) %>% 
  mutate(Cruise = str_replace(Cruise, "Aloha", "ALOHA")) %>%
  arrange(Cruise, desc(`Number of trophic predictions`)) %>%
  write_csv("g1G2G3/speciesInCruiseSupplementaryTable.csv")


