library(tidyverse)
library(readxl)

setwd("~/Dropbox/grad/research/g3/")

merged <- read_csv("G3PA.contig_dat_all_speciesWithTrophicPredictions_summary_updatedMarferret_marmicroDb2023_noOutliers_fall2023_28.csv")

merged_notOfInt <- read_csv("G3PA.contig_dat_all_otherThanSpeciesWithTrophicPredictions_summary_updatedMarferret_marmicroDb2023_noOutliers_fall2023_28.csv")

#downloaded from https://docs.google.com/spreadsheets/d/1lKwvlNjUV8dXHKO5HkMNa2m70QQYydrTf-Fb6PRQ9oE/edit#gid=0
sample <- read_csv("../metaT sample log - metaT.csv")

sample <- sample %>% filter(Cruise == "Gradients3" | Cruise == "Gradients 3")

sample <- sample %>% filter(Type == "polyA")

merged <- merged %>% mutate(Sample.ID = str_extract(sample, "UW[0-9]{1,}"))

merged_notOfInt <- merged_notOfInt %>% mutate(Sample.ID = str_extract(sample, "UW[0-9]{1,}"))

merged <- merged %>% mutate(Sample.ID = ifelse(is.na(Sample.ID), sample, Sample.ID))

merged_notOfInt <- merged_notOfInt %>% mutate(Sample.ID = ifelse(is.na(Sample.ID), sample, Sample.ID))

sample <- sample %>% mutate(Sample.ID = ifelse(Exp == "Diel", str_c("S",Station,"C",as.character(Cast)), Sample.ID))

merged %>% anti_join(sample, by = c("Sample.ID"))
merged_notOfInt %>% anti_join(sample, by = c("Sample.ID"))

sample <- sample %>% select(Cruise:Exp)

nrow(merged)
merged <- merged %>% left_join(sample %>% select(-c(Date, Time, Date.Time)) %>% distinct(), by = c("Sample.ID"))
nrow(merged)

nrow(merged_notOfInt)
merged_notOfInt <- merged_notOfInt %>% left_join(sample %>% select(-c(Date, Time, Date.Time)) %>% distinct(), by = c("Sample.ID"))
nrow(merged_notOfInt)

stan <- read_csv("G3PA_underway_dawn_norm_factors.txt", col_names = TRUE)

stan <- stan %>% mutate(sample = str_replace(sample_name, "\\.flash", ""))
stan <- stan %>% mutate(sample = str_replace(sample, "G3PA\\.diel\\.", ""))
stan <- stan %>% mutate(sample = str_replace(sample, "\\.[A-Z]$", ""))
stan <- stan %>% select(sample, NORM_FACTOR)

merged %>% anti_join(stan, by = c("sample"))
merged_notOfInt %>% anti_join(stan, by = c("sample"))

nrow(merged)
merged <- merged %>% left_join(stan, by = c("sample"))
nrow(merged)

nrow(merged_notOfInt)
merged_notOfInt <- merged_notOfInt %>% left_join(stan, by = c("sample"))
nrow(merged_notOfInt)

str(merged)
str(merged_notOfInt)

merged <- merged %>% mutate(absoluteCounts = est_counts*NORM_FACTOR)
merged_notOfInt <- merged_notOfInt %>% mutate(absoluteCounts = est_counts*NORM_FACTOR)

speciesPred <- read_csv("G3_trophicPredictions_cleanedUp_updatedMarferret_marmicroDb2023_noOutliers.csv")

speciesPred %>% anti_join(merged, by = c("taxa" = "tax_name"))

phot <- read_csv("../g1G2G3/phototrophicSpeciesStrict_dielIncubations.csv") %>% mutate(type = "Phototrophic")
het <- read_csv("../g1G2G3/heterotrophicSpeciesStrict_dielIncubations.csv") %>% mutate(type = "Heterotrophic")
mix <- read_csv("../g1G2G3/mixedSpeciesStrict_dielIncubations.csv") %>% mutate(type = "Mixed predictions")

labs = bind_rows(phot, het, mix)

labs <- labs %>% distinct(taxa, type)

merged %>% anti_join(labs, by = c("tax_name" = "taxa"))

nrow(merged)
merged <- merged %>% left_join(labs %>% select(1:2), by = c("tax_name" = "taxa"))
nrow(merged)

merged %>% group_by(Latitude, type) %>%
  summarize(sd = sd(absoluteCounts), absoluteCounts = mean(absoluteCounts), n = n()) %>% 
  write_csv("absoluteCountsByTrophicMode.csv")

merged %>% 
  write_csv("absoluteCountsByTrophicMode_full.csv")

both <- bind_rows(merged %>% mutate(taxaGroup = "speciesWithTrophicMode"), merged_notOfInt %>% mutate(taxaGroup = "Rest"))

both <- both %>% filter(!is.na(tax_name))

both <- both %>% 
  filter(str_detect(tax_name, " ") & !str_detect(tax_name, "unclassified") & !str_detect(tax_name, "uncultured") & 
           !(tax_name %in% c("Chlorella clade", "core chlorophytes", "CS clade", "Elliptochloris clade", "cellular organisms",
                             "Marine_stramenopiles_MASTs group", "Micromonas <green algae>", "Nitzschia <diatoms>", "PX clade", "Rhodomonas <cryptomonads>",
                             "Calanus finmarchicus", "Hydra vulgaris", "cellular organisms", "core chlorophytes",
                             "uncultured eukaryote", "Micromonas <green algae>", "Rhodomonas <cryptomonads>", "PX clade",
                             "CS clade", "Chlorella clade", "Elliptochloris clade", "Nitzschia <diatoms>", "Marine_stramenopiles_MASTs group",
                             "Neocalanus flemingeri", "Oikopleura dioica", "Paracyclopina nana", "Vannella robusta", "Pleuromamma xiphias")))

both <- both %>% mutate(genus = str_extract(tax_name, "[A-Za-z0-9]{1,}"))
both %>% filter(is.na(genus))

exclude <- c("Apostichopus", "Caulerpa", "Acartia", "Adineta", "Amphimedon", "Anisakis", "Aplysia", "Bugula", "Calanus", "Capitella", "Chara", "Chlorokybus", "Chondrus", 
             "Cladosiphon", "Compsopogon", "Crassostrea", "Daphnia", "Debaryomyces", "Dibothriocephalus", "Ectocarpus", "Eucyclops", 
             "Eurytemora", "Gracilariopsis", "Helobdella", "Hippoglossus", "Homarus", "Klebsormidium", "Labidocera", "Lepeophtheirus", 
             "Lingula", "Nemacystus", "Nematostella", "Neopyropia", "Octopus", "Opisthorchis", "Porphyra", "Saccharina", "Strongylocentrotus", 
             "Tigriopus", "Trichuris", "Tursiops", "Ulva", "Undaria", "Vaucheria")

both <- both %>% filter(!(genus %in% exclude))

withTrophMode <- both %>% group_by(taxaGroup) %>% summarize(absoluteCounts = sum(absoluteCounts)) %>% 
  filter(taxaGroup == "speciesWithTrophicMode") %>% select(absoluteCounts)

rest <- both %>% group_by(taxaGroup) %>% summarize(absoluteCounts = sum(absoluteCounts)) %>% 
  filter(taxaGroup == "Rest") %>% select(absoluteCounts)

#73.29 28 species
round(100*withTrophMode$absoluteCounts/(withTrophMode$absoluteCounts+rest$absoluteCounts),2)




