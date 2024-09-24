library(tidyverse)
library(readxl)

setwd("~/Dropbox/grad/research/g3Diel/")

dat <- read_csv("g3_diel_allSamples_processed_updatedMarferret_marmicroDb2023_noOutliers_speciesWithTrophicPredictions_fall2023_28.csv")

notOfInt <- read_csv("g3_diel_allSamples_processed_updatedMarferret_marmicroDb2023_noOutliers_otherThanSpeciesWithTrophicPredictions_fall2023_28.csv")

stan <- read_csv("G3PA_diel_norm_factors.csv") %>% select(sample_name, NORM_FACTOR)

stan <- stan %>% mutate(sample = str_extract(sample_name, "S[0-9]{1,}C[0-9]{1,}\\.[ABC]"))

dat <- dat %>% mutate(sample = str_replace(sample, "G3PA.diel.", ""))
notOfInt <- notOfInt %>% mutate(sample = str_replace(sample, "G3PA.diel.", ""))

dat %>% anti_join(stan, by = c("sample"))
notOfInt %>% anti_join(stan, by = c("sample"))

head(stan)
colnames(stan)[2] <- "convFactor"

nrow(dat)
dat <- dat %>% left_join(stan, by = c("sample"))
nrow(dat)

nrow(notOfInt)
notOfInt <- notOfInt %>% left_join(stan, by = c("sample"))
nrow(notOfInt)

str(dat)
str(notOfInt)

dat <- dat %>% mutate(absoluteCounts = est_counts*convFactor)
notOfInt <- notOfInt %>% mutate(absoluteCounts = est_counts*convFactor)

notOfInt <- notOfInt %>% filter(!is.na(tax_name))

dat %>% select(tax_name, sample, absoluteCounts) %>% write_csv("speciesWithTrophicModePredictionsAbundance.csv")


both <- bind_rows(dat %>% mutate(taxaGroup = "speciesWithTrophicMode"), notOfInt %>% mutate(taxaGroup = "Rest"))

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

#63.68 28 species
round(100*withTrophMode$absoluteCounts/(withTrophMode$absoluteCounts+rest$absoluteCounts),2)

