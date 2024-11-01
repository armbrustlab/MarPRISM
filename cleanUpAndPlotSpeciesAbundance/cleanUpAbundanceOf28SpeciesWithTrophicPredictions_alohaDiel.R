library(tidyverse)
library(readxl)

setwd("~/Dropbox/grad/research/alohaDiel/")

dat <- read_csv("G1_diel_allSamples_processed_updatedMarferret_marmicroDb2023_noOutliers_speciesWithTrophicPredictions_fall2023_28.csv")

notOfInt <- read_csv("G1_diel_allSamples_processed_updatedMarferret_marmicroDb2023_noOutliers_otherThanSpeciesWithTrophicPredictions_fall2023_28.csv")

stan <- read_csv("diel1_SUM_norm_factors.csv")

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

type1 <- read_excel("~/Dropbox/grad/research/g1Surface/typeOfTaxaInPredictions_noOutliers.xlsx")
type2 <- read_excel("~/Dropbox/grad/research/g2Surface/typeOfTaxaInPredictions_noOutliers.xlsx")
type3 <- read_excel("~/Dropbox/grad/research/g3/typeOfTaxaInPredictions_noOutliers.xlsx")

type <- bind_rows(type1, type2, type3)

type <- type %>% mutate(genus = str_extract(taxa, "[A-Za-z]{1,}"))
type <- type %>% distinct(genus, group)
type <- type %>% filter(!is.na(genus))

nrow(both)
both <- both %>% left_join(type %>% distinct(genus, group), by = c("genus"))
nrow(both)

group <- read_csv("~/Dropbox/grad/research/g2Surface/cleanedUpIFCB.csv")

nrow(both)
both <- both %>% left_join(group %>% distinct(var, group.x), by = c("genus" = "var"))
nrow(both)

both <- both %>% mutate(group = ifelse(is.na(group), group.x, group))
both <- both %>% select(-group.x)

taxInfo <- read_csv("../tax_v11.csv")

taxInfo <- taxInfo %>% distinct(genus, class_)
exclude <- taxInfo %>% group_by(genus) %>% distinct(class_) %>% summarize(n = n()) %>% filter(n > 1) %>% distinct(genus)

taxInfo <- taxInfo %>% filter(!(genus %in% exclude$genus))

nrow(both)
both <- both %>% left_join(taxInfo %>% distinct(genus, class_), by = c("genus"))
nrow(both)

both <- both %>% mutate(group = ifelse((is.na(group) & class_ == "Dinophyceae"), "Dinoflagellate", group))

both <- both %>% mutate(class_ = ifelse(is.na(class_), genus, class_))

#"Abedinium", "Apocalathium", "Breviolum", "Fugacium" are the only dinoflagellates in this list
both %>% filter(is.na(group)) %>% distinct(class_)

dinos <- c("Abedinium", "Apocalathium", "Breviolum", "Fugacium")

both <- both %>% mutate(group = ifelse(class_ %in% dinos, "Dinoflagellate", group))

both <- both %>% mutate(group = ifelse(is.na(group), "Nondinoflagellate", group))

both <- both %>% mutate(group = ifelse(group != "Dinoflagellate", "Nondinoflagellate", group))

both %>% distinct(group)

both <- both %>% mutate(absoluteCounts_dinoCorr = ifelse(group == "Dinoflagellate", absoluteCounts/6.4, absoluteCounts))

both %>% distinct(tax_name, group) %>% write_csv("../dinoflagellateSpecies.csv")

withTrophMode <- both %>% group_by(taxaGroup) %>% summarize(absoluteCounts_dinoCorr = sum(absoluteCounts_dinoCorr)) %>% 
  filter(taxaGroup == "speciesWithTrophicMode") %>% select(absoluteCounts_dinoCorr)

rest <- both %>% group_by(taxaGroup) %>% summarize(absoluteCounts_dinoCorr = sum(absoluteCounts_dinoCorr)) %>% 
  filter(taxaGroup == "Rest") %>% select(absoluteCounts_dinoCorr)

#47.43 28 species
round(100*withTrophMode$absoluteCounts_dinoCorr/(withTrophMode$absoluteCounts_dinoCorr+rest$absoluteCounts_dinoCorr),2)







