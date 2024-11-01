library(tidyverse)
library(readxl)

setwd("~/Dropbox/grad/research/g2Surface/")

dat <- read_csv("G2PA.contig_dat_all_speciesWithTrophicPredictions_summary_updatedMarferret_marmicroDb2023_noOutliers_fall2023_28.csv")

notOfInt <- read_csv("G2PA.contig_dat_all_otherThanSpeciesWithTrophicPredictions_summary_updatedMarferret_marmicroDb2023_noOutliers_fall2023_28.csv")

dat <- dat %>% mutate(size = str_extract(var, "[0-9]um"))

notOfInt <- notOfInt %>% mutate(size = str_extract(var, "[0-9]um"))

dat %>% distinct(size)

sample <- read_excel("Gradients2_discrete_samples.xlsx", sheet = 2)

head(sample)
colnames(sample) <- sample[1,]
sample <- sample[-1,]
head(sample)

sample <- sample %>% select(Date:Notes)

dat <- dat %>% mutate(station = str_extract(var, "S[0-9]{1,}C[0-9]{1,}"))
notOfInt <- notOfInt %>% mutate(station = str_extract(var, "S[0-9]{1,}C[0-9]{1,}"))

dat <- dat %>% mutate(Station = parse_number(station))
notOfInt <- notOfInt %>% mutate(Station = parse_number(station))

dat <- dat %>% mutate(Cast = str_extract(var, "C[0-9]{1,}"))
dat <- dat %>% mutate(Cast = parse_number(Cast))

notOfInt <- notOfInt %>% mutate(Cast = str_extract(var, "C[0-9]{1,}"))
notOfInt <- notOfInt %>% mutate(Cast = parse_number(Cast))

sample$Station <- as.numeric(sample$Station)
notOfInt$Station <- as.numeric(notOfInt$Station)

sample$Cast <- as.numeric(sample$Cast)
notOfInt$Cast <- as.numeric(notOfInt$Cast)

dat <- dat %>% mutate(Size = ifelse(size == "2um", "0.2um", size))
notOfInt <- notOfInt %>% mutate(Size = ifelse(size == "2um", "0.2um", size))
dat %>% distinct(size)

sample <- sample %>% distinct(Station, Cast, `Depth (m)`, `Time Start (HST)`, Latitude, Size, `Sample ID`)

dat <- dat %>% mutate(`Sample ID` = str_extract(var, "[ABC]$"))
notOfInt <- notOfInt %>% mutate(`Sample ID` = str_extract(var, "[ABC]$"))

dat %>% anti_join(sample, by = c("Station", "Cast", "Size", "Sample ID"))
notOfInt %>% anti_join(sample, by = c("Station", "Cast", "Size", "Sample ID"))

nrow(dat)
dat <- dat %>% left_join(sample, by = c("Station", "Cast", "Size", "Sample ID"))
nrow(dat)

nrow(notOfInt)
notOfInt <- notOfInt %>% left_join(sample, by = c("Station", "Cast", "Size", "Sample ID"))
nrow(notOfInt)

dat %>% select(Latitude) %>% str()

dat$Latitude <- as.numeric(dat$Latitude)
notOfInt$Latitude <- as.numeric(notOfInt$Latitude)

stan <- read_csv("G2PA.standard_counts.NORM_FACTORS2.csv", col_names = FALSE)

stan <- stan %>% mutate(Station = str_extract(X1, "S[0-9]{1,}C[0-9]{1,}"))
stan %>% filter(is.na(Station))

stan <- stan %>% mutate(Size = ifelse(str_detect(X1, "0_2um"), "0.2um", NA))
stan <- stan %>% mutate(Size = ifelse(str_detect(X1, "3um"), "3um", Size))
stan %>% filter(is.na(Size))

stan <- stan %>% mutate(Replicate = str_extract(X1, "A|B|C$"))
stan %>% distinct(Replicate)

dat %>% anti_join(stan, by = c("station" = "Station", "Size", "Sample ID" = "Replicate"))

notOfInt %>% anti_join(stan, by = c("station" = "Station", "Size", "Sample ID" = "Replicate"))

head(stan) 
colnames(stan)[2] <- "convFactor"

nrow(dat)
dat <- dat %>% left_join(stan %>% select(-X1), by = c("station" = "Station", "Size", "Sample ID" = "Replicate"))
nrow(dat)

nrow(notOfInt)
notOfInt <- notOfInt %>% left_join(stan %>% select(-X1), by = c("station" = "Station", "Size", "Sample ID" = "Replicate"))
nrow(notOfInt)

dat <- dat %>% mutate(absoluteCounts = est_counts*convFactor)
notOfInt <- notOfInt %>% mutate(absoluteCounts = est_counts*convFactor)

taxaPred <- read_csv("G2_surface_trophicPredictions_cleanedUp_updatedMarferret_marmicroDb2023_noOutliers.csv")

taxaPred %>% anti_join(dat, by = c("taxa" = "tax_name"))

phot <- read_csv("../g1G2G3/phototrophicSpeciesStrict_dielIncubations.csv") %>% mutate(type = "Phototrophic")
het <- read_csv("../g1G2G3/heterotrophicSpeciesStrict_dielIncubations.csv") %>% mutate(type = "Heterotrophic")
mix <- read_csv("../g1G2G3/mixedSpeciesStrict_dielIncubations.csv") %>% mutate(type = "Mixed predictions")

labs <- bind_rows(phot, het, mix)

labs <- labs %>% distinct(taxa, type)

dat %>% anti_join(labs, by = c("tax_name" = "taxa"))

nrow(dat)
dat <- dat %>% left_join(labs %>% select(1:2), by = c("tax_name" = "taxa"))
nrow(dat)

dat %>% group_by(Latitude, type) %>%
  summarize(sd = sd(absoluteCounts), absoluteCounts = mean(absoluteCounts), n = n()) %>% 
  write_csv("absoluteCountsByTrophicMode.csv")

dat %>% 
  write_csv("absoluteCountsByTrophicMode_full.csv")

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

dinos <- read_csv("../dinoflagellateSpecies.csv")

both %>% anti_join(dinos, by = c("tax_name"))

nrow(both)
both <- both %>% left_join(dinos, by = c("tax_name"))
nrow(both)

both <- both %>% mutate(absoluteCounts_dinoCorr = ifelse(group == "Dinoflagellate", absoluteCounts/6.4, absoluteCounts))

withTrophMode <- both %>% group_by(taxaGroup) %>% summarize(absoluteCounts_dinoCorr = sum(absoluteCounts_dinoCorr)) %>% 
  filter(taxaGroup == "speciesWithTrophicMode") %>% select(absoluteCounts_dinoCorr)

rest <- both %>% group_by(taxaGroup) %>% summarize(absoluteCounts_dinoCorr = sum(absoluteCounts_dinoCorr)) %>% 
  filter(taxaGroup == "Rest") %>% select(absoluteCounts_dinoCorr)


#66.61 28 species
round(100*withTrophMode$absoluteCounts_dinoCorr/(withTrophMode$absoluteCounts_dinoCorr+rest$absoluteCounts_dinoCorr),2)




