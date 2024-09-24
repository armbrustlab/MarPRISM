library(tidyverse)

setwd("~/Dropbox/grad/research/g1Surface/")

dat <- read_csv("G1PA.contig_dat_all_speciesWithTrophicPredictions_summary_updatedMarferret_marmicroDb2023_noOutliers_fall2023_28.csv")

notOfInt <- read_csv("G1PA.contig_dat_all_otherThanSpeciesWithTrophicPredictions_summary_updatedMarferret_marmicroDb2023_noOutliers_fall2023_28.csv")

dat <- dat %>% mutate(size = str_extract(var, "[0-9\\.]{1,}um"))

notOfInt <- notOfInt %>% mutate(size = str_extract(var, "[0-9\\.]{1,}um"))

dat %>% distinct(size)

sample <- read_csv("ctd_g1_g2_15m_source_Stephen.csv")

sample <- sample %>% select(1:11)

sample <- sample %>% mutate(SAMPLE_ID = str_replace(SAMPLE_ID, "0_2", "0.2"))

dat %>% anti_join(sample, by = c("var" = "SAMPLE_ID"))
notOfInt %>% anti_join(sample, by = c("var" = "SAMPLE_ID"))

nrow(dat)
dat <- dat %>% left_join(sample, by = c("var" = "SAMPLE_ID"))
nrow(dat)

nrow(notOfInt)
notOfInt <- notOfInt %>% left_join(sample, by = c("var" = "SAMPLE_ID"))
nrow(notOfInt)

stan <- read_csv("gradients1.norm_factor_SUMS.csv", col_names = FALSE)

stan <- stan %>% mutate(X1 = str_replace(X1, "_2", ".2"))

dat %>% anti_join(stan, by = c("var" = "X1"))

notOfInt %>% anti_join(stan, by = c("var" = "X1"))

head(stan)
colnames(stan)[2] <- "convFactor"

nrow(dat)
dat <- dat %>% left_join(stan, by = c("var" = "X1"))
nrow(dat)

nrow(notOfInt)
notOfInt <- notOfInt %>% left_join(stan, by = c("var" = "X1"))
nrow(notOfInt)

dat %>% select(est_counts, convFactor) %>% str()
dat <- dat %>% mutate(absoluteCounts = est_counts*convFactor)

notOfInt <- notOfInt %>% mutate(absoluteCounts = est_counts*convFactor)

taxaPred <- read_csv("G1_surface_trophicPredictions_marFerret2023_cleanedUp_noOutliers.csv")

taxaPred %>% anti_join(dat, by = c("tax_name"))

phot <- read_csv("../g1G2G3/phototrophicSpeciesStrict_dielIncubations.csv") %>% mutate(type = "Phototrophic")
het <- read_csv("../g1G2G3/heterotrophicSpeciesStrict_dielIncubations.csv") %>% mutate(type = "Heterotrophic")
mix <- read_csv("../g1G2G3/mixedSpeciesStrict_dielIncubations.csv") %>% mutate(type = "Mixed predictions")

labs <- bind_rows(phot, het, mix)

labs <- labs %>% distinct(taxa, type)

dat %>% anti_join(labs, by = c("tax_name" = "taxa"))

nrow(dat)
dat <- dat %>% left_join(labs %>% select(1:2), by = c("tax_name" = "taxa"))
nrow(dat)

dat %>% group_by(LATITUDE, type) %>%
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

withTrophMode <- both %>% group_by(taxaGroup) %>% summarize(absoluteCounts = sum(absoluteCounts)) %>% 
  filter(taxaGroup == "speciesWithTrophicMode") %>% select(absoluteCounts)

rest <- both %>% group_by(taxaGroup) %>% summarize(absoluteCounts = sum(absoluteCounts)) %>% 
  filter(taxaGroup == "Rest") %>% select(absoluteCounts)

#76.44 28 species
round(100*withTrophMode$absoluteCounts/(withTrophMode$absoluteCounts+rest$absoluteCounts),2)




       