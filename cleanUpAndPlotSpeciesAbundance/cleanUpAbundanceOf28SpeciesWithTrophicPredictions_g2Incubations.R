library(tidyverse)

setwd("~/Dropbox/grad/research/g2Incubation//")

dat <- read_csv("g2Inc_ofInt.csv")

notOfInt <- read_csv("g2Inc_notOfInt.csv")

stan <- read_csv("G2PA_REXP_norm_factors.csv")

dat <- dat %>% mutate(sample = str_replace(sample, "G2PA.", ""))
dat <- dat %>% mutate(sample = str_replace(sample, ".abundance.tsv", ""))

dat %>% anti_join(stan, by = c("sample" = "sample_id"))

nrow(dat)
dat <- dat %>% left_join(stan %>% select(sample_id, NORM_FACTOR), by = c("sample" = "sample_id"))
nrow(dat)

notOfInt <- notOfInt %>% mutate(sample = str_replace(sample, "G2PA.", ""))
notOfInt <- notOfInt %>% mutate(sample = str_replace(sample, ".abundance.tsv", ""))

notOfInt %>% anti_join(stan, by = c("sample" = "sample_id"))

nrow(notOfInt)
notOfInt <- notOfInt %>% left_join(stan %>% select(sample_id, NORM_FACTOR), by = c("sample" = "sample_id"))
nrow(notOfInt)

str(dat)
str(notOfInt)

dat <- dat %>% mutate(absoluteCounts = est_counts*NORM_FACTOR)
notOfInt <- notOfInt %>% mutate(absoluteCounts = est_counts*NORM_FACTOR)

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

both %>% filter(taxaGroup == "speciesWithTrophicMode") %>% select(tax_name, sample, absoluteCounts) %>% 
  write_csv("absoluteCountsByIncubation.csv")

withTrophMode <- both %>% group_by(taxaGroup) %>% summarize(absoluteCounts = sum(absoluteCounts)) %>% 
  filter(taxaGroup == "speciesWithTrophicMode") %>% select(absoluteCounts)

rest <- both %>% group_by(taxaGroup) %>% summarize(absoluteCounts = sum(absoluteCounts)) %>% 
  filter(taxaGroup == "Rest") %>% select(absoluteCounts)

#56.82 28 species 
round(100*withTrophMode$absoluteCounts/(withTrophMode$absoluteCounts+rest$absoluteCounts),2)



       