library(tidyverse)
library(readxl)
library(ggpattern)
library(ggpattern)

setwd("~/Dropbox/grad/research/g1G2G3/")

g1 <- read_csv("../g1Surface/absoluteCountsByTrophicMode_full.csv")
g2 <- read_csv("../g2Surface/absoluteCountsByTrophicMode_full.csv")
g3 <- read_csv("../g3/absoluteCountsByTrophicMode_full.csv")

g1 <- g1 %>% mutate(Cruise = "G1 surface") %>% select(Cruise, tax_name, LATITUDE, DEPTH, size, absoluteCounts, REPLICATE)
g2 <- g2 %>% mutate(Cruise = "G2 surface") %>% select(Cruise, tax_name, Latitude, `Depth (m)`, size, absoluteCounts, `Sample ID`)
g3 <- g3 %>% mutate(Cruise = "G3 surface") %>% select(Cruise, tax_name, Latitude, Depth, Filter, absoluteCounts, Replicate)

colnames(g1) <- c("Cruise", "Species", "Latitude", "Depth (m)", "Filter (¬µm)", "Transcripts per L", "Replicate")
colnames(g2) <- c("Cruise", "Species", "Latitude", "Depth (m)", "Filter (¬µm)", "Transcripts per L", "Replicate")
colnames(g3) <- c("Cruise", "Species", "Latitude", "Depth (m)", "Filter (¬µm)", "Transcripts per L", "Replicate")

g1$`Filter (¬µm)` <- parse_number(g1$`Filter (¬µm)`)

g2$`Filter (¬µm)` <- parse_number(g2$`Filter (¬µm)`)

g <- bind_rows(g1, g2, g3)

g %>% distinct(`Filter (¬µm)`)
g <- g %>% mutate(`Filter (¬µm)` = ifelse(`Filter (¬µm)` == 2, 0.2, `Filter (¬µm)`))

aloha <- read_csv("../alohaDiel/speciesWithTrophicModePredictionsAbundance.csv")
meta <- read_excel("../alohaDiel/SCOPE_cruise_sample_log.xlsx")

head(meta)
colnames(meta) <- meta[5,]
meta <- meta[-c(1:5),]

head(meta)
colnames(meta)[20:22] <- c("date", "cast", "time")
meta <- meta %>% select(20:22)

meta <- meta %>% filter(!is.na(date))

meta <- meta %>% mutate(stationCast = str_replace_all(cast, " ", "_"))
meta <- meta %>% mutate(stationCast = str_replace(stationCast, ":", ""))
meta <- meta %>% mutate(stationCast = str_replace(stationCast, "S6", "S06"))
meta <- meta %>% mutate(stationCast = str_replace(stationCast, "S7", "S07"))
meta <- meta %>% mutate(stationCast = str_replace(stationCast, "S8", "S08"))

meta <- meta %>% mutate(Replicate = str_extract(cast, "[ ABC]"))
meta <- meta %>% mutate(Replicate = str_replace(Replicate, " ", ""))

aloha %>% anti_join(meta, by = c("sample" = "stationCast"))

nrow(aloha)
aloha <- aloha %>% left_join(meta %>% select(stationCast, date, Replicate), by = c("sample" = "stationCast"))
nrow(aloha)

aloha <- aloha %>% mutate(Time = str_extract(sample, "[0-9]{1,}$"))
aloha$Time <- as.numeric(aloha$Time)
aloha %>% filter(is.na(Time))

aloha$date <- as.character(aloha$date)

aloha %>% filter(date == "2015-07-26")

aloha <- aloha %>% mutate(Time = ifelse(date == "2015-07-27", Time+2400, Time))
aloha <- aloha %>% mutate(Time = ifelse(date == "2015-07-28", Time+2400+2400, Time))
aloha <- aloha %>% mutate(Time = ifelse(date == "2015-07-29", Time+2400+2400+2400, Time))
aloha <- aloha %>% mutate(Time = ifelse(date == "2015-07-30", Time+2400+2400+2400+2400, Time))

aloha <- aloha %>% mutate(Latitude = 22.75)

aloha <- aloha %>% select(-date)

aloha <- aloha %>% mutate(Cruise = "Aloha diel")

aloha <- aloha %>% mutate(`Filter (¬µm)` = 0.2)

aloha <- aloha %>% mutate(`Depth (m)` = 15)

aloha <- aloha %>% select(Cruise, tax_name, Latitude, Time, `Depth (m)`, `Filter (¬µm)`, absoluteCounts, Replicate)

colnames(aloha)[2]
colnames(aloha)[2] <- "Species"
colnames(aloha)[7]
colnames(aloha)[7] <- "Transcripts per L"

g <- bind_rows(g, aloha)

g3Diel <- read_csv("../g3Diel/speciesWithTrophicModePredictionsAbundance.csv")

g3Diel <- g3Diel %>% mutate(Cruise = "G3 diel")

g3Diel <- g3Diel %>% mutate(`Filter (¬µm)` = 0.2)

meta <- read_csv("~/Dropbox/grad/research/metaT sample log - metaT.csv")

meta <- meta %>% filter(Cruise == "Gradients3" | Cruise == "Gradients 3")

meta <- meta %>% filter(Type == "polyA")

meta <- meta %>% mutate(sample = str_c("S",Station, "C", as.character(Cast), ".", Replicate))

meta <- meta %>% select(sample, Time, Date, Replicate, Latitude)

meta <- meta %>% mutate(Time = str_replace(Time, "\\:", ""))
meta$Time <- as.numeric(meta$Time)

g3Diel %>% anti_join(meta, by = c("sample"))

nrow(g3Diel)
g3Diel <- g3Diel %>% left_join(meta, by = c("sample"))
nrow(g3Diel)

g3Diel <- g3Diel %>% mutate(Time = ifelse(Date == "4/17/19", Time+2400, Time))
g3Diel <- g3Diel %>% mutate(Time = ifelse(Date == "4/18/19", Time+2400+2400, Time))

g3Diel <- g3Diel %>% mutate(`Depth (m)` = 15)
g3Diel <- g3Diel %>% mutate(Cruise = "G3 diel") %>% select(Cruise, tax_name, Latitude, `Depth (m)`, `Filter (¬µm)`, Time, absoluteCounts, Replicate)

colnames(g3Diel)[2]
colnames(g3Diel)[2] <- "Species"
colnames(g3Diel)[7]
colnames(g3Diel)[7] <- "Transcripts per L"

g <- bind_rows(g, g3Diel)

g3Depth <- read_csv("../g3Depth/speciesWithTrophicModePredictionsAbundance.csv")

g3Depth <- g3Depth %>% mutate(Depth = str_extract(`Add'l ID`, "DCM|15|75|125"))

g3Depth$`Filter (¬µm)` <- 0.2

g3Depth <- g3Depth %>% mutate(station = str_extract(`Add'l ID`, "S[0-9]{1,}"))

#1 S8      41   
#2 S6      130  
#3 S4      50   
#4 S5      55   
g3Depth <- g3Depth %>% mutate(Depth = ifelse(Depth == "DCM" & station == "S8", 41, Depth))
g3Depth <- g3Depth %>% mutate(Depth = ifelse(Depth == "DCM" & station == "S6", 130, Depth))
g3Depth <- g3Depth %>% mutate(Depth = ifelse(Depth == "DCM" & station == "S4", 50, Depth))
g3Depth <- g3Depth %>% mutate(Depth = ifelse(Depth == "DCM" & station == "S5", 55, Depth))

meta <- read_csv("~/Dropbox/grad/research/g3Depth/Gradients3_discrete_samples - RNA_cast.csv")

meta <- meta %>% mutate(station = str_extract(`Sample ID`, "S[0-9]{1,}C[0-9]{1,}"))

meta <- meta %>% mutate(Replicate = str_extract(`Sample ID`, " [ABC]"))
meta <- meta %>% mutate(Replicate = str_replace(Replicate, " ", ""))
meta <- meta %>% distinct(station, `Lat (dec)`, Replicate)

g3Depth <- g3Depth %>% mutate(station = str_extract(`Add'l ID`, "S[0-9]{1,}C[0-9]{1,}"))

g3Depth <- g3Depth %>% mutate(Replicate = str_extract(`Add'l ID`, " [ABC]"))
g3Depth <- g3Depth %>% mutate(Replicate = str_replace(Replicate, " ", ""))

g3Depth %>% anti_join(meta, by = c("station", "Replicate"))

nrow(g3Depth)
g3Depth <- g3Depth %>% left_join(meta, by = c("station", "Replicate"))
nrow(g3Depth)

g3Depth <- g3Depth %>% mutate(Cruise = "G3 depth")

g3Depth <- g3Depth %>% select(Cruise, tax_name, `Lat (dec)`, Depth, `Filter (¬µm)`, absoluteCounts, sample, Replicate)

colnames(g3Depth) <- c("Cruise", "Species", "Latitude", "Depth (m)", "Filter (¬µm)", "Transcripts per L", "sample", "Replicate")

g3Depth$`Depth (m)` <- as.numeric(g3Depth$`Depth (m)`)

g3Depth %>% distinct(`Depth (m)`)

g <- bind_rows(g, g3Depth %>% select(-sample))

g <- g %>% select(Cruise, Species, Latitude, `Depth (m)`, Time, `Filter (¬µm)`, Replicate, `Transcripts per L`)

g %>% filter(is.na(Cruise))
g %>% filter(is.na(Species))
g %>% filter(is.na(Latitude))
g %>% filter(is.na(`Depth (m)`)) %>% distinct(Cruise)
g %>% filter(!is.na(Time)) %>% distinct(Cruise)
g %>% filter(is.na(`Filter (¬µm)`))
g %>% filter(is.na(`Transcripts per L`))

g <- g %>% arrange(Cruise, Species, Latitude, `Depth (m)`, Time, `Filter (¬µm)`, Replicate)

g %>% group_by(Cruise) %>% distinct(Species) %>% summarize(n = n())

g2Inc <- read_csv("../g2Incubation/absoluteCountsByIncubation.csv")

g2Inc <- g2Inc %>% mutate(Cruise = "G2 incubations")

g2Inc <- g2Inc %>% mutate(`Depth (m)` = 15)

g2IncMeta <- read_csv("../g2Incubation/G2.RR_exp.metadata.csv")

g2IncMeta <- g2IncMeta %>% distinct(`Sample name`, `Size Fraction`, Expt, Timepoint, Treatment, Rep)

g2Inc %>% anti_join(g2IncMeta, by = c("sample" = "Sample name"))

nrow(g2Inc)
g2Inc <- g2Inc %>% left_join(g2IncMeta, by = c("sample" = "Sample name"))
nrow(g2Inc)

g2Inc <- g2Inc %>% select(-sample)

g2Inc <- g2Inc %>% mutate(latitude = ifelse(Expt == "REXP1", 41.42, NA))
g2Inc <- g2Inc %>% mutate(latitude = ifelse(Expt == "REXP2", 37.00, latitude))
g2Inc <- g2Inc %>% mutate(latitude = ifelse(Expt == "REXP3", 32.93, latitude))

g2Inc <- g2Inc %>% select(-Expt)

g2Inc <- g2Inc %>% mutate(`Size Fraction` = parse_number(`Size Fraction`))

colnames(g2Inc)
colnames(g2Inc) <- c("Species", "Transcripts per L", "Cruise", "Depth (m)", "Filter (¬µm)", "Timepoint (hr)", "Treatment", "Replicate", "Latitude")

g <- bind_rows(g, g2Inc)

g <- g %>% arrange(Cruise)

g %>% distinct(Species)

dinos <- c("Azadinium spinosum", "Brandtodinium nutricula", "Dinophysis acuminata", "Gymnodinium catenatum GC744", 
           "Karenia brevis", "Karlodinium veneficum", "Pelagodinium beii", "Prorocentrum minimum", "Scrippsiella trochoidea", 
           "Tripos fusus")

g <- g %>% mutate(`Transcripts per L with dinoflagellate correction` = ifelse(Species %in% dinos, `Transcripts per L`/6.4, `Transcripts per L`))

g %>% 
  mutate(Cruise = str_replace(Cruise, "Aloha", "ALOHA")) %>% 
  mutate(Cruise = str_replace(Cruise, "G", "Gradients")) %>% 
  select(Cruise, Species, Latitude, `Depth (m)`, Time, Treatment, `Timepoint (hr)`, 6, Replicate, 8, `Transcripts per L with dinoflagellate correction`) %>%
  write_csv("~/Dropbox/grad/research/speciesWithTrophicModePredictionsAbundance_supplementaryTable.csv")


###

g1 <- read_csv("../g1Surface/absoluteCountsByTrophicMode_full.csv")
g2 <- read_csv("../g2Surface/absoluteCountsByTrophicMode_full.csv")
g3 <- read_csv("../g3/absoluteCountsByTrophicMode_full.csv")

g1 <- g1 %>% select(tax_name, LATITUDE, type, absoluteCounts, var, REPLICATE)
g2 <- g2 %>% select(tax_name, Latitude, type, absoluteCounts, var, `Sample ID`)
g3 <- g3 %>% select(tax_name, Latitude, type, absoluteCounts, sample, Replicate)

colnames(g1)[2]
colnames(g1)[2] <- "Latitude"

colnames(g1)[5]
colnames(g1)[5] <- "sample"

colnames(g2)[5]
colnames(g2)[5] <- "sample"

colnames(g1)[6]
colnames(g1)[6] <- "Replicate"

colnames(g2)[6]
colnames(g2)[6] <- "Replicate"

g <- bind_rows(g1 %>% mutate(cruise = "Gradients1: 2016"), g2 %>% mutate(cruise = "Gradients2: 2017"), g3 %>% mutate(cruise = "Gradients3: 2019"))

g %>% distinct(Replicate)

g3Depth_surf <- g3Depth %>% filter(`Depth (m)` == 15)

g3Depth_surf <- g3Depth_surf %>% mutate(size = "0.2")

g3Depth_surf <- g3Depth_surf %>% select(Species, Latitude, `Transcripts per L`, sample, Replicate, size)

colnames(g3Depth_surf)
colnames(g3Depth_surf) <- c("tax_name", "Latitude", "absoluteCounts", "sample", "Replicate", "size")

#from plotAllTrophicModePredictions.R
mixedPreds <- read_csv("~/Dropbox/grad/research/g1G2G3/mixedSpeciesStrict_dielIncubations.csv")
photPreds <- read_csv("~/Dropbox/grad/research/g1G2G3/phototrophicSpeciesStrict_dielIncubations.csv")
hetPreds <- read_csv("~/Dropbox/grad/research/g1G2G3/heterotrophicSpeciesStrict_dielIncubations.csv")

types <- bind_rows(mixedPreds %>% mutate(type = "Mixed predictions"), photPreds %>% mutate(type = "Phototrophic"), 
                   hetPreds %>% mutate(type = "Heterotrophic"))

nrow(g3Depth_surf)
g3Depth_surf <- g3Depth_surf %>% left_join(types, by = c("tax_name" = "taxa"))
nrow(g3Depth_surf)

g3Depth_surf %>% filter(is.na(type))

g3Depth_surf <- g3Depth_surf %>% mutate(cruise = "Gradients3: 2019", cruise2 = "g3Depth")

g3Depth_surf$sample <- as.character(g3Depth_surf$sample)

merged <- g %>% 
  bind_rows(g3Depth_surf)

merged <- merged %>% mutate(size = ifelse(is.na(size), str_extract(sample, "0.2|3um"), size))
merged %>% filter(is.na(size)) %>% distinct(cruise, sample)

sample <- read_csv("../metaT sample log - metaT.csv")
sample <- sample %>% filter(Cruise == "Gradients3" | Cruise == "Gradients 3")
sample <- sample %>% filter(Type == "polyA")

sample <- sample %>% filter(str_detect(Sample.ID, "UW")) %>% 
  mutate(Sample.ID = str_c("G3PA.", Sample.ID))

sample <- sample %>% select(Sample.ID, Filter)

nrow(merged)
merged <- merged %>% left_join(sample, by = c("sample" = "Sample.ID"))
nrow(merged)

merged <- merged %>% mutate(size = ifelse(is.na(size), Filter, size)) %>% select(-Filter)

merged %>% distinct(size)
merged <- merged %>% mutate(size = ifelse(size == "0_2", "0.2", size))
merged <- merged %>% mutate(size = ifelse(size == "3um", "3", size))
merged %>% distinct(size)

merged <- merged %>% select(-sample) 

merged %>% group_by(tax_name, Latitude, type, Replicate, cruise, cruise2) %>% summarize(n = n()) %>% 
  filter(n > 2)

merged <- merged %>% spread(key = size, value = absoluteCounts, fill = 0)

merged <- merged %>% mutate(absoluteCounts = `0.2`+`3`)


nutr_g1 <- read_csv("~/Dropbox/grad/research/g1Surface/Cruise KOK1606 Gradients 1 Organic and Inorganic Nutrients/Cruise KOK1606 Gradients 1 Organic and Inorganic Nutrients.csv")
nutr_g1 <- nutr_g1 %>% filter(depth < 24) %>% dplyr::select(lat, lon, SiO4:NH4)

nutr_g1 <- nutr_g1 %>% gather(SiO4:NH4, key = "var", value = "value")
nutr_g1 <- nutr_g1 %>% distinct()

nutr_g2 <- read_csv("~/Dropbox/grad/research/mixotrophyIFCB/Cruise MGL1704 Gradients 2 Organic and Inorganic Nutrients.csv")
nutr_g2 <- nutr_g2 %>% filter(depth < 24) %>% dplyr::select(lat, lon, SiO4:TOC)
nutr_g2 <- nutr_g2 %>% gather(SiO4:TOC, key = "var", value = "value")
nutr_g2 <- nutr_g2 %>% distinct() 

nutr_g3 <- read_csv("~/Dropbox/grad/research/g3/Gradients 3 KM1906 Organic and Inorganic Nutrients/Gradients 3 KM1906 Organic and Inorganic Nutrients.csv")
nutr_g3 <- nutr_g3 %>% filter(depth < 24) %>% dplyr::select(lat, lon, SiO4:NO2)
nutr_g3 <- nutr_g3 %>% gather(SiO4:NO2, key = "var", value = "value")
nutr_g3 <- nutr_g3 %>% distinct() 

nut <- bind_rows(nutr_g1 %>% mutate(cruise = "Gradients1: 2016"), nutr_g2 %>% mutate(cruise = "Gradients2: 2017"), nutr_g3 %>% mutate(cruise = "Gradients3: 2019"))

nut <- nut %>% filter(var %in% c("NO3_NO2", "NO3_plus_NO2"))
nut <- nut %>% mutate(var = "NO3_NO2")
nut <- nut %>% filter(!is.na(value))

#load modis satellite light data
g1_modis <- read_excel("../g1_PAR.xlsx")
g2_modis <- read_excel("../g2_PAR.xlsx")
g3_modis <- read_excel("../g3_PAR.xlsx")

#modis days are one behind what they should be because of time format
#they are behind by 7 hours
modis <- bind_rows(g1_modis %>% mutate(cruise = "Gradients1: 2016"), g2_modis %>% mutate(cruise = "Gradients2: 2017"), g3_modis %>% mutate(cruise = "Gradients3: 2019"))

modis <- modis %>% filter(!is.na(PAR))

head(modis)

#fixes time by adding 7 hours 
modis <- modis %>% mutate(time = time + hm("7:00"))

latDay <- read_csv("~/Dropbox/grad/research/trophicModePredictionsSupplementaryTable.csv")
latDay <- latDay %>% select(Cruise, Date, Time, Latitude, `Depth (m)`)

latDay <- latDay %>% filter(Cruise %in% c("Gradients1", "Gradients2", "Gradients3", "Gradients3 depth"))

latDay <- latDay %>% filter(!(Cruise == "Gradients3 depth" & `Depth (m)` != 15))

latDay <- latDay %>% mutate(Cruise = ifelse(Cruise == "Gradients1", "Gradients1: 2016", Cruise))
latDay <- latDay %>% mutate(Cruise = ifelse(Cruise == "Gradients2", "Gradients2: 2017", Cruise))
latDay <- latDay %>% mutate(Cruise = ifelse(Cruise == "Gradients3" | Cruise == "Gradients3 depth", "Gradients3: 2019", Cruise))

modis <- modis %>% mutate(Date = format(time, "%m-%d-%Y"))

latDay$Date <- strptime(latDay$Date, format="%m-%d-%y")
latDay$Date <- strftime(latDay$Date, format="%m-%d-%Y")

colnames(latDay)[1]
colnames(latDay)[1] <- "cruise"

latDay %>% anti_join(modis, by = c("Date", "cruise")) %>% distinct(Date, cruise)

modis %>% group_by(cruise) %>% summarize(min = min(Date), max = max(Date))
latDay %>% group_by(cruise) %>% summarize(min = min(Date), max = max(Date))

modis <- modis %>% filter(Date != "04-19-2016")
modis <- modis %>% filter(Date != "05-28-2017")
modis <- modis %>% filter(Date != "04-10-2019")

modis %>% group_by(cruise) %>% summarize(min = min(Date), max = max(Date))
latDay %>% group_by(cruise) %>% summarize(min = min(Date), max = max(Date))

merged %>% group_by(cruise) %>% distinct(tax_name) %>% summarize(n = n())

merged %>% distinct(cruise)
merged %>% distinct(type)

merged <- merged %>% mutate(absoluteCounts_dinoCorr = ifelse(tax_name %in% dinos, absoluteCounts/6.4, absoluteCounts))

toPlot <- merged %>% 
  group_by(cruise, Latitude, type, tax_name, cruise2) %>% 
  summarize(absoluteCounts_dinoCorr = mean(absoluteCounts_dinoCorr), n = n()) %>% 
  mutate(tax_name = ifelse(!(tax_name %in% c("Pelagomonas calceolata", "Oxytricha trifallax", "Triparma pacifica", "Pelagodinium beii", "Karlodinium veneficum", "Chloroparvula japonica")), "Other species with trophic predictions", tax_name)) %>%
  mutate(tax_name = factor(tax_name, levels = c("Oxytricha trifallax", "Karlodinium veneficum", "Pelagodinium beii", "Chloroparvula japonica", "Pelagomonas calceolata", "Triparma pacifica", "Other species with trophic predictions"))) %>%
  mutate(type = ifelse(type == "Heterotrophic", "Species with only\nheterotrophy predictions", type)) %>%
  mutate(type = ifelse(type == "Mixed predictions", "Species with mixed\ntrophic predictions", type)) %>%
  mutate(type = ifelse(type == "Phototrophic", "Species with only\nphototrophy predictions", type)) %>%
  mutate(type = factor(type, levels = c("Species with only\nheterotrophy predictions", "Species with mixed\ntrophic predictions", "Species with only\nphototrophy predictions")))

toPlot %>%
  ggplot(aes(x = Latitude, y = absoluteCounts_dinoCorr/1E9, color = type)) +
  geom_jitter(data = toPlot %>% 
                filter(tax_name == "Other species with trophic predictions"), width = .5, alpha = .2, size = 2) +
  geom_jitter(data = toPlot %>% 
                filter(tax_name != "Other species with trophic predictions"), aes(shape = tax_name), width = .5, size = 3.5, stroke = 1) +
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
  labs(y = "Billion transcripts per liter", x = "Latitude (°N)", shape = "") + 
  geom_vline(data=filter(g, cruise=="Gradients1: 2016"), aes(xintercept=32.15), colour="red", linetype="dashed") + 
  geom_vline(data=filter(g, cruise=="Gradients1: 2016"), aes(xintercept=33), colour="black", linetype="dashed") + 
  
  geom_vline(data=filter(g %>% mutate(cruise = "Gradients2: 2017"), cruise=="Gradients2: 2017"), aes(xintercept=32.5), colour="red", linetype="dashed") + 
  geom_vline(data=filter(g %>% mutate(cruise = "Gradients2: 2017"), cruise=="Gradients2: 2017"), aes(xintercept=36.2), colour="black", linetype="dashed") + 
  
  geom_vline(data=filter(g %>% mutate(cruise = "Gradients3: 2019"), cruise=="Gradients3: 2019"), aes(xintercept=32.45), colour="red", linetype="dashed") + 
  geom_vline(data=filter(g %>% mutate(cruise = "Gradients3: 2019"), cruise=="Gradients3: 2019"), aes(xintercept=35), colour="black", linetype="dashed") + 
  facet_wrap(~cruise) + labs(color = "") + 
  scale_shape_manual(values=c(4,0,2,3,5,6)) + 
  scale_color_manual(values = c("orange", "purple", "forestgreen")) +
  guides(color = FALSE)  + 
  scale_x_continuous(breaks = c(25,30,35,40), limits = c(23, 44))

ggsave("abundanceOverG1G2G3surface_byTrophicMode_taxa_avgAcrossReplicates.png", height = 7, width = 20)


merged %>%
  ungroup() %>%
  group_by(cruise, cruise2, Latitude, type, Replicate) %>%
  summarize(absoluteCounts_dinoCorr = sum(absoluteCounts_dinoCorr), n = n()) %>%
  ungroup() %>% 
  group_by(cruise, Latitude, type) %>%
  summarize(absoluteCounts_dinoCorr = mean(absoluteCounts_dinoCorr), n = n()) %>%
  ungroup() %>%
  group_by(Latitude, cruise) %>%
  arrange(desc(absoluteCounts_dinoCorr)) %>%
  slice(1) %>%
  ggplot(aes(x = Latitude, color = type)) +
  geom_point(y = 1, size = 6) +
  geom_point(y = 1,pch=21, color = "black", size = 6) + 
  facet_wrap(~cruise) +
  scale_color_manual(values = c("orange", "purple", "forestgreen")) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text.x = element_text(size = 26, color = 'black'),
    strip.text.y = element_text(size = 26, color = 'black'),
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.line.x.bottom = element_blank(),   # Remove bottom x-axis line
    axis.line.x.top = element_line(size = 0.5, color = "black"),  # Add top border
    axis.line.y = element_line(size = 0.5, color = "black"),  # Add left border
    axis.text.y = element_text(size = 22, color = 'black'),
    axis.title.y = element_text(size = 26, color = 'black'),
    legend.position = "none",  # Remove legend
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(y = "", x = "", shape = "", fill = "") + 
  scale_x_continuous(breaks = c(25,30,35,40), limits = c(23, 44))

ggsave("dominantTrophicModeByAbundance.png", height = 4, width = 20)


merged %>%
  ungroup() %>%
  group_by(cruise, cruise2, Latitude, type, Replicate) %>%
  summarize(absoluteCounts_dinoCorr = sum(absoluteCounts_dinoCorr)) %>%
  group_by(cruise, Latitude, type) %>%
  summarize(absoluteCounts_dinoCorr = mean(absoluteCounts_dinoCorr)) %>%
  ungroup() %>%
  group_by(Latitude, cruise) %>%
  arrange(desc(absoluteCounts_dinoCorr)) %>%
  slice(2) %>%
  ggplot(aes(x = Latitude, color = type)) +
  geom_point(y = 1, size = 6) +
  geom_point(y = 1,pch=21, color = "black", size = 6) + 
  facet_wrap(~cruise) +
  scale_color_manual(values = c("orange", "purple", "forestgreen")) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text.x = element_text(size = 26, color = 'black'),
    strip.text.y = element_text(size = 26, color = 'black'),
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.line.x.bottom = element_blank(),   # Remove bottom x-axis line
    axis.line.x.top = element_line(size = 0.5, color = "black"),  # Add top border
    axis.line.y = element_line(size = 0.5, color = "black"),  # Add left border
    axis.text.y = element_text(size = 22, color = 'black'),
    axis.title.y = element_text(size = 26, color = 'black'),
    legend.position = "none",  # Remove legend
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(y = "", x = "", shape = "", fill = "") + 
  scale_x_continuous(breaks = c(25,30,35,40), limits = c(23, 44))

ggsave("secondDominantTrophicModeByAbundance.png", height = 4, width = 20)


merged %>% 
  mutate(type = ifelse(type == "Heterotrophic", "Heterotrophic species bins", type)) %>% 
  mutate(type = ifelse(type == "Phototrophic", "Phototrophic species bins", type)) %>% 
  mutate(type = ifelse(type == "Mixed predictions", "Species bins with\nmixotrophic capabilities", type)) %>%
  mutate(type = factor(type, levels = c("Heterotrophic species bins", "Species bins with\nmixotrophic capabilities", "Phototrophic species bins"))) %>% 
  ungroup() %>% 
  group_by(cruise, cruise2, Latitude, type, Replicate) %>% 
  summarize(absoluteCounts_dinoCorr = sum(absoluteCounts_dinoCorr), n = n()) %>%
  group_by(cruise, Latitude, type) %>% 
  summarize(absoluteCounts_dinoCorr = mean(absoluteCounts_dinoCorr), n = n()) %>%
  ggplot(aes(x = Latitude, y = absoluteCounts_dinoCorr/1E9, fill = type)) +
  geom_area() +
  scale_fill_manual(values = c("orange", "purple", "forestgreen")) +
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
  labs(y = "Billion transcripts per liter", x = "Latitude (°N)", shape = "", fill = "") + 
  geom_vline(data=filter(g, cruise=="Gradients1: 2016"), aes(xintercept=32.15), colour="red", linetype="dashed") + 
  geom_vline(data=filter(g, cruise=="Gradients1: 2016"), aes(xintercept=33), colour="black", linetype="dashed") + 
  
  geom_vline(data=filter(g %>% mutate(cruise = "Gradients2: 2017"), cruise=="Gradients2: 2017"), aes(xintercept=32.5), colour="red", linetype="dashed") + 
  geom_vline(data=filter(g %>% mutate(cruise = "Gradients2: 2017"), cruise=="Gradients2: 2017"), aes(xintercept=36.2), colour="black", linetype="dashed") + 
  
  geom_vline(data=filter(g %>% mutate(cruise = "Gradients3: 2019"), cruise=="Gradients3: 2019"), aes(xintercept=32.45), colour="red", linetype="dashed") + 
  geom_vline(data=filter(g %>% mutate(cruise = "Gradients3: 2019"), cruise=="Gradients3: 2019"), aes(xintercept=35), colour="black", linetype="dashed") + 
  facet_wrap(~cruise) + labs(color = "") +
  scale_x_continuous(breaks = c(25,30,35,40), limits = c(23, 44)) + 
  geom_point(data = nut %>% mutate(type = "Phototrophic species bins") %>% filter(lat > 23.49) %>% filter(lat < 42.335), aes(x = lat, y = value*5), size = .7, color = 'blue') +
  geom_line(data = nut %>% mutate(type = "Phototrophic species bins") %>% filter(lat > 23.49) %>% filter(lat < 42.335), aes(x = lat, y = value*5), color = 'blue') +
  
  geom_point(data = modis %>% mutate(type = "Phototrophic species bins") %>% group_by(cruise, type, lat) %>% summarize(PAR = mean(PAR)) %>% filter(lat < 42.335), aes(x = lat, y = PAR/2), size = .7, color = 'red') +
  geom_line(data = modis %>% mutate(type = "Phototrophic species bins") %>% group_by(cruise, type, lat) %>% summarize(PAR = mean(PAR)) %>% filter(lat < 42.335), aes(x = lat, y = PAR/2), color = 'red') +
  
  scale_y_continuous(limits = c(0, 40), name = "Billion transcripts per liter",
                     sec.axis = sec_axis(trans=~./5, name = expression(NO[3] * "_" * NO[2]~(μmol/L)))) + theme(axis.text.y.right = element_text(colour = "blue"), axis.title.y.right = element_text(colour = "blue"))

ggsave("abundanceOverG1G2G3surface_byTrophicMode_sum_area.png", height = 7, width = 20)


ggplot() + 
  geom_point(data = modis %>% mutate(type = "Phototrophic species bins") %>% group_by(cruise, type, lat) %>% summarize(PAR = mean(PAR)), aes(x = lat, y = PAR/2), size = .7, color = 'red') +
  geom_line(data = modis %>% mutate(type = "Phototrophic species bins") %>% group_by(cruise, type, lat) %>% summarize(PAR = mean(PAR)), aes(x = lat, y = PAR/2), color = 'red') +
  geom_point(data = nut %>% mutate(type = "Phototrophic species bins"), aes(x = lat, y = value*5), size = .7, color = 'blue') +
  geom_line(data = nut %>% mutate(type = "Phototrophic species bins"), aes(x = lat, y = value*5), color = 'blue') + 
  facet_wrap(~cruise) + labs(color = "") +
  scale_x_continuous(breaks = c(25,30,35,40)) + theme_classic() +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size = 26, color = 'black')) + 
  theme(strip.text.y = element_text(size = 26, color = 'black')) + 
  theme(axis.text.x = element_text(size = 22, color = 'black'))  + 
  theme(axis.text.y = element_text(size = 22, color = 'blue')) + 
  theme(axis.title.y = element_text(size = 20, color = 'blue')) + 
  theme(legend.text = element_text(size = 20, color = 'black')) +
  theme(legend.title = element_text(size = 22, color = 'black')) +
  theme(axis.title.x = element_text(size = 26, color = 'black')) +  
  labs(y = expression(NO[3] * "_" * NO[2]~(μmol/L)), x = "", shape = "", fill = "") + 
  scale_y_continuous(limits = c(0, 40), sec.axis = sec_axis(trans=~.*2, name=expression("PAR (" * "mol quanta " * m^-2 * " " *day^-1 * ")"))) + theme(axis.text.y.right = element_text(colour = "red"), axis.title.y.right = element_text(size = 20, colour = "red"))

ggsave("nitPARagainstLat.png", height = 3.7, width = 20)

merged %>% 
  mutate(type = ifelse(type == "Heterotrophic", "Heterotrophic species bins", type)) %>% 
  mutate(type = ifelse(type == "Phototrophic", "Phototrophic species bins", type)) %>% 
  mutate(type = ifelse(type == "Mixed predictions", "Species bins with\nmixotrophic capabilities", type)) %>%
  mutate(type = factor(type, levels = c("Heterotrophic species bins", "Species bins with\nmixotrophic capabilities", "Phototrophic species bins"))) %>% 
  ungroup() %>% 
  group_by(cruise, cruise2, Latitude, type, Replicate) %>% 
  summarize(absoluteCounts_dinoCorr = sum(absoluteCounts_dinoCorr)) %>%
  group_by(cruise, Latitude, type) %>% 
  summarize(absoluteCounts_dinoCorr = mean(absoluteCounts_dinoCorr)) %>%
  ggplot(aes(x = Latitude, y = absoluteCounts_dinoCorr/1E9, fill = type)) +
  geom_area() +
  scale_fill_manual(values = c("orange", "purple", "forestgreen")) +
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
  labs(y = "Billion transcripts per liter", x = "Latitude (°N)", shape = "", fill = "") + 
  geom_vline(data=filter(g, cruise=="Gradients1: 2016"), aes(xintercept=32.15), colour="red", linetype="dashed") + 
  geom_vline(data=filter(g, cruise=="Gradients1: 2016"), aes(xintercept=33), colour="black", linetype="dashed") + 
  
  geom_vline(data=filter(g %>% mutate(cruise = "Gradients2: 2017"), cruise=="Gradients2: 2017"), aes(xintercept=32.5), colour="red", linetype="dashed") + 
  geom_vline(data=filter(g %>% mutate(cruise = "Gradients2: 2017"), cruise=="Gradients2: 2017"), aes(xintercept=36.2), colour="black", linetype="dashed") + 
  
  geom_vline(data=filter(g %>% mutate(cruise = "Gradients3: 2019"), cruise=="Gradients3: 2019"), aes(xintercept=32.45), colour="red", linetype="dashed") + 
  geom_vline(data=filter(g %>% mutate(cruise = "Gradients3: 2019"), cruise=="Gradients3: 2019"), aes(xintercept=35), colour="black", linetype="dashed") + 
  facet_wrap(~cruise) + labs(color = "") +
  scale_x_continuous(breaks = c(25,30,35,40)) + 
  geom_point(data = nut %>% mutate(type = "Phototrophic species bins"), aes(x = lat, y = value*5), size = .7, color = 'blue') +
  geom_line(data = nut %>% mutate(type = "Phototrophic species bins"), aes(x = lat, y = value*5), color = 'blue') +
  
  geom_point(data = modis %>% mutate(type = "Phototrophic species bins") %>% group_by(cruise, type, lat) %>% summarize(PAR = mean(PAR)), aes(x = lat, y = PAR/2), size = .7, color = 'red') +
  geom_line(data = modis %>% mutate(type = "Phototrophic species bins") %>% group_by(cruise, type, lat) %>% summarize(PAR = mean(PAR)), aes(x = lat, y = PAR/2), color = 'red') +
  
  
  scale_y_continuous(limits = c(0, 40), name = "Billion transcripts per liter",
                     sec.axis = sec_axis(trans=~.*2, name=expression("PAR (" * "mol quanta " * m^-2 * " " *day^-1 * ")"))) + theme(axis.text.y.right = element_text(colour = "red"), axis.title.y.right = element_text(colour = "red"))

ggsave("abundanceOverG1G2G3surface_byTrophicMode_sum_area_par.png", height = 7, width = 20)



