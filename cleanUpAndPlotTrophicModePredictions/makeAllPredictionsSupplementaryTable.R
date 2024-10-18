library(tidyverse)
library(readxl)

setwd("~/Dropbox/grad/research/")

#from plotAllTrophicModePredictions.R
all <- read_csv("g1G2G3/allSpeciesTrophicModePredictions_depthDielIncubations.csv")

all <- all %>% select(sample:probability, size, depth, time, date, station, cruise, group)

##g3 diel
g3Diel <- read_csv("g3Diel/G3_trophicPredictions_cleanedUp_updatedMarferret_marmicroDb2023_diel_noOutliers.csv")

nrow(all)
all <- all %>% left_join(g3Diel %>% distinct(Sample.ID, Latitude, Longitude), by = c("sample" = "Sample.ID"))
nrow(all)

all %>% filter(!is.na(Latitude)) %>% distinct(cruise)

##g1 and g2

g1 <- read_csv("~/Dropbox/grad/research/g1Surface/ctd_g1_g2_15m_source_Stephen.csv")

g1 <- g1 %>% select(SAMPLE_ID, LATITUDE, LONGITUDE, DATETIME)

head(g1)
g1 <- g1 %>% mutate(date = str_extract(DATETIME, "[0-9\\/]{1,}"))
g1 <- g1 %>% mutate(date = str_replace_all(date, "\\/", "-"))
g1 %>% filter(is.na(date))

g1 <- g1 %>% mutate(time = str_extract(DATETIME, "[0-9]{1,}\\:[0-9]{1,}"))
g1 <- g1 %>% mutate(time = str_replace(time, "\\:", ""))
g1$time <- as.numeric(g1$time)
g1 %>% distinct(time)

g1 <- g1 %>% mutate(time = time-1000)

g1 %>% filter(time < 0)

g1 <- g1 %>% mutate(date = ifelse(time < 0 & date == "4-29-16", "4-28-16", date))
g1 <- g1 %>% mutate(date = ifelse(time < 0 & date == "6-10-17", "6-9-17", date))

g1 <- g1 %>% mutate(time = ifelse(time < 0, time+2400, time))

g1 <- g1 %>% select(-DATETIME)

colnames(g1)
colnames(g1) <- c("sample", "latitude", "longitude", "date", "time")

g1 <- g1 %>% mutate(sample = str_replace(sample, "0_2", "0.2"))

all %>% filter(cruise == "g1") %>% anti_join(g1, by = c("sample"))

head(g1)
g1 <- g1 %>% filter(!str_detect(date, "-17$"))

nrow(all)
all <- all %>% left_join(g1, by = c("sample"))
nrow(all)

all %>% filter(!is.na(latitude)) %>% distinct(cruise)
all <- all %>% mutate(time.x = ifelse(!is.na(time.y), time.y, time.x))
all <- all %>% mutate(date.x = ifelse(!is.na(date.y), date.y, date.x))
all <- all %>% select(-date.y, -time.y)

colnames(all)[7:8]
colnames(all)[7:8] <- c("time", "date")

##g3 

#downloaded from https://docs.google.com/spreadsheets/d/1lKwvlNjUV8dXHKO5HkMNa2m70QQYydrTf-Fb6PRQ9oE/edit#gid=0
g3 <- read_csv("~/Dropbox/grad/research/metaT sample log - metaT.csv")

g3 <- g3 %>% filter(Cruise == "Gradients3" | Cruise == "Gradients 3")

g3 <- g3 %>% filter(Type == "polyA")

g3 <- g3 %>% mutate(sample = str_c("G3PA.", Sample.ID))

all %>% filter(cruise == "g3") %>% anti_join(g3, by = c("sample"))

g3 <- g3 %>% select(sample, Date, Time, Latitude, Longitude)

head(g3)
g3 <- g3 %>% mutate(Time = str_replace(Time, "\\:", ""))
g3$Time <- as.numeric(g3$Time)
g3 %>% distinct(Time)

all <- all %>% mutate(latitude = ifelse(is.na(latitude), Latitude, latitude))
all <- all %>% mutate(longitude = ifelse(is.na(longitude), Longitude, longitude))

all <- all %>% select(-Latitude, -Longitude)

nrow(all)
all <- all %>% left_join(g3, by = c("sample"))
nrow(all)

all <- all %>% mutate(latitude = ifelse(!is.na(Latitude), Latitude, latitude))
all <- all %>% mutate(longitude = ifelse(!is.na(Longitude), Longitude, longitude))
all <- all %>% mutate(date = ifelse(!is.na(Date), Date, date))
all <- all %>% mutate(time = ifelse(!is.na(Time), Time, time))

all <- all %>% select(-Latitude, - Longitude, -Date, -Time)

#g3 depth 
g3D <- read_csv("~/Dropbox/grad/research/g3Depth/Gradients3_discrete_samples - RNA_cast.csv")

g3D <- g3D %>% filter(Type == "DEPTH")

head(g3D)

g3D <- g3D %>% mutate(stationCast = str_extract(`Sample ID`, "S[0-9]{1,}C[0-9]{1,}"))
g3D <- g3D %>% mutate(Depth = str_extract(`Sample ID`, " [0-9mDCM]{1,}"))
g3D <- g3D %>% mutate(Depth = str_replace(Depth, " ", ""))
g3D <- g3D %>% mutate(replicate = str_extract(`Sample ID`, " [AB]"))
g3D <- g3D %>% mutate(replicate = str_replace(replicate, " ", ""))

g3D <- g3D %>% mutate(sample = str_c("G3PA.depth.", stationCast, ".", Depth, ".", replicate))
g3D %>% distinct(sample)
g3D %>% filter(is.na(sample))

g3D <- g3D %>% select(sample, Date, `Time Sampled`, `Lat (dec)`, `Lon (dec)`)

all %>% filter(cruise == "g3Depth") %>% anti_join(g3D, by = c("sample"))

g3D$`Time Sampled` <- as.character(g3D$`Time Sampled`)
g3D <- g3D %>% mutate(`Time Sampled` = str_replace(`Time Sampled`, "\\:00$", ""))
g3D <- g3D %>% mutate(`Time Sampled` = str_replace(`Time Sampled`, "\\:", ""))
g3D$`Time Sampled` <- as.numeric(g3D$`Time Sampled`)

g3D <- g3D %>% mutate(Date = str_replace_all(Date, "\\/", "-"))

nrow(all)
all <- all %>% left_join(g3D, by = c("sample"))
nrow(all)

all <- all %>% mutate(latitude = ifelse(!is.na(`Lat (dec)`), `Lat (dec)`, latitude))
all <- all %>% mutate(longitude = ifelse(!is.na(`Lon (dec)`), `Lon (dec)`, longitude))
all <- all %>% mutate(time = ifelse(!is.na(`Time Sampled`), `Time Sampled`, time))
all <- all %>% mutate(date = ifelse(!is.na(Date), Date, date))

colnames(all)
all <- all %>% select(-c(Date:`Lon (dec)`))

#g2 
g2 <- read_csv("~/Dropbox/grad/research/metaT sample log - metaT.csv")

g2 <- g2 %>% filter(Cruise == "Gradients2")

g2 <- g2 %>% mutate(sample = str_c("G2PA.S", as.character(Station), "C", as.character(Cast), ".15m.", as.character(Filter), "um.", Replicate))

g2 <- g2 %>% mutate(sample = str_replace(sample, "S2C", "S02C"))
g2 <- g2 %>% mutate(sample = str_replace(sample, "S5C", "S05C"))
g2 <- g2 %>% mutate(sample = str_replace(sample, "S6C", "S06C"))
g2 <- g2 %>% mutate(sample = str_replace(sample, "S7C", "S07C"))
g2 <- g2 %>% mutate(sample = str_replace(sample, "S9C", "S09C"))
g2 <- g2 %>% mutate(sample = str_replace(sample, "0.2um", "0_2um"))

all %>% filter(cruise == "g2") %>% anti_join(g2, by = c("sample")) %>% select(sample)

g2 <- g2 %>% mutate(Date = str_extract(Date.Time, "[0-9\\/]{1,}"))
g2 <- g2 %>% mutate(Date = str_replace_all(Date, "\\/", "-"))

g2 <- g2 %>% mutate(Time = str_extract(Date.Time, "[0-9]{1,}:[0-9]{1,}$"))
g2 <- g2 %>% mutate(Time = str_replace(Time, ":", ""))
g2$Time <- parse_number(g2$Time)

g2 <- g2 %>% select(sample, Date, Time, Latitude, Longitude)

all %>% semi_join(g2, by = c("sample")) %>% distinct(cruise)

g2 <- g2 %>% distinct()

g2 <- g2 %>% mutate(Time = Time-1000)

nrow(all)
all <- all %>% left_join(g2, by = c("sample"))
nrow(all)

all %>% filter(Time < 0) %>% distinct(Date, Time)
all <- all %>% mutate(Date = ifelse(Time < 0, "6-9-17", Date))
all <- all %>% mutate(Time = ifelse(Time < 0, Time+2400, Time))

all <- all %>% mutate(latitude = ifelse(cruise == "g2", Latitude, latitude))
all <- all %>% mutate(longitude = ifelse(cruise == "g2", Longitude, longitude))
all <- all %>% mutate(date = ifelse(cruise == "g2", Date, date))
all <- all %>% mutate(time = ifelse(cruise == "g2", Time, time))

all <- all %>% select(-Latitude, -Longitude, -Date, -Time)

all <- all %>% mutate(size = ifelse(cruise == "g3Depth", 0.2, size))
all <- all %>% mutate(size = ifelse(cruise == "alohaDiel", 0.2, size))

all <- all %>% mutate(depth = ifelse(cruise == "alohaDiel", "15", depth))

all <- all %>% mutate(latitude = ifelse(cruise == "alohaDiel", 22.75, latitude))

all <- all %>% mutate(longitude = ifelse(cruise == "alohaDiel", -158, longitude))

all %>% filter(latitude < 0)
all %>% filter(longitude > 0)
all <- all %>% mutate(longitude = ifelse(longitude > 0, -1*longitude, longitude))
all %>% summarize(min = min(longitude, na.rm = TRUE), max = max(longitude, na.rm = TRUE))
all %>% summarize(min = min(latitude, na.rm = TRUE), max = max(latitude, na.rm = TRUE))

#fix date formats 
all <- all %>% mutate(date = str_replace_all(date, "\\/", "-"))
all <- all %>% mutate(date = str_replace_all(date, "\\-06\\-", "-6-"))
all <- all %>% mutate(date = str_replace_all(date, "\\-07\\-", "-7-"))
all <- all %>% mutate(date = ifelse(str_detect(date, "^2017"), str_c(date, -17), date))
all <- all %>% mutate(date = str_replace(date, "^2017\\-", ""))
all <- all %>% mutate(date = ifelse(str_detect(date, "^2015"), str_c(date, -15), date))
all <- all %>% mutate(date = str_replace(date, "^2015\\-", ""))

all <- all %>% mutate(date = ifelse(date == "6-09-17", "6-9-17", date))
all <- all %>% mutate(date = ifelse(date == "6-06-17", "6-6-17", date))

all <- all %>% mutate(time = ifelse(cruise == "alohaDiel" & date == "7-27-15", time - 2400, time))
all <- all %>% mutate(time = ifelse(cruise == "alohaDiel" & date == "7-28-15", time - 4800, time))
all <- all %>% mutate(time = ifelse(cruise == "alohaDiel" & date == "7-29-15", time - 7200, time))
all <- all %>% mutate(time = ifelse(cruise == "alohaDiel" & date == "7-30-15", time - 9600, time))

#29.58° (95 m), 32.00° (87 m), 34.00° (56 m), and 37.00° (50 m) 
all <- all %>% mutate(depth = ifelse(latitude == 32 & cruise == "g2DCM", "87", depth))
all <- all %>% mutate(depth = ifelse(latitude == 34 & cruise == "g2DCM", "56", depth))
all <- all %>% mutate(depth = ifelse(latitude == 37 & cruise == "g2DCM", "50", depth))
all <- all %>% mutate(depth = ifelse(latitude == 29.58 & cruise == "g2DCM", "95", depth))

colnames(all)
all <- all %>% select(-c(sample, station))

all <- all %>% 
  mutate(cruise = ifelse(cruise == "alohaDiel", "ALOHA diel", cruise)) %>% 
  mutate(cruise = ifelse(cruise == "g1", "Gradients1", cruise)) %>% 
  mutate(cruise = ifelse(cruise == "g2", "Gradients2", cruise)) %>% 
  mutate(cruise = ifelse(cruise == "g3", "Gradients3", cruise)) %>% 
  mutate(cruise = ifelse(cruise == "g3Diel", "Gradients3 diel", cruise)) %>% 
  mutate(cruise = ifelse(cruise == "g3Depth", "Gradients3 depth", cruise)) %>% 
  mutate(cruise = ifelse(cruise == "G2 incubations", "Gradients2 incubations", cruise))

all %>% distinct(cruise)

all %>% filter(is.na(latitude)) %>% distinct(cruise)

colnames(all)
colnames(all) <- c("Species", "Trophic mode prediction", "Trophic mode prediction probability", "Filter (µm)", "Depth (m)",  
                   "Time", "Date", "Cruise", "Taxonomic group", "Latitude", "Longitude")

all <- all %>% select(Cruise, Species, `Taxonomic group`, Date, Time, Latitude, Longitude, `Depth (m)`,
                      `Filter (µm)`,
                      `Trophic mode prediction`, `Trophic mode prediction probability`)


all <- all %>% arrange(Cruise, Species, `Taxonomic group`, Date, Time, Latitude, Longitude, `Depth (m)`, `Filter (µm)`)

all <- all %>% filter(Cruise != "Gradients2 incubations")

g2Inc <- read_csv("~/Dropbox/grad/research/g2Incubation/g2Inc_trophicPredictions_cleanedUp_updatedMarferret_marmicroDb2023_noOutliers.csv")

g2Inc <- g2Inc %>% select(taxa, xg_pred, probability, `Size Fraction`, Expt, Timepoint, Treatment)

type1 <- read_excel("~/Dropbox/grad/research/g1Surface/typeOfTaxaInPredictions_noOutliers.xlsx")
type1 <- type1 %>% mutate(species = str_extract(tax_name, "[A-Z]{1,}[ a-z\\.0-9\\-A-Z']{1,}$"))
type2 <- read_excel("~/Dropbox/grad/research/g2Surface/typeOfTaxaInPredictions_noOutliers.xlsx")
type2 <- type2 %>% mutate(species = str_extract(taxa, "[A-Z]{1,}[ a-z\\.0-9\\-A-Z]{1,}$"))
type3 <- read_excel("~/Dropbox/grad/research/g3/typeOfTaxaInPredictions_noOutliers.xlsx")
type3 <- type3 %>% mutate(species = str_extract(taxa, "[A-Z]{1,}[ a-z\\.0-9\\-A-Z]{1,}$"))

colnames(type2)[1] <- "tax_name"
colnames(type3)[1] <- "tax_name"

type <- bind_rows(type1, type2, type3)

nrow(g2Inc)
g2Inc <- g2Inc %>% left_join(type %>% distinct(species, group), by = c("taxa" = "species"))
nrow(g2Inc)

g2Inc %>% filter(is.na(group))

colnames(all)
head(g2Inc)

g2Inc <- g2Inc %>% mutate(Cruise = "Gradients2 incubations")

g2Inc <- g2Inc %>% mutate(Latitude = ifelse(Expt == "REXP1", 41.42, NA))
g2Inc <- g2Inc %>% mutate(Latitude = ifelse(Expt == "REXP2", 37.00, Latitude))
g2Inc <- g2Inc %>% mutate(Latitude = ifelse(Expt == "REXP3", 32.93, Latitude))

g2Inc <- g2Inc %>% select(-Expt)

colnames(g2Inc)
colnames(g2Inc) <- c("Species", "Trophic mode prediction", "Trophic mode prediction probability", "Filter (µm)", "Incubation time", 
                     "Incubation treatment", "Taxonomic group", "Cruise", "Latitude")

g2Inc <- g2Inc %>% mutate(Longitude = -158)

g2Inc$`Filter (µm)` <- parse_number(g2Inc$`Filter (µm)`)

all <- all %>% bind_rows(g2Inc)

all %>% distinct(`Trophic mode prediction`)
all$`Trophic mode prediction` <- str_replace(all$`Trophic mode prediction`, "ic", "y")

all <- all %>% select(Cruise, Species, `Taxonomic group`, Date, Time, Latitude, Longitude, `Depth (m)`, `Filter (µm)`, 
                      `Incubation treatment`, `Incubation time`, `Trophic mode prediction`, `Trophic mode prediction probability`)

all <- all %>% arrange(Cruise, Species, Date, Time)

all %>% filter(Cruise == "Gradients2 incubations")
all <- all %>% mutate(`Depth (m)` = ifelse(Cruise == "Gradients2 incubations", 15, `Depth (m)`))

g2IncMeta <- read_csv("g2Incubation/Gradients2_discrete_samples - RNA_exp.csv")

g2IncMeta <- g2IncMeta %>% filter(str_detect(Experiment, "REXP"))

g2IncMeta <- g2IncMeta %>% mutate(treatment = str_extract(`Sample ID`, "[A-Za-z]{2,}"))

g2IncMeta <- g2IncMeta %>% select(Experiment, Size, Date, `Time Start`, treatment) %>% 
  group_by(Experiment, Size, Date, treatment) %>% arrange(`Time Start`) %>% slice(1)

g2IncMeta <- g2IncMeta %>% mutate(timepoint = str_extract(Experiment, "[096]{1,}"))

g2IncMeta <- g2IncMeta %>% mutate(treatment = ifelse(timepoint == "0", "Ctrl", treatment))

g2IncMeta <- g2IncMeta %>% mutate(latitude = ifelse(str_detect(Experiment, "REXP1"), 41.42, NA))
g2IncMeta <- g2IncMeta %>% mutate(latitude = ifelse(str_detect(Experiment, "REXP2"), 37.00, latitude))
g2IncMeta <- g2IncMeta %>% mutate(latitude = ifelse(str_detect(Experiment, "REXP3"), 32.93, latitude))

all %>% filter(Cruise == "Gradients2 incubations") %>% head()
head(g2IncMeta)

g2IncMeta <- g2IncMeta %>% mutate(timepoint = as.numeric(timepoint))

g2IncMeta$Size <- parse_number(g2IncMeta$Size)

all %>% anti_join(g2IncMeta, by = c("Latitude" = "latitude", "Incubation treatment" = "treatment", "Incubation time" = "timepoint", "Filter (µm)" = "Size")) %>% 
  distinct(Cruise)

nrow(all)
all <- all %>% left_join(g2IncMeta, by = c("Latitude" = "latitude", "Incubation treatment" = "treatment", "Incubation time" = "timepoint", "Filter (µm)" = "Size"))
nrow(all)

all <- all %>% mutate(Date.x = ifelse(Cruise == "Gradients2 incubations", Date.y, Date.x))
all <- all %>% mutate(Time = ifelse(Cruise == "Gradients2 incubations", `Time Start`, Time))

all <- all %>% select(-c(Experiment, Date.y, `Time Start`))

colnames(all)[4]
colnames(all)[4] <- "Date"

all <- all %>% mutate(Date = str_replace_all(Date, "\\/", "-"))
all <- all %>% mutate(Date = str_replace(Date, "\\-2917$", "-17"))
all <- all %>% mutate(Date = str_replace(Date, "\\-2017$", "-17"))

colnames(all)
all %>% filter(is.na(Cruise))
all %>% filter(is.na(`Taxonomic group`))
all %>% filter(is.na(Date))
all %>% filter(is.na(Time))
all %>% filter(is.na(Latitude))
all %>% filter(is.na(Longitude))
all %>% filter(is.na(`Depth (m)`))
all %>% distinct(`Filter (µm)`)
all %>% filter(is.na(`Trophic mode prediction`))
all %>% filter(is.na(`Trophic mode prediction probability`))

all %>%
 write_csv("~/Dropbox/grad/research/trophicModePredictionsSupplementaryTable.csv")

all %>% mutate(Time = round(Time, -2)) %>% filter(!str_detect(Cruise, "diel")) %>% 
  group_by(Cruise, Time) %>% summarize(n = n()) %>% arrange(Cruise, desc(n))

