library(tidyverse)
library(readxl)

setwd("~/Dropbox/grad/research/")

dat <- read_csv("g3Depth/speciesWithTrophicModePredictionsAbundance.csv")

dat <- dat %>% mutate(station = parse_number(`Add'l ID`))

dat <- dat %>% mutate(Latitude = ifelse(station == 8, 42.33, NA))
dat <- dat %>% mutate(Latitude = ifelse(station == 4, 41.67, Latitude))
dat <- dat %>% mutate(Latitude = ifelse(station == 5, 37, Latitude))
dat <- dat %>% mutate(Latitude = ifelse(station == 6, 32.92, Latitude))

dat <- dat %>% mutate(depth = ifelse(str_detect(`Add'l ID`, "15m"), "Surface", NA))
dat <- dat %>% mutate(depth = ifelse(str_detect(`Add'l ID`, "75m"), "75", depth))
dat <- dat %>% mutate(depth = ifelse(str_detect(`Add'l ID`, "125m"), "125", depth))
dat <- dat %>% mutate(depth = ifelse(str_detect(`Add'l ID`, "DCM"), "DCM", depth))

dat <- dat %>% select(tax_name, absoluteCounts, Latitude, depth, sample)
dat$sample <- as.character(dat$sample)

dat <- dat %>% mutate(cruise = "Gradients3: 2019")

mixed <- read_csv("g1G2G3/mixedSpeciesStrict_dielIncubations.csv") %>% mutate(type = "mixed")
phot <- read_csv("g1G2G3/phototrophicSpeciesStrict_dielIncubations.csv") %>% mutate(type = "phot")
het <- read_csv("g1G2G3/heterotrophicSpeciesStrict_dielIncubations.csv") %>% mutate(type = "het")

type <- bind_rows(mixed, phot, het)

type <- type %>% select(taxa, type)

dat %>% anti_join(type, by = c("tax_name" = "taxa"))

nrow(dat)
dat <- dat %>% left_join(type, by = c("tax_name" = "taxa"))
nrow(dat)

dat %>% distinct(depth)
dat <- dat %>% mutate(depthNum = ifelse(depth == "Surface", 15, NA))
dat <- dat %>% mutate(depthNum = ifelse(depth == "75", 75, depthNum))
dat <- dat %>% mutate(depthNum = ifelse(depth == "125", 125, depthNum))

dat %>% distinct(depth, Latitude, cruise) %>% filter(depth == "DCM") %>% filter(cruise == "Gradients3: 2019")

#DCM was at 130m at 32.9 °N, 55m at 37 °N, 50m at 41.7 °N, and 41m at 42.3 °N
dat <- dat %>% mutate(depthNum = ifelse(depth == "DCM" & cruise == "Gradients3: 2019" & Latitude == 32.92, 130, depthNum))
dat <- dat %>% mutate(depthNum = ifelse(depth == "DCM" & cruise == "Gradients3: 2019" & Latitude == 37, 55, depthNum))
dat <- dat %>% mutate(depthNum = ifelse(depth == "DCM" & cruise == "Gradients3: 2019" & Latitude == 41.67, 50, depthNum))
dat <- dat %>% mutate(depthNum = ifelse(depth == "DCM" & cruise == "Gradients3: 2019" & Latitude == 42.33, 41, depthNum))

bySpecies <- dat

top <- bySpecies %>% group_by(Latitude, type, depthNum) %>% arrange(desc(absoluteCounts)) %>% 
  slice(1:2) %>% ungroup() %>% distinct(tax_name)

bySpecies <- bySpecies %>% group_by(cruise, Latitude, type, depthNum, sample, tax_name) %>% summarize(absoluteCounts = sum(absoluteCounts))

bySpecies <- bySpecies %>% ungroup() %>% ungroup() %>% group_by(cruise, Latitude, type, depthNum, tax_name) %>% 
  summarize(sd = sd(absoluteCounts), absoluteCounts = mean(absoluteCounts), numSamples = n())

bySpecies <- bySpecies %>% ungroup()

bySpecies <- bySpecies %>% mutate(se = sd/sqrt(numSamples))

bySpecies <- bySpecies %>% mutate(depth_tax_name = str_c(depthNum, "  ", tax_name))

bySpecies %>% filter(Latitude > 40) %>% filter(type == "phot") %>% select(Latitude, depthNum, tax_name, absoluteCounts) %>% 
  mutate(absoluteCounts = absoluteCounts/1e8) %>%
  spread(key = depthNum, value = absoluteCounts) %>% filter(`50`<`15`, `75`<`50`)

bySpecies %>% filter(Latitude > 40) %>% filter(type == "phot") %>% select(Latitude, depthNum, tax_name, absoluteCounts) %>% 
  mutate(absoluteCounts = absoluteCounts/1e8) %>%
  spread(key = depthNum, value = absoluteCounts) %>% filter(`41`<`15`, `75`<`41`)

bySpecies %>% filter(Latitude > 40) %>% filter(type == "het") %>% select(Latitude, depthNum, tax_name, absoluteCounts) %>% 
  mutate(absoluteCounts = absoluteCounts/1e8) %>%
  spread(key = depthNum, value = absoluteCounts) %>% filter(`50`>`15`)

bySpecies %>% filter(Latitude > 40) %>% filter(type == "het") %>% select(Latitude, depthNum, tax_name, absoluteCounts) %>% 
  mutate(absoluteCounts = absoluteCounts/1e8) %>%
  spread(key = depthNum, value = absoluteCounts) %>% filter(`41`>`15`)

bySpecies %>% filter(Latitude > 40) %>% filter(type == "mixed") %>% select(Latitude, depthNum, tax_name, absoluteCounts) %>% 
  mutate(absoluteCounts = absoluteCounts/1e8) %>%
  spread(key = depthNum, value = absoluteCounts) %>% filter(`50`>`15`)

bySpecies %>% filter(Latitude > 40) %>% filter(type == "mixed") %>% select(Latitude, depthNum, tax_name, absoluteCounts) %>% 
  mutate(absoluteCounts = absoluteCounts/1e8) %>%
  spread(key = depthNum, value = absoluteCounts) %>% filter(`41`>`15`)

order <- bySpecies %>% arrange(desc(depthNum), desc(type), desc(tax_name)) %>% distinct(depth_tax_name)

bySpecies %>% mutate(absoluteCounts = absoluteCounts/1e9, se = se/1e9) %>% arrange(desc(absoluteCounts))

bySpecies %>% filter(cruise == "Gradients3: 2019") %>% 
  mutate(type = factor(type, levels = c("phot", "mixed", "het"))) %>%
  mutate(Latitude = str_c(as.character(Latitude), " °N")) %>%
  mutate(depth_tax_name = factor(depth_tax_name, levels = order$depth_tax_name)) %>%
  ggplot(aes(x = depth_tax_name, y = absoluteCounts/1e9, fill = type)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  geom_errorbar(aes(ymin=absoluteCounts/1e9-se/1e9, ymax=absoluteCounts/1e9+se/1e9), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~Latitude, nrow = 1, scales = "free_y") + 
  coord_flip() + 
  theme_classic() +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size = 26, color = 'black')) + 
  theme(strip.text.y = element_text(size = 26, color = 'black')) + 
  theme(axis.text.x = element_text(size = 22, color = 'black'))  + 
  theme(axis.text.y = element_text(size = 12, color = 'black')) + 
  theme(axis.title.y = element_text(size = 26, color = 'black')) + 
  theme(legend.text = element_text(size = 20, color = 'black')) +
  theme(legend.title = element_text(size = 22, color = 'black')) +
  theme(axis.title.x = element_text(size = 26, color = 'black')) +  
  labs(y = "Billion transcripts per liter", x = "Depth (m)") + 
  scale_fill_manual(values = c("forestgreen", "purple", "orange"))

ggsave("g1G2G3/abundanceByDepthTrophicMode_species.png", dpi = 300, height = 14, width = 30)

g3_par <- read_csv("2019 SCOPE Gradients Downcast CTD Data/KM1906_Gradients3_CTD.csv")

g3_nut <- read_csv("g3/Gradients 3 KM1906 Organic and Inorganic Nutrients/Gradients 3 KM1906 Organic and Inorganic Nutrients.csv")

g3_fluor <- g3_par %>% select(time,lat,lon,depth,Fluor)
g3_par <- g3_par %>% select(time,lat,lon,depth,PAR)
g3_nut <- g3_nut %>% select(time, lat,lon,depth,NO2)

dat <- dat %>% group_by(cruise, Latitude, type, depthNum, sample) %>% summarize(absoluteCounts = sum(absoluteCounts))

dat <- dat %>% ungroup() %>% ungroup() %>% group_by(cruise, Latitude, type, depthNum) %>% 
  summarize(sd = sd(absoluteCounts), absoluteCounts = mean(absoluteCounts), numSamples = n())

dat <- dat %>% ungroup()

dat <- dat %>% mutate(se = sd/sqrt(numSamples))

g3_nut <- g3_nut %>% filter(lat %in% c(32.93300, 36.99900, 41.68100, 42.33400))

g3_nut <- g3_nut %>% mutate(lat = ifelse(lat == 32.93300, 32.92, lat))
g3_nut <- g3_nut %>% mutate(lat = ifelse(lat == 36.99900, 37.00, lat))
g3_nut <- g3_nut %>% mutate(lat = ifelse(lat == 41.68100, 41.67, lat))
g3_nut <- g3_nut %>% mutate(lat = ifelse(lat == 42.33400, 42.33, lat))

g3_nut <- g3_nut %>% filter(depth <= 150)

colnames(g3_nut)[2]
colnames(g3_nut)[2] <- "Latitude"

g3_par <- g3_par %>% mutate(lat = round(lat, 5)) %>% filter(lat %in% c(31.09583, 36.91750, 41.66917, 42.33361))
g3_fluor <- g3_fluor %>% mutate(lat = round(lat, 5)) %>% filter(lat %in% c(31.09583, 36.91750, 41.66917, 42.33361))

g3_par <- g3_par %>% mutate(lat = ifelse(lat == 31.09583, 32.92, lat))
g3_par <- g3_par %>% mutate(lat = ifelse(lat == 36.91750, 37.00, lat))
g3_par <- g3_par %>% mutate(lat = ifelse(lat == 41.66917, 41.67, lat))
g3_par <- g3_par %>% mutate(lat = ifelse(lat == 42.33361, 42.33, lat))

g3_fluor <- g3_fluor %>% mutate(lat = ifelse(lat == 31.09583, 32.92, lat))
g3_fluor <- g3_fluor %>% mutate(lat = ifelse(lat == 36.91750, 37.00, lat))
g3_fluor <- g3_fluor %>% mutate(lat = ifelse(lat == 41.66917, 41.67, lat))
g3_fluor <- g3_fluor %>% mutate(lat = ifelse(lat == 42.33361, 42.33, lat))

g3_par <- g3_par %>% filter(as.character(time) %in% c("2019-04-13 03:36:52", "2019-04-15 02:17:25", "2019-04-17 03:19:27", "2019-04-25 01:32:40"))
g3_fluor <- g3_fluor %>% filter(as.character(time) %in% c("2019-04-13 03:36:52", "2019-04-15 02:17:25", "2019-04-17 03:19:27", "2019-04-25 01:32:40"))

g3_par <- g3_par %>% filter(depth <= 150)
g3_fluor <- g3_fluor %>% filter(depth <= 150)

colnames(g3_par)[2]
colnames(g3_par)[2] <- "Latitude"

colnames(g3_fluor)[2]
colnames(g3_fluor)[2] <- "Latitude"

surfacePar <- g3_par %>% group_by(Latitude) %>% arrange(depth) %>% slice(1) %>% ungroup()

g3_par <- g3_par %>% ungroup()

nrow(g3_par)
g3_par <- g3_par %>% left_join(surfacePar %>% select(-depth), by = c("time", "Latitude", "lon"))
nrow(g3_par)

g3_par <- g3_par %>% mutate(percPAR = PAR.x/PAR.y)
g3_par <- g3_par %>% mutate(percPAR = percPAR*100)

dat %>% filter(cruise == "Gradients3: 2019") %>% 
  mutate(type = factor(type, levels = c("phot", "mixed", "het"))) %>%
  mutate(Latitude = str_c(as.character(Latitude), " °N")) %>%
  ggplot(aes(x = depthNum, y = absoluteCounts/1e9, fill = type)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  geom_errorbar(aes(ymin=absoluteCounts/1e9-se/1e9, ymax=absoluteCounts/1e9+se/1e9), width=5,
                position=position_dodge(22)) +
  coord_flip() + 
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
  labs(y = "Billion transcripts per liter", x = "Depth (m)") +
  scale_fill_manual(values = c("forestgreen", "purple", "orange")) + 
  geom_point(data = g3_nut %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = NO2*10), color = 'blue', se = FALSE, size = .6)+ 
  geom_point(data = g3_par %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = percPAR/5), color = 'red', se = FALSE, size = .6) +
  geom_line(data = g3_nut %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = NO2*10), color = 'blue', se = FALSE)+ 
  geom_line(data = g3_par %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = percPAR/5), color = 'red', se = FALSE) +
  
  geom_point(data = g3_fluor %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = Fluor*5), color = 'lightgreen', se = FALSE, size = .6) +
  geom_line(data = g3_fluor %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = Fluor*5), color = 'lightgreen', se = FALSE) +
  scale_x_reverse(limits = c(150,0)) +
  facet_wrap(~Latitude, nrow = 1)

ggsave("g1G2G3/abundanceBySpeciesType_depth.png", dpi = 600, height = 5, width = 14)

dat %>% filter(cruise == "Gradients3: 2019") %>% 
  mutate(type = factor(type, levels = c("phot", "mixed", "het"))) %>%
  mutate(Latitude = str_c(as.character(Latitude), " °N")) %>%
  ggplot(aes(x = depthNum, y = absoluteCounts/1e9, fill = type)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  geom_errorbar(aes(ymin=absoluteCounts/1e9-se/1e9, ymax=absoluteCounts/1e9+se/1e9), width=5,
                position=position_dodge(22)) +
  coord_flip() + 
  theme_classic() +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size = 26, color = 'black')) + 
  theme(strip.text.y = element_text(size = 26, color = 'black')) + 
  theme(axis.text.x = element_text(size = 14, color = 'black'))  + 
  theme(axis.text.y = element_text(size = 22, color = 'black')) + 
  theme(axis.title.y = element_text(size = 26, color = 'black')) + 
  theme(legend.text = element_text(size = 20, color = 'black')) +
  theme(legend.title = element_text(size = 22, color = 'black')) +
  theme(axis.title.x = element_text(size = 26, color = 'black')) +  
  labs(y = "Billion transcripts per liter", x = "Depth (m)") +
  scale_fill_manual(values = c("forestgreen", "purple", "orange")) + 
  geom_point(data = g3_nut %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = NO2*10), color = 'blue', se = FALSE, size = .6)+ 
  geom_point(data = g3_par %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = percPAR/5), color = 'red', se = FALSE, size = .6) +
  geom_line(data = g3_nut %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = NO2*10), color = 'blue', se = FALSE)+ 
  geom_line(data = g3_par %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = percPAR/5), color = 'red', se = FALSE) +
  
  geom_point(data = g3_fluor %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = Fluor*5), color = 'lightgreen', se = FALSE, size = .6) +
  geom_line(data = g3_fluor %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = Fluor*5), color = 'lightgreen', se = FALSE) +
  
  scale_x_reverse(limits = c(150,0)) +
  facet_wrap(~Latitude, nrow = 1) + 
  scale_y_continuous(sec.axis = sec_axis(trans=~./10, name = expression(NO[3] * "_" * NO[2]~(μmol/L))))

ggsave("g1G2G3/abundanceBySpeciesType_depth_nut.svg", dpi = 600, height = 6, width = 14)

dat %>% filter(cruise == "Gradients3: 2019") %>% 
  mutate(type = ifelse(is.na(type), "NA", type)) %>%
  mutate(type = factor(type, levels = c("phot", "mixed", "het"))) %>%
  mutate(Latitude = str_c(as.character(Latitude), " °N")) %>%
  ggplot(aes(x = depthNum, y = absoluteCounts/1e9, fill = type)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  geom_errorbar(aes(ymin=absoluteCounts/1e9-se/1e9, ymax=absoluteCounts/1e9+se/1e9), width=5,
                position=position_dodge(22)) +
  coord_flip() + 
  theme_classic() +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size = 26, color = 'black')) + 
  theme(strip.text.y = element_text(size = 26, color = 'black')) + 
  theme(axis.text.x = element_text(size = 14, color = 'black'))  + 
  theme(axis.text.y = element_text(size = 22, color = 'black')) + 
  theme(axis.title.y = element_text(size = 26, color = 'black')) + 
  theme(legend.text = element_text(size = 20, color = 'black')) +
  theme(legend.title = element_text(size = 22, color = 'black')) +
  theme(axis.title.x = element_text(size = 26, color = 'black')) +  
  labs(y = "Billion transcripts per liter", x = "Depth (m)") +
  scale_fill_manual(values = c("forestgreen", "purple", "orange")) + 
  geom_point(data = g3_nut %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = NO2*10), color = 'blue', se = FALSE, size = .6)+ 
  geom_point(data = g3_par %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = percPAR/5), color = 'red', se = FALSE, size = .6) +
  geom_line(data = g3_nut %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = NO2*10), color = 'blue', se = FALSE)+ 
  geom_line(data = g3_par %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = percPAR/5), color = 'red', se = FALSE) +
  
  geom_point(data = g3_fluor %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = Fluor*5), color = 'lightgreen', se = FALSE, size = .6) +
  geom_line(data = g3_fluor %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = Fluor*5), color = 'lightgreen', se = FALSE) +
  
  scale_x_reverse(limits = c(150,0)) +
  facet_wrap(~Latitude, nrow = 1) + 
  scale_y_continuous(sec.axis = sec_axis(trans=~./5, name="Fluorescence (mg/m^3)"))

ggsave("g1G2G3/abundanceBySpeciesType_depth_fluor.svg", dpi = 600, height = 5, width = 14)

dat %>% filter(cruise == "Gradients3: 2019") %>% 
  mutate(type = factor(type, levels = c("phot", "mixed", "het"))) %>%
  mutate(Latitude = str_c(as.character(Latitude), " °N")) %>%
  ggplot(aes(x = depthNum, y = absoluteCounts/1e9, fill = type)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  geom_errorbar(aes(ymin=absoluteCounts/1e9-se/1e9, ymax=absoluteCounts/1e9+se/1e9), width=5,
                position=position_dodge(22)) +
  coord_flip() + 
  theme_classic() +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size = 26, color = 'black')) + 
  theme(strip.text.y = element_text(size = 26, color = 'black')) + 
  theme(axis.text.x = element_text(size = 12, color = 'black'))  + 
  theme(axis.text.y = element_text(size = 22, color = 'black')) + 
  theme(axis.title.y = element_text(size = 26, color = 'black')) + 
  theme(legend.text = element_text(size = 20, color = 'black')) +
  theme(legend.title = element_text(size = 22, color = 'black')) +
  theme(axis.title.x = element_text(size = 26, color = 'black')) +  
  labs(y = "Billion transcripts per liter", x = "Depth (m)") +
  scale_fill_manual(values = c("forestgreen", "purple", "orange")) + 
  geom_point(data = g3_nut %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = NO2*10), color = 'blue', se = FALSE, size = .6)+ 
  geom_point(data = g3_par %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = percPAR/5), color = 'red', se = FALSE, size = .6) +
  geom_line(data = g3_nut %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = NO2*10), color = 'blue', se = FALSE)+ 
  geom_line(data = g3_par %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = percPAR/5), color = 'red', se = FALSE) +
  
  geom_point(data = g3_fluor %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = Fluor*5), color = 'lightgreen', se = FALSE, size = .6) +
  geom_line(data = g3_fluor %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = Fluor*5), color = 'lightgreen', se = FALSE) +
  
  scale_x_reverse(limits = c(150,0)) +
  facet_wrap(~Latitude, nrow = 1) + 
  scale_y_continuous(sec.axis = sec_axis(trans=~.*5, name="% surface PAR"))
                     
ggsave("g1G2G3/abundanceBySpeciesType_depth_par.svg", dpi = 600, height = 5, width = 14)


g3_nut %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot") %>% 
  ggplot() +
  geom_point(aes(x = depth, y = NO2), color = 'blue', se = FALSE, size = .6)+ 
  geom_line(aes(x = depth, y = NO2), color = 'blue', se = FALSE)+ 
  coord_flip() + 
  theme_classic() +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size = 14, color = 'black')) + 
  theme(strip.text.y = element_text(size = 14, color = 'black')) + 
  theme(axis.text.x = element_text(size = 6, color = 'black'))  + 
  theme(axis.text.y = element_text(size = 16, color = 'black')) + 
  theme(axis.title.y = element_text(size = 18, color = 'black')) + 
  theme(legend.text = element_text(size = 18, color = 'black')) +
  theme(legend.title = element_text(size = 18, color = 'black')) +
  theme(axis.title.x = element_text(size = 18, color = 'black')) +  
  labs(y = "NO3_NO2 (μmol/L)", x = "Depth (m)") +

  geom_point(data = g3_par %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = percPAR/70), color = 'red', se = FALSE, size = .6) +
  geom_line(data = g3_par %>% mutate(Latitude = str_c(as.character(Latitude), " °N")) %>% mutate(type = "phot"), aes(x = depth, y = percPAR/70), color = 'red', se = FALSE) +
    scale_x_reverse(limits = c(130,0), breaks = c(125,100,75,50,25,0)) +
  scale_y_continuous(breaks = c(0,5,10,15,20,25), limits = c(0,25)) +
  facet_wrap(~Latitude, nrow = 1) + 
  scale_y_continuous(limits = c(0,1.5), sec.axis = sec_axis(trans=~.*70, name="% surface PAR"))

ggsave("g1G2G3/abundanceBySpeciesType_depth_parNut.svg", dpi = 600, height = 9, width = 4)

