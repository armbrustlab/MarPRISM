library(tidyverse)
library(readxl)

setwd("~/Dropbox/grad/research/")

g1 <- read_csv("g1Surface/excluded.csv")
g2 <- read_csv("g2Surface/excluded.csv")
g3 <- read_csv("g3/excluded.csv")

g2Inc <- read_csv("g2Incubation/excluded.csv")
aloha <- read_csv("alohaDiel/excluded.csv")

g3Depth <- read_csv("g3Depth/excluded.csv")

head(g1)
g1 <- g1 %>% select(tax_name:probability, Size, DEPTH, LATITUDE)
head(g2)
g2 <- g2 %>% select(taxa:probability, Size, `Depth (m)`, Latitude)
head(g3)
g3 <- g3 %>% select(taxa:probability, Filter, Depth, Latitude)

head(g2Inc)
g2Inc <- g2Inc %>% select(taxa:probability, `Size Fraction`, Expt:Treatment)

head(g3Depth)
g3Depth <- g3Depth %>% select(taxa:probability, depth, latitude)

g2Inc <- g2Inc %>% mutate(depth = 15)

g2Inc <- g2Inc %>% mutate(latitude = ifelse(Expt == "REXP1", 41.42, NA))
g2Inc <- g2Inc %>% mutate(latitude = ifelse(Expt == "REXP2", 37.00, latitude))
g2Inc <- g2Inc %>% mutate(latitude = ifelse(Expt == "REXP3", 32.93, latitude))

head(aloha)
aloha <- aloha %>% select(taxa:probability, time.x, date)
aloha <- aloha %>% mutate(depth = 15, size = 0.2, latitude = 24.4)

##add taxonomic group 
##add cruise variable 

taxaGroup <- read_csv("trophicModePredictionsSupplementaryTable.csv") %>% select(2,3)

head(g1)
colnames(g1) <- c("Species", "Trophic mode prediction", "Trophic mode prediction probability", "Filter (µm)", "Depth (m)", "Latitude")
head(g2)
colnames(g2) <- c("Species", "Trophic mode prediction", "Trophic mode prediction probability", "Filter (µm)", "Depth (m)", "Latitude")
head(g3)
colnames(g3) <- c("Species", "Trophic mode prediction", "Trophic mode prediction probability", "Filter (µm)", "Depth (m)", "Latitude")

head(aloha)
colnames(aloha) <- c("Species", "Trophic mode prediction", "Trophic mode prediction probability", "Sampling time (HST)", 
                     "Date", "Depth (m)", "Filter (µm)", "Latitude")

head(aloha)
aloha %>% distinct(Date)
#aloha <- aloha %>% mutate(`Sampling time (HST)` = `Sampling time (HST)` - 2400)

aloha <- aloha %>% mutate(Date = "2015-07-27")

head(g2Inc)
g2Inc <- g2Inc %>% select(-Expt)
colnames(g2Inc) <- c("Species", "Trophic mode prediction", "Trophic mode prediction probability", "Filter (µm)", 
                     "Incubation time (hr)", "Incubation treatment", "Depth (m)", "Latitude")


head(g3Depth)
colnames(g3Depth) <- c("Species", "Trophic mode prediction", "Trophic mode prediction probability", "Depth (m)", "Latitude")

g3Depth <- g3Depth %>% mutate(`Filter (µm)` = 0.2)

head(g1)
g1 <- g1 %>% mutate(`Filter (µm)` = parse_number(`Filter (µm)`))
head(g2)
g2 <- g2 %>% mutate(`Filter (µm)` = parse_number(`Filter (µm)`))
head(g3)
head(aloha)
head(g2Inc)
g2Inc <- g2Inc %>% mutate(`Filter (µm)` = parse_number(`Filter (µm)`))

aloha$`Trophic mode prediction probability` <- as.numeric(aloha$`Trophic mode prediction probability`)

dat <- bind_rows(g1 %>% mutate(Cruise = "Gradients1"), g2 %>% mutate(Cruise = "Gradients2"), g3 %>% mutate(Cruise = "Gradients3"), 
                 aloha %>% mutate(Cruise = "ALOHA diel"), g2Inc %>% mutate(Cruise = "Gradients2 incubations"), 
                 g3Depth %>% mutate(Cruise = "Gradients3 depth"))

nrow(dat)
dat <- dat %>% left_join(taxaGroup %>% distinct(Species, `Taxonomic group`), by = c("Species"))
nrow(dat)

head(dat)
dat %>% distinct(`Filter (µm)`)
dat %>% distinct(`Depth (m)`)

head(dat)
dat <- dat %>% select(Cruise, Species, `Taxonomic group`, Date, `Sampling time (HST)`, Latitude, `Depth (m)`, `Filter (µm)`, `Incubation treatment`, 
                      `Incubation time (hr)`, `Trophic mode prediction`, `Trophic mode prediction probability`)

dat <- dat %>% arrange(Cruise, Species, Date, `Sampling time (HST)`)

head(dat)

dat %>% write_csv("g1G2G3/excludedPredictionsSupplementaryTable.csv")

dat %>% group_by(`Taxonomic group`) %>% summarize(n = n()) %>% arrange(desc(n))
dat %>% group_by(Species) %>% summarize(n = n()) %>% arrange(desc(n))

read_csv("trophicModePredictionsSupplementaryTable.csv") %>% nrow()

dat %>% group_by(Species) %>% summarize(n = n()) %>% arrange(desc(n)) %>% 
  left_join(read_csv("trophicModePredictionsSupplementaryTable.csv") %>% group_by(Species) %>% summarize(total = n()), 
                  by = c("Species")) %>% 
  mutate(prop = n/(n+total)*100) %>% 
  arrange(desc(prop)) %>% select(Species, prop)


dat %>% group_by(`Taxonomic group`) %>% summarize(n = n()) %>% arrange(desc(n)) %>% 
  left_join(read_csv("trophicModePredictionsSupplementaryTable.csv") %>% group_by(`Taxonomic group`) %>% summarize(total = n()), 
            by = c("Taxonomic group")) %>% 
  mutate(prop = n/(n+total)*100) %>% 
  arrange(desc(prop))
