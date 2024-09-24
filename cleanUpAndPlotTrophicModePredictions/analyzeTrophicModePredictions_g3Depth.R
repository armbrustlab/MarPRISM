library(tidyverse)
library(readxl)

setwd("~/Dropbox/grad/research/g3Depth")

trop <- read_csv("g3Depth_trophicModePredictions_updatedMarferret_marmicroDb2023_noOutliers_fall2023")
meta <- read_csv("g3Depth_tpm_updatedMarferret_marmicroDb2023_sampleTaxa_noOutliers_fixedTPM_fall2023.csv")
prob <- read_csv("g3Depth_trophicModePredictionsProbabilities_updatedMarferret_marmicroDb2023_noOutliers_fall2023")

merged <- cbind(meta, trop, prob)

#Mix, Het, Phot
#1, 0, 2,
merged <- merged %>% mutate(xg_pred = ifelse(xg_pred == 1, "Mix", ifelse(xg_pred == 0, "Het", "Phot")))

colnames(merged)[2]
colnames(merged)[2] <- "taxa"

merged <- merged %>% mutate(xg_pred = str_c(xg_pred, "otrophic"))
merged <- merged %>% mutate(xg_pred = str_replace(xg_pred, "Hetotrophic", "Heterotrophic"))
merged %>% distinct(xg_pred)

#downloaded from https://docs.google.com/spreadsheets/d/1lKwvlNjUV8dXHKO5HkMNa2m70QQYydrTf-Fb6PRQ9oE/edit#gid=0
sample <- read_excel("lookup_armbrust_grc_rnaseq_6.xlsx")

sample <- sample %>% filter(str_detect(`Investigator Sample Name`, "D"))
sample <- sample %>% select(2,3,4)

colnames(sample)[1] <- "sample"

sample <- sample %>% mutate(stationCast = str_extract(`Add'l ID`, "S[0-9]{1,}C[0-9]{1,}"))
sample <- sample %>% mutate(Depth = str_extract(`Add'l ID`, " [0-9mDCM]{1,}"))
sample <- sample %>% mutate(Depth = str_replace(Depth, " ", ""))
sample <- sample %>% mutate(replicate = str_extract(`Add'l ID`, " [AB]"))
sample <- sample %>% mutate(replicate = str_replace(replicate, " ", ""))

sample <- sample %>% mutate(sample = str_c("G3PA.depth.", stationCast, ".", Depth, ".", replicate))
sample %>% distinct(sample)
merged %>% filter(is.na(sample))

merged %>% anti_join(sample, by = c("sample"))

nrow(merged)
merged <- merged %>% left_join(sample, by = c("sample"))
nrow(merged)

merged <- merged %>% mutate(station = str_extract(`Add'l ID`, "S[0-9]{1,}"))
merged %>% distinct(station)

merged %>% distinct(Depth)

merged <- merged %>% mutate(depth = parse_number(Depth))
merged %>% distinct(depth)
merged <- merged %>% mutate(depth = ifelse(Depth == "DCM" & station == "S4", 50, depth))
merged <- merged %>% mutate(depth = ifelse(Depth == "DCM" & station == "S5", 55, depth))
merged <- merged %>% mutate(depth = ifelse(Depth == "DCM" & station == "S6", 130, depth))
merged <- merged %>% mutate(depth = ifelse(Depth == "DCM" & station == "S8", 41, depth))

merged <- merged %>% mutate(latitude = ifelse(station == "S8", 42.33, NA))
merged <- merged %>% mutate(latitude = ifelse(station == "S4", 41.67, latitude))
merged <- merged %>% mutate(latitude = ifelse(station == "S5", 37, latitude))
merged <- merged %>% mutate(latitude = ifelse(station == "S6", 32.92, latitude))

merged %>% filter(is.na(depth))
merged %>% filter(is.na(latitude))

#exclude nonspecies and nonprotists
merged <- merged %>% filter(!taxa %in% c("cellular organisms", "Micromonas <green algae>", 
                                         "Calanus finmarchicus", "Eucyclops serrulatus", 
                                         "Hydra vulgaris", "Oikopleura dioica", "Paracyclopina nana", 
                                         "unclassified Acanthoeca", "unclassified Micromonas", 
                                         "unclassified Phaeocystis", "unclassified Chrysochromulina") & 
                              str_detect(taxa, " "))

totalSample <- merged %>% group_by(taxa, latitude, depth) %>% summarize(total = n()) 

exclude <- merged %>% group_by(taxa, xg_pred, latitude, depth) %>% summarize(numSamples = n()) %>% 
  left_join(totalSample, by = c("taxa", "latitude", "depth")) %>% mutate(prop = numSamples/total*100) %>% 
  filter(xg_pred %in% c("Heterotrophic", "Phototrophic")) %>% ungroup() %>% 
  filter(prop > 25) %>%
  group_by(taxa, latitude, depth) %>% distinct(xg_pred) %>% summarize(n = n()) %>% 
  filter(n > 1) %>% ungroup() %>% distinct(taxa, latitude, depth)

merged %>% semi_join(exclude, by = c("taxa", "latitude", "depth")) %>% nrow()

merged %>% semi_join(exclude, by = c("taxa", "latitude", "depth")) %>% 
  write_csv("excluded.csv")

merged <- merged %>% anti_join(exclude, by = c("latitude", "taxa", "depth"))

#write cleaned up trophic mode predictions to csv file
merged %>% write_csv("g3Depth_trophicPredictions_cleanedUp_updatedMarferret_marmicroDb2023_noOutliers.csv")


##plot along with G3 surface transect

#from analyzeTrophicModePredictions_g3Surface.R
g3Surf <- read_csv("../g3/G3_trophicPredictions_cleanedUp_updatedMarferret_marmicroDb2023_noOutliers.csv")

g3Surf <- g3Surf %>% semi_join(merged, by = c("taxa"))

merged %>% distinct(depth, Depth)

merged <- merged %>% mutate(Depth = ifelse(Depth == "15m", "Surface", Depth))

g3Surf <- g3Surf %>% mutate(Depth = "Surface")

colnames(g3Surf)[13]
colnames(g3Surf)[13] <- "latitude"

merged %>% distinct(Depth)

all <- bind_rows(merged, g3Surf)

all %>% distinct(depth, Depth)

mixedPreds <- read_csv("../g1G2G3/mixedSpeciesStrict_dielIncubations.csv")

all %>% semi_join(mixedPreds, by = c("taxa")) %>% distinct(taxa)
all %>% anti_join(mixedPreds, by = c("taxa")) %>% distinct(taxa)

all %>% distinct(depth, Depth)

all <- all %>% mutate(depth = ifelse(Depth == "Surface", "Surface", depth))

all %>% distinct(depth, Depth)

str(all$depth)

all <- all %>% mutate(depth = ifelse(depth != "Surface", str_c(depth, "m"), depth))

all <- all %>% mutate(xg_pred = str_replace(xg_pred, "ic", "y"))

all %>% distinct(depth, Depth)

all %>% group_by(latitude, depth, taxa, xg_pred) %>%
  summarize(n = n()) %>% arrange(desc(n))

all %>% semi_join(mixedPreds, by = c("taxa")) %>% 
  filter(taxa == "Karlodinium veneficum") %>%
  mutate(xg_pred = str_replace(xg_pred, "ic", "y")) %>%
  group_by(latitude, depth, taxa, xg_pred) %>%
  summarize(n = n()) %>%
  ungroup() %>% 
  bind_rows(data.frame(latitude = 25.9, depth = "Surface", taxa = "Karlodinium veneficum", xg_pred = "Phototrophy", n = 0)) %>%
  mutate(depth = factor(depth, levels = c("Surface", "41m", "50m", "55m", "75m", "125m", "130m"))) %>%
  ggplot(aes(x = latitude, y = xg_pred, color = xg_pred)) + geom_point(aes(size = n)) +  
  geom_point(aes(size = n),pch=21, color = "black") + 
  facet_grid(rows = vars(depth), cols = vars(taxa)) +
  theme_bw() + 
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.y = element_text(size = 20, color = 'black'), strip.text.x = element_text(size = 20, color = 'black'), 
        axis.text.y = element_text(size = 18, color = 'black'), axis.title = element_text(size = 26, color = 'black'), 
        axis.text.x = element_text(size = 18, color = 'black')) + 
  theme(legend.title=element_text(size=22), legend.text=element_text(size=20)) +
  labs(x = "Latitude (°N)", y = "", color = "", size = "Number of predictions") + 
  geom_vline(xintercept = 32.45, color = 'red', linetype="dashed") + 
  geom_vline(xintercept = 35, color = 'black', linetype="dashed") + 
  scale_color_manual(values = c("Phototrophy" = "deepskyblue2", "Heterotrophy" = "red", "Mixotrophy" = "black")) + 
  scale_size_continuous(limits = c(1,8), breaks = c(2,4,6,8)) + 
  scale_x_continuous(limits = c(25,43), breaks = c(25,30,35,40)) + 
  guides(color = FALSE)

ggsave("g3DepthAndSurface_allOrganismsTrophicPredictions_dotPlot_updatedMarferret_marmicroDb2023_notGroupedBySize_noOutliers_depths.png", height = 9, width = 9)

