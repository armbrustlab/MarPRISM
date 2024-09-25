library(tidyverse)
library(readxl)
library(scales)

setwd("~/Dropbox/grad/research/")

g1 <- read_csv("g1Surface/G1_surface_trophicPredictions_marFerret2023_cleanedUp_noOutliers.csv")

g1 <- g1 %>% distinct(tax_name)

g2 <- read_csv("g2Surface/G2_surface_trophicPredictions_cleanedUp_updatedMarferret_marmicroDb2023_noOutliers.csv")

g2 <- g2 %>% distinct(taxa)

g2DCM <- read_csv("g2DCM/G2_DCM_trophicPredictions_cleanedUp_updatedMarferret_marmicroDb2023_noOutliers.csv")

g3 <- read_csv("g3/G3_trophicPredictions_cleanedUp_updatedMarferret_marmicroDb2023_noOutliers.csv")

g3 <- g3 %>% distinct(taxa)

colnames(g1) <- "taxa"

g <- bind_rows(g1 %>% mutate(gradients = "G1"), g2 %>% mutate(gradients = "G2"), g3 %>% mutate(gradients = "G3"))

g <- g %>% distinct(taxa, gradients) %>% arrange(taxa) 

g %>% write_csv("asvAnalysis/gradientsSpeciesWithTrophicModePredictions_noOutliers.csv")

#from Becca, updated version
asv <- read_excel("asvAnalysis/elaina_g123_18s_master 2.xlsx")

asv %>% filter(kingdom == "Haptophyta") %>% nrow()
asv %>% filter(kingdom == "Haptophyta") %>% distinct(taxa) %>% nrow()

##get haptophyte sequences for phylogenetic placement
hapt <- asv %>% filter(kingdom == "Haptophyta") %>% distinct(taxa) 

id <- data.frame(id = str_c(">contig", 1:nrow(hapt)))

hapt <- hapt %>% mutate(id = str_c(">contig", 1:nrow(hapt)))

idToContig <- hapt

hapt <- bind_rows(id, hapt)

hapt <- hapt %>% arrange(id)

hapt <- hapt %>% mutate(taxa = ifelse(is.na(taxa), id, taxa))

hapt %>% select(taxa) %>% 
  write.table("asvAnalysis/haptophytes.fasta", col.names = FALSE, row.names = FALSE, quote = FALSE)

haptTax <- read_csv("asvAnalysis/haptophytePhylogenyASVs.csv")

haptTax %>% filter(is.na(contig)) %>% nrow()
haptTax %>% filter(contig == "") %>% nrow()

idToContig <- idToContig %>% mutate(id = str_replace(id, ">", ""))

haptTax %>% anti_join(idToContig, by = c("contig" = "id"))

nrow(haptTax)
haptTax <- haptTax %>% left_join(idToContig, by = c("contig" = "id"))
nrow(haptTax)

haptTax <- haptTax %>% select(taxa.y, species)

colnames(haptTax)[1] <- "taxa"

haptTax <- haptTax %>% mutate(species = str_replace(species, "\\-", "_"))

nrow(asv)
asv <- asv %>% left_join(haptTax, by = c("taxa"))
nrow(asv)

colnames(asv)
asv <- asv %>% select(-c(domain:family))

asv <- asv %>% mutate(species.x = ifelse(!is.na(species.y), species.y, species.x))

asv <- asv %>% mutate(genus = ifelse(!is.na(species.y), str_extract(species.x, "[A-Z][a-z]{1,}"), genus))

asv <- asv %>% select(-species.y)

colnames(asv)[4] <- "species"

asv %>% distinct(depth) %>% arrange(depth)

asv %>% filter(species == "NA")
asv %>% filter(is.na(species))

asv <- asv %>% mutate(species = ifelse(str_detect(species, "_sp.$"), "NA", species))

asv %>% filter(str_detect(species, "sp$")) %>% distinct(species)

asv <- asv %>% mutate(species = ifelse(str_detect(species, "_sp$"), "NA", species))

asv <- asv %>% mutate(species = ifelse(species == "NA", genus, species))

#get rid of asvs without species or genus identification
asv <- asv %>% filter(species != "NA")

asv <- asv %>% mutate(genus = str_replace(species, "_[a-zA-Z0-9_\\.]{1,}", ""))

asv <- asv %>% mutate(genus = ifelse(str_detect(genus, "MAST-4"), "Stramenopiles", genus))

asv <- asv %>% mutate(species = ifelse(str_detect(species, "MAST-4"), "Stramenopiles_MAST-4", species))

#same asv sequence between Gephyrocapsa_oceanica and Emiliania_huxleyi
asv <- asv %>% mutate(species = ifelse(species == "Gephyrocapsa_oceanica", "Emiliania_huxleyi", species))
asv <- asv %>% mutate(genus = ifelse(genus == "Gephyrocapsa", "Emiliania", genus))

distinctSpecies <- asv %>% distinct(genus, species) %>% arrange(genus, species)

g <- g %>% mutate(genus = str_replace(taxa, " [a-zA-Z0-9_\\.\\- ]{1,}", ""))

g %>% filter(str_detect(genus, "\\'lucimarinus\\'")) %>% distinct(genus)
g <- g %>% mutate(genus = str_replace(genus, "\\'lucimarinus\\'", ""))

#https://roscoff-culture-collection.org/pelagomonadaceaecladecsp
g <- g %>% mutate(genus = ifelse(taxa == "Pelagophyceae sp. RCC1024", "Pelagomonadaceae", genus))
g <- g %>% mutate(taxa = ifelse(taxa == "Pelagophyceae sp. RCC1024", "Pelagomonadaceae_clade_C_sp.", taxa))

asv <- asv %>% select(sampleID, cruise, filter, latitude, longitude, depth, featureID, genus:counts)

g <- g %>% mutate(species = str_replace(taxa, " ", "_"))

#https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1333877
g <- g %>% mutate(species = ifelse(taxa == "Brandtodinium nutricula", "Zooxanthella_nutricula", species))
g <- g %>% mutate(genus= ifelse(taxa == "Brandtodinium nutricula", "Zooxanthella", genus))

g <- g %>% mutate(species = ifelse(species == "Ostreococcus_sp. 'lucimarinus'", "Ostreococcus_lucimarinus", species))

g <- g %>% mutate(species = ifelse(species == "Phaeocystis_sp. RCC851", "Phaeocystis_globosa", species))

g <- g %>% mutate(species = ifelse(species == "Prasinoderma_singulare", "Prasinoderma_singularis", species))

g <- g %>% mutate(species = ifelse(taxa == "Stramenopiles sp. TOSAG23-2", "Stramenopiles_MAST-4", species))
g <- g %>% mutate(species = ifelse(taxa == "Stramenopiles sp. TOSAG23-3", "Stramenopiles_MAST-4", species))

g <- g %>% mutate(genus = ifelse(genus == "Bolidomonas", "Triparma", genus))

#I checked all of these
g %>% anti_join(distinctSpecies, by = c("species")) %>% semi_join(distinctSpecies, by = c("genus")) %>% distinct(species)

#Chloroparvula_cf. pacifica RCC999 = Chloroparvula_pacifica
#Chloroparvula_sp. RCC696
#Prasinoderma_singulare Prasinoderma_singularis
#Stramenopiles_MAST-4 = MAST-4A, MAST-4B, MAST-4C, MAST-4D, MAST-4E
#TOSAG23-2
#TOSAG23-3

#Bolidomonas sp. 1657 = Parmales sp. RCC1657	= Triparma sp. 1657 
#bolidomonas = triparma which are genera
#Parmales is the order so not all parmales are bolidomonas/triparma

#Dinophysis acuminata 
#Kryptoperidinium foliaceum = Glenodinium foliaceum (not there either)
#Oxyrrhis marina
#Oxytricha trifallax = Sterkiella histriomuscorum (not there either)
g %>% anti_join(distinctSpecies, by = c("species")) %>% anti_join(distinctSpecies, by = c("genus")) %>% distinct(taxa)

asv_1 <- asv %>% semi_join(g, by = c("species"))

miss <- g %>% anti_join(asv, by = c("species"))

asv_2 <- asv %>% semi_join(miss, by = c("genus")) %>% 
  semi_join(g, by = c("species" = "genus"))

asv <- bind_rows(asv_1 %>% mutate(level = "species"), asv_2 %>% mutate(level = "genus"))

asv %>% filter(cruise == "G1") %>% distinct(genus, species) %>% arrange(genus)

asv %>% distinct(species) %>% nrow()
g %>% distinct(species) %>% nrow()

asv <- asv %>% mutate(oldSpecies = ifelse(species == "Ostreococcus_lucimarinus", "Ostreococcus sp. 'lucimarinus'", species))
asv <- asv %>% mutate(oldSpecies = ifelse(species == "Pelagomonadaceae_clade_C_sp.", "Pelagophyceae sp. RCC1024", oldSpecies))
asv <- asv %>% mutate(oldSpecies = ifelse(species == "Zooxanthella_nutricula", "Brandtodinium nutricula", oldSpecies))
asv <- asv %>% mutate(oldSpecies = ifelse(species == "Prasinoderma_singularis", "Prasinoderma_singulare", oldSpecies))

asv <- asv %>% mutate(oldSpecies = str_replace(oldSpecies, "_", " "))

asv %>% distinct(cruise, depth)

asv %>% filter(is.na(depth))

asv <- asv %>% arrange(species, featureID) 

distinctASV <- asv %>% distinct(species, featureID)

distinctASV <- distinctASV %>% mutate(asvNum = ifelse(species != lag(species), 1, NA))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID == "022891d41c76e21fa0e755e5fe758dad", 1, asvNum))

distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(featureID != lag(featureID) & species == lag(species), lag(asvNum)+1, asvNum))

distinctASV %>% filter(is.na(asvNum)) %>% nrow()
distinctASV <- distinctASV %>% mutate(asvNum = ifelse(species == "Aureococcus_anophagefferens", 1:87, asvNum))
distinctASV %>% filter(is.na(asvNum))

distinctASV$asvNum <- as.character(distinctASV$asvNum)

asv %>% distinct(cruise)

all3 <- asv

nrow(all3)
all3 <- all3 %>% left_join(distinctASV, by = c("species", "featureID"))
nrow(all3)

mixed <- read_csv("g1G2G3/mixedSpeciesStrict_dielIncubations.csv")

all3 %>% semi_join(mixed, by = c("oldSpecies" = "taxa")) %>% filter(level == "species") %>% distinct(oldSpecies) 
all3 %>% anti_join(mixed, by = c("oldSpecies" = "taxa")) %>% filter(level == "species") %>%  distinct(oldSpecies)

dashed <- all3 %>% filter(depth <= 15) %>%
  semi_join(mixed, by = c("oldSpecies" = "taxa")) %>% 
  mutate(latitude = round(latitude)) %>%
  group_by(latitude, oldSpecies, featureID, asvNum, cruise) %>%
  summarize(sd = sd(counts), counts = sum(counts), n = n()) %>%
  mutate(oldSpecies = str_replace(oldSpecies, " ", "\n")) %>%
  mutate(cruise = ifelse(cruise == "G1", "Gradients1: 2016", cruise)) %>%
  mutate(cruise = ifelse(cruise == "G2", "Gradients2: 2017", cruise)) %>%
  mutate(cruise = ifelse(cruise == "G3", "Gradients3: 2019", cruise)) %>%
  mutate(oldSpecies = factor(oldSpecies, levels = c("Karlodinium\nveneficum", "Pelagodinium\nbeii", "Tripos\nfusus")))

dashed <- dashed %>% ungroup() %>% select(latitude, oldSpecies, cruise) %>% mutate(n = 0)
dashed <- dashed %>% distinct() %>% spread(key = cruise, value = n, fill = 0)
dashed <- dashed %>% gather(`Gradients1: 2016`:`Gradients3: 2019`, key = "cruise", value = "n")

all3 %>% filter(depth <= 15) %>%
  semi_join(mixed, by = c("oldSpecies" = "taxa")) %>% 
  mutate(latitude = round(latitude)) %>%
  group_by(latitude, oldSpecies, featureID, asvNum, cruise) %>%
  summarize(sd = sd(counts), counts = sum(counts), n = n()) %>%
  mutate(oldSpecies = str_replace(oldSpecies, " ", "\n")) %>%
  
  mutate(cruise = ifelse(cruise == "G1", "Gradients1: 2016", cruise)) %>%
  mutate(cruise = ifelse(cruise == "G2", "Gradients2: 2017", cruise)) %>%
  mutate(cruise = ifelse(cruise == "G3", "Gradients3: 2019", cruise)) %>%
  
  mutate(oldSpecies = factor(oldSpecies, levels = c("Karlodinium\nveneficum", "Pelagodinium\nbeii", "Tripos\nfusus"))) %>%
  
  ggplot(aes(x = latitude, y = counts, fill = asvNum)) +
  geom_bar(stat = 'identity', width = .75, color = 'black') +
  theme_classic() +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size = 26, color = 'black')) + 
  theme(strip.text.y = element_text(size = 18, color = 'black')) + 
  theme(axis.text.x = element_text(size = 22, color = 'black'))  + 
  theme(axis.text.y = element_text(size = 18, color = 'black')) + 
  theme(axis.title.y = element_text(size = 26, color = 'black')) + 
  theme(legend.text = element_text(size = 20, color = 'black')) +
  theme(legend.title = element_text(size = 22, color = 'black')) +
  theme(axis.title.x = element_text(size = 26, color = 'black')) + 
  labs(color = "Taxa", y = "Relative abundance", x = "Latitude (°N)") +
  facet_grid(rows = vars(oldSpecies), cols = vars(cruise), scales = 'free_y') + 
  theme(legend.position = 'none') +
  scale_y_continuous(breaks=pretty_breaks(n = 3)) + 
  scale_x_continuous(breaks = c(20,25,30,35,40)) +
  geom_vline(data=filter(dashed, cruise=="Gradients1: 2016"), aes(xintercept=32.15), colour="red", linetype="dashed") + 
  geom_vline(data=filter(dashed, cruise=="Gradients1: 2016"), aes(xintercept=33), colour="black", linetype="dashed") + 
  
  geom_vline(data=filter(dashed %>% mutate(cruise = "Gradients2: 2017"), cruise=="Gradients2: 2017"), aes(xintercept=32.5), colour="red", linetype="dashed") + 
  geom_vline(data=filter(dashed %>% mutate(cruise = "Gradients2: 2017"), cruise=="Gradients2: 2017"), aes(xintercept=36.2), colour="black", linetype="dashed") + 
  
  geom_vline(data=filter(dashed %>% mutate(cruise = "Gradients3: 2019"), cruise=="Gradients3: 2019"), aes(xintercept=32.45), colour="red", linetype="dashed") + 
  geom_vline(data=filter(dashed %>% mutate(cruise = "Gradients3: 2019"), cruise=="Gradients3: 2019"), aes(xintercept=35), colour="black", linetype="dashed")

ggsave("asvAnalysis/speciesWithTrophicModePredictions_ASVAbundanceOverG1G2G3_surface_differentTrophicModes_noOutliers.png", dpi = 600, width = 15, height = 8)
