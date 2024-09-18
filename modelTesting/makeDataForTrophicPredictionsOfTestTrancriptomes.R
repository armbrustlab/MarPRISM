library(tidyverse)
library(readxl)

#from cleanUpTestTranscriptomes.R
test <- read_csv("~/Dropbox/grad/research/forElaina/cleanedUpTestTranscriptomes.csv")

merg_ex_lab <- read_csv("~/Dropbox/grad/research/trophicModePrediction/Field_training_labels_noOutliers_newContaminationMetric_excluded.csv")
merg_ex_dat <- read_csv("~/Dropbox/grad/research/trophicModePrediction/Field_training_data_noOutliers_newContaminationMetric_excluded.csv")
meta <- read_excel("~/Dropbox/grad/papers/chryso/theDynamicTrophicArchitectureOfOpenOceanProtistCommunitiesRevealedThroughMachineGuidedMetatranscriptomics/dataset.xlsx")

colnames(meta) <- meta[13,]
meta <- meta[-c(1:13),]

merg_ex <- cbind(merg_ex_lab, merg_ex_dat)

nrow(merg_ex)
merg_ex <- merg_ex %>% left_join(meta %>% select(-`Trophic Mode`), by = c("M_id" = "MMETSP ID"))
nrow(merg_ex)

merg_ex <- merg_ex %>% filter(`Trophic mode` != "Un")

merg_ex <- merg_ex %>% gather(PF00001:PF17125, key = "pfam", value = "TPM")

colnames(merg_ex)[1] <- "sample"

merg_ex <- merg_ex %>% mutate(sample = str_c(sample, " ", `Organism Name`))
merg_ex <- merg_ex %>% select(-`Organism Name`)

dat <- bind_rows(test, merg_ex)

dat %>% group_by(Name, sample) %>% summarize(n = n()) %>% filter(n > 1)

dat <- dat %>% select(sample, `Trophic mode`, `Other Factors`, Reference, pfam, TPM)

dat$TPM <- as.numeric(dat$TPM)

features <- read_csv("~/Dropbox/grad/research/trophicModePrediction/Extracted_Pfams_noOutliers_newContaminationMetric_xg.csv")

supp <- dat %>% filter(!is.na(pfam))

supp <- supp %>% group_by(sample, pfam) %>% summarize(TPM = sum(TPM))

supp <- supp %>% ungroup()

supp <- supp %>% mutate(sample = str_replace(sample, "AM", "day"))
supp <- supp %>% mutate(sample = str_replace(sample, "PM", "night"))

supp <- supp %>% mutate(sample = str_replace(sample, "altemp", "altemp12"))
supp <- supp %>% mutate(sample = str_replace(sample, "kbha", "kbha01"))

supp <- supp %>% mutate(sample = str_replace(sample, "nitz", "nitzschia"))

supp <- supp %>% mutate(sample = str_replace(sample, "O1393", "CCMP1393"))

supp <- supp %>% mutate(sample = str_replace(sample, "micpol", "micromonas"))

supp <- supp %>% mutate(sample = str_replace(sample, "Pyr", "pyramimonas"))

supp <- supp %>% mutate(sample = str_replace(sample, "_hn_", "_highNutrient_"))
supp <- supp %>% mutate(sample = str_replace(sample, "_ln_", "_lowNutrient_"))

supp <- supp %>% mutate(sample = str_replace(sample, "HN_", "_highNutrient_"))
supp <- supp %>% mutate(sample = str_replace(sample, "LN_", "_lowNutrient_"))

supp <- supp %>% mutate(sample = str_replace(sample, "CCMP1393D", "CCMP1393_bacteriaDark"))
supp <- supp %>% mutate(sample = str_replace(sample, "CCMP1393L_", "CCMP1393_axenic_"))
supp <- supp %>% mutate(sample = str_replace(sample, "CCMP1393LB_", "CCMP1393_bacteriaLight_"))

supp <- supp %>% mutate(sample = str_replace(sample, "OchrD_", "BG1_bacteriaDark_"))
supp <- supp %>% mutate(sample = str_replace(sample, "OchrL_", "BG1_depletedBacteria_"))
supp <- supp %>% mutate(sample = str_replace(sample, "OchrLB_", "BG1_bacteriaLight_"))

supp <- supp %>% mutate(sample = str_replace(sample, "_r", "_rep"))
supp <- supp %>% mutate(sample = str_replace(sample, "_1", "_rep1"))
supp <- supp %>% mutate(sample = str_replace(sample, "_2", "_rep2"))
supp <- supp %>% mutate(sample = str_replace(sample, "_3", "_rep3"))
supp <- supp %>% mutate(sample = str_replace(sample, "_4", "_rep4"))

union <- bind_rows(features, read_csv("~/Dropbox/grad/research/trophicModePrediction/Extracted_Pfams.csv") %>% 
                     mutate(pfam = str_replace(Pfam, "\\.[0-9]{1,}", "")) %>% select(pfam))
union <- union %>% distinct(pfam)

supp <- supp %>% semi_join(union, by = c("pfam"))

toBind <- union %>% anti_join(supp, by = c("pfam")) %>% mutate(TPM = 0, sample = "CCMP1393_bacteriaDark_rep1")

supp <- bind_rows(supp, toBind)

supp <- supp %>% spread(key = pfam, value = TPM, fill = 0)

nonMM <- supp %>% distinct(sample) %>% filter(!str_detect(sample, "MMET")) %>% as.list()
nonMM <- nonMM$sample

mm <- supp %>% distinct(sample) %>% filter(str_detect(sample, "MMET")) %>% as.list()
mm <- mm$sample

supp$sample <- factor(supp$sample, levels = c(nonMM, mm))

supp <- supp %>% arrange(sample)

supp %>% write_csv("/Users/elainathomas/Library/CloudStorage/Dropbox/grad/research/ForElaina/testTranscriptomesSupplementaryTable.csv")


##make dataframe to run in testTranscriptomes_xgModel_xgFeatures_contaminationRemoved.ipynb
dat_update <- dat %>% semi_join(features, by = c("pfam"))

dat_update <- dat_update %>% group_by(sample, pfam) %>% summarize(TPM = sum(TPM))

miss <- features %>% anti_join(dat, by = c("pfam"))

dat_update <- bind_rows(dat_update, miss %>% mutate(sample = "O1393D_r1"))

dat_update <- dat_update %>% spread(key = pfam, value = TPM, fill = 0)

dat_update %>% ungroup() %>% select(-sample) %>% write_csv("~/Dropbox/grad/research/ForElaina/testTranscriptomes_updated.csv")


##make dataframe to run in testTranscriptomes_lambertModel.ipynb

dat_old <- dat %>% filter(!str_detect(sample, "MMET"))

features <- read_csv("~/Dropbox/grad/research/trophicModePrediction/Extracted_Pfams.csv")

colnames(features) <- "pfam"

features <- features %>% mutate(pfam = str_replace(pfam, "\\.[0-9]{1,}", ""))

dat_old <- dat_old %>% semi_join(features, by = c("pfam"))

dat_old <- dat_old %>% group_by(sample, pfam) %>% summarize(TPM = sum(TPM))

miss <- features %>% anti_join(dat_old, by = c("pfam"))

dat_old <- bind_rows(dat_old, miss %>% mutate(sample = "O1393D_r1"))

dat_old <- dat_old %>% spread(key = pfam, value = TPM, fill = 0)

dat_old %>% ungroup() %>% select(-sample) %>% write_csv("~/Dropbox/grad/research/ForElaina/testTranscriptomes_old.csv")


##analyze trophic predictions

#from testTranscriptomes_xgModel_xgFeatures_contaminationRemoved.ipynb
preds <- read_csv("~/Dropbox/grad/research/ForElaina/testTranscriptomes_predictions")

preds <- dat_update %>% ungroup() %>% select(sample) %>% 
  cbind(preds)

preds <- preds %>% mutate(xg_pred = ifelse(xg_pred == 1, "Mix", ifelse(xg_pred == 0, "Het", "Phot")))

#from testTranscriptomes_lambertModel.ipynb
preds_old <- read_csv("~/Dropbox/grad/research/ForElaina/testTranscriptomes_predictions_old")

preds_old <- dat_old %>% ungroup() %>% select(sample) %>% 
  cbind(preds_old)

colnames(preds_old)[2] <- "xg_pred_old"

nrow(preds)
preds <- preds %>% left_join(preds_old, by = c("sample"))
nrow(preds)

nrow(preds)
preds <- preds %>% left_join(merg_ex %>% distinct(sample, `Trophic mode`, `Other Factors`, Reference), by = c("sample"))
nrow(preds)

nonMM <- preds %>% distinct(sample) %>% filter(!str_detect(sample, "MMET")) %>% as.list()
nonMM <- nonMM$sample

mm <- preds %>% distinct(sample) %>% filter(str_detect(sample, "MMET")) %>% as.list()
mm <- mm$sample

preds$sample <- factor(preds$sample, levels = c(nonMM, mm))

preds <- preds %>% arrange(sample)

preds %>% filter(str_detect(sample, "MMET")) %>% filter(xg_pred == `Trophic mode`) %>% nrow()
preds %>% filter(str_detect(sample, "MMET")) %>% nrow()

preds <- preds %>% select(sample, `Other Factors`, Reference, `Trophic mode`, xg_pred)

preds %>% filter(str_detect(sample, "MMET")) %>% write_csv("~/Dropbox/grad/research/ForElaina/contaminatedTranscriptomesPredictions.csv")
