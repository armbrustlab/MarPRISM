library(tidyverse)

##nitz 

setwd("~/Dropbox/grad/research/ForElaina/Nitz/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("Trinity.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("Trinity.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

nitz <- quant


##altemp am

setwd("~/Dropbox/grad/research/ForElaina/ALT/am/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("Trinity_AL_AM.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("Trinity_AL_AM.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

altemp_AM <- quant


##altemp pm 

setwd("~/Dropbox/grad/research/ForElaina/ALT/pm/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("Trinity_AL_PM.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("Trinity_AL_PM.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

altemp_PM <- quant


##kbha am

setwd("~/Dropbox/grad/research/ForElaina/KAHB/am/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("Trinity_KBHA_AM.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("Trinity_KBHA_AM.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

kbha_AM <- quant


##kbha pm

setwd("~/Dropbox/grad/research/ForElaina/KAHB/pm/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("Trinity.Trinity.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("Trinity.Trinity.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

kbha_PM <- quant


##mic pol hn r1

setwd("~/Dropbox/grad/research/ForElaina/MicPolHN/r1")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

micpol_hn_1 <- quant


##mic pol hn r2

setwd("~/Dropbox/grad/research/ForElaina/MicPolHN/r2")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

micpol_hn_2 <- quant



##mic pol hn r3

setwd("~/Dropbox/grad/research/ForElaina/MicPolHN/r3")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

micpol_hn_3 <- quant


##mic pol hn r4

setwd("~/Dropbox/grad/research/ForElaina/MicPolHN/r4")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

micpol_hn_4 <- quant


##mic pol ln r1

setwd("~/Dropbox/grad/research/ForElaina/MicPolLN/r1/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

micpol_ln_1 <- quant


##mic pol ln r2

setwd("~/Dropbox/grad/research/ForElaina/MicPolLN/r2/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

micpol_ln_2 <- quant


##mic pol ln r3

setwd("~/Dropbox/grad/research/ForElaina/MicPolLN/r3/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

micpol_ln_3 <- quant

##mic pol ln r4

setwd("~/Dropbox/grad/research/ForElaina/MicPolLN/r4/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

micpol_ln_4 <- quant


##O1393D r1

setwd("~/Dropbox/grad/research/ForElaina/O1393D/r1/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

O1393D_r1 <- quant


##O1393D r2

setwd("~/Dropbox/grad/research/ForElaina/O1393D/r2/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

O1393D_r2 <- quant


##O1393D r3

setwd("~/Dropbox/grad/research/ForElaina/O1393D/r3/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

O1393D_r3 <- quant


##O1393D r4

setwd("~/Dropbox/grad/research/ForElaina/O1393D/r4/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

O1393D_r4 <- quant


##O1393L r1

setwd("~/Dropbox/grad/research/ForElaina/O1393L/r1/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

O1393L_r1 <- quant


##O1393L r2

setwd("~/Dropbox/grad/research/ForElaina/O1393L/r2/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

O1393L_r2 <- quant


##O1393L r3

setwd("~/Dropbox/grad/research/ForElaina/O1393L/r3/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

O1393L_r3 <- quant


##O1393L r4

setwd("~/Dropbox/grad/research/ForElaina/O1393L/r4/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

O1393L_r4 <- quant


##O1393LB r1

setwd("~/Dropbox/grad/research/ForElaina/O1393LB/r1/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

O1393LB_r1 <- quant


##O1393LB r2

setwd("~/Dropbox/grad/research/ForElaina/O1393LB/r2/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

O1393LB_r2 <- quant


##O1393LB r3

setwd("~/Dropbox/grad/research/ForElaina/O1393LB/r3/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

O1393LB_r3 <- quant


##O1393LB r4

setwd("~/Dropbox/grad/research/ForElaina/O1393LB/r4/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

O1393LB_r4 <- quant


##OchrD r1

setwd("~/Dropbox/grad/research/ForElaina/OchrD/r1/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

OchrD_r1 <- quant


##OchrD r2

setwd("~/Dropbox/grad/research/ForElaina/OchrD/r2/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

OchrD_r2 <- quant


##OchrD r3

setwd("~/Dropbox/grad/research/ForElaina/OchrD/r3/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

OchrD_r3 <- quant


##OchrL r1

setwd("~/Dropbox/grad/research/ForElaina/OchrL/r1/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

OchrL_r1 <- quant


##OchrL r2

setwd("~/Dropbox/grad/research/ForElaina/OchrL/r2/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

OchrL_r2 <- quant


##OchrL r3

setwd("~/Dropbox/grad/research/ForElaina/OchrL/r3/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

OchrL_r3 <- quant


##OchrLB r1

setwd("~/Dropbox/grad/research/ForElaina/OchrLB/r1/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

OchrLB_r1 <- quant


##OchrLB r2

setwd("~/Dropbox/grad/research/ForElaina/OchrLB/r2/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

OchrLB_r2 <- quant


##OchrLB r3

setwd("~/Dropbox/grad/research/ForElaina/OchrLB/r3/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

OchrLB_r3 <- quant


##PyrHN r1

setwd("~/Dropbox/grad/research/ForElaina/PyrHN/r1/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

PyrHN_r1 <- quant


##PyrHN r2

setwd("~/Dropbox/grad/research/ForElaina/PyrHN/r2/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

PyrHN_r2 <- quant


##PyrHN r3

setwd("~/Dropbox/grad/research/ForElaina/PyrHN/r3/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

PyrHN_r3 <- quant


##PyrHN r4

setwd("~/Dropbox/grad/research/ForElaina/PyrHN/r4/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

PyrHN_r4 <- quant


##PyrLN r1

setwd("~/Dropbox/grad/research/ForElaina/PyrLN/r1/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

PyrLN_r1 <- quant


##PyrLN r2

setwd("~/Dropbox/grad/research/ForElaina/PyrLN/r2/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

PyrLN_r2 <- quant


##PyrLN r3

setwd("~/Dropbox/grad/research/ForElaina/PyrLN/r3/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

PyrLN_r3 <- quant


##PyrLN r4

setwd("~/Dropbox/grad/research/ForElaina/PyrLN/r4/")

quant <- read.table("quant.sf")

colnames(quant) <- quant[1,]
quant <- quant[-1,]

name <- read_csv("sal.fasta.dammit.namemap.csv")

name <- name %>% mutate(Name = str_replace(original, " .*", ""))

quant %>% anti_join(name, by = c("Name"))

nrow(quant)
quant <- quant %>% left_join(name %>% select(Name, renamed), by = c("Name"))
nrow(quant)

gff <- read_csv("sal.fasta.dammit.gff3")

colnames(gff) <- "var1"

gff <- gff %>% filter(str_detect(var1, "HMMER"))

gff <- gff %>% mutate(renamed = str_extract(var1, "Transcript_[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_extract(var1, "Pfam:PF[0-9]{1,}"))

gff <- gff %>% mutate(pfam = str_replace(pfam, "Pfam\\:", ""))

gff <- gff %>% mutate(evalue = str_extract(var1, "[0-9]{1,}\\.[0-9]{1,}e[\\-+][0-9]{1,}"))

gff %>% filter(is.na(evalue))

gff$evalue <- as.numeric(gff$evalue)

gff <- gff %>% filter(evalue <1e-05)

gff <- gff %>% group_by(renamed) %>% arrange(evalue) %>% slice(1)

nrow(quant)
quant <- quant %>% left_join(gff %>% distinct(renamed, pfam), by = c("renamed"))
nrow(quant)

PyrLN_r4 <- quant

setwd("~/Dropbox/grad/research/cristatum/")

files <- list.files(pattern = "*Salmon.txt")

df_list <- c()

for (i in 1:length(files)) {
  df <- read.table(files[[i]])
  colnames(df) <- df[1,]
  df <- df[-1,]
  df <- df %>% mutate(sample = files[[i]])
  df_list[[i]] <- df
}


##cristatum 

setwd("~/Dropbox/grad/research/cristatum/")

files <- list.files(pattern = "*Salmon.txt")

df_list <- c()

for (i in 1:length(files)) {
  df <- read.table(files[[i]])
  colnames(df) <- df[1,]
  df <- df[-1,]
  df <- df %>% mutate(sample = files[[i]])
  df_list[[i]] <- df
}

cris <- bind_rows(df_list)

pfam <- read_csv("cristatum_pfamAnnotations_bestPfam.tab")

cris <- cris %>% mutate(nt_id = str_extract(Name, "NODE_[0-9]{1,}"))
pfam <- pfam %>% mutate(nt_id = str_extract(nt_id, "NODE_[0-9]{1,}"))

pfam %>% anti_join(cris, by = c("nt_id"))

nrow(cris)
cris <- cris %>% left_join(pfam, by = c("nt_id"))
nrow(cris)

cris$TPM <- as.numeric(cris$TPM)

cris <- cris %>% group_by(sample, V5) %>% summarize(TPM = sum(TPM))

cris <- cris %>% mutate(pfam = str_replace(V5, "\\.[0-9]{1,}", ""))

cris <- cris %>% select(-V5)

cris %>% distinct(sample)

cris <- cris %>% mutate(sample = str_replace(sample, "_.*", ""))

cris %>% distinct(sample)

cris <- cris %>% mutate(sample = ifelse(sample == "GSM8292005", "cristatum_f/2_day11_r1", sample)) 
cris <- cris %>% mutate(sample = ifelse(sample == "GSM8292006",	"cristatum_f/2_day11_r2", sample)) 
cris <- cris %>% mutate(sample = ifelse(sample == "GSM8292007", "cristatum_f/2_day11_r3", sample)) 
cris <- cris %>% mutate(sample = ifelse(sample == "GSM8292008", "cristatum_f/2_day11_r4", sample)) 
cris <- cris %>% mutate(sample = ifelse(sample == "GSM8292009", "cristatum_f/2_day11_r5", sample)) 
cris <- cris %>% mutate(sample = ifelse(sample == "GSM8292010", "cristatum_f/20_day11_r1", sample)) 
cris <- cris %>% mutate(sample = ifelse(sample == "GSM8292011", "cristatum_f/20_day11_r2", sample))
cris <- cris %>% mutate(sample = ifelse(sample == "GSM8292012", "cristatum_f/20_day11_r3", sample)) 
cris <- cris %>% mutate(sample = ifelse(sample == "GSM8292013", "cristatum_f/20_day11_r4", sample)) 
cris <- cris %>% mutate(sample = ifelse(sample == "GSM8292014", "cristatum_f/20_day11_r5", sample)) 
cris <- cris %>% mutate(sample = ifelse(sample == "GSM8292015", "cristatum_f/20_day16_r1", sample)) 
cris <- cris %>% mutate(sample = ifelse(sample == "GSM8292016", "cristatum_f/20_day16_r2", sample)) 
cris <- cris %>% mutate(sample = ifelse(sample == "GSM8292017", "cristatum_f/20_day16_r3", sample)) 
cris <- cris %>% mutate(sample = ifelse(sample == "GSM8292018", "cristatum_f/20_day16_r4", sample)) 
cris <- cris %>% mutate(sample = ifelse(sample == "GSM8292019", "cristatum_f/20_day16_r5", sample)) 

cris %>% distinct(sample)

cris$TPM <- as.character(cris$TPM)


##shiri diatoms 

setwd("~/Dropbox/grad/research/diatomShiri/")

amp <- read_csv("GSE217467_amphora_TPM.csv") %>% mutate(taxa = "amphora")
chat <- read_csv("GSE217467_chaetoceros_TPM.csv") %>% mutate(taxa = "chaetoceros")
cylin <- read_csv("GSE217467_cylindrotheca_TPM.csv") %>% mutate(taxa = "cylindrotheca")

amp <- amp %>% gather(H2O2.1:no_Fe.3, key = "sample", value = "TPM")

chat <- chat %>% gather(cont.1:H2O2.1, key = "sample", value = "TPM")

cylin <- cylin %>% gather(cont.1:H2O2.3, key = "sample", value = "TPM")

diatom <- bind_rows(amp, chat, cylin)

pfam_amp <- read_csv("amphora_pfamAnnotations_bestPfam.tab") %>% mutate(taxa = "amphora")
pfam_chat <- read_csv("chaetoceros_pfamAnnotations_bestPfam.tab") %>% mutate(taxa = "chaetoceros")
pfam_cylin <- read_csv("cylindrotheca_pfamAnnotations_bestPfam.tab") %>% mutate(taxa = "cylindrotheca")

pfam <- bind_rows(pfam_amp, pfam_chat, pfam_cylin)

pfam <- pfam %>% select(nt_id, V5, taxa)

pfam %>% anti_join(diatom, by = c("nt_id" = "contig_id", "taxa"))

nrow(diatom)
diatom <- diatom %>% left_join(pfam, by = c("taxa", "contig_id" = "nt_id"))
nrow(diatom)

diatom <- diatom %>% mutate(sample = str_c(taxa, "_", sample))

diatom <- diatom %>% group_by(sample, V5) %>% summarize(TPM = sum(TPM))

diatom <- diatom %>% mutate(pfam = str_replace(V5, "\\.[0-9]{1,}", ""))

diatom <- diatom %>% select(-V5)

diatom %>% distinct(sample)

diatom <- diatom %>% mutate(sample = str_replace(sample, "\\.", "_"))

diatom %>% distinct(sample)

diatom$TPM <- as.character(diatom$TPM)


##merge all datasets

dat <- bind_rows(altemp_AM %>% mutate(sample = "altemp_AM"), 
                 altemp_PM %>% mutate(sample = "altemp_PM"), 
                 kbha_AM %>% mutate(sample = "kbha_AM"), 
                 kbha_PM %>% mutate(sample = "kbha_PM"), 
                 micpol_hn_1 %>% mutate(sample = "micpol_hn_1"), 
                 micpol_hn_2 %>% mutate(sample = "micpol_hn_2"), 
                 micpol_hn_3 %>% mutate(sample = "micpol_hn_3"), 
                 micpol_hn_4 %>% mutate(sample = "micpol_hn_4"), 
                 
                 micpol_ln_1 %>% mutate(sample = "micpol_ln_1"), 
                 micpol_ln_2 %>% mutate(sample = "micpol_ln_2"), 
                 micpol_ln_3 %>% mutate(sample = "micpol_ln_3"), 
                 micpol_ln_4 %>% mutate(sample = "micpol_ln_4"), 
                 
                 nitz %>% mutate(sample = "nitz"), 
                 
                 
                 O1393D_r1 %>% mutate(sample = "O1393D_r1"), 
                 O1393D_r2 %>% mutate(sample = "O1393D_r2"), 
                 O1393D_r3 %>% mutate(sample = "O1393D_r3"), 
                 O1393D_r4 %>% mutate(sample = "O1393D_r4"), 
                 
                 
                 O1393L_r1 %>% mutate(sample = "O1393L_r1"), 
                 O1393L_r2 %>% mutate(sample = "O1393L_r2"), 
                 O1393L_r3 %>% mutate(sample = "O1393L_r3"), 
                 O1393L_r4 %>% mutate(sample = "O1393L_r4"), 
                 
                 
                 O1393LB_r1 %>% mutate(sample = "O1393LB_r1"), 
                 O1393LB_r2 %>% mutate(sample = "O1393LB_r2"), 
                 O1393LB_r3 %>% mutate(sample = "O1393LB_r3"), 
                 O1393LB_r4 %>% mutate(sample = "O1393LB_r4"), 
                 
                 OchrD_r1 %>% mutate(sample = "OchrD_r1"), 
                 OchrD_r2 %>% mutate(sample = "OchrD_r2"), 
                 OchrD_r3 %>% mutate(sample = "OchrD_r3"), 
                 
                 OchrL_r1 %>% mutate(sample = "OchrL_r1"), 
                 OchrL_r2 %>% mutate(sample = "OchrL_r2"), 
                 OchrL_r3 %>% mutate(sample = "OchrL_r3"), 
                 
                 OchrLB_r1 %>% mutate(sample = "OchrLB_r1"), 
                 OchrLB_r2 %>% mutate(sample = "OchrLB_r2"), 
                 OchrLB_r3 %>% mutate(sample = "OchrLB_r3"), 
                 
                 PyrHN_r1 %>% mutate(sample = "PyrHN_r1"), 
                 PyrHN_r2 %>% mutate(sample = "PyrHN_r2"), 
                 PyrHN_r3 %>% mutate(sample = "PyrHN_r3"), 
                 PyrHN_r4 %>% mutate(sample = "PyrHN_r4"), 
                 
                 PyrLN_r1 %>% mutate(sample = "PyrLN_r1"), 
                 PyrLN_r2 %>% mutate(sample = "PyrLN_r2"), 
                 PyrLN_r3 %>% mutate(sample = "PyrLN_r3"), 
                 PyrLN_r4 %>% mutate(sample = "PyrLN_r4"), 
                 
                 cris, 
                 diatom)

dat %>% write_csv("../cleanedUpTestTranscriptomes.csv")
