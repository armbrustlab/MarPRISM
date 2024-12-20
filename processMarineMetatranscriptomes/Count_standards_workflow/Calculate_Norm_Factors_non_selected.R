#by Sacha Coesel, Jan 6, 2024
#aim, calculate normalization factors from spiked-in standards, bacterial only.
# note, the Norm factor is a multiplication factor!

library(dplyr)
setwd("/Users/sachacoesel/Documents/Gradients/Standards/G3NS_underway_standards")
list.files()

df <- read.csv("G3NS_underway_norm_factors_temp.csv")
head(df)
colnames(df)
df <- df[,1:16]


#calculate mean of the commercial BMS prok standards. Do not use custom euk standards
#basically ignore smallest size class of both custome and commercial standards
df$MEAN <- rowMeans(df[,c(12:15)])
df$added <- 13580968026
df$ratio <- df$MEAN/df$added

#Bryn transformed ratio into multiplication factor. I will be consistent and do the same
df$xFactor <- 1/df$ratio
#correct for total amount of seawater filtered
df$NORM_FACTOR <-df$xFactor /df$Volume_L
write.csv(df, file = "G3NS_underway_norm_factors.csv")
