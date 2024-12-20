#by Sacha Coesel, April 29, 2020
#aim, calculate normalization factors from spiked-in standards.
# note, the Norm factor is a multiplication factor!

library(dplyr)

df <- read.csv("standard_counts.csv")

#For non-selected RNA: calculate mean of the commercial BMS prok standards (exclude the smallest size standards).
#For poly-A selected RNA: calculate mean of small to large Euk standards (exclude extra small standards).
#basically ignore smallest size class of both custome and commercial standards, they are not properly recovered in our pipeline.

df$MEAN <- rowMeans(df[,c(4:9)]) #select the appropriate range, depending on NS or Poly-A
df$ratio <- df$MEAN/df$added
#Bryn transformed ratio into multiplication factor. I will be consistent and do the same
df$xFactor <- 1/df$ratio
#correct for total amount of seawater filtered
df$NORM_FACTOR <-df$xFactor /df$volume
write.csv(df, file = "norm_factors.csv")
#