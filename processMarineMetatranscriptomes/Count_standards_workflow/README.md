## **Scripts to calculate absolute transcripts for polyA metatranscriptomes**  

We used these scripts to calculate transcripts per liter. 
   - Transcripts per liter were not used as input to MarPRISM. 
   - Rather, transcripts per liter were used to assess the abundance of protists. 

G1-G3 surface, ALOHA diel, and G3 diel metatranscriptome absolute counts were taken from [North Pacific Eukaryotic Gene Catalog](https://www.nature.com/articles/s41597-024-04005-5), but did use these scripts. 

We used the scripts in this directory to coculate transcripts per liter for the G2 incubation and G3 depthp profile samples. 
   - Resulting normalization factors can be found in 'Normalization_factors'.

**Bowtie2_map_count_standards_pipe.txt**: Bowtie2 is used to map trimmed reads to the custom standard sequences

**CustomStandardSequences.fasta**: Custom standard sequences that were spiked in

**count_custom_standards.py**: Counts the number of reads mapped to the custom standard sequences

**Calculate_Norm_Factors.R**: Calculates normalization factors based on mapping of reads to custom standard sequences. Normalization factors can then be used to convert from number of transcripts mapped to transcripts per liter. 

