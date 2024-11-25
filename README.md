# **MarPRISM: Marine PRotist In-Situ Trophic Mode Predictor**  

## **Overview**  
To examine the **in situ** activity of protists, [Lambert et al., 2022](https://www.pnas.org/doi/abs/10.1073/pnas.2100916119) developed a machine learning model to predict the trophic mode of marine protist species based on gene expression from metatranscriptomes.  
- The Lambert model code can be found [here](https://github.com/armbrustlab/trophic-mode-ml).  

Recent studies ([Groussman et al., 2023](https://www.nature.com/articles/s41597-024-04005-5); [Lasek-Nesselquist and Johnson, 2019](https://academic.oup.com/gbe/article/11/11/3218/5610072); [Van Vlierberghe et al., 2021](https://link.springer.com/article/10.1186/s13104-021-05717-2)) identified that a number of the Marine Microbial Eukaryote Transcriptome Sequencing Project ([MMETSP](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001889)) transcriptomes used for training the model have a low number of sequences and/or high contamination.  

## **Improvements to the Lambert model**  
We set out to improve the Lambert model as we expected the inclusion of low-sequence and contaminated transcriptomes used for model training likely added noise to the selection of feature Pfams and reduced the accuracy of trophic predictions.  

### Key changes  
- **Updated model software**  
- **Removed contaminated and low-sequence transcriptomes**  

These changes did not increase the accuracy of trophic predictions. However:  
- The set of feature Pfams needed for reliable predictions was reduced from **1046 to 183 feature Pfams**.  

## **Model Performance**  

### **Cross-Validation**  
The performance of MarPRISM was estimated using cross-validation:  
- The model was trained on **83% of the training data** and tested on the remaining data.  
- Performance was evaluated using **F1 score** (Mean F1 score ± standard error)
  - **Overall mean**: 0.944 ± 0.0154  
  - **Heterotrophy mean**: 0.958 ± 0.0271  
  - **Mixotrophy mean**: 0.888 ± 0.0254  
  - **Phototrophy mean**: 0.961 ± 0.0124  

**Mixotrophy** was the most difficult trophic mode to predict, likely due to overlapping Pfams with both phototrophy and heterotrophy.  

### **Testing on cultured protist transcriptomes**  
We further quantified MarPRISM's performance by testing its ability to make trophic predictions for cultured protist transcriptomes not present in the training data:  
- **21/27 (77.78%) protist cultures** were correctly predicted across all replicate transcriptomes.  
- **60/76 (78.95%) transcriptomes** were correctly predicted when replicate transcriptomes were considered individually

How to run model on marine metatranscriptomes
Collect poly(A)-selected metatranscriptomes
Trim, quality control and de novo assemble RNA sequences
Map transcripts to de novo assemblies with kallisto
Functionally annotate transcripts with [Pfam database](https://www.ebi.ac.uk/interpro/download/pfam/) with hmmsearch (E-value < 1e-05) 
Taxonomically annotate assemblies with Diamond last common ancestor (Buchfink et al., 2015), using the [Marine Functional EukaRyotic Reference Taxa (MarFERReT) reference sequence library](https://www.nature.com/articles/s41597-023-02842-4) (E-value < 1e-05)
Sum the estimated number of reads mapped to each contig (outputted by kallisto) by taxonomic annotation and sample (metatranscriptome) 
To normalize sequence reads to TPM, the estimated number of reads mapped to each contig (outputted by kallisto) was divided by the contig nucleotide length (kilobases) of the contig to generate reads per kilobase (RPK). The RPK for each contig was summed by species and sample (metatranscriptome) and divided by one million to generate a conversion factor; RPK divided by the conversion factor is the TPM per contig. TPMs were summed by Pfam for each species and each sample (metatranscriptome). 
For each sample (metatranscriptome), get only the species bins that have at least 70% of eukaryotic core transcribed genes (CTGs) expressed ([MarFERReT.v1.core_genes.csv](https://zenodo.org/records/10278540), filter for lineage Eukaryota) 
Make a dataframe (filling in missing Pfams in a species/sample with 0
	Pfam1	Pfam2	Pfam3	etc.
Species1_sample1	TPM	TPM	TPM	TPM
Species2_sample1	TPM	TPM	TPM	TPM
Species1_sample2	TPM	TPM	TPM	TPM
etc.	TPM	TPM	TPM	TPM![image](https://github.com/user-attachments/assets/9127def6-d338-4e20-887c-d70ac23a42e4)


Create conda environment for MarPRISM
conda env create -f MarPRISM_environment.mlk.yml                               

Activate conda environment for MarPRISM
conda activate MarPRISM

Run your datafile with TPM by Pfam by species bins that pass 70% CTGs detected



