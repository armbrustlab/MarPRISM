# **MarPRISM: Marine PRotist In Situ Trophic Mode Predictor**

<img src="https://github.com/user-attachments/assets/dbeac577-6127-4633-8e3e-dfdd6355532e" alt="marPRISMLogo" width="300"/>

## **Overview**  
To examine the **in situ** activity of protists, [Lambert et al., 2022](https://www.pnas.org/doi/abs/10.1073/pnas.2100916119) developed a machine learning model to predict the trophic mode of marine protist species based on gene expression from metatranscriptomes.  
- The Lambert model code can be found [here](https://github.com/armbrustlab/trophic-mode-ml).  

Recent studies ([Groussman et al., 2023](https://www.nature.com/articles/s41597-024-04005-5); [Lasek-Nesselquist and Johnson, 2019](https://academic.oup.com/gbe/article/11/11/3218/5610072); [Van Vlierberghe et al., 2021](https://link.springer.com/article/10.1186/s13104-021-05717-2)) identified that a number of the Marine Microbial Eukaryote Transcriptome Sequencing Project ([MMETSP](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001889)) transcriptomes used for training the model have a low number of sequences and/or high contamination.  

## **Improvements to the Lambert model**  
We set out to improve the Lambert model as we expected the inclusion of low-sequence and contaminated transcriptomes used for model training likely added noise to the selection of feature Pfams and reduced the accuracy of trophic predictions.  

### Key changes  
- **Updated model software**  
- **Removed contaminated and low-sequence transcriptomes.**

These changes did not increase the accuracy of trophic predictions. However:  
- The set of feature Pfams needed for reliable predictions was reduced from **1046 to 183 feature Pfams**.
- More info on 183 feature Pfams can be found in `furtherInfo/featurePfams.xlsx`.

## **Model development**  

After removing contaminated and low-sequence transcriptomes, we conducted feature selection using mean decrease in accuracy to identify the feature Pfams essential for model performance.  

1. **Balancing training data**  
   - Training data was imbalanced:  
     - 44 heterotrophic entries  
     - 85 mixotrophic entries  
     - 258 phototrophic entries  
   - To address this, phototrophic transcriptomes were randomly undersampled to create five balanced training datasets:  
     - Number of phototrophic transcriptomes = 50, 80, 100, 120, 140  
     - Mixotrophic and heterotrophic transcriptomes were included in full.  

2. **Feature selection**  
   - Conducted on the balanced training datasets to identify essential Pfams.  
   - Script: `modelDevelopmentTesting/mda.py`  

3. **Hyperparameter optimization**  
   - A grid search was performed to optimize model parameters.  
   - Used one dataset with 100 phototrophic transcriptomes and all mixotrophic and heterotrophic transcriptomes.  
   - Script: `modelDevelopmentTesting/parameter_gridsearch.py`  

## **Model performance**  

### **Cross-validation**  
The performance of MarPRISM was estimated using cross-validation:  
- Cross-validation can be run in `modelDevelopmentTesting/marPRISM_crossValidation.ipynb`.
- The model was trained on **83% of the training data** and tested on the remaining data.  
- Performance was evaluated using **F1 score** (Mean F1 score ± standard error)
  - **Overall mean**: 0.944 ± 0.0154  
  - **Heterotrophy mean**: 0.958 ± 0.0271  
  - **Mixotrophy mean**: 0.888 ± 0.0254  
  - **Phototrophy mean**: 0.961 ± 0.0124  

**Mixotrophy** was the most difficult trophic mode to predict, likely due to overlapping Pfams with both phototrophy and heterotrophy.  

### **Testing on cultured protist transcriptomes**  
We further quantified MarPRISM's performance by testing its ability to make trophic predictions for cultured protist transcriptomes not present in the training data. 
- `modelDevelopmentTesting/testTranscriptomes.csv` can be run through `MarPRISM.ipynb` to recreate trophic predictions.
- **21/27 (77.78%) protist cultures** were correctly predicted across all replicate transcriptomes.  
- **60/76 (78.95%) transcriptomes** were correctly predicted when replicate transcriptomes were considered individually


## **How to run the model on marine metatranscriptomes**  
Details of how we process metatranscriptomes can be found [here](https://www.nature.com/articles/s41597-024-04005-5).

1. **Collect poly(A)-selected metatranscriptomes.**  
2. **Trim, quality control, and de novo assemble RNA sequences.**  
3. **Map transcripts to de novo assemblies.** We used `kallisto` estimated counts (est_counts). 
4. **Functionally annotate transcripts** using the [Pfam database](https://www.ebi.ac.uk/interpro/download/pfam/) with `hmmsearch` (E-value < 1e-05).  
5. **Taxonomically annotate assemblies.** We used `Diamond last common ancestor`, using the [Marine Functional EukaRyotic Reference Taxa (MarFERReT) reference sequence library](https://www.nature.com/articles/s41597-023-02842-4) (E-value < 1e-05).  
6. **Sum the estimated number of reads mapped to each contig** by taxonomic annotation and sample (metatranscriptome).  
7. **Normalize sequence reads to TPM:**  
   - Divide the estimated number of reads mapped to each contig by its nucleotide length (in kilobases) to generate reads per kilobase (RPK).  
   - Sum the RPK by species and sample, then divide by one million to generate a conversion factor.  
   - Divide the RPK by the conversion factor to calculate TPM per contig.  
   - Sum TPMs by Pfam for each species and sample.  
8. **Filter species bins:** Retain only species bins with at least 70% of eukaryotic core transcribed genes (CTGs) expressed ([MarFERReT.v1.core_genes.csv](https://zenodo.org/records/10554340), filter for lineage Eukaryota). Retain only species bins identified as protists. 
9. **Create a data frame:**  
   Fill in missing Pfams for a species, sample pair with `0`.  

   |                  | Pfam1 | Pfam2 | Pfam3 | ...  |
   |------------------|-------|-------|-------|------|
   | Species1_sample1 | TPM   | TPM   | TPM   | TPM  |
   | Species2_sample1 | TPM   | TPM   | TPM   | TPM  |
   | Species1_sample2 | TPM   | TPM   | TPM   | TPM  |

10. **Create the Conda environment for MarPRISM**  
    ```bash
    conda env create -f MarPRISM_environment.mlk.yml
    ```

11. **Activate the Conda environment for MarPRISM**  
    ```bash
    conda activate MarPRISM
    ```
    
12. **Run your dataframe through Jupyter Notebook `MarPRISM.ipynb`:**  
    - Open the Jupyter Notebook `MarPRISM.ipynb` and follow the instructions within to process the DataFrame and make predictions.
    ```bash
    jupyter-notebook MarPRISM.ipynb
    ```
    - **Files:**  
      - Training Data: `trainingDataMarPRISM.csv`
      - Features: `MarPRISM_featurePfams.csv`
      - Example DataFrame: `exampleDataset.csv`
     
     - **Output:**
       - To check that your Jupyter Notebook is working correctly, if `exampleDataset.csv` is used as input, your output `exampleDataset_trophicPredictions.csv` should match `exampleDataset_trophicPredictions_toCompare.csv`.
     
13. **Deactivate the Conda environment for MarPRISM**  
    ```bash
    conda deactivate
    ```

14. **Filter predictions based on replicate consistency:**  
    - Exclude phototrophy and heterotrophy predictions that are evenly split between replicate metatranscriptomes for the same species bin.  
    - For non-diel samples, we excluded instances where ≥25% of trophic predictions across replicates for one species bin fall into phototrophy and heterotrophy categories.  

15. **Prioritize replicate-supported predictions:**  
    - When interpreting results, put more trust in trophic predictions when multiple replicates give you the same trophic mode prediction.
