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

2. **Feature selection**  
   - Conducted on the balanced training datasets to identify essential Pfams.
   - - To address this, phototrophic transcriptomes were randomly undersampled to create five balanced training datasets:  
     - Number of phototrophic transcriptomes = 50, 80, 100, 120, 140  
     - Mixotrophic and heterotrophic transcriptomes were included in full.
     - Balanced training datasets can be found [here](https://zenodo.org/uploads/14518902).
   - Script: `modelDevelopmentTesting/mda.py`
   - Using 

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
- Cross-validation was conducted using the script located at: `modelDevelopmentTesting/crossValidation.ipynb`.
- We tested other versions of the model in the same script, including a Random Forest model and different sets of feature Pfams.

#### To run cross-validation:

```bash
conda env create -f MarPRISM_environment.mlk.yml
conda activate MarPRISM
cd modelDevelopmentTesting
jupyter-notebook crossValidation.ipynb
conda deactivate
```
For each model version, this will output two csv files: one for the mean and standard error of the F1 score for the overall model, and one for the mean and standard error of the F1 score treating each trophic mode separately. For MarpRISM, a csv file will be outputted with the mean and standard error of the F1 score by percentage of training data used. 

- To run cross-validation on the previous version of the model (Lambert et al. 2022):
```bash
cd modelDevelopmentTesting
conda env create -f environment.mlk.yml
conda activate environment
jupyter-notebook lambertModel_crossValidation.ipynb
conda deactivate
```
This will output two csv files for the previous version of the model (Lambert et al. 2022): one for the mean and standard error of the F1 score for the overall model, and one for the mean and standard error of the F1 score treating each trophic mode separately. 

**Mixotrophy** was the most difficult trophic mode to predict, likely due to overlapping Pfams with both phototrophy and heterotrophy.  

### **Testing on cultured protist transcriptomes**  
We further quantified MarPRISM's performance by testing its ability to make trophic predictions for cultured protist transcriptomes not present in the training data. 
- `modelDevelopmentTesting/testTranscriptomes.csv` can be run through `MarPRISM.ipynb` to recreate trophic predictions.
- **21/27 (77.78%) protist cultures** were correctly predicted across all replicate transcriptomes.  
- **60/76 (78.95%) transcriptomes** were correctly predicted when replicate transcriptomes were considered individually
- Transcript per million counts for the test transcriptomes are located [here](https://zenodo.org/uploads/14518902)
- The transcriptomes used for testing are from publicly available sources: accession IDs for transcriptomes used to test MarPRISM are available in [testTranscriptomes.xlsx](https://zenodo.org/uploads/14518902) 
- To run the test transcriptomes through MarPRISM:
 ```bash
 conda env create -f MarPRISM_environment.mlk.yml
 conda activate MarPRISM
 jupyter-notebook MarPRISM.ipynb
 substitute exampleDataset.csv for testTranscritomes.csv (https://zenodo.org/uploads/14518902) 
 conda deactivate
 ```
The trophic mode predictions for the test transcriptomes can be compared to their metadata: [testTranscriptomes.xlsx](https://zenodo.org/uploads/14518902) 

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
8. **Filter species bins:** Retain only transcript bins identified as protists and at the species level. 
10. **Create a data frame:**  
   Fill in missing Pfams for a species, sample pair with `0`.  

   |                  | Pfam1 | Pfam2 | Pfam3 | ...  |
   |------------------|-------|-------|-------|------|
   | Species1_sample1 | TPM   | TPM   | TPM   | TPM  |
   | Species2_sample1 | TPM   | TPM   | TPM   | TPM  |
   | Species1_sample2 | TPM   | TPM   | TPM   | TPM  |

11. **Create the Conda environment for MarPRISM**  
    ```bash
    conda env create -f MarPRISM_environment.mlk.yml
    ```

12. **Activate the Conda environment for MarPRISM**  
    ```bash
    conda activate MarPRISM
    ```
    
13. **Run your dataframe through Jupyter Notebook `MarPRISM.ipynb`:**  
    - Open the Jupyter Notebook `MarPRISM.ipynb` and follow the instructions within to process the DataFrame and make predictions.
    ```bash
    jupyter-notebook MarPRISM.ipynb
    ```
    - MarPRISM will output a warning and not make predictions for species bin sample pairs with <70% of eukaryote core transcribed genes (CTGs) expressed.
      - Trophic predictions are not reliable for species bins with low coverage.
      - The eukaryote CTGs `MarFERReT.v1.core_genes_eukaryota.csv` are from [MarFERReT.v1.core_genes.csv](https://zenodo.org/records/10554340), filtering for lineage Eukaryota.

    - **Files:**  
      - Training Data: `trainingDataMarPRISM.csv`
      - Features: `MarPRISM_featurePfams.csv`
      - Example DataFrame: `exampleDataset.csv`
      - Eukaryote CTGs: `MarFERReT.v1.core_genes_eukaryota.csv`

     - **Output:**
       - To check that your Jupyter Notebook is working correctly, if `exampleDataset.csv` is used as input, your output `exampleDataset_trophicPredictions.csv` should match    
         `exampleDataset_trophicPredictions_toCompare.csv`.
     
13. **Deactivate the Conda environment for MarPRISM**  
    ```bash
    conda deactivate
    ```

14. **Filter predictions based on replicate consistency:**  
    - Exclude phototrophy and heterotrophy predictions that are evenly split between replicate metatranscriptomes for the same species bin.  
    - For non-diel samples, we excluded instances where ≥25% of trophic predictions across replicates for one species bin fall into phototrophy and heterotrophy categories.  

15. **Prioritize replicate-supported predictions:**  
    - When interpreting results, put more trust in trophic predictions when multiple replicates give you the same trophic mode prediction.
   
### Example use of MarPRISM:

We ran MarPRISM on TPM counts collected across the North Pacific Ocean, from surface and euphotic depths, across the diel cycle, and onboard nutrient amendment incubations.  
The TPM counts for these samples can be found on [10.5281/zenodo.14519070](https://www.nature.com/articles/s41597-024-04005-5).

You can substitute `exampleDataset.csv` with the following datasets from the Gradients (G) cruises to generate trophic mode predictions. 

Descriptions for **G1-G3 surface**, **ALOHA diel**, and **G3 diel** samples are provided by [North Pacific Ocean study](https://www.nature.com/articles/s41597-024-04005-5).

#### Available Datasets:
- **G1PA.tpm_counts.csv.gz** (G1 surface)
- **G2PA.tpm_counts.csv.gz** (G2 surface)
- **G3PA.tpm_counts.csv.gz** (G3 surface)
- **G3PA_depth.tpm_counts.csv.gz** (G3 depth profiles)
- **D1PA.tpm_counts.csv.gz** (ALOHA diel)
- **G3PA_diel.tpm_counts.csv.gz** (G3 diel)
- **G2PA_incubations.tpm_counts.csv.gz** (G2 onboard nutrient amendment incubations)

#### Additional Resources:
- **G2 Incubations Read Processing and Mapping Scripts, and Metadata**:  
  [GitHub - G2 Read Processing](https://github.com/armbrustlab/armbrust-metat/tree/main/gradients2/g2_dcm_rr_pa_metat)

- **G3 Depth Read Processing, Assembly, and Mapping Scripts, and Metadata**:  
  [GitHub - G3 Depth Read Processing](https://github.com/armbrustlab/armbrust-metat/tree/main/gradients3/g3_depth_pa_metat)

- **Code to Generate ALOHA Diel Assemblies**:  
  [GitHub - ALOHA Diel Assemblies and North Pacific Eukaryotic Gene Catalog](https://github.com/armbrustlab/NPac_euk_gene_catalog)
