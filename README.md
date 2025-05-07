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

1. **Feature selection**  
   - Training dataset is unbalanced, more phototrophic transcriptomes than heterotrophic and mixotrophic transcriptomes.
   - To address this, phototrophic transcriptomes were randomly undersampled to create four more balanced datasets:  
     - Number of phototrophic transcriptomes = 50, 80, 100, 120
     - Mixotrophic and heterotrophic transcriptomes were included in full.
     - Datasets with undersampled phototrophic transcriptomes can be found [here](https://zenodo.org/uploads/14518902).
     - Script: `modelDevelopmentTesting/mda.py`.

3. **Hyperparameter optimization**  
   - A grid search was performed to optimize model parameters.  
   - Used one dataset with 100 phototrophic transcriptomes and all mixotrophic and heterotrophic transcriptomes, located [here](https://zenodo.org/uploads/14518902).
   - Script: `modelDevelopmentTesting/parameter_gridsearch.py`.

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

### **Testing on Cultured Protist Transcriptomes**  
To further validate, MarPRISM was tested on cultured protist transcriptomes excluded from training:  
- **Accuracy (correct prediction across all replicates):** 21/27 (77.78%)  
- **Accuracy (by replicate):** 60/76 (78.95%)

## **How to run MarPRISM on marine metatranscriptomes to generate in situ trophic mode predictions**  
Details of how we process metatranscriptomes can be found [here](https://www.nature.com/articles/s41597-024-04005-5).

1. **Collect poly(A)-selected metatranscriptomes.**  
2. **Trim, quality control, and de novo assemble RNA sequences.**  
3. **Map transcripts to de novo assemblies.** We used `kallisto` estimated counts (est_counts). 
4. **Functionally annotate transcripts** using the [Pfam database](https://www.ebi.ac.uk/interpro/download/pfam/) with `hmmsearch` (E-value < 1e-05).  
5. **Taxonomically annotate assemblies.** We used `Diamond last common ancestor`, using the [Marine Functional EukaRyotic Reference Taxa (MarFERReT) reference sequence library](https://www.nature.com/articles/s41597-023-02842-4) (E-value < 1e-05).  
6. **Sum the estimated number of reads mapped to each contig** by taxonomic annotation and sample (metatranscriptome).  
7. **Normalize sequence reads to transcripts per million (TPM):**  
   - Divide the estimated number of reads mapped to each contig by its nucleotide length (in kilobases) to generate reads per kilobase (RPK).  
   - Sum the RPK by species and sample, then divide by one million to generate a conversion factor.  
   - Divide the RPK by the conversion factor to calculate transcripts per million per contig.  
   - Sum transcripts per million by Pfam for each species and sample.  
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
    - For non-diel samples, we excluded instances where >25% of trophic predictions across replicates for one species bin fell into phototrophy and heterotrophy categories.  
    - For diel samples, we excluded instances where >25% of trophic predictions across samples and timepoints in a day for one species bin fell into phototrophy and heterotrophy             categories.  

15. **Prioritize replicate-supported predictions:**  
    - When interpreting results, put more trust in trophic predictions when multiple replicates give you the same trophic mode prediction.
   
## **Example application of MarPRISM to marine metatranscriptomes**  

We ran MarPRISM on processed metatranscriptomes collected across the North Pacific Ocean: 
   - Surface
   - Euphotic depths
   - Diel cycle
   - Onboard nutrient amendment incubations

Most of these metatranscriptomes came from the Gradients (G) cruises.

G1-G3 surface, ALOHA diel, and G3 diel metatranscriptomes were processed for the [North Pacific Eukaryotic Gene Catalog](https://www.nature.com/articles/s41597-024-04005-5).

Scripts for proceasing G2 incubations and G3 depth profiles are available in `processMarineMetatranscriptomes`.

Estimated counts outputted by `kallisto` were converted to transcripts per million as outlined in previous section.  
   - Transcripts per million for each set of samples can be found [here](https://zenodo.org/uploads/14519070).

#### Datasets:
- G1PA.tpm_counts.csv (G1 surface transect)
- G2PA.tpm_counts.csv (G2 surface transect)
- G3PA.tpm_counts.csv (G3 surface transect)
- G3PA_depth.tpm_counts.csv (G3 depth profiles)
- D1PA.tpm_counts.csv (ALOHA diel study)
- G3PA_diel.tpm_counts.csv (G3 diel study)
- G2PA_incubations.tpm_counts.csv (G2 onboard nutrient amendment incubations)

#### To run MarPRISM on these transcriptomes:

```bash
conda env create -f MarPRISM_environment.mlk.yml
conda activate MarPRISM
jupyter-notebook MarPRISM.ipynb
#substitute exampleDataset.csv in MarPRISM.ipynb for one of the above datasets
conda deactivate
```
#### Output:
   - `exampleDataset_trophicPredictions.csv` which contains the trophic predictions for each taxonomic bin and sample pair that have ≥70% eukaryote core transcribed genes expressed. 

Then we **retained only trophic predictions for transcript bins identified as protists and at the species level**. 

We **excluded predictions split between phototrophy and heterotrophy**:
   - For non-diel samples, exclude instances where >25% of trophic predictions across replicates for one species bin were in both phototrophy and heterotrophy.  
   - For diel samples, exclude instances where >25% of trophic predictions across samples and timepoints in a day for one species bin were in both phototrophy and heterotrophy.  

## **How we ran feature selection**  

#### Datasets:
   - Phototrophic transcriptomes were undersampled to select features
      - Datasets used for feature selection can be found [here](https://zenodo.org/uploads/14518902).
      - `trainingData_contamLowSeqsRemoved_50phot`
      - `trainingData_contamLowSeqsRemoved_80phot`
      - `trainingData_contamLowSeqsRemoved_100phot`
      - `trainingData_contamLowSeqsRemoved_120phot`

#### To run feature selection for XGBoost model:

```bash
conda env create -f MarPRISM_environment.mlk.yml
conda activate MarPRISM
cd modelDevelopmentTesting
#find feature Pfams for xgboost model
python mda.py -model xgboost -data trainingData_contamLowSeqsRemoved_50phot -o features_contamLowSeqsRemoved_50phot_xg
python mda.py -model xgboost -data trainingData_contamLowSeqsRemoved_80phot -o features_contamLowSeqsRemoved_80phot_xg
python mda.py -model xgboost -data trainingData_contamLowSeqsRemoved_100phot -o features_contamLowSeqsRemoved_100phot_xg
python mda.py -model xgboost -data trainingData_contamLowSeqsRemoved_120phot -o features_contamLowSeqsRemoved_120phot_xg
conda deactivate
```

We then took the **union of Pfams** in `features_contamLowSeqsRemoved_50phot_xg`, `features_contamLowSeqsRemoved_80phot_xg`, `features_contamLowSeqsRemoved_100phot_xg`, and `features_contamLowSeqsRemoved_120phot_xg` that had an **importance score greater than 0**. 

## **How we ran hyperparameter optimization**  

#### Datasets:
   - Phototrophic transcriptomes were undersampled to find best performing hyperparameters
         - Dataset with 100 phototrophic transcriptomes and all mixotrophic and heterotrophic transcriptomes, located [here](https://zenodo.org/uploads/14518902).
         - `trainingData_contamLowSeqsRemoved_100phot`

#### To run hyperparameter search for XGBoost model:

```bash
conda env create -f MarPRISM_environment.mlk.yml
conda activate MarPRISM
cd modelDevelopmentTesting
#find best performing hyperparameters for xgboost model
python parameter_gridsearch.py -model xgboost -data trainingData_contamLowSeqsRemoved_100phot/data_contamLowSeqsRemoved_phot100.csv -labels trainingData_contamLowSeqsRemoved_100phot/labels_contamLowSeqsRemoved_phot100.csv -o hyperparameters_contamLowSeqsRemoved_100phot_xg
conda deactivate
```

#### **How we ran cross-validation**

We used cross-validation to **quantify the performance of MarPRISM** as well as other model versions. 

#### To run cross-validation for MarPRISM and other model versions:

```bash
conda env create -f MarPRISM_environment.mlk.yml
conda activate MarPRISM
cd modelDevelopmentTesting
jupyter-notebook crossValidation.ipynb
conda deactivate
```

#### Output:
   - 'model_overall_f1_scores.csv': the overall mean and standard error of the F1 score for each model tested. 
   - 'models_byClass_f1_scores.csv': the mean and standard error of the F1 score by trophic mode for each model tested. 
   - 'marPRISM_k_train_size_vs_f1_score_by_class.csv': For MarPRISM, the mean and standard error of the F1 score by percentage of training data used.
   - 'marPRISM_cumulative_confusion_matrix.csv': For MarPRISM, data for confusion matrix, summarizing performance across folds of cross-validation.

#### To run cross-validation on the previous version of the model (Lambert et al. 2022):

```bash
cd modelDevelopmentTesting
conda env create -f environment.mlk.yml
conda activate environment
jupyter-notebook lambertModel_crossValidation.ipynb
conda deactivate
```
#### Output:
   - 'lambert_model_overall_f1_score.csv': the overall mean and standard error of the F1 score for the previous version of the model (Lambert et al. 2022). 
   - 'lambert_model_overall_f1_score_byClass.csv': the mean and standard error of the F1 score by trophic mode for the previous version of the model (Lambert et al. 2022). 
   
## **How we ran MarPRISM on test transcriptomes**  

#### Dataset:
   - `testTranscriptomes.csv.gz` located [here](https://zenodo.org/uploads/14518902)
   - Has transcript per million counts for the test transcriptomes.
   - Transcriptomes used for testing are from publicly available sources: accession IDs for transcriptomes used to test MarPRISM are available in `testTranscriptomes.xlsx` located [here](https://zenodo.org/uploads/14518902).

#### To run MarPRISM on test transcriptomes:

 ```bash
 conda env create -f MarPRISM_environment.mlk.yml
 conda activate MarPRISM
 jupyter-notebook MarPRISM.ipynb
 #substitute exampleDataset.csv for testTranscriptomes.csv
 conda deactivate
 ```
#### Output:
   - `exampleDataset_trophicPredictions.csv` has the trophic mode predictions for the test transcriptomes. 

The trophic mode predictions can then be compared to their expected trophic mode: [testTranscriptomes.xlsx](https://zenodo.org/uploads/14518902).
