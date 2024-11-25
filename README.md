MarPRISM: Marine PRotist In-Situ trophic Mode predictor

To examine the in situ activity of protists, [Lambert et al., 2022](https://www.pnas.org/doi/abs/10.1073/pnas.2100916119) developed a machine learning model to predict the trophic mode of marine protist species based on gene expression from metatranscriptomes. The Lambert model code can be found [here](https://github.com/armbrustlab/trophic-mode-ml).

Recent studies ([Groussman et al., 2023](https://www.nature.com/articles/s41597-024-04005-5); [Lasek-Nesselquist and Johnson, 2019](https://academic.oup.com/gbe/article/11/11/3218/5610072); [Van Vlierberghe et al., 2021](https://link.springer.com/article/10.1186/s13104-021-05717-2)) identified that a number of the Marine Microbial Eukaryote Transcriptome Sequencing Project ([MMETSP](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001889)) transcriptomes used for training the model have a low number of sequences and/or high contamination. 

We set out to improve the Lambert model as we expected the inclusion of low-sequence and contaminated transcriptomes used for model training likely added noise to the selection of feature Pfams and reduced the accuracy of trophic predictions.

Updating model software and removing contaminated and low-sequence transcriptomes from the training data did not increase the accuracy of trophic predictions. However, these changes did greatly reduce the set of feature Pfams needed to make predictions: from 1046 to 183 feature Pfams. 
