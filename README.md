
# Single-cell atlas of the first intra-mammalian developmental stage of the human parasite Schistosoma mansoni

This repository contains the various scripts used in the original analysis of the paper published in NatComms: Single-cell atlas of the first intra-mammalian developmental stage of the human parasite Schistosoma mansoni.

The main purpose of this code is to document the details of the analysis methodology used in the paper.

1) #Mapping: Chromium 10X single-cell RNAseq is mapped using CellRanger. Please note that this part of the analysis was done by the core facility at the Sanger Institute. However, the code for CellRanger should be accessible from the 10X website. We encourange the reader to check this with 10X. 

2) #Seurat: The analysis of this dataset is carried out by Seurat. 

3) #Random Forest Data prep: We are providing a script that was used to update the raw matrix with the planaria orthologues.

4) #Random Forest analysis: We compared planarian and schistosomes datasets as described in the methods

5) #TopGo: We also include a figure in the paper where we analyse the Go terms for each of the clusters.





