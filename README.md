
# Single-cell atlas of the first intra-mammalian developmental stage of the human parasite Schistosoma mansoni

This repository contains the various scripts used in the original analysis of the paper published in NatComms: Single-cell atlas of the first intra-mammalian developmental stage of the human parasite Schistosoma mansoni.

The main purpose of this code is to document the details of the analysis methodology used in the paper. The files are the following:

1) #Paper_Somules_MainSeurat.Rmd: The analysis of this dataset using Seurat

2) #Paper_Somules_orthomcl_to_matrix.py: We are providing a script that was used to update the raw matrix with the planaria orthologues

3) #Paper_Somules_RandomForest.Rmd: We compared planarian and schistosomes datasets as described in the methods

4) #Paper_Somules_TopGO.R: We also include a figure in the paper where we analyse the Go terms for each of the clusters

5) #Planaria_Seurat_annot_CD: Annotation used for the planaria dataset (Join_Id7)

#Mapping Notes: Chromium 10X single-cell RNAseq is mapped using CellRanger. Please note that this part of the analysis was done by the core facility at the Sanger Institute. However, the code for CellRanger should be accessible from the 10X website. We encourange the reader to check this with 10X. 


