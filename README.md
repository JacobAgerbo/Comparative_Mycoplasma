# Comparative_Mycoplasma
Genome-resolved metagenomics suggests a mutualistic 2 relationship between Mycoplasma and salmonid hosts

Metagenomic pipeline were performed each individual host seperately, inclduing: 
-|Host|No. individual|Seq Chemistry|-
-|Rainbow Trout|n =5|  PE150|-
-|Atlantic Salmon|n = 8|PE100|-
• European Whitefish (n = 1, but three experimental replicates)
Host filtered metagenomic sequences can be found on ENA archive PRJEB40990. 

## Binning
After assembly of host filtered reads were contigs:
• Processed in Anvio, following their guidelines: http://merenlab.org/2016/06/22/anvio-tutorial-v2/
• Binned automatically using CONCOCT in Anvio.
• Bin and MAGs where curated, following guideline from Veronika Kivenson: http://merenlab.org/2017/05/11/anvi-refine-by-veronika/
• Contigs were annotated, using Kaiju
• SCG were annotated based HMMs, using HMMER
• Functional Annotations were carried out, using COG, PFAM, KEGG (GhostKoala)

## Metabolic Reconstruction, using RAST

