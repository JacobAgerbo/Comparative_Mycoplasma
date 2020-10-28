# Genome-resolved metagenomics suggests a mutualistic relationship between Mycoplasma and salmonid hosts
_Rasmussen et al. 2020b_

Metagenomic pipeline were performed each individual host seperately, inclduing: 

Host|No. individual|Seq Chemistry
--- | --- | --- 
Rainbow Trout|n =5|  PE150
Atlantic Salmon|n = 8|PE100
European Whitefish|n = 1 (three experimental replicates)|PE100

Host filtered metagenomic sequences can be found on ENA archive PRJEB40990.

_Please find information of parameters and modules used in Metagenomic_pipeline.sh_

## Binning and MAG curation
After assembly of host filtered reads were contigs processed in Anvio. 
The subsequent workflow is outlined at http://merenlab.org/2016/06/22/anvio-tutorial-v2/. Briefly; 
* anvi’o was used to profile the scaffolds using Prodigal/v2.6.338 with default parameters to identify genes and HMMER/v.3.339 to identify genes matching to archaeal, protista (based on http://merenlab.org/delmont-euk-scgs), and bacterial single-copy core gene collections. Also, ribosomal RNA based HMMs were identified (based on https://github.com/tseemann/barrnap). The HMMs were used to determine completeness of  metagenome assembled genomes (MAGs); 
* Kaiju 41 was used with NCBI’s non-redundant protein database ‘nr’ to infer the taxonomy of genes (as described in http://merenlab.org/2016/06/18/importing-taxonomy/); 
* we mapped short reads from the metagenomic set to the scaffolds using BWA/v0.7.1596 (minimum identity of 95%) and stored the recruited reads as BAM files using samtools 
* anvi'o profiled each BAM file to estimate the coverage and detection statistics of each scaffold, and combined mapping profiles into a merged profile database for each metagenomic set. Contigs were binned automatically, using CONCOCT, by constraining the number of clusters per metagenomic set to 10.
* Bin and MAGs where curated, following guideline from Veronika Kivenson: http://merenlab.org/2017/05/11/anvi-refine-by-veronika/

_Please find information of parameters and modules used in Anvio_pipeline.sh_

_MAGs will be publicaly available on: __figshare__ upon acceptance_

## Phylogenomic of _Mycoplasma_
Will be uploaded ASAP, co-author is unavailable for a week.

## Comparative Genomics of _Mycoplasma_
* Curated MAGs and external genomes were analysed, using comparative genomics: http://merenlab.org/2016/11/08/pangenomics-v2/
* Information, Genbank IDs, about external genome can be found in the pangenome repository
* Similarities of each amino acid sequence in every genome were calculated against every other amino acid sequence across all genomes, using blastp. minbit heuristics were implemented to eliminate weak matches between two amino acid sequences
* The MCL algorithms were used to identify gene clusters in amino acid sequence similarity search results
* Euclidean distance and ward linkage were used to organise gene clusters and genomes. Average Nucleotide Identity (ANI) was calculated, using PyANI

## Metabolic Reconstruction, using RAST
* Functional Annotations were carried out, using COG, PFAM, KEGG (GhostKoala)
* MAGs were recovered from Anvio DB and uploaded to RAST: https://rast.nmpdr.org/rast.cgi for subsystem annotation
* MAGs and external genomes were downloaded from RAST and concatenated, using pandas (see RAST_merge.py)
* RAST systems of all MAGs and genomes were imported to Gephi
* Genomes and MAGs were classified according to have been reported in intestinal environment previously or not
* In Gephi were only subsystem classified to Amino Acid metabolism kept for this analysis
* MAGs, Genomes, and RAST system were analysed, using Force Atlas 2 algorithm
* To confirm our analysis we did enrichment analysis, using anvi’o pipeline http://merenlab.org/2016/11/08/pangenomics-v2/
  * Summary of this is available on: __figshare__ upon acceptance
