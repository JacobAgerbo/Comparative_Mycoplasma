# Comparative_Mycoplasma
Genome-resolved metagenomics suggests a mutualistic 2 relationship between Mycoplasma and salmonid hosts

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

## Comparative Genomics of _Mycoplasma_
* Curated MAGs and external genomes were analysed, using comparative genomics: http://merenlab.org/2016/11/08/pangenomics-v2/
* 

## Metabolic Reconstruction, using RAST
* Functional Annotations were carried out, using COG, PFAM, KEGG (GhostKoala)
* MAGs were recovered from Anvio DB and uploaded to RAST: https://rast.nmpdr.org/rast.cgi for subsystem annotation
* MAGs and external genomes were downloaded from RAST and concatenated, using pandas (see RAST_merge.py)
* RAST systems of all MAGs and genomes were imported to Gephi
* Genomes and MAGs were classified according to have been reported in intestinal environment previously or not
* In Gephi were only subsystem classified to Amino Acid metabolism kept for this analysis
* MAGs, Genomes, and RAST system were analysed, using Force Atlas 2 algorithm
* To confirm our analysis we did enrichment analysis, using anvi’o pipeline http://merenlab.org/2016/11/08/pangenomics-v2/
