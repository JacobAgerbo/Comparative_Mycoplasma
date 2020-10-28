#-----------------------------------------------------------------------------------------------#
#      ---------------Bin and curate salmonid related Mycoplasma MAGs---------------            #
#-----------------------------------------------------------------------------------------------#
#
# It should be noted that this is not a full functional script, but more a list of parameters and modules used to obtain this.
# Great guidelines have been made thanks to amazing work of Meren Lab: http://merenlab.org/
#
# If needed then you can install anvio, using conda environment by following this guideline: http://merenlab.org/2016/06/26/installation-v2/
## for mac:
# first get a copy of the following file (if you get a "command not found"
# error for wget, you can either install it first, or manually download the
# file from that URL:
wget http://merenlab.org/files/anvio-conda-environments/anvio-environment-6.2.yml

# just to make sure there is not a v6.2 environment already: 
conda env remove --name anvio-6.2

# create a new v6.2 environment using the file you just downloaded:
conda env create -f anvio-environment-6.2.yml

# activate that environment
conda activate anvio-6.2

# Remember to run self-test first time: anvi-self-test --suite mini

### Get your contig.fa files straight!
CON='path/to/contigs'
find $CON/*.contigs.fa > temp
sed 's/.contigs.fa//g' temp > temp2
uniq temp2 > contigs_list.txt
rm -f temp*
contig_list=$(cat contigs_list.txt)

# Prior to pipeline set up some depencies
anvi-setup-pfams --pfam-data-dir /PATH/TO/Pfam


for a in $contig_list
  do
    anvi-script-reformat-fasta $CON/"$a".contigs.fa -o $CON/"$a".contigs-fixed.fa -l 1000 --simplify-names
    anvi-gen-contigs-database -f $CON/"$a".contigs-fixed.fa -o $CON/"$a".CONTIGS.db -n 'CONTIG DB of Salmonid Samples'
    bwa index $CON/"$a".contigs-fixed.fa
    echo -e "yay! "$a" DB is generated! Ready to annotate HMMs, COGs, and Pfams for "$a" contig database"
    anvi-run-hmms -c $CON/"$a".CONTIGS.db --num-threads 10
    anvi-run-ncbi-cogs -c $CON/"$a".CONTIGS.db --num-threads 10
    anvi-run-pfams -c $CON/"$a".CONTIGS.db --pfam-data-dir /PATH/TO/Pfam -T 10
    echo -e "yay! annotations are done for "$a", please proceed to profiling"
done


### Remap with fixed contigs and start profiling
module load htslib/v1.9 samtools/v1.9 bedtools/v2.26.0 bwa
# define paths for contigs, output, and metagenomic fastq files
FQ='path/to/fastq'
BAM='path/to/profiledbam'
CON='path/to/contigs'

find $FQ/*fq.gz > temp
sed 's/.fq.gz//g' temp > temp2
uniq temp2 > sample_list.txt
rm -f temp*
sample_list=$(cat sample_list.txt)

for a in $contig_list
  do
    for b in $sample_list
      do
        bwa mem -t 20 $CON/"$a".contigs-fixed.fa $FQ/"$b"_1.fastq.gz $FQ/"$b"_2.fastq.gz > $BAM/"$b"_aln_pe.sam
        samtools view -bS $CON/"$b"_aln_pe.sam > $CON/"$b"_aln_pe.bam
        rm $BAM/"$b"_aln_pe.sam #to much waste of space
        anvi-init-bam "$b"_aln_pe.bam -o "$b"_out.bam
        rm "$b"_aln_pe.bam #to much waste of space
        anvi-profile -i $CON/"$b"_out.bam -c $CON/"$a".CONTIGS.db -o $CON/"$a"/"$b"/ -T 10 --cluster-contigs --profile-SCVs -M 1000 --write-buffer-size 100
      done
done


### Merge profiles into a single profile database for each respective contig database
anvi-merge <profile_1.db> <profile_2.db> <profile_3.db> <profile_4.db> <profile_5.db> -o $CON/"$a"/MERGED_Profile --enforce-hierarchical-clustering -c $CON/"$a".CONTIGS.db # please see http://merenlab.org/2016/06/22/anvio-tutorial-v2/ to minimise confusion

for a in $contig_list
  do
    anvi-cluster-contigs -p $CON/"$a"/MERGED_Profile/PROFILE.db -c $CON/"$a".CONTIGS.db -C CONCOCT --driver concoct --just-do-it -T 40 --clusters 10
done


# import of taxonomy
#!/bin/sh
### KAIJU on HPC
module load kaiju/v1.5.0
makeDB.sh -e -t 20 # only do this once, it takes forever

for a in $contig_list
  do
    anvi-get-sequences-for-gene-calls -c $CON/"$a".CONTIGS.db -o $CON/"$a".gene_calls.fa
done

for a in $contig_list
  do
    kaiju -t /path/to/kaiju/nodes.dmp -f /path/to/kaiju/nr_euk/kaiju_db_nr_euk.fmi -i $CON/"$a".gene_calls.fa -o $CON/"$a".gene_calls_nr.out -z 8 -v
    addTaxonNames -t /path/to/kaiju/nodes.dmp -n /path/to/kaiju/names.dmp -i $CON/"$a".gene_calls_nr.out -o $CON/"$a".gene_calls_nr.names -r superkingdom,phylum,order,class,family,genus,species
    anvi-import-taxonomy-for-genes -i $CON/"$a".gene_calls_nr.names -c $CON/"$a".CONTIGS.db -p kaiju --just-do-it
done
#

## continue curation of mycoplasma MAGs, using anvio interface following guidelines: http://merenlab.org/2017/05/11/anvi-refine-by-veronika/
