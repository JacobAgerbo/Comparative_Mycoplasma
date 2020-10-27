
###-----------------------------------------------------------------###
###						     HappyFish Metagenomic Pipeline				            ###
###-----------------------------------------------------------------###

# Credits: Jacob Agerbo Rasmussen
# Contact: Genomicsisawesome@gmail.com
# Version: 1.0.5

-----------------------------------------------#
#							Initial QC							  #
#-----------------------------------------------------------------#
module load java/1.8.0  fastqc/0.11.8
mkdir FastQC_initial
cd FastQC_initial

find *.fq.gz > list
sed 's/.fq.gz//g' list > list2
uniq list2 > sample_list
rm -f list*
sample_list=$(cat sample_list)
echo ${sample_list[@]}


module load java/1.8.0  fastqc/0.11.8
for a in $sample_list
  do
    fastqc -o out/ "$a"
done    

### Create short report for all samples
unzip '*fastqc.zip'
FASTA_DIR='<path/to/FASTQ/>'
WORK_DIR='<path/to/WD/>'

cd $FASTA_DIR/
find *_fastqc.zip > temp
sed 's/_fastqc.zip//g' temp > temp2
uniq temp2 > sample_list.txt
rm -f temp*
sample_list=$(cat sample_list.txt)
cd $WORK_DIR
for a in $sample_list
  do
  cd "$a"_fastqc
    total_seqs=`cat fastqc_data.txt | grep 'Total Sequences' | cut -f 2`
    gc_percent=`cat fastqc_data.txt | grep '%GC' | cut -f 2`
    seq_length=`cat fastqc_data.txt | grep -A1 '#Length' | tail -n +2 | cut -f 1`
    seq_qual=`cat fastqc_data.txt | awk '/>>Per base sequence quality/,/>>END_MODULE/' | tail -n +3 | head -n -1 | awk '{total+=$2} END {print total/NR}'`
    n_count=`cat fastqc_data.txt | awk '/>>Per base N content/,/>>END_MODULE/' | tail -n +3 | head -n -1 | awk '{total+=$2} END {print total/NR}'`
    echo -e "File Name:\t"$a"\nNumber of Sequences:\t${total_seqs}\nGC%:\t${gc_percent}\nSequence Length:\t${seq_length}\nAverage per base sequence quality:\t${seq_qual}\nN%\t${n_count}" > ../"$a"_short.txt
    cd ..
    rm -r "$a"_fastqc
done
cat *.txt > initial_report.txt
rm *_short.txt

#-----------------------------------------------------------------#
#			          	Pre-processing, without collapse		       		  #
#-----------------------------------------------------------------#
### Trim adapters
mkdir 1-trimmed
cd 1-trimmed

# Load all required modules for the job
module load adapterremoval/2.2.4
FASTA_DIR='<path/to/FASTQ/>'
WORK_DIR='<path/to/WD/>'

cd $FASTA_DIR/
for a in $sample_list
  do
    AdapterRemoval --file1 ../"$a"_read_1.fq.gz --file2 ../"$a"_read_2.fq.gz --basename "$a" --output1 "$a"_filtered_1.fq.gz --output2 "$a"_filtered_2.fq.gz --adapter1 AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA --adapter2 GAACGACATGGCTACGATCCGACTT --qualitybase 30 --minlength 50 --threads 20 --gzip
done

#-----------------------------------------------------------------#
#														Remove duplicates										  #
#-----------------------------------------------------------------#

module load pigz/2.3.4 seqkit/0.7.1
FASTA_DIR='<path/to/FASTQ/>'
WORK_DIR='<path/to/WD/>'
cd $FASTA_DIR/
find *.fq > temp
sed 's/_filtered.fq//g' temp > temp2
uniq temp2 > sample_list.txt
rm -f temp*
sample_list=$(cat sample_list.txt)
cd $WORK_DIR

for a in $sample_list
do
cat $FASTA_DIR/"$a"_filtered_1.fq.gz| seqkit rmdup -j 20 -s -o $WORK_DIR//"$a"_1_clean.fastq
cat $FASTA_DIR/"$a"_filtered_1.fq.gz| seqkit rmdup -j 20 -s -o $WORK_DIR//"$a"_1_clean.fastq
done


#-----------------------------------------------------------------#
#		         		Re-pair fastq	after rm duplicates    	         	  #
#-----------------------------------------------------------------#

module load jre/1.8.0 bbmap/38.35
FASTA_DIR='<path/to/FASTQ/>'
WORK_DIR='<path/to/WD/>'
cd $FASTA_DIR/
find *.fq > temp
sed 's/_[1-2]_clean.fq//g' temp > temp2
uniq temp2 > sample_list.txt
rm -f temp*
sample_list=$(cat sample_list.txt)
cd $WORK_DIR

for a in $sample_list
do
repair.sh in=$FASTA_DIR/"$a"_1_clean.fq in2=$FASTA_DIR/"$a"_2_clean.fq out=$WORK_DIR/"$a"_1_repaired.fq out2=$WORK_DIR/"$a"_2_repaired.fq outsingle=$WORK_DIR/Singletons/"$a"singletons.fq overwrite=t
done

#-----------------------------------------------------------------#
#		         		Filtering of Host DNA, using BWA	   		      	  #
#-----------------------------------------------------------------#

#	Filter against:
#	phiX174     Remove phiX174, which were spiked in during sequencing
#	Human       HG19, i might have been a messy head
# Host DNA, using Omyk_1.0 reference genome:		https://www.ncbi.nlm.nih.gov/assembly/GCF_002163495.1/

mkdir 1-Host_Removal
cd 1-Host_removal

# 1) Build BWA Host_DB
module load samtools/1.9 bedtools/2.28.0 bwa/0.7.15
REF_DIR='<path/to/REF/>'
FASTA_DIR='<path/to/FASTQ/>'
WORK_DIR='<path/to/WD/>'

cd $FASTA_DIR/
find *.fq > temp
sed 's/_[1-2]_repaired.fq//g' temp > temp2
uniq temp2 > sample_list.txt
rm -f temp*
sample_list=$(cat sample_list.txt)
cd $WORK_DIR

#phiX174 filtering
Ref='phiX174'
for a in $sample_list
  do
    bwa mem -t 28 $REF_DIR/$Ref $FASTA_DIR/"$a"_1_repaired.fq $FASTA_DIR/"$a"_2_repaired.fq | samtools view -f 12 -T $REF_DIR/phiX174 -b > $WORK_DIR/noPhiX/"$a"_noPhiX.bam
    bedtools bamtofastq -i $WORK_DIR/noPhiX/"$a"_noPhiX.bam -fq $WORK_DIR/noPhiX/"$a"_noPhiX_1.fq -fq2 $WORK_DIR/noPhiX/"$a"_noPhiX_2.fq
done

####  Human filtering
REF_DIR='<path/to/REF/>'
FASTA_DIR='<path/to/FASTQ/>'
WORK_DIR='<path/to/WD/>'
cd $FASTA_DIR/
find *.fq > temp
sed 's/_noPhiX_[1-2].fq//g' temp > temp2
uniq temp2 > sample_list.txt
rm -f temp*
sample_list=$(cat sample_list.txt)
cd $WORK_DIR

Ref='human.fna' # human genome HG19
for a in $sample_list
  do
    bwa mem -t 20 $REF_DIR/human.fna $FASTA_DIR/"$a"_noPhiX_1.fq $FASTA_DIR/"$a"_noPhiX_2.fq | samtools view -f 12 -T $REF_DIR/human.fna -b > $WORK_DIR/noHuman/"$a"_noHuman.bam
    bedtools bamtofastq -i $WORK_DIR/noHuman/"$a"_noHuman.bam -fq $WORK_DIR/noHuman/"$a"_noHuman_1.fq -fq2 $WORK_DIR/noHuman/"$a"_noHuman_2.fq
done

##### Filter trout
REF_DIR='<path/to/REF/>'
FASTA_DIR='<path/to/FASTQ/>'
WORK_DIR='<path/to/WD/>'

cd $FASTA_DIR/
find *.fq > temp
sed 's/_noPhiX_[1-2].fq//g' temp > temp2
uniq temp2 > sample_list.txt
rm -f temp*
sample_list=$(cat sample_list.txt)
cd $WORK_DIR

Ref='Omyk_1.0_genome.fna' # change to other

for a in $sample_list
  do
    bwa mem -t 28 $REF_DIR/$Ref $FASTA_DIR/"$a"_noPhiX_1.fq $FASTA_DIR/"$a"_noPhiX_2.fq | samtools view -f 13 -T $REF_DIR/human.fna -S -b > $WORK_DIR/noHost/"$a"_noHost.bam
    bedtools bamtofastq -i $WORK_DIR/noHost/"$a"_noHost.bam -fq $WORK_DIR/noHost/"$a"_noHost_1.fq -fq2 $WORK_DIR/noHost/"$a"_noHost_2.fq
done

#-------------------------------------------------------------------#
#			               		Co-Assembly, using MegaHit			      			#
#-------------------------------------------------------------------#

# I'm using MegaHit with metagenomic presets. Keep the minimum length at 1000, else we get fucked up binning.
module load megahit/1.1.1 anaconda2/4.4.0
FASTA_DIR='<path/to/FASTQ/>'
WORK_DIR='<path/to/WD/>'
FASTA=$(ls $FASTA_DIR/*.fq | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])')
# Check the FASTA list
echo $FASTA
# run megahit
megahit -r $FASTA --min-contig-len 1000 -t 28 --presets meta-sensitive -o $WORK_DIR

#-------------------------------------------------------------------#
#														Assembly QC															#
#-------------------------------------------------------------------#

# Load all required modules for the job
module load anaconda2/4.4.0 quast/5.0.2
for a in final
  do
    quast.py -o Quast_Report/ $WORK_DIR/"$a".contigs.fa
done
