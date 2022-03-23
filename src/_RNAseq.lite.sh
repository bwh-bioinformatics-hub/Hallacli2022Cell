###########################################
# bash script for running paired-end RNAseq
# author: Xianjun Dong
# email: xdong@rics.bwh.harvard.edu
# date: 9/16/2013
# version: 2.0
# Note: call this script in the folder of fastq file
## TODO:
#1. use PINSEQ for QC
#2. use STAR for mapping
###########################################
#!/bin/bash

if [ $# -ne 2 ]
then
  echo "Usage: $0 HC_BN10-39_2.R1.fastq.gz HC_BN10-39_2.R2.fastq.gz"
  exit
fi

###########################################
echo "["`date`"] STEP 1. Configuring"
###########################################

modulename=`basename $0`
set +o posix  #  enables the execution of process substitution e.g. http://www.linuxjournal.com/content/shell-process-redirection

R1=$1  # filename of R1 
R2=$2  # filename of R2 (for paired-end reads)
samplename=${R1/[.|_]R1*/}

## $pipeline_path and $CONFIG_FILE already exported in the main pipeline.sh
config_file=./config.txt
[ -e $config_file ] || config_file=$HOME/neurogen/pipeline/RNAseq/config.txt
echo "Using configuration file at:" $config_file;
source $config_file

#strandoption="--library-type fr-unstranded"; split="-nosplit"; 
#[[ $samplename == *stranded* ]] && strandoption="--library-type fr-secondstrand"  # by default we use Illumina SMARTer stranded RNA-Seq kit
#[[ $samplename == *stranded* ]] && split="-split"

inputdir=$PWD
outputdir=$inputdir/../run_output
[ -d $outputdiri/$samplename ] || mkdir -p $outputdir/$samplename

# IMPORTANT: copy config.txt to the output file, in order to record the software/assembly/annotation version etc.
cp $config_file $outputdir/$samplename

###########################################
echo "["`date`"] STEP 2. quality filter: adaptor removal/clip"
###########################################
## TODO: Reads with PolyA/T tail were not trimmed. Using PRINSEQ instead for removal

##### adaptor removal
[ -d $inputdir/../filtered ] || mkdir $inputdir/../filtered

[ -e $inputdir/../filtered/adaptor.fa ] || echo -e ">Truseq_complementary_part_3p\nAGATCGGAAGAGC" > $inputdir/../filtered/adaptor.fa

# remove adaptor with fastq-mcf (https://code.google.com/archive/p/ea-utils/wikis/FastqMcf.wiki)
[ ! -f $outputdir/$samplename/.status.$modulename.adaptorremoval ] && \
fastq-mcf -o $inputdir/../filtered/$R1 -o $inputdir/../filtered/$R2 -t 0 -x 10 -l 15 -w 4 -q 10 -u $inputdir/../filtered/adaptor.fa <(zcat $R1) <(zcat $R2) && \
touch $outputdir/$samplename/.status.$modulename.adaptorremoval 

cd $inputdir/../filtered

#############################################
echo "["`date`"] STEP 3. QC"
############################################

[ ! -f $outputdir/$samplename/.status.$modulename.fastqc ] && \
fastqc --outdir $outputdir/$samplename --extract -t 2 $R1 $R2 && \
rm $outputdir/$samplename/*fastqc.zip && \
touch $outputdir/$samplename/.status.$modulename.fastqc

############################################
echo "["`date`"] STEP 4. mapping"
############################################
cd $inputdir/../filtered

[ ! -f $outputdir/$samplename/.status.$modulename.mapping ] && \
star  --genomeDir /data/bioinformatics/referenceGenome/hg38sacCer3_combined/starIndex \
      --runMode alignReads \
      --twopassMode Basic \
      --outFileNamePrefix $outputdir/$samplename/ \
      --readFilesCommand zcat \
      --readFilesIn $R1 $R2 \
      --outSAMtype BAM SortedByCoordinate \
      --outFilterType BySJout \
      --outFilterMultimapNmax $MAX_HIT \
      --outFilterMismatchNmax 999 \
      --outFilterMismatchNoverLmax 0.1 \
      --alignIntronMax 1000000 \
      --alignMatesGapMax 1000000 \
      --alignSJoverhangMin 8 \
      --alignSJDBoverhangMin 1 \
      --chimOutType Junctions \
      --chimSegmentMin 10 \
      --chimJunctionOverhangMin 15 \
      --runThreadN $CPU \
      --outSAMstrandField intronMotif \
      --outSAMunmapped Within && \
touch $outputdir/$samplename/.status.$modulename.mapping

cd $outputdir/$samplename
ln -fs Aligned.sortedByCoord.out.bam accepted_hits.bam

###########################################
echo "["`date`"] STEP 5. post-processing, format converting"
###########################################

cd $outputdir/$samplename

[ ! -f .status.$modulename.sam2bam ] && \
samtools index accepted_hits.bam && \
touch .status.$modulename.sam2bam

[ ! -f .status.$modulename.bam2stat ] && \
echo `samtools view -cF 0x100 accepted_hits.bam` "primary alignments (from samtools view -cF 0x100)" > accepted_hits.bam.stat && \
samtools flagstat accepted_hits.bam >> accepted_hits.bam.stat && \
touch .status.$modulename.bam2stat


# [ ! -f $outputdir/$samplename/.status.$modulename.bam2annotation ] && \
# _bam2annotation.sh accepted_hits.bam > accepted_hits.bam.bam2annotation && \
# Rscript $pipeline_path/modules/_bam2annotation.r accepted_hits.bam.bam2annotation accepted_hits.bam.bam2annotation.pdf && \
# touch $outputdir/$samplename/.status.$modulename.bam2annotation
# 
# [ ! -f $outputdir/$samplename/.status.$modulename.sam2bw ] && \
# echo "## normalization using total reads mapped to non_rRNA_mt" && \
# total_mapped_reads=`grep -w total_non_rRNA_mt accepted_hits.bam.bam2annotation | cut -f2 -d' '` && \
# bam2bigwig.sh accepted_hits.bam $split $total_mapped_reads && \
# touch $outputdir/$samplename/.status.$modulename.sam2bw

#rm accepted_hits.bed accepted_hits.*bedGraph

###########################################
echo "["`date`"] STEP 7. assembly and quantification"
###########################################

cd $outputdir/$samplename

echo "## run htseq for reads count"
[ ! -f .status.$modulename.htseqcount ] && \
htseq-count -m intersection-strict -t exon -i gene_id -s no -q -f bam accepted_hits.bam $ANNOTATION_GTF > hgseqcount.by.gene.tab 2> hgseqcount.by.gene.tab.stderr && \
touch .status.$modulename.htseqcount

echo "## run cufflinks to get FPKM"
# Note: "-b" option (for bias correction) can lead to segementation fault error.
[ ! -f .status.$modulename.cufflinks ] && \
cufflinks --no-update-check $strandoption -o ./ -p $CPU -G $ANNOTATION_GTF -M $MASK_GTF --compatible-hits-norm --multi-read-correct accepted_hits.bam && \
touch .status.$modulename.cufflinks

############################################
echo "["`date`"] STEP 9. prepare for tracks files to display on UCSC / IGV"
############################################
#
[ -d $inputdir/../for_display ] || mkdir $inputdir/../for_display
cd $inputdir/../for_display

### others
[ ! -f $outputdir/$samplename/.status.$modulename.makelinks ] && \
# make soft link
for i in $outputdir/$samplename/{accepted_hits.bam,accepted_hits.bam.bai,*tracking,hgseqcount.by.gene.tab}; do ii=${i/.*\//}; ln -fs $i $samplename.$ii; done && \

## QC
ln -fs $outputdir/$samplename/$samplename.R1_fastqc $samplename.R1_fastqc && \
ln -fs $outputdir/$samplename/$samplename.R2_fastqc $samplename.R2_fastqc && \

touch $outputdir/$samplename/.status.$modulename.makelinks

echo "["`date`"] DONE: $modulename job for sample $samplename is done !!"
