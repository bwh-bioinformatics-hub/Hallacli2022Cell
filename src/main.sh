####################################
# Pipeline for RNAseq data Analysis
# Authors: Xianjun Dong
# Email: xdong@rics.bwh.harvard.edu
# Date: 10/6/2014
# Version: 2.0
####################################
#!/bin/bash

########################
## 0. setting
########################
# project folders
input_dir=/data/bioinformatics/projects/vik2021/rawfiles

config_file=$input_dir/config.txt
[ -e "$config_file" ] || cp $HOME/neurogen/pipeline/RNAseq/config.txt $config_file
echo "Using configuration file at:" $config_file;
source $config_file

# create the subfolders (e.g. filtered, run_output, for_display, results)
filtered_dir=$input_dir/../filtered
[ -d $filtered_dir ] || mkdir $filtered_dir

output_dir=$input_dir/../run_output
[ -d $output_dir ] || mkdir $output_dir

fordisplay_dir=$input_dir/../for_display
[ -d $fordisplay_dir ] || mkdir $fordisplay_dir

result_dir=$input_dir/../results
[ -d $result_dir ] || mkdir $result_dir

####################################
# data clean-up
####################################
cd $input_dir;
wget --no-check-certificate -qO - "https://docs.google.com/spreadsheets/d/e/2PACX-1vQUv4qshWcbOyqC28RfRhSXtj_2OOeoL3HOkf5sw1WB7o24WcQc2DkyrWnxLaa21qHlhtJL2JR-rK4n/pub?gid=2020811866&single=true&output=tsv" | cut -f1,8 | dos2unix | awk 'NR>1{print "ln -fs",$2, $1".R1.fastq.gz"; sub("R1_001","R2_001",$2); print "ln -fs",$2, $1".R2.fastq.gz";}' | bash

####################################
# run pipeline per sample
####################################
cd $input_dir;
for i in *PFFyes*.R1.fastq.gz;
do
    R1=$i
    R2=${i/R1/R2};
    samplename=${R1/.R1*/}

    [ -e $i ] || continue
    
    # run the QC/mapping/assembly/quantification for RNAseq
    bsub -J $samplename -oo $output_dir/$samplename/_RNAseq.log -eo $output_dir/$samplename/_RNAseq.log -q $QUEUE -n $CPU -M $MEMORY -R rusage[mem=$MEMORY] -R "select[hname!=cmu066]" -u $EMAIL -N ~/projects/vik2021/src/_RNAseq.lite.sh $R1 $R2;

done


########################
## merge all samples to get big matrix for expression (e.g. one row per gene/Tx, one col per sample)
########################
[ -d $result_dir/merged ] || mkdir $result_dir/merged
cd $result_dir/merged

# uniq mapper
Rscript $HOME/neurogen/pipeline/RNAseq/modules/_mergeSamples.R `ls ../../run_output/*/genes.fpkm_tracking` genes.fpkm.cufflinks.allSamples.xls
Rscript $HOME/neurogen/pipeline/RNAseq/modules/_mergeSamples_htseq.R `ls ../../run_output/*/hgseqcount.by.gene.tab` genes.htseqcount.cufflinks.allSamples.xls

#########################
### QC: detect possible outlier
#########################
cd $result_dir/merged
Rscript ~/projects/vik2021/src/_normQC.R 

# ## debug: map only to yeast to see if geting a simiarl # of reads
# cd /data/bioinformatics/projects/vik2021/filtered
# mkdir ../run_output/SNCA4_PFFyes_h6_Rep3/yeastOnly/; 
# module load star/2.7.3
# bsub -q big-multi -n 4 -M 4000 -R rusage[mem=4000] -J test_yeast_only star --genomeDir /data/bioinformatics/referenceGenome/Saccharomyces_cerevisiae/Ensembl/release103/starIndex --runMode alignReads --twopassMode Basic --outFileNamePrefix ../run_output/SNCA4_PFFyes_h6_Rep3/yeastOnly/ --readFilesCommand zcat --readFilesIn SNCA4_PFFyes_h6_Rep3.R1.fastq.gz SNCA4_PFFyes_h6_Rep3.R2.fastq.gz --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --outFilterMultimapNmax 100 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --chimOutType Junctions --chimSegmentMin 10  --chimJunctionOverhangMin 15 --runThreadN 4 --outSAMstrandField intronMotif --outSAMunmapped Within

####################################
# DE analysis
####################################
cd $result_dir/merged
Rscript ~/projects/vik2021/src/_DE2.R