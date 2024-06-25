#!/bin/bash

##############################################################################
## This script align RNA-seq reads to annotated transposable elements (TE) 
## and get the expression of different TEs 
#############################################################################

# folder with all folders of fastq files
top_level_folder="/data/guang/raw_data"

######### Parameters related to STAR ###############
genome_STAR_TE_index_Dir="/home/guang/mouse_genome_index/RNA_STAR_GRCm38.101_with_ERCC92_Overhang100_for_TE"
####################################################

######### Parameter for StringTie #############
TE_annotation="/home/guang/mouse_genome_index/RNA_STAR_GRCm38.101_with_ERCC92_Overhang100_for_TE/GRCm38_rmsk_TE_plus_ERCC92.gtf"
###############################################
##******************* Parameter Definition Finished Here *********************##

##*************** The actual alignment script starts here. ********************##
##*************** Use pre-trimmed data; No need to do trimming ****************##
################################################################################

## Global variables used by all script pieces.
## The global variables will be re-defined in each part. Although this re-definition
## is not necessary for execution of this merged script, re-definition of these
## global variables will make each part still a complete script and can be copied
## out to run independently.
## Re-define of the global variables to make  this part independent.
## top_level_folder="/data/guang/raw_data"
## Initialize the top_level folder. Must use full path!!
## This folder should contain the sample folders with single-end fatq.gz files.
## This folder will be used by all script pieces below.
cd $top_level_folder
sample_folder_names="$(ls -l | grep "^d" | awk '{print $NF}')"
## Get names of sample folders:
## The command on the right site will first use 'ls -l' to check all files and folders
## then 'grep "^d" ' will select the ones will a 'd' property at the beginning of 
## 'ls -l' command which are folders
## finally 'awk '{print $NF}' will print out the last column of 'ls -l' which are 
## the actual folder names. 
## NOTE: There should be no space in the folder names. If space exists in folder names
##       only the last word of folder name will be print out. This is not what we want!
## sample_folder_names will be used by all script pieces below.



## Generate index for TEtranscripts; Run only ONCE!
STAR --sjdbOverhang 100 --sjdbGTFfile $TE_annotation \
     --runMode genomeGenerate --genomeDir $genome_STAR_TE_index_Dir \
     --genomeFastaFiles /home/guang/mouse_genome_index/RNA_STAR_GRCm38.101_with_ERCC92_Overhang100_for_TE/Mus_musculus.GRCm38.101.dna.primary_assembly_plus_ERCC92.fa \
     --runThreadN 6


cd $top_level_folder
sample_folder_names="$(ls -l | grep "^d" | awk '{print $NF}')"

for sample_folder in $sample_folder_names
do
	cd $top_level_folder/$sample_folder
	#get into sample_folder
	mkdir -p $top_level_folder/$sample_folder/$sample_folder\_TE_RNA_STAR_Results
	
	trimmed_fastq_gz_files="$(ls -d $PWD/$sample_folder\_Trimmomatic\_trimmed/*.fq.gz)"
	#$PWD is current path
	#This command get the full path of trimmed fast.gz files
	
	num_trimmed_fastq="$(echo "$trimmed_fastq_gz_files" | wc -l)"
	if [ $num_trimmed_fastq -ge 2 ] 
	# Check numer of trimmed fastq.gz files. If great than or equals to 2, they must come from 
	# paired-end read files.
	# Treat as trimmed fastq files from paired-end fastq files
	then
	paired_trimmed_reads="$(ls -d $PWD/$sample_folder\_Trimmomatic_trimmed/*_paired.*)"
	# use paired_reads to store the filenames of the 2 paired-end read files
	unpaired_trimmed_reads="$(ls -d $PWD/$sample_folder\_Trimmomatic_trimmed/*_unpaired.*)"
	# use unpaired_reads to store the filenames of the 2 unpaired read files trimmed out by trimmomatic
	
	# Only align pair-end part
	STAR --runThreadN 6 \
    --genomeDir $genome_STAR_TE_index_Dir \
    --sjdbGTFfile $TE_annotation \
    --sjdbOverhang 100 \
    --readFilesIn $paired_trimmed_reads \
    --readFilesCommand zcat \
    --outSAMtype BAM Unsorted \
    --winAnchorMultimapNmax 200 --outFilterMultimapNmax 100 \
    --outReadsUnmapped Fastx \
    --outFileNamePrefix ./$sample_folder\_TE_RNA_STAR_Results/STAR_alignment\_$sample_folder\_paired_
	
	
	# index the big merged and sorted bam file
	# samtools sort STAR_alignment\_$sample_folder\_paired_Aligned.sortedByCoord.out.bam
	samtools sort -@ 6 -m 4G -o ./$sample_folder\_TE_RNA_STAR_Results/STAR_alignment\_$sample_folder\_paired_Aligned.out.sorted.bam \
                     ./$sample_folder\_TE_RNA_STAR_Results/STAR_alignment\_$sample_folder\_paired_Aligned.out.bam
	fi
done
##****************************Part2 Finished here *******************************##

##******** Part 3.1 Using featureCounts to get TE expression table from .bam and .gtf files **********##

## Re-define of the global variables to make Part3.1 independent script.
## top_level_folder="/data/guang/raw_data"
## Initialize the top_level folder. Must use full path!!
## This folder should contain the sample folders with single-end fatq.gz files.
cd $top_level_folder
sample_folder_names="$(ls -l | grep "^d" | awk '{print $NF}')"
## Get names of sample folders
## The command on the right site will first use 'ls -l' to check all files and folders
## then 'grep "^d" ' will select the ones will a 'd' property at the beginning of 
## 'ls -l' command which are folders
## finally 'awk '{print $NF}' will print out the last column of 'ls -l' which are 
## the actual folder names. 
## NOTE: There should be no space in the folder names. If space exists in folder names
##       only the last word of folder name will be print out. This is not what we want!

for sample_folder in $sample_folder_names
do
	cd $top_level_folder/$sample_folder
	#get into sample_folder
	mkdir -p $top_level_folder/$sample_folder/$sample_folder\_TE_featureCounts_Results
	featureCounts -M --fraction -F GTF -T 6 -s 0 -p -a $TE_annotation \
	-o $top_level_folder/$sample_folder/$sample_folder\_TE_featureCounts_Results/$sample_folder\_outfeatureCounts.txt \
	./$sample_folder\_TE_RNA_STAR_Results/STAR_alignment\_$sample_folder\_paired_Aligned.out.sorted.bam

done
