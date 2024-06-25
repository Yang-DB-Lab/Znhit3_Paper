# Download primary assembly of mouse genome sequence
wget ftp://ftp.ensembl.org/pub/release-101/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz -O Mus_musculus.GRCm38.101.dna.primary_assembly.fa.gz
# unzip mouse genome sequence
gunzip Mus_musculus.GRCm38.101.dna.primary_assembly.fa.gz
# move the ERCC92.zip file to the folder for RNA_STAR index
# unzip the ERCC92.zip file to get ERCC92.fa and ERCC92.gtf

# add the ERCC92.fa to the end of mouse chromosome sequence
cat ERCC92.fa >> Mus_musculus.GRCm38.101.dna.primary_assembly.fa
# Change the name of the merged .fa file(reference genome plus ERCC)
mv Mus_musculus.GRCm38.101.dna.primary_assembly.fa Mus_musculus.GRCm38.101.dna.primary_assembly_plus_ERCC92.fa


# Download gtf file
wget ftp://ftp.ensembl.org/pub/release-101/gtf/mus_musculus/Mus_musculus.GRCm38.101.gtf.gz
gunzip Mus_musculus.GRCm38.101.gtf.gz

# Add ERCC92 annotation file to the end of mouse genome annotation gtf file
cat ERCC92.gtf >> Mus_musculus.GRCm38.101.gtf
# rename the merged gtf file
mv Mus_musculus.GRCm38.101.gtf Mus_musculus.GRCm38.101_plus_ERCC92.gtf

# Build the index for RNA-STAR 
# (For read >= 100bp, use sjdbOverhang 100 is good enough)
# Otherwise use SequencingLength-1 as sdjbOverhang
STAR --runMode genomeGenerate --runThreadN 8 \
--genomeDir /home/guang/mouse_genome_index/RNA_STAR_GRCm38.101_Overhang100_with_ERCC92 \
--genomeFastaFiles /home/guang/mouse_genome_index/RNA_STAR_GRCm38.101_Overhang100_with_ERCC92/Mus_musculus.GRCm38.101.dna.primary_assembly_plus_ERCC92.fa \
--sjdbGTFfile /home/guang/mouse_genome_index/RNA_STAR_GRCm38.101_Overhang100_with_ERCC92/Mus_musculus.GRCm38.101_plus_ERCC92.gtf \\
--sjdbOverhang 100
