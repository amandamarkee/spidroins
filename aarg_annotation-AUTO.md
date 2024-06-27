# Automated _spidroin_ annotation for _Argiope argentata_ 

Here, I describe an automated pipeline for pulling and annotating _spider fibroin_ genes. Contents in this repository were developed by myself, 
or modified from existing scripts created by Dr. Paul Frandsen and Joe Arguelles. The goal of this repository is to provide guidence and scripts
for the methods used in my first dissertation chapter, assessing the intraspecific allelic diversity in the fibroin gene family.

The goals for this pipeline are to automate the process of:
1) Identifying spidroins in long-read whole genome haplotype assemblies
2) Extracting all putitive spidroins per haplotype
3) Running _ab initio_ and trained gene prediction to determine intron/exon boundaries
4) Determining full length (bp), CDS length, and # of exons for spidroins from each haplotype across multiple populations.
  - 5 California: AargB1, AargC1, AargH1, AargF1 and AargG2 
  - 1 Texas: AargTX
  - 1 Florida: AargVK1


## Step 1) Identifying and Sorting _Spidroins_ with BlastX/BlastN & Excel


## Step 2) Extracting Putitive _Spidroins_ with Samtools and BioPython


## Step 3) _Ab initio_ or Trained Gene Prediction with Augustus 


## Step 4) Measuring Allelic Diversity with Geneious and Excel


