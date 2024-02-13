# _Spidroin_ annotation in the silver garden spider _Argiope argentata_ 

## Cluster Setup and File Pathways | (October 11, 2023)
For this workflow, I'll be working in my own huxley workspace in the following directory:
```
/home/amarkee/nas4/aargentata_genome
```

I will use the following pre-existing files for _A. argentata_ to work with:
```
Aarg_ON_FCI_Polca.fa # California assembled whole-genome from Oxford Nanopore sequencing
AargTX_ON_FlyeCorrI.fa # Texas assembled whole-genome from Oxford Nanopore sequencing

AargTX_hifiasm_I.fasta # Texas assembled whole-genome from PacBio sequencing
AargTX_hifiasm_hap1.fasta # Texas assembled haplotig-1 from PacBio sequencing
AargTX_hifiasm_hap2.fasta # Texas assembled haplotig-2 from PacBio sequencing 

Orbic_NandCtermini_proII_db.fasta # Termini sequence database for blastx
```

## Notes Before Starting

1) Rick dropped the five genomes, and the termini database in my huxley working directory
   - Note: Paul Frandsen BYU sequenced five additional PacBio genomes (with three others unassembled, recieved 2/9/2024 BYU) to be included in this analysis at a later date
   - Completed assemblies include: _A. lobata, A.australis, A.argentata H1, A. argentata F1,_ and _A.argentata C1._
   - Unassembled genomes include: _A. aurantia PAR3_, _A. argentata G2_ and _A. argentata B1_
3) The database file contains only terminal regions (150bp for n-termini, and 100bp for c-termini). No repetitive domains
4) Because of the repetitive nature of spidroins, be sure to increase the allowed number of hits (~200 hits) for a first pass
5) Once the blastx results come back, we will sort, assess and assign best seq hit per gene in Excel
6) The final list will contain close to 60 hits per genome (30 n-termini, 30 c-termini, for ~30 spidroins, masp and aciniforms)
7) Pair up best hits for matching n-termini / c-termini, whould be little overlap with approx. 25 Kbp in between each termini
8) Extract each gene out of the assembly (~30) after we identify where they are, including flanking regions (~5000bp up/downstream
   -  Note: This step will be done in Geneious
9) Use minimap or hisat2 (with splice unaware mapping enabled) to map RNAseq data back to genes. This determines intron/exon borders


## Step 1) BLAST Argiope Genome x Termini Database
Here, I'll take the whole assembled genomes and blast against Rick's termini database file using blastx. The goal is to identify complete spidroins present 
in each genome, and use the conserved terminal regions to find them.

Navigate to the working directory
```
cd /home/amarkee/nas4/aargentata_genome/blastx/aarg_pbTX # directory for blastx results for PacBio Texas A.arg genome
module load ncbi-blast-2.12.0+
```

Create a protein database using the termini database file Rick sent:
```
makeblastdb -in <database>.fa -dbtype prot
makeblastdb -in /home/amarkee/nas4/aargentata_genome/Orbic_NandCtermini_proII_db.fasta -dbtype prot

##########################

Output:
Building a new DB, current time: 10/11/2023 13:18:42
New DB name:   /home/amarkee/nas4/aargentata_genome/Orbic_NandCtermini_proII_db.fasta
New DB title:  /home/amarkee/nas4/aargentata_genome/Orbic_NandCtermini_proII_db.fasta
Sequence type: Protein
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 617 sequences in 0.118198 seconds.

```

Run blastx command to blast Texas PacBio genome against the generated db:
```
blastx -query <seqs.fasta> -db /apps/blast/databases/nr -outfmt 5 -max_target_seqs 200 -evalue 1e-5 -out xml/results.xml
blastx -query /home/amarkee/nas4/aargentata_genome/AargTX_hifiasm_I.fasta -db /home/amarkee/nas4/aargentata_genome/Orbic_NandCtermini_proII_db.fasta -outfmt 5 -max_target_seqs 200 -evalue 1e-5 -out /home/amarkee/nas4/aargentata_genome/blastx/aarg_pbTX/results.xml
```

Note: Because I'll be running this for the five genomes (2 full-length Oxford Nanopore, 2 full-length PacBio Haplotypes, and 1 full-length PacBio Primary), I edited the 
following script for each input genome fasta file:
```
#PBS -V
#PBS -q batch 
#PBS -l select=1:ncpus=25
#PBS -o /home/amarkee/nas4/aargentata_genome/blastx/aarg_pbTX
#PBS -e /home/amarkee/nas4/aargentata_genome/blastx/aarg_pbTX
#PBS -M amarkee@amnh.org
#PBS -m abe
#PBS -N blastx_aarg
#PBS -l walltime=999:00:00

module load ncbi-blast-2.12.0+

blastx -query /home/amarkee/nas4/aargentata_genome/AargTX_hifiasm_hap2.fasta -db /home/amarkee/nas4/aargentata_genome/Orbic_NandCtermini_proII_db.fasta -outfmt 6 -max_target_seqs 200 -evalue 1e-5 -out /home/amarkee/nas4/aargentata_genome/blastx/aarg_pbTX/hap2_results.txt
```

I tried to convert this below as a script 'universal_blastx.sh' that uses environmental variables, but have not tested it yet since I already submit all 5 scripts:
```
qsub universal_blastx.sh <input.fa> <protein_db.fa> <output.file>
```

```
#!/bin/bash
#PBS -V
#PBS -q batch 
#PBS -l select=1:ncpus=25
#PBS -o $PBS_OUT_DIR #manually change out dir
#PBS -e $PBS_ERR_DIR #manually change error dir
#PBS -M $PBS_EMAIL
#PBS -m abe
#PBS -N blastx_aarg
#PBS -l walltime=9999:00:00

module load ncbi-blast-2.12.0+

inputfasta=${1}
blastdb=${2}
outputfile=${3}


blastx -query ${inputfasta} \
-db ${blastdb} \
-outfmt 6 \
-max_target_seqs 200 \
-evalue 1e-5 \
-out ${outputfile}
```

## Step 2) Complete first pass: Determine best hits for each gene

After running blastx, output format -6 will produce a textfile that contains the hits in the following format:
```
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
# qseqid = query, or contig id
# seqid = db hit id
# qstart/qend = coordinates in genome
# evalue = significance threshold for match happening by chance. similar to p-value
```

There will be maximum 200 hits for each unique qstart position, and now we must manually assign the best hits for
both N and C termini for each gene. We determine the "gene" by the approximate length between chunks of hits, and 
typically use the following clues to determine if the paired N/C hits are a gene:
  1) the N/C terms are ~15,000 - 90,000 bp apart. we will look closer in Geneious to confirm in second pass
  2) the database ID (e.g. the positive hit ID) will be from the same/closely related species
  3) the database ID (e.g. the positive hit ID) will have the same gene ID, such as MaSp, Flag, MiSp, etc.
  4) the database ID (e.g. the positive hit ID) will end with the matching _C or _N in the naming convention
    - Note: there is no particular order in which the N/C term will appear, because they can be in any direction
     When assigning best hit, just take the first best hit regardless of term identity, and look for it's matching N/C

Here is an example of results from blastx, uncurated but sorted by qseqid (A-Z), qstart (smallest - largest), then e-value (smallest-largest).
This ensures that most of the time, the best hit and most accurate hit will be first in the group of similar qstart values, and that the hits
are organized by genes occuring across a contig, by the ascending start position:

<img width="982" alt="Screenshot 2024-02-13 at 1 44 22 PM" src="https://github.com/amandamarkee/spidroins/assets/56971761/fb258996-9e5e-406b-b063-01b653a98ec6">

Note: In the above example, the Aaur_TuSp_1238.2_N and "_C pair are the best hit due to them having the lowest e-value, belonging to the same genus in 
the sister species _Argiope aurantia_ and reasonable gene length from each terminal (Aaur_TuSp_1238.2_N at 4,715,052 bp, and Aaur_TuSp_1238.1_C at 4,752,200 bp)
making the gene ~37,148 bp long and within reasonable range. 





