# Workflow for pulling major ampullate (MaSp) spidroins from _Argiope argentata_ (October 11, 2023)

## Cluster Setup and File Pathways
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
1) Rick dropped the five genomes, and the termini database in my working directory.
2) The database file contains only terminal regions (150bp for n-termini, and 100bp for c-termini). No repetitive domains
3) Because of the repetitive nature of spidroins, be sure to increase the allowed number of hits (~200 hits) for a first pass
4) Once the blastx results come back, we will sort, assess and assign best seq hit per MaSp gene in Excel
5) The final list will contain close to 60 hits (30 n-termini, 30 c-termini, for ~30 spidroins, masp and aciniforms)
6) Pair up which n-termini go with each c-termini, whould be little overlap with approx. 25 Kbp in between each terminus
7) Extract each MaSp gene out of the assembly (~30) after we identify where they are. (Note: will be done in Geneious)


## Step 1) BLAST Argiope Genome x Termini Database
Here, I'll take the whole assembled genomes and blast against Rick's termini database file. The goal is to identify complete spidroins present 
in each genome, and use the conserved terminal regions to find them.

Navigate to the working directory
```
cd /home/amarkee/nas4/aargentata_genome/blastx/aarg_pbTX
module load blast-2.2.26
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
