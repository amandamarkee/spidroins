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
- Run initial BLAST on each haplotype genome vs. _Argiope argentata_ termini database
- Sorting steps from Joe in Excel
- Formula sheet for determining spidroin full length, notes on MiSp(rev) and MaSp2a-e, and including flanking region in list file.
- Safe as xlsx, and create a list text file called "coordinates.txt" with the following information in the following order:
    - contig, gene start coordinate, gene end coordinate, gene name

## Step 2) Extracting Putitive _Spidroins_ with Samtools and BioPython
- Save list file from Step 1 to working directory
- Run extract_sequences.py to pull coordinates and rename with gene name appended
- Saves to a concattonated spidroin 
Save the following script as "extract_sequences.py"   
```
import sys
from Bio import SeqIO

def read_coordinates(file):
    coordinates = []
    with open(file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) != 4:
                continue
            contig, start, end, gene = parts
            coordinates.append((contig, int(start), int(end), gene))
    return coordinates

def fetch_sequences(reference_file, coordinates):
    sequences = SeqIO.to_dict(SeqIO.parse(reference_file, "fasta"))
    fetched_sequences = []

    for contig, start, end, gene in coordinates:
        if contig in sequences:
            seq_record = sequences[contig]
            subseq = seq_record.seq[start-1:end]  # Adjust for 0-based index
            subseq_id = "{}_{}_{}_{}".format(seq_record.id, start, end, gene)
            subseq_record = seq_record[start-1:end]
            subseq_record.id = subseq_id
            subseq_record.description = ""
            fetched_sequences.append(subseq_record)
        else:
            print("Contig {} not found in reference genome.".format(contig))

    return fetched_sequences

def write_sequences(sequences, output_file):
    with open(output_file, "w") as f:
        SeqIO.write(sequences, f, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_sequences.py <coordinates_file> <reference_genome_file> <output_file>")
        sys.exit(1)

    coordinates_file = sys.argv[1]
    reference_genome_file = sys.argv[2]
    output_file = sys.argv[3]

    coordinates = read_coordinates(coordinates_file)
    sequences = fetch_sequences(reference_genome_file, coordinates)
    write_sequences(sequences, output_file)
    print("Sequences have been written to {}".format(output_file))
```

Bash script for running extract_sequences.py. Save this script as "pull.sh"
```
#PBS -q batch
#PBS -m abe
#PBS -M amarkee@amnh.org
#PBS -l select=1:ncpus=30
#PBS -l walltime=4:00:00
#PBS -o /home/amarkee/nas4/aargentata_genome/automated_pull/b1h1/bits.olog
#PBS -e /home/amarkee/nas4/aargentata_genome/automated_pull/b1h1/bits.elog
#PBS -N AargB1H1_automated

# Load any necessary modules (e.g., Python), Adjust according to your environment
conda activate spidroins
module load python-3.6.1

# Set the file paths
COORDINATES_FILE="/home/amarkee/nas4/aargentata_genome/automated_pull/b1h1/coordinates.txt"
REFERENCE_GENOME_FILE="/home/amarkee/nas4/aargentata_genome/assemblies/haplotype_asm/aarg_b1_assembly.asm.bp.hap1.p_ctg.fa"
OUTPUT_FILE="/home/amarkee/nas4/aargentata_genome/automated_pull/b1h1/aargb1h1_bits.fa"

# Execute the Python script with the specified file paths
python /home/amarkee/nas4/aargentata_genome/automated_pull/extract_sequences.py $COORDINATES_FILE $REFERENCE_GENOME_FILE $OUTPUT_FILE
```

## Step 3) _Ab initio_ or Trained Gene Prediction with Augustus 

_Ab initio_ gene prediction step:
```
#!/bin/bash
#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4G   # memory per CPU core
#SBATCH -J "augustus"   # job name
#SBATCH --mail-user=amarkee@amnh.org   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

species=$1
gene="${species}_gene.fa"
output="${species}.gff"

module load miniconda3/4.12-pws-472
conda activate busco

export AUGUSTUS_CONFIG_PATH=~/fsl_groups/fslg_pws472/apps/miniconda3/envs/busco/config/
augustus --strand=both --singlestrand=true \
--extrinsicCfgFile=./extrinsic.cfg \
--alternatives-from-evidence=true \
--gff3=on \
--uniqueGeneId=true \
--UTR=off \
--species=fly \
$gene > \
$output
```

```
sbatch augustus.sh aarg_g2
```

Trained gene prediction (RNA) – to be tested
```
#!/bin/bash
#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4G   # memory per CPU core
#SBATCH -J "augustus"   # job name
#SBATCH --mail-user=amarkee@amnh.org   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

species=$1
gene="${species}_gene.fa"
rna_seq_dir="/home/amarkee/nas4/aargentata_genome/rna_seq"  # Directory containing RNA-seq fasta files
output="${species}.gff"

module load miniconda3/4.12-pws-472
conda activate busco

export AUGUSTUS_CONFIG_PATH=~/fsl_groups/fslg_pws472/apps/miniconda3/envs/busco/config/

# Concatenate all RNA-seq fasta files into one file
cat ${rna_seq_dir}/*.fa > combined_rna_seq.fa

# Train Augustus using the combined RNA-seq data
autoAug.pl --species=$species --useexisting --singleCPU --pasa --genome=$gene --cdna=combined_rna_seq.fa

# Run Augustus with the new trained parameters
augustus --strand=both --singlestrand=true \
--extrinsicCfgFile=./extrinsic.cfg \
--alternatives-from-evidence=true \
--gff3=on \
--uniqueGeneId=true \
--UTR=off \
--species=$species \
$gene > \
$output
```

```
sbatch augustus.sh aarg_g2
```

## Step 4) Measuring Allelic Diversity with Geneious and Excel

First pass annotation for finding intron exon boundaries.

Example of MaSp2.2c output from _ab initio_ gene prediction with Augustus. 

<img width="1727" alt="Screenshot 2024-06-28 at 10 35 17 AM" src="https://github.com/amandamarkee/spidroins/assets/56971761/c43bde70-51e1-4ad3-bd3c-71fca360b94d">

Above, we import both the concatonated spidroin file (FASTA of all 35-40 spidroins per haplotype assembly), and the .gtf annotation that gets output from Augustus. Just checking by eye, prior to confirming with RNAseq mapped data, the structure of MaSp with homogenized introns looks to be in line with what we expect the gene structure to look like.


Example of Flag output from _ab initio_ gene prediction with Augustus.

<img width="1728" alt="Screenshot 2024-06-28 at 10 34 54 AM" src="https://github.com/amandamarkee/spidroins/assets/56971761/c0f3a34e-7ab7-4704-a1ef-020299c00ff3">

Above, we again import the same files, and check the structure to confirm the intro/exon structure is in line with what has been seen in annotations done by hand. 

It's important to note that Augustus _ab initio_ gene prediction should serve as a first step, and has been known to sometimes include/exclude additional introns that may not be present. It will be important to prune the _ab initio_ annotations carefully when measuring allelic diversity through gene length, CDS region length, and number of exons. 

To try and improve this first pass annotation, we will also be trying to train the Augustus annotation with RNAseq data, or a concatonated CDS fasta for all known spidroins.
