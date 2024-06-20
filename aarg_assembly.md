# Whole genome assembly of _Argiope argentata_ (PacBio HiFi)

## **Background**

The goal of this workflow is to take note of assembly details for the five California _Argiope argentata_ genomes, one Florida _Argiope argentata_ genome, and one Texas _Argiope argentata_ genome. Once whole genomes are assembled, each pair of haplotype assemblies will be used to assess intraspecific (allelic) variation of silk genes within the spidroin gene family. I will be trying to run the genome assembler (hifiasm) on both American Museum of Natural History high-performance computing systems (Mendel and Huxley), and Brigham Young University's computing system (Fulton).
- California samples (assembled by PBF): Aarg_C1, Aarg_F1, and Aarg_H1
- California samples (assembled by AM): Aarg_B1 and Aarg_G2
- Texas sample (assembled by AM): AargTX
- Florida sample (assembled by AM): AargVK1

## **Project Set-Up and File Conversion**

Prior to running preliminary QC and genome assembly, I will convert my Revio output (.bam) into usable formats (.fastq) for both programs. 
```
# raw data directory on huxley
/home/amarkee/nas5/Aarg_B1/hifi_reads
/home/amarkee/nas5/Aarg_G2/hifi_reads

# raw data directory on fulton
/home/fslcollab384/compute/genomics_workshop/aarg_g2
/home/fslcollab384/compute/genomics_workshop/aarg_vk1

# working directory for all haplotype assemblies on huxley
/home/amarkee/nas4/aargentata_genome/assemblies/haplotype_asm

# working directory for all primary assemblies on huxley
/home/amarkee/nas4/aargentata_genome/assemblies/primary_asm
```

I use the following script on Huxley to convert .bam to .fastq files using bedtools, changing file path locations for each of the 5 genomes. We do this because the PacBio assembler hifiasm requires a fastq file for input. Note: in the future, I'll transform this into a cleaner environmental variable script:
```
#!/bin/bash
#PBS -V
#PBS -q batch
#PBS -l select=1:ncpus=25
#PBS -o /home/amarkee/nas5/Aarg_B1/hifi_reads                
#PBS -e /home/amarkee/nas5/Aarg_B1/hifi_reads                
#PBS -M amarkee@amnh.org
#PBS -m abe
#PBS -N aarg_bedtools
#PBS -l walltime=99:00:00

module load bedtools-2.30.0

bedtools bamtofastq -i /home/amarkee/nas5/Aarg_B1/hifi_reads/m84100_240126_204718_s4.hifi_reads.bc1046.bam  -fq /home/amarkee/nas5/Aarg_B1/hifi_reads/m84100_240126_204718_s4.hifi_reads.bc1046.fastq
```

These are the resulting paths to the extracted fastq files:
```
/home/amarkee/nas5/Aarg_B1/hifi_reads/m84100_240126_204718_s4.hifi_reads.bc1046.fastq
/home/amarkee/nas5/Aarg_G2/hifi_reads/m84100_240128_044323_s3.hifi_reads.bc1045.fastq
```

## **Raw Read Quality Assessment with FastQC**

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a quality assessment tool for assessing raw read quality of next generation sequencing. 

Copy the following script into the directory containing HiFi reads (SLURM for Mendel; PBS for Huxley)
```
#!/bin/bash
#PBS -V
#PBS -q batch
#PBS -l select=1:ncpus=25
#PBS -o /home/amarkee/nas4/aargentata_genome/fastqc/Aarg_B1                
#PBS -e /home/amarkee/nas4/aargentata_genome/fastqc/Aarg_B1               
#PBS -M amarkee@amnh.org
#PBS -m abe
#PBS -N aarg_fastqc
#PBS -l walltime=40:00:00

#conda init
source ~/.bash_profile
conda activate spidroins

module load fastqc-0.11.9

## run fastqc on the data 
fastqc /home/amarkee/nas5/Aarg_B1/hifi_reads/m84100_240126_204718_s4.hifi_reads.bc1046.fastq

```

The FastQC output files include an .html file, which will contain visualizations for the following results:
- Basic Statistics
- Per base sequence quality
- Per sequence quality scores (Phred scores)
- Per base sequence content
- Per sequence GC content
- Per base N content
- Sequence Length Distribution
- Sequence Duplication Levels
- Overrepresented sequences
- Adapter content


## **Genome Assembly With Hifiasm**

[hifiasm](https://hifiasm.readthedocs.io/en/latest/) is a fast and easy haplotype-resolved de novo assembly software for PacBio HiFi reads
 - hifiasm documentation explaining input parameters: https://hifiasm.readthedocs.io/en/latest/pa-assembly.html
 - hifiasm documentation explaining output files: https://hifiasm.readthedocs.io/en/latest/interpreting-output.html

Below is the script for assembling the genome with standard duplicate purging enambled (option -l 2)
```
#!/bin/bash
#PBS -V
#PBS -q batch 
#PBS -l select=1:ncpus=25
#PBS -o /home/amarkee/nas4/aargentata_genome/assemblies/primary_asm/Aarg_B1.out
#PBS -e /home/amarkee/nas4/aargentata_genome/assemblies/primary_asm/Aarg_B1.err
#PBS -M amarkee@amnh.org
#PBS -m abe
#PBS -N aarg-b1_hifiasm
#PBS -l walltime=99:00:00

#conda init
source ~/.bash_profile
conda activate spidroins

hifiasm -o /home/amarkee/nas4/aargentata_genome/assemblies/primary_asm/Aarg_B1/aarg_b1_assembly.asm -l 2 -t 32 /home/amarkee/nas5/Aarg_B1/hifi_reads/m84100_240126_204718_s4.hifi_reads.bc1046.fastq

```

Because of extensive resource use between Huxley and Mendel, I'm running the Aarg_B1 assembly on Huxley, and the Aarg_G2 assembly on Mendel. Below is the SLURM submission required to run the assembly on Mendel:
```
#!/bin/bash
#SBATCH --job-name aargb1_hifiasm
#SBATCH --nodes=16
#SBATCH --mem=100gb
#SBATCH --tasks-per-node=1 # Number of cores per node
#SBATCH --time=90:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amarkee@amnh.org
#SBATCH --output=slurm-%j-%x.out

#conda init
source ~/.bash_profile
conda activate spidroins

hifiasm -o /home/amarkee/mendel-nas1/aarg_pbCA/assemblies/aarg_b1_assembly.asm -l 2 -t 32 /home/amarkee/mendel-nas1/raw_data/aarg_raw/Aarg_B1/m84100_240126_204718_s4.hifi_reads.bc1046.fq.gz
```


<br />

## **Genome Assembly QC With assemblystats.py**

- After assembly with hifiasm, we can assess assembly quality using the [assemblystats.py script](https://github.com/MikeTrizna/assembly_stats/tree/0.1.4) created by Mike Trizna.
- The version of assemblystats.py used here was modified by Paul Frandsen (Brigham Young University).

First, I copied this script into my working directory, and called it assemblystats.py

```
#!/usr/bin/env python

import numpy as np
from itertools import groupby
import json
import sys


def fasta_iter(fasta_file):
    """Takes a FASTA file, and produces a generator of Header and Sequences.
    This is a memory-efficient way of analyzing a FASTA files -- without
    reading the entire file into memory.

    Parameters
    ----------
    fasta_file : str
        The file location of the FASTA file

    Returns
    -------
    header: str
        The string contained in the header portion of the sequence record
        (everything after the '>')
    seq: str
        The sequence portion of the sequence record
    """

    fh = open(fasta_file)
    fa_iter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in fa_iter:
        # drop the ">"
        header = next(header)[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.upper().strip() for s in next(fa_iter))
        yield header, seq


def read_genome(fasta_file):
    """Takes a FASTA file, and produces 2 lists of sequence lengths. It also
    calculates the GC Content, since this is the only statistic that is not
    calculated based on sequence lengths.

    Parameters
    ----------
    fasta_file : str
        The file location of the FASTA file

    Returns
    -------
    contig_lens: list
        A list of lengths of all contigs in the genome.
    scaffold_lens: list
        A list of lengths of all scaffolds in the genome.
    gc_cont: float
        The percentage of total basepairs in the genome that are either G or C.
    """

    gc = 0
    total_len = 0
    contig_lens = []
    scaffold_lens = []
    for _, seq in fasta_iter(fasta_file):
        scaffold_lens.append(len(seq))
        if "NN" in seq:
            contig_list = seq.split("NN")
        else:
            contig_list = [seq]
        for contig in contig_list:
            if len(contig):
                gc += contig.count('G') + contig.count('C')
                total_len += len(contig)
                contig_lens.append(len(contig))
    gc_cont = (gc / total_len) * 100
    return contig_lens, scaffold_lens, gc_cont


def calculate_stats(seq_lens, gc_cont):
    stats = {}
    seq_array = np.array(seq_lens)
    stats['sequence_count'] = seq_array.size
    stats['gc_content'] = gc_cont
    sorted_lens = seq_array[np.argsort(-seq_array)]
    stats['longest'] = int(sorted_lens[0])
    stats['shortest'] = int(sorted_lens[-1])
    stats['median'] = np.median(sorted_lens)
    stats['mean'] = np.mean(sorted_lens)
    stats['total_bps'] = int(np.sum(sorted_lens))
    csum = np.cumsum(sorted_lens)
    for level in [10, 20, 30, 40, 50]:
        nx = int(stats['total_bps'] * (level / 100))
        csumn = min(csum[csum >= nx])
        l_level = int(np.where(csum == csumn)[0])
        n_level = int(sorted_lens[l_level])

        stats['L' + str(level)] = l_level
        stats['N' + str(level)] = n_level
    return stats


if __name__ == "__main__":
    infilename = sys.argv[1]
    contig_lens, scaffold_lens, gc_cont = read_genome(infilename)
    contig_stats = calculate_stats(contig_lens, gc_cont)
    scaffold_stats = calculate_stats(scaffold_lens, gc_cont)
    stat_output = {'Contig Stats': contig_stats,
                   'Scaffold Stats': scaffold_stats}
    print(json.dumps(stat_output, indent=2, sort_keys=True))
```

Next, I changed permissions as follows to allow execution permissions.
```
chmod +x assemblystats.py
```

Then, I produced a FASTA file from the initial GFA output files from the hifiasm assembly output for haplotype 1, haplotype 2, and total. I used the primary contig file, as indicated with asm.bp.p_ctg.fa (p_ctg meaning primary contig)
```
awk '/^S/{print ">"$2;print $3}' aarg_b1_assembly.asm.bp.hap1.p_ctg.gfa > aarg_b1_assembly.asm.bp.hap1.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' aarg_b1_assembly.asm.bp.hap2.p_ctg.gfa > aarg_b1_assembly.asm.bp.hap2.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' aarg_b1_assembly.asm.bp.p_ctg.gfa > aarg_b1_assembly.asm.bp.p_ctg.fa
```

Lastly, I ran the assemblystats.py script on the newly generated fasta file of the fga in the format of scriptfilepath/scirptname.py nameofassembly.fa and save as a txt file 
```
./assemblystats.py aarg_b1_assembly.asm.bp.hap1.p_ctg.fa >> aarg_b1_assembly.hap1.p.ctg.stats.txt
./assemblystats.py aarg_b1_assembly.asm.bp.hap2.p_ctg.fa >> aarg_b1_assembly.hap2.p.ctg.stats.txt
./assemblystats.py aarg_b1_assembly.asm.bp.p_ctg.fa >> aarg_b1_assembly.p.ctg.stats.txt
```

## **Genome Assembly QC With Quast**
On the BYU Fulton cluster, I use a conda environment from [the SPIN grant genomics workshop](https://github.com/pbfrandsen/SPIN_workshop/tree/main) to run Quast to assess quality statistics similar to the previous script.

```
cd 
module load miniconda3/4.12-pws-472
conda activate quast
```

Then we can simply run quast on the assembly file. It is pretty fast so we can run it interactively and won't need to put it into a job file.

```
quast <assembly_name>.p_ctg.fasta
```

Now it will run on your genome and the results will be added to `quast_results/latest`. Once your run is complete, you can navigate to that folder. There should be a file called `report.pdf`. Go ahead and download that file with `scp` and examine it's contents. Below is a an example of the output:

<img width="701" alt="Screenshot 2024-06-20 at 9 36 45 AM" src="https://github.com/amandamarkee/spidroins/assets/56971761/ed45c0eb-da2d-453a-aff3-a3ab88eba164">

## **Checking Assembly Completeness with BUSCO** 

[BUSCO](https://busco.ezlab.org/busco_userguide.html) is a program that estimates genome completeness based on evolutionarily-informed expectations of gene content of near-universal single-copy orthologs.

Since the worker nodes of the AMNH's computational clusters don't have access to the internet, it is necesarry to install BUSCO locally:
```
#clone repository
git clone https://gitlab.com/ezlab/busco.git
cd busco/

#pip install busco
python -m pip install .
```

Script for running BUSCO using the arachnida_odb10 database on Mendel:
```
#!/bin/sh
#SBATCH --job-name busco_aargb1
#SBATCH --nodes=1
#SBATCH --mem=40gb
#SBATCH --time=144:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amarkee@amnh.org

source ~/.bash_profile
conda activate spidroins

export BUSCO_CONFIG_FILE=/home/amarkee/mendel-nas1/aarg_pbCA/busco/config/config.ini
echo $BUSCO_CONFIG_FILE

busco -i /home/amarkee/mendel-nas1/aarg_pbCA/assemblies/aarg_b1_assembly.asm.bp.p_ctg.fa \
-o BUSCO_aargentata_b1 -l arachnida_odb10 \
-m genome
```

Script for running BUSCO using the arachnida_odb10 database on Fulton:
```
#!/bin/bash
#SBATCH --time=22:00:00   # walltime
#SBATCH --ntasks=24   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=10240M   # memory per CPU core
#SBATCH -J "aarg_busco"   # job name
#SBATCH --mail-user=amarkee@amnh.org   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Set the max number of threads to use for programs using OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load miniconda3/4.12-pws-472
conda activate busco
busco -o busco_arachnida -i aarg_hifiasm.asm.bp.p_ctg.fasta -l arachnida_odb10 -c 24 -m genome --offline

```
