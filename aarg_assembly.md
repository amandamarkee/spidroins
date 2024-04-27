## **Background**

The goal of this workflow is to take note of assembly details for the five California _Argiope argentata_ genomes. Once genomes are assembled, each pair of haplotype assemblies will be used to assess interspecific (allelic) variation in silk genes
within the spidroin gene family. I will be trying to run the genome assembler (hifiasm) on both American Museum of Natural History high-performance computing systems (Mendel and Huxley)
- California samples (assembled by PBF): Aarg_C1, Aarg_F1, and Aarg_H1
- California samples (not assembled): Aarg_B1 and Aarg_G2

## **Project Set-Up and File Conversion**

Prior to running preliminary QC and genome assembly, I will convert my Revio output (.bam) into usable formats (.fastq) for both programs. 
```
# raw data directory on huxley
/home/amarkee/nas5/Aarg_B1/hifi_reads
/home/amarkee/nas5/Aarg_G2/hifi_reads

# working directory for assemblies
/home/amarkee/nas4/aargentata_genome/assemblies/primary_asm
```

I use the following script on Huxley to convert files using bedtools, changing file path locations for each of the 5 genomes. Note: in the future, I'll transform this into a cleaner environmental variable script:
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
awk '/^S/{print ">"$2;print $3}' XXX_assembly.asm.bp.hap1.p_ctg.gfa > XXX_assembly.hap1.p.ctg.fa
awk '/^S/{print ">"$2;print $3}' XXX_assembly.asm.bp.hap2.p_ctg.gfa > XXX_assembly.hap2.p.ctg.fa
awk '/^S/{print ">"$2;print $3}' XXX_assembly.asm.bp.p_ctg.gfa > XXX_assembly.p.ctg.fa
```

Lastly, I ran the assemblystats.py script on the newly generated fasta file of the fga in the format of scriptfilepath/scirptname.py nameofassembly.fa and save as a txt file 
```
./assemblystats.py XXX_assembly.hap1.p.ctg.fa >> XXX_assembly.hap1.p.ctg.stats.txt
./assemblystats.py XXX_assembly.hap2.p.ctg.fa >> XXX_assembly.hap2.p.ctg.stats.txt
./assemblystats.py XXX_assembly.p.ctg.fa >> XXX_assembly.p.ctg.stats.txt
```

## **Checking Assembly Completeness with BUSCO** 

[BUSCO](https://busco.ezlab.org/busco_userguide.html) is a program that estimates genome completeness based on evolutionarily-informed expectations of gene content of near-universal single-copy orthologs.

Since the worker nodes of the AMNH's computational clusters don't have access to the internet, it is necesarry to install BUSCO locally:
```
#clone repository
git clone https://gitlab.com/ezlab/busco.git
cd busco/

#pip install
python -m pip install .
```

Make sure you have all the [dependencies](https://busco.ezlab.org/busco_userguide.html#editing-busco-run-configuration) installed for the type of BUSCO run you are planning on running.

