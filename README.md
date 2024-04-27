# Genomic characterization of spider silk genes in the silver garden spider _Argiope argentata_

The silver garden spider (_Argiope argentata_) is an orb-weaving spider in the family Araneidae, with seven distinct types of silk: major ampullate (MaSp), minor ampullate (MiSp), flagelliform (FlAg), aggregate (AgSp), aciniform (AcSp), tubuliform (TuSp) and pyriform (PySp). Each silk type is characterized based on the gland they originate from, and have variation in the properties of silk the produce. Spider fibroins, further referred to as _spidroins_, are a family of silk proteins produced by long, highly repetitive, glycine-rich genes.

</br>
<img width="981" alt="Screenshot 2024-04-26 at 10 06 26 PM" src="https://github.com/amandamarkee/spidroins/assets/56971761/5df66f95-9789-46d8-9c9b-b9acee02dd9e">
<br/><br/>


This repository will serve as documenting the bioinformatic workflow required for full assembly and annotation of six _Argiope argentata_ genomes from disperate populations in California (N=5) and Texas (N=1), as part of the first chapter of my dissertation. The ultimate goal is to assess interspecfic variation (allelic variation) of different spidroins to glean insight into the diversification and molecular evolution of this unique gene family. 

## Workflow

1) Raw Read Quality Control with FastQC
2) Whole Genome Assembly with Hifiasm
3) Assembly Quality Control with BUSCO and assemblystats.py
4) Putitive Spidroin Identification with BLASTx
5) Spidroin Extraction with Geneious
6) Annotation of Intron/Exon Boundaries with IGV and Geneious
