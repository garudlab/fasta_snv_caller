# Fasta SNV Caller

This repository provides a simple way to call SNPs and annotate them as 
synonymous or non-synonymous from a set of high-quality fasta files (i.e. with
large, overlapping contigs). It builds up from code in Almeida et al. 2021 (see github:
[snv_analysis_almeida2019](https://github.com/zjshi/snv_analysis_almeida2019.git)
), and some code to annotate "haplotype" dataframes with genes and syn/non-syn 
variants from the 
[Garud lab's codebase](https://github.com/garudlab/UHGG.git) to work with UHGG 
data. 

## A more detailed description

The scripts included attempt to align fasta files one-by-one from a 
number of different fastas to a reference fasta using the nucmer aligner from 
MUMmer. Single nucleotide polymorphisms (SNPs) are then called from the 
concatenated set of SNPs with respect to the reference. A gff file for the 
reference is then used to pull coding sequences (CDS) and assign individual SNPs 
to genes and annotate them as synonymous (no change to amino acid sequence) or 
non-synonymous (amino acid change). An amino acid fasta (*.faa) derived from the 
gff's CDS can be used as a check for the reference amino acid state. Finally, I 
include a script to download RefSeq data by passing in accession numbers for a 
Bioproject (e.g. a set of assemblies from an individual study) and a reference 
accession. The default is to download *Neisseria gonorrhoeae* sequences from 
Grad et al. 2014 and a well-studied *Neisseria gonorrhoeae* reference isolate
FA1090. 

## A quick tutorial (for *Neisseria gonorrhoeae*, ~20 min):

0. (`conda env create -f fasta_snv_caller.yml`)
1. `bash make_snps/download_command.sh`
2. `bash make_snps/run_all_snv_catalog.sh`
3. `bash annotate_haps/run_annotations.sh`

## Variables to change 

In `download_command.sh`: species, species_snvs_dir, (genome_accession), (ref_accession)

In `run_all_snv_catalog.sh`: SPECIES, ref_name, ref_fasta, species_snvs_dir

In `run_annotations.sh`: SPECIES, species_snvs_dir, hap_dir

Other considerations: 
* Currently, the `annotation_utils.py` code depends on 
having a `.faa` file to check the reference state of the amino acid for that SNP. 
If required, you may need to manually make this optional in the script.
* This repo uses default settings for `nucmer` and MUMmer-associated code. 
But incorporating user knowledge on the lengths and general alignability of contigs, 
etc. in the `nucmer` settings may improve your SNP calling :) 

## What's required

Environment dependencies: 
* Python (3.8.5) (for env see fasta_snv_caller.yaml)
* MUMmer (I'm using [4.0.0rc1)](https://github.com/mummer4/mummer/releases) 
* (optional) [NCBI's datasets tool](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/) 
(if planning on downloading RefSeq data straight from NCBI)

Data dependencies:

* If using NCBI datasets:
    * reference accession number
    * BioProject (or genome) accession number
* If using own data:
    * Reference DNA fasta file (`.fna`)
    * Reference general feature format file (`.gff`)
    * (optional) Reference amino acid fasta file (`.faa`)
    * Non-reference DNA fasta files


## Output: 

In species_snvs_dir: 

* haplotypes/
    * [contig]_haplotypes.csv (haplotype df w genes and syn/nonsyn annotation)

* snp_annotations/
    * [contig]_annotations.csv (snps df w genes, syn/nonsyn, and amino acid annotation)

* snps/
    * [ref_name].catalog.final.tsv (matches 
[snv_analysis_almeida2019](https://github.com/zjshi/snv_analysis_almeida2019.git) 
format)

(and other intermediate files)



## References

> Almeida, Alexandre, et al. "A unified catalog of 204,938 reference genomes from the human gut microbiome." Nature biotechnology 39.1 (2021): 105-114.

> Grad, Yonatan H., et al. "Genomic epidemiology of Neisseria gonorrhoeae with reduced susceptibility to cefixime in the USA: a retrospective observational study." The Lancet infectious diseases 14.3 (2014): 220-226.

Major thanks to Ricky Wolff for doing some heavy lifting with the annotations code, 
and first suggesting the Almeida pipeline as a way to call SNVs. 

Feel free to submit a pull request or raise an issue on GitHub!


