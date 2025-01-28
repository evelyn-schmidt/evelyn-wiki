# Note that some of these paths may not be correct since I orignally ran this code locally

## install rust
`curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`

I selected 1) Proceed with standard installation (default - just press enter)

## download the zipped source code and unzip it

```
If you need to compile from source, install Rust, then type cargo build --release from within the directory containing the scHLAcount source code. The executable will appear at target/release/sc_hla_count. As usual it's important to use a release build to get good performance.
```

I had to do this:
https://users.rust-lang.org/t/error-e0310-im-a-complete-beginner/114775
`cargo update -p rustc-serialize`

Before this command that builds the executable would work:
`cargo build --release`

## Get reference files
- install samtools locally (could have also used a docker)
- get files from here https://github.com/ANHIG/IMGTHLA

```
/usr/local/bin/samtools faidx hla_nuc.fasta
/usr/local/bin/samtools faidx hla_gen.fasta
```

Make sure these are in your `scHLAcount-0.2.0` directory

## Create a file of the known genotypes

Questions:
Do I include both samples in the file?
	yes i think because that will genotype both?
Do I have to have them listed in this `DQA1*01:02:01:01` specfic format?
What if I include duplicated?
We strongly recommend that if genotypes are unknown for any of the genes, you put the reference genome allele for those genes in the known genotypes file. Alleles represented in the GRCh38 primary assembly are listed below:?
```
A*11
A*30
B*27
B*49
C*02
C*07
DRB1*07:01
DRB4*01
DQB1*03:02
DQB1*02:02
DPB1*04:01
DPB1*02:01
A*29:02
A*11:01
B*27:05
B*35:01
C*04:01
C*02:02
DRB1*01:01
DRB1*04:04
DRB4*01
DQB1*03:02
DQB1*05:01
```

## Setting up BAM

### Make Sure the bam file is in the correct format
Make chromosomes in the correct format :melting_face:

The contig name must just be a number, no `chr` or any other characters. **This might be the case for files produced by newer Cellranger** runs.

To view the header of your bam files:
`/usr/local/bin/samtools view -h possorted_genome_bam_numeric.bam`

```
SN:3
```

```
bsub -Is -n 8 -G compute-oncology -q oncology-interactive -M 32G -R 'select[mem>32G] span[hosts=1] rusage[mem=32G]' -a 'docker(quay.io/biocontainers/samtools:1.11--h6270b1f_0)' /bin/bash

samtools reheader <(samtools view -H possorted_genome_bam.bam | sed 's/GRCh38_//g') possorted_genome_bam.bam > possorted_genome_bam_numeric.bam`

/usr/local/bin/samtools index possorted_genome_bam_numeric.bam
```

### Get all the bardcodes from the bam

`bsub -n 8 -G compute-oncology -q oncology -eo err_barcodes.log -oo out_barcodes.log -M 32G -R 'select[mem>32G] span[hosts=1] rusage[mem=32G]' -a 'docker(quay.io/biocontainers/samtools:1.11--h6270b1f_0)' /bin/bash get_cell_barcodes.sh`


`/usr/local/bin/samtools view possorted_genome_bam.bam | cut -f 12- | tr "\t" "\n" | grep "CB:Z:" | uniq > cell_barcodes_all.txt `


2. subset bam

`/usr/local/bin/samtools view -hb --threads 8 possorted_genome_bam_numeric.bam 6  > chr6.bam`

`/usr/local/bin/samtools index chr6.bam`

## Run prepare_reference.sh

`/Users/evelynschmidt/Bioinformatics_tools/scHLAcount-0.2.0/prepare_reference.sh genotype.txt `

gave me this error
```
[W::fai_get_val] Reference HLA:HLA00404 not found in FASTA file, returning empty sequence
[W::fai_get_val] Reference HLA:HLA00404 not found in FASTA file, returning empty sequence
[faidx] Failed to fetch sequence in HLA:HLA00404
```

## Run Commnad
```
/Users/evelynschmidt/Bioinformatics_tools/scHLAcount-0.2.0/target/release/sc_hla_count -f cds.fasta -g gen.fasta -b /Volumes/jennifer.a.foltz/Active/Evelyn/possorted_genome_bam.bam -c /Volumes/jennifer.a.foltz/Active/Evelyn/cell_barcodes_test.txt 
```
