# Running Soup Or Cell

[Paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7617080)
[GitHub](https://github.com/wheaton5/souporcell)
[Docker](https://hub.docker.com/r/cumulusprod/souporcell)


## Interactive Job Docker 

```
bsub -Is -n 8 -G compute-oncology -q oncology-interactive -M 32G -R 'select[mem>32G] span[hosts=1] rusage[mem=32G]' -a 'docker(cumulusprod/souporcell:2.5)' /bin/bash
```

## Submitting the whole pipeline

I did this step but only the first 'renamer' step succeeded, switching to a step by step approach

bsub -n 8 -eo err.log -oo out.log -G compute-oncology -q general -M 32G -R 'select[mem>32G] span[hosts=1] rusage[mem=32G]' -a 'docker(cumulusprod/souporcell:2.5)' python3 /opt/souporcell/souporcell_pipeline.py  -i /storage1/fs1/jennifer.a.foltz/Active/Evelyn/genotyping_nk_cells/scHLAcount/data/possorted_genome_bam.bam -f /storage1/fs1/bga/Active/gmsroot/gc2560/core/reference_sequences/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa -t 8 -b /storage1/fs1/jennifer.a.foltz/Active/Evelyn/genotyping_nk_cells/souporcell/barcodes.tsv.gz -o outdir_bjobs -k 8


## 1. Remapping

Run the renamery.py script, seemed to have succeed and outputs are here `outdir_bjobs/tmp.fq`

```
bsub -n 4 -eo err.log -oo out.log -G compute-oncology -q general -M 32G -R 'select[mem>32G] span[hosts=1] rusage[mem=32G]' -a 'docker(cumulusprod/souporcell:2.5)' python3 /opt/souporcell/renamer.py --bam possorted_genome_bam.bam --barcodes barcodes.tsv --out fq.fq
```

Now remap using minimap, have to put this into bash scirpt

```
bsub -n 8 -eo err_minimap.log -oo out_minimap.log -G compute-oncology -q oncology -M 64G -R 'select[mem>64G] span[hosts=1] rusage[mem=32G]' -a 'docker(cumulusprod/souporcell:2.5)' /bin/bash ./scripts/minimap_cmd.sh
```

Now we must retag the reads with their cell barcodes and UMIs

```
bsub -n 4 -eo err.log -oo out.log -G compute-oncology -q general -M 32G -R 'select[mem>32G] span[hosts=1] rusage[mem=32G]' -a 'docker(cumulusprod/souporcell:2.5)' python3 /opt/souporcell/retag.py --sam minimap.sam --out minitagged.bam
```

Then we must sort and index our bam. Requires samtools

```
bsub -n 8 -eo err_sort.log -oo out_sort.log -G compute-oncology -q oncology -M 32G -R 'select[mem>32G] span[hosts=1] rusage[mem=32G]' -a 'docker(cumulusprod/souporcell:2.5)' /opt/samtools/samtools sort minitagged.bam -o minitagged_sorted.bam --threads 8
```

```
bsub -n 4 -eo err_index.log -oo out_index.log -G compute-oncology -q oncology -M 32G -R 'select[mem>32G] span[hosts=1] rusage[mem=32G]' -a 'docker(cumulusprod/souporcell:2.5)' /opt/samtools/samtools index minitagged_sorted.bam -@ 4
```

## 2. Calling candidate variants

```
bsub -n 1 -eo err_freebayes.log -oo out_freebayes.log -G compute-oncology -q oncology -M 64G -R 'select[mem>64G] span[hosts=1] rusage[mem=64G]' -a 'docker(cumulusprod/souporcell:2.5)'   \ 
-f /storage1/fs1/bga/Active/gmsroot/gc2560/core/reference_sequences/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa 
-iXu # Remove indel observations from input,Remove MNP observations from input, Remove complex allele observations from input
-C 2 # min-alternate-count
-q 20 # min-base-quality
-n 3 # Evaluate only the best N SNP alleles, ranked by sum of supporting quality scores.  (Set to 0 to use all; default: all)
-E 1 # max-complex-gap
-m 30 # min-mapping-quality
--min-coverage 6 
--limit-coverage 100000 minitagged_sorted.bam -v freebayes_var.vcf
```

From the Freebayes manual on -iXu
```
These flags are meant for testing.
  They are not meant for filtering the output.
  They actually filter the input to the algorithm by throwing away alignments.
  This hurts performance by hiding information from the Bayesian model.
  Do not use them unless you can validate that they improve results!
```

Maybe this to get more variants?
```
bsub -n 1 -eo err_freebayes.log -oo out_freebayes.log -G compute-oncology -q oncology -M 64G -R 'select[mem>64G] span[hosts=1] rusage[mem=64G]' -a 'docker(cumulusprod/souporcell:2.5)'  /opt/freebayes -f /storage1/fs1/bga/Active/gmsroot/gc2560/core/reference_sequences/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa minitagged_sorted.bam -v freebayes_defaultparams.vcf
```

## 2. Cell allele counting 

```
bsub -n 8 -eo err_vartrix.log -oo out_vartrix.log -G compute-oncology -q oncology -M 32G -R 'select[mem>32G] span[hosts=1] rusage[mem=32G]' -a 'docker(cumulusprod/souporcell:2.5)' /opt/vartrix --umi --mapq 30 -b minitagged_sorted.bam -c barcodes.tsv.gz --scoring-method coverage --threads 8 --ref-matrix ref.mtx --out-matrix alt.mtx -v freebayes_var.vcf --fasta /storage1/fs1/bga/Active/gmsroot/gc2560/core/reference_sequences/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa
```

^^ Did i use the right BAM here?


## 4. Clustering cells by genotype (souporcell) 

```
bsub -n 8 -eo err_souporcell.log -oo out_souporcell.log -G compute-oncology -q oncology -M 32G -R 'select[mem>32G] span[hosts=1] rusage[mem=32G]' -a 'docker(cumulusprod/souporcell:2.5)' /bin/bash ./scripts/souporcell_cmd.sh
```

Cluster but the number of invoduals that you have

## 5. Calling doublets

```
/opt/souporcell/troublet/target/release/troublet -a alt.mtx -r ref.mtx --singlet_threshold 0.7 --clusters clusters_tmp.tsv > clusters.tsv
```

```
1260740 loaded 0 counts, is this a problem? # with five clusters
```

```
78856 loaded 0 counts, is this a problem? # with two clusters 
```
## 6. 

bsub -n 8 -eo err_consensus.log -oo out_consensus.log -G compute-oncology -q general -M 32G -R 'select[mem>32G] span[hosts=1] rusage[mem=32G]' -a 'docker(cumulusprod/souporcell:2.5)' python3 /opt/souporcell/consensus.py -c clusters.tsv -a alt.mtx -r ref.mtx --soup_out soup.txt -v freebayes_var.vcf --vcf_out cluster_genotypes.vcf --output_dir .

Recieve this error
```
INFO:pystan:COMPILING THE C++ CODE FOR MODEL anon_model_ca32e407e94c33afe8f72cdc7357f09f NOW.
Traceback (most recent call last):
  File "/opt/souporcell/consensus.py", line 199, in <module>
    cluster = int(tokens[2])
ValueError: invalid literal for int() with base 10: '-50.924984'
```
Manually open the cluster tsv in excel and converted the values to numeric...

Okay I just figured out how to run the doublet finder and got this message:

```
INFO:pystan:COMPILING THE C++ CODE FOR MODEL anon_model_ca32e407e94c33afe8f72cdc7357f09f NOW.
63763 excluded for potential RNA editing
1897 doublets excluded from genotype and ambient RNA estimation
304 not used for soup calculation due to possible RNA edit
Initial log joint probability = -13461.8
    Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
       5      -12445.2   0.000449259    0.00553004           1           1        6   
Optimization terminated normally: 
  Convergence detected: relative gradient magnitude is below tolerance
```

And with two clusters:

```
INFO:pystan:COMPILING THE C++ CODE FOR MODEL anon_model_ca32e407e94c33afe8f72cdc7357f09f NOW.
63763 excluded for potential RNA editing
1761 doublets excluded from genotype and ambient RNA estimation
311 not used for soup calculation due to possible RNA edit
Initial log joint probability = -6841.93
    Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
       4      -5856.82   0.000156917   1.03373e-05           1           1        8   
Optimization terminated normally: 
  Convergence detected: relative gradient magnitude is below tolerance

```

```ambient RNA estimated as 58.92420228793033%```


# Creating a Validation Dataset

## Figuring out away to see if this works

bsub -n 8 -Is -G compute/ -g /evelyn/default -q oncology-interactive -M 16G -R 'rusage[mem=16G]' -a 'docker(chrisamiller/docker-genomic-analysis)' /bin/bash

- maybe we can combine two bams

I got A BAM from Weihang

### 1. subset to smaller bam (like I did with scHLAcount)

`samtools view -hb --threads 8 weiheng_test_alignments.bam chr6  > weiheng_test_alignments_chr6.bam`

### 2. Add a new tag to bam

```
/usr/bin/java -Xmx16g -jar /usr/picard/picard.jar AddOrReplaceReadGroups       I=../../scHLAcount/data/chr6.bam       O=chr6_labeled.bam       RGID=1     RGLB=lib1       RGPL=illumina       RGPU=unit1       RGSM="TestSample"
```

```
/usr/bin/java -Xmx16g -jar /usr/picard/picard.jar AddOrReplaceReadGroups       I=weiheng_test_alignments_chr6.bam       O=weiheng_test_alignments_chr6_labeled.bam       RGID=2     RGLB=lib1       RGPL=illumina       RGPU=unit1       RGSM="WeihangSample"
```

3. append one bam onto another

```
samtools merge -f --threads 8 merged_chr6_labled.bam chr6_labeled.bam weiheng_test_alignments_chr6_labeled.bam
```

4. See if it worked

```
samtools view -H merged_chr6_labled.bam | grep "@RG"
```

```
@RG	ID:1	LB:lib1	PL:illumina	SM:TestSample	PU:unit1
@RG	ID:2	LB:lib1	PL:illumina	SM:WeihangSample	PU:unit1
```
SLAYYYYYY

5. Sort and index?

```
samtools sort -@ 8 merged_chr6_labled.bam -o merged_chr6_labled_sort.bam
samtools index merged_chr6_labled_sort.bam
```

6. rerun souporcell make sure that sample name is perserved!!
-- So I think if I edit the renamer to add the RG ID to the cell name that might be a way to fix this, but I would have to create a docker with VIM inside it....

FINE I WILL DO IT MYSELF

```
bsub -G compute-oncology -q oncology-interactive -Is -a 'docker_build(evelyns2000/souporcell)' -- --tag evelyns2000/souporcell .


bsub -Is -n 8 -G compute-oncology -q oncology-interactive -M 32G -R 'select[mem>32G] span[hosts=1] rusage[mem=32G]' -a 'docker(evelyns2000/souporcell)' /bin/bash
```

SOOOO after the renamer step the ID is in the fastq read name SLAY

Also had to fix the retag.py step but it was really easy

```
    else:
        print(tokens)
        print(len(tokens))
        assert len(tokens) == 4
        read.set_tag(CELL_TAG, tokens[-3])
        read.set_tag(UMI_TAG, tokens[-2])
        read.set_tag("RG", tokens[-1])
    bamout.write(read)
```

#### Creating a new header for the freebayes bam

```
samtools view -h minitagged_sorted.bam | head -200 | grep '^@' > bam_header.txt
cat readgroups.txt >> bam_header.txt
samtools reheader bam_header.txt minitagged_sorted.bam > new_minitagged_sorted.bam
samtools index new_minitagged_sorted.bam
```

readgroups.txt contains

```
@RG ID:1    SM:sample1  LB:lib1 PL:ILLUMINA
@RG ID:2    SM:sample2  LB:lib2 PL:ILLUMINA
```

```
bsub -n 1 -eo err_freebayes.log -oo out_freebayes.log -G compute-oncology -q oncology -M 64G -R 'select[mem>64G] span[hosts=1] rusage[mem=64G]' -a 'docker(cumulusprod/souporcell:2.5)'  /opt/freebayes -f /storage1/fs1/bga/Active/gmsroot/gc2560/core/reference_sequences/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa -iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 6 --limit-coverage 100000 new_minitagged_sorted.bam -v freebayes_var.vcf

```

# Further exploration of if it works

the other file was FLEX data... not good, so we try another!!

BEWARE!
- This sample has mouse genes, which shouldn't effect it because we are just taking chr6

```
bsub -n 8 -Is -G compute/ -g /evelyn/default -q oncology-interactive -M 16G -R 'rusage[mem=16G]' -a 'docker(chrisamiller/docker-genomic-analysis)' /bin/bash

cd /storage1/fs1/jennifer.a.foltz/Active/Evelyn/genotyping_nk_cells/souporcell/test_case/test_part2_CorrectData/bam_prep

samtools view -hb --threads 8 TCIML7_Donor_possorted_genome_bam.bam 6  > TCIML7_Donor_possorted_genome_chr6.bam

/usr/bin/java -Xmx16g -jar /usr/picard/picard.jar AddOrReplaceReadGroups       I=TCIML7_Donor_possorted_genome_chr6.bam       O=TCIML7_Donor_possorted_genome_chr6_labeled.bam       RGID=2     RGLB=lib1       RGPL=illumina       RGPU=unit1       RGSM="TCIML7_Donor"

samtools merge -f --threads 8 merged_chr6_labled.bam ../../test_part1_verysilly/bam_prep/chr6_labeled.bam TCIML7_Donor_possorted_genome_chr6.bam

samtools sort -@ 8 merged_chr6_labled.bam -o merged_chr6_labled_sort.bam
samtools index merged_chr6_labled_sort.bam

bsub -n 8 -G compute-oncology -q oncology -eo err_barcodes.log -oo out_barcodes.log -M 32G -R 'select[mem>32G] span[hosts=1] rusage[mem=32G]' -a 'docker(quay.io/biocontainers/samtools:1.11--h6270b1f_0)' /bin/bash ../scripts/get_cell_barcodes.sh
```

1. Editing Renamer
```
bsub -Is -n 8 -G compute-oncology -q oncology-interactive -M 32G -R 'select[mem>32G] span[hosts=1] rusage[mem=32G]' -a 'docker(evelyns2000/souporcell)' /bin/bash

```

```
        tag = read.get_tag("RG")
        if read.has_tag(CELL_TAG) and read.get_tag(CELL_TAG) in cell_barcodes:
            if args.no_umi:
                fastq.write("@"+read.qname+";"+cell_barcode+";"+tag+"\n")
            else:
                fastq.write("@"+read.qname+";"+cell_barcode+";"+UMI+";"+tag+"\n")
            fastq.write(read.seq+"\n")
            fastq.write("+\n")
            fastq.write(read.qual+"\n")
```

```
python3 /opt/souporcell/renamer.py --bam bam_prep/merged_chr6_labled_sort.bam --barcodes bam_prep/cell_barcodes.txt --out output/fq.fq
```

2. Minimap
```
bsub -n 8 -eo logs/err_minimap.log -oo logs/out_minimap.log -G compute-oncology -q oncology -M 64G -R 'select[mem>64G] span[hosts=1] rusage[mem=32G]' -a 'docker(cumulusprod/souporcell:2.5)' /bin/bash ./scripts/minimap_cmd.sh
```

3. Edit retag 
```
    else:
        print(tokens)
        print(len(tokens))
        assert len(tokens) == 4
        read.set_tag(CELL_TAG, tokens[-3])
        read.set_tag(UMI_TAG, tokens[-2])
        read.set_tag("RG", tokens[-1])
    bamout.write(read)
```

```
python3 /opt/souporcell/retag.py --sam output/minimap.sam --out output/minitagged.bam

/opt/samtools/samtools sort minitagged.bam -o minitagged_sorted.bam --threads 8
/opt/samtools/samtools index minitagged_sorted.bam -@ 4
```

4. Freebayes
```
samtools view -h minitagged_sorted.bam | head -200 | grep '^@' > bam_header.txt
cat readgroups.txt >> bam_header.txt
samtools reheader bam_header.txt minitagged_sorted.bam > new_minitagged_sorted.bam
samtools index new_minitagged_sorted.bam
```

readgroups.txt contains

```
@RG ID:1    SM:sample1  LB:lib1 PL:ILLUMINA
@RG ID:2    SM:TCIML7_Donor LB:lib2 PL:ILLUMINA
```

Not that I change the min coverage from 6 to 3
```
bsub -n 1 -eo logs/err_freebayes.log -oo logs/out_freebayes.log -G compute-oncology -q oncology -M 64G -R 'select[mem>64G] span[hosts=1] rusage[mem=64G]' -a 'docker(cumulusprod/souporcell:2.5)'  /opt/freebayes -f /storage1/fs1/bga/Active/gmsroot/gc2560/core/reference_sequences/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa -iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 3 --limit-coverage 100000 output/new_minitagged_sorted.bam -v output/freebayes_var.vcf

```

5. Vartrix

```
bsub -n 8 -eo err_vartrix.log -oo out_vartrix.log -G compute-oncology -q oncology -M 32G -R 'select[mem>32G] span[hosts=1] rusage[mem=32G]' -a 'docker(cumulusprod/souporcell:2.5)' /opt/vartrix --umi --mapq 30 -b output/minitagged_sorted.bam -c bam_prep/cell_barcodes.txt --scoring-method coverage --threads 8 --ref-matrix output/ref.mtx --out-matrix outpu/alt.mtx -v output/freebayes_var.vcf --fasta /storage1/fs1/bga/Active/gmsroot/gc2560/core/reference_sequences/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa
```

6. Clustering by Genotype
```
bsub -n 8 -eo err_souporcell.log -oo out_souporcell.log -G compute-oncology -q oncology -M 32G -R 'select[mem>32G] span[hosts=1] rusage[mem=32G]' -a 'docker(cumulusprod/souporcell:2.5)' /bin/bash ./scripts/souporcell_cmd.sh
```

7. Doublet Calling

```
/opt/souporcell/troublet/target/release/troublet -a output/alt.mtx -r output/ref.mtx --singlet_threshold 0.7 --clusters output/clusters_tmp.tsv > output/clusters.tsv
```

Most cells say 'error' which i think means unasigned... 

8. Consensus

```
python3 /opt/souporcell/consensus.py -c output/clusters.tsv -a output/alt.mtx -r output/ref.mtx --soup_out output/soup.txt -v output/freebayes_var.vcf --vcf_out output/cluster_genotypes.vcf --output_dir output
```


#### Get the variants with different genotypes

```
[evelyn@compute1-exec-280:test_part2_CorrectData]$ grep -v '^#' output/cluster_genotypes.vcf | cut -f1 | wc -l
46634
[evelyn@compute1-exec-280:test_part2_CorrectData]$ grep -v '^#' output/cluster_genotypes.vcf | cut -f1 | grep "6" | wc -l
45060
[evelyn@compute1-exec-280:test_part2_CorrectData]$ grep -v '^#' output/cluster_genotypes.vcf | awk -F'\t' '{split($10,a,":"); split($11,b,":"); if(a[1]!=b[1]) print $0}' | wc -l
35116
```

487 variants with different genotypes founf in other chromosomes beside 6

So there are 34631 variants that seem like they might be real... that is pretty good

```
grep "CHROM" output/freebayes_var.vcf > freebayes_variants.tsv

grep -v '^#' output/freebayes_var.vcf >> freebayes_variants.tsv

grep "CHROM" output/cluster_genotypes.vcf > souporcell_variants.tsv

grep -v '^#' output/cluster_genotypes.vcf >> souporcell_variants.tsv
```

```
grep -v '^#' output/cluster_genotypes.vcf | awk -F'\t' '{split($10,a,":"); split($11,b,":"); if(a[1]!=b[1]) print $0}' > different_genotypes.vcf
```

### So I am kinda suspicious that the genotype info is used in freebayes

```
cd /storage1/fs1/jennifer.a.foltz/Active/Evelyn/genotyping_nk_cells/souporcell/test_case/test_part3_NoReadGroups

```

```
/usr/bin/java -Xmx16g -jar /usr/picard/picard.jar AddOrReplaceReadGroups       I=new_minitagged_sorted.bam      O=new_minitagged_sorted_unifiedRG.bam       RGID=1     RGLB=lib1       RGPL=illumina       RGPU=unit1       RGSM="Test"

samtools sort -@ 8 new_minitagged_sorted_unifiedRG.bam -o new_minitagged_sorted_unifiedRG.bam

samtools index new_minitagged_sorted_unifiedRG.bam

/opt/freebayes -f /storage1/fs1/bga/Active/gmsroot/gc2560/core/reference_sequences/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa -iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 3 --limit-coverage 100000 new_minitagged_sorted_unifiedRG.bam -v freebayes_var.vcf
```