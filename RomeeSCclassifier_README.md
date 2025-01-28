Evelyn Schmidt
02 2024

# Downloading files from SRA
https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/

## 1. obtain accession list
Basically obtain list of SRRs to download

- you start with the [GSE](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198141)
- Then you got to the [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA813854&o=acc_s%3Aa)
- Obtaining these identifiers wiill be how you track all of this down
		SRA: SRP363984
		BioProject: PRJNA813854
		GEO: GSE198141
- Searching the BioPorject number will get you to an overview of the project. There seems to be three sections on this page: Common Fields, Select, Found X Items. To get the aasertion list got to Selct and in the Total row click Accession List which will automically download a txt file
		SRR18324462
		SRR18324463
		SRR18324464
		SRR18324465
		SRR18324466
		SRR18324467
		SRR18324468
		SRR18324469
		SRR18324470
		SRR18324471
		SRR18324472
		SRR18324473
		SRR18324474
		SRR18324475
		SRR18324476
		SRR18324477


## 2.  Make sure the toolkit  can fine the preconfigured settings 
`mkdir $HOME/.ncbi`

Check if its there:
`bsub -n 1 -Is -G compute/ -g /evelyn/default -q general-interactive -M 16G -R 'rusage[mem=16G]' -a 'docker(ncbi/sra-tools)' vdb-config -o n NCBI_SETTINGS`

## 3. Prefetch?
`bsub -n 1 -Is -G compute/ -g /evelyn/default -q general-interactive -M 16G -R 'rusage[mem=16G]' -a 'docker(ncbi/sra-tools)' prefetch --option-file ../sra_acc_list.txt `

You will have to specfy this max size option for any realtively big fastq so might as well so it for all of them
`bsub -n 4  -G compute-oncology  -M 64G -R 'rusage[mem=64G]' -a 'docker(ncbi/sra-tools)' prefetch --max-size 100000000 SRR18324471`

## 4. Convert to fastq
`bsub -n 4  -G compute-oncology  -M 64G -R 'rusage[mem=64G]' -a 'docker(ncbi/sra-tools)' fasterq-dump --split-files SRR18324475/SRR18324475.sra`

## 5. Zip Fastqs
This takes exactly TOO LONG -- run overnight/over the weekend if you cannn

`gzip SRR18324464_S1_L001_R2_001.fastq`

```
#!/bin/bash

for file in $(find . -type f -name "*.fastq"); do 
	echo "Zipping $file"
	gzip "$file"
done
```

`bsub -J 'zippingFastqs' -oo out.log -M 200000000 -q oncology -G compute-oncology -R "select[mem>200000] span[hosts=1] rusage[mem=200000]" -n 4 -a 'docker(chrisamiller/docker-genomic-analysis)' /bin/bash zip_fastq.sh`

 ### So the package pigz seemed to zip with mutiple cores and is therefore faster!

```
#!/bin/bash

for file in $(find . -type f -name "*.fastq"); do 
	echo "Zipping $file"
	/usr/local/bin/pigz/pigz "$file"
done
```

`bsub -J 'pigzFastqs' -oo pigz_out.log -M 64GB -q oncology -G compute-oncology -R "select[mem>64GB] span[hosts=8] rusage[mem=200000]" -n 8 -a 'docker(evelyns2000/foltz_tools:test)' /bin/bash pigz_fastqs.sh`

# Cellranger
**FASTQS MUST BE NAMED IN THIS FORMAT: sample_S1_L001_R1_001.fastq**

`bsub -J 'SRR18324478' -oo out.log -M 200000000 -q oncology -G compute-oncology -R "select[mem>200000] span[hosts=1] rusage[mem=200000]" -n 4 -a 'docker(registry.gsc.wustl.edu/alex.paul/cellranger:6.1.2)' /apps/cellranger-6.1.2/cellranger count \
--localmem=200 \
--localcores=4 \
--id="SRR18324478" \
--transcriptome=/storage1/fs1/bga/Active/gmsroot/gc2560/core/reference_sequences/refdata-cellranger-GRCh38-3.0.0 \
--fastqs=/storage1/fs1/tfehnige/Active/scRNA-seq/ML_NK/Mechanisms_paper/Romee-JCI/GSE198141/GSE198141_RAW/sra-data-repo/SRR18324478 \
--expect-cells 10000`

Make sure the files are named with the suffix `_S1_L001_R1_001` -- I did this manually which I know, I know is not the CS way of doing it but it was 16 files so i didn't feel like the risk of messing them all up was worth the 10 minutes I spent manually renaming them..

# Creating Seurat Objects
https://www.jci.org/articles/view/154334

## Identifying NK Cell populations:

### Looking At Jennifer's Previous Identifications (do not trust)
Tcells 				-- CD3D, CD3E, TRAC, TCF7 is 
Bcells 			    -- CD79A, ENTPD1 
NK Bcells			-- GNLY, FGFBP2, KLRC2, KIR2DL4, KIT,  ~TCF7 
plasma 				-- CD34 

### Paper
 T cells 			-- CD3E, CD4, CD8A (CD3E and CD8A can also be present in NK cells)
 naive T cells 		-- IL7R, CCR7 
 B cells 			-- CD19, CD20, SDC1??, MS4A1==CD20 
 CD14+ monocytes 	-- CD14, LYZ
 NK cells 			-- NKG7, GNLY
 platelets 			-- PPBP

### After Some Research:
NK 					--- KLRC1, KLRC2, GNLY, KIR2DL4 (lacking CD3?) ~FGFBP2, ~KIT, ~ENTPD1
T  					--- GZMK, CD3D, TRAC, CD3E, TCF7, TRDV2, FOXP3  (also CD5)
B  					--- CD79A, MZB1, IGKC, ~ENTPD1, PLD4(??)

myloid 				-- CD68
Progenitor Cell 	-- CD34, KIT
LYZ, CST3 			-- some immune cell, monocytes
ENTPD1, JCHAIN 		-- Plasma
FGFBP2

**Only CD3E, LYZ, and GNLY are used by Jennnifer and Paper**


## Progress

- Patient 1 ============================
	- step 1 						DONE
	- step 3 						DONE
	- step 4a						DONE
	- step 4b						DONE
	- step 5 						DONE
	- Export R matricies			DONE
	- run classifer					DONE
- Patient 2 ============================
	- step 1 						DONE
	- step 3 						DONE
	- step 4a						DONE
	- step 4b						DONE
	- step 5 						DONE
	- Export R matricies			DONE
	- Create Query					DONE
- Patient 3 ============================
	- step 1 						DONE
	- step 3 						DONE
	- step 4a						DONE
	- step 4b						DONE
	- step 5 						DONE
	- Export R matricies			DONE
	- Create Query					DONE
- Patient 4 ============================
	- step 1 						DONE
	- step 3 						DONE
	- step 4a						DONE
	- step 4b						DONE 	** incorrect number of clusters when remerging sub object (rerun 05/16 -- DONE)	
	- step 5                        DONE
	- Export R matricies			DONE
	- Create Query					DONE
- Patient 5 ============================
	- step 1 						DONE					
	- step 3 						DONE					
	- step 4a						DONE
	- step 4b						DONE
	- step 5 						DONE
	- Export R matricies			DONE
	- Create Query					DONE
- Patient 6 ============================
	- step 1 						DONE		
	- step 3 						DONE			
	- step 4a						DONE		    	
	- step 4b						DONE 	** incorrect number of clusters when remerging sub object(rerun 05/16 -- )
	- step 5 						DONE
	- Export R matricies			DONE
	- Create Query					DONE
	


## Running Model on All Time Points of Data to See if eML1 and eML2 are detected at infusion

```R
library("Seurat")
library("dplyr")

object <- readRDS(file =  "56546_9792_overnightplusinvitroSeuratfile122222.rds")

Assays(object) # make sure ADT data is in object

mtxC <- GetAssayData(object, assay = "ADT", slot = "counts")
write.csv(mtxC, sprintf("ADTcounts.csv"))


object <- RunUMAP(object, dims = 1:15) # not sure about these dims
umap <- object[["umap"]]@cell.embeddings
umap <- as.data.frame(umap)
write.csv(umap, "umapcoordinates.csv")
```






