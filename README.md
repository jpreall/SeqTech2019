# SeqTech2019 
## CSHL 2019 Advanced Sequencing Technologies and Applications Course Materials

### Mapping your reads with Cellranger mkfastq
```bash
cellranger mkfastq \
	--localcores=12 \
	--run=$BCL \
	--samplesheet=$BASEDIR/SampleSheet.csv \
```

-------

These are the files produced by Cellranger mkfastq from a NextSeq500 sequencing run:

```bash
$ ls /path/to/fastqs/

SeqCourse2018-10XGEX-LPLard_S1_L001_I1_001.fastq.gz  SeqCourse2018-10XGEX-LPLard_S1_L003_I1_001.fastq.gz
SeqCourse2018-10XGEX-LPLard_S1_L001_R1_001.fastq.gz  SeqCourse2018-10XGEX-LPLard_S1_L003_R1_001.fastq.gz
SeqCourse2018-10XGEX-LPLard_S1_L001_R2_001.fastq.gz  SeqCourse2018-10XGEX-LPLard_S1_L003_R2_001.fastq.gz
SeqCourse2018-10XGEX-LPLard_S1_L002_I1_001.fastq.gz  SeqCourse2018-10XGEX-LPLard_S1_L004_I1_001.fastq.gz
SeqCourse2018-10XGEX-LPLard_S1_L002_R1_001.fastq.gz  SeqCourse2018-10XGEX-LPLard_S1_L004_R1_001.fastq.gz
SeqCourse2018-10XGEX-LPLard_S1_L002_R2_001.fastq.gz  SeqCourse2018-10XGEX-LPLard_S1_L004_R2_001.fastq.gz

```

-------


**...I1... = 8bp Sample index read:**
```bash
@NB551387:259:H25HLBGXC:1:11101:7081:1841 1:N:0:GTTGCAGC
GTTGCAGC
+
AAAAAEEE

```
These files contain the Sample Index read information.  Note that the Cellranger mkfastq pipeline also writes the error-corrected sample index into the name line of the corresponding R1 and R2 reads.  See below.  

-------

**...R1... = Illumina Read 1.  This contains the cell barcode and UMI:**
```bash
@NB551387:259:H25HLBGXC:1:11101:7081:1841 1:N:0:GTTGCAGC
CATCAAGGTCAGATAAGGTCGATCCGTT
+
AAAAAEEEEEEEEEEEEEEEE/EEEEE<

```
In the most recent version of the chemistry, this will be a 28bp read. 
  * Bases 1-16: Cell Barcode
  * Bases 17-28: UMI

```bash
CATCAAGGTCAGATAAGGTCGATCCGTT
< Cell Barcode >
                <    UMI   >
```

The cell barcodes are only accepted if they are a close match to the curated list of possible barcodes generated by the synthesis chemistry. This list is buried within the Cellranger directory:

```bash
$ head -n 4 chromium_whitelist_3M_2018.txt 
AAACCCAAGAAACACT
AAACCCAAGAAACCAT
AAACCCAAGAAACCCA
AAACCCAAGAAACCCG
```

-------



**...R2... = Illumina Read 2.  This contains the gene body read (or Cell Hash / CITE-seq tags):**

Our current sequencing method of choice is the NextSeq500 HighOutput SE75 flow cell.  These kits allow us to do 56 base pair reads for the gene body.  This is because the 75 cycle kit actually ships enough extra reagent to sequence the Illumina indices, plus a little more to spare.  In reality, we can squeeze out 92bp from a SE75 kit. 

** 56 (gene body) + 28 (cell barcode / UMI) + 8 (Sample Index) = 92 bases
 
```bash
@NB551387:259:H25HLBGXC:1:11101:7081:1841 2:N:0:GTTGCAGC
CTGCAAACCATCTCCTGTGCAGGGTCCTGCTGGCACCATGGTCTCACAGCCACCCG
+
A//AAEEE///EA<////</E/</EEAA/AA//EEE<<A<//EAA/6A/<EEEA</
```
Officially, 10X recommends quite long reads to map the gene body:

|Read	|Read 1	|i7 Index	|i5 Index	|Read 2|
|------------ |:------------:| -------:| -------:| --------:|
|Purpose	|Cell barcode & UMI	|Sample Index	|N/A	|Insert|
|Length	|28	|8	|0	|91|

**Note from 10X:**
* *Shorter transcript reads may lead to reduced transcriptome alignment rates. Cell barcode, UMI and Sample index reads must not be shorter than indicated. Any read can be longer than recommended. Additional bases in Sample index reads must be trimmed using cellranger mkfastq or Illumina's bcl2fastq prior to further analysis. Additional bases in Cell barcode or UMI reads will automatically be ignored by Cell Ranger 1.3 or later.* *

-------

### Primary data analysis with Cellranger Count 

For each sample in your experiment, you'll need to run Cellranger count.  This will map the reads to the reference genome that you specify and count digital gene expression according to the transcriptome model that was used during building of the reference.  In most cases, we use the pre-built references for the [human and mouse genomes provided by 10X:](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)


```bash
#!/bin/sh

SAMPLE=My_Sample
TRANSCRIPTOME=/path/to/transcriptome/folder/


cellranger count \
  --id=$SAMPLE
  --jobmode=local
  --localcores=12
  --transcriptome=$TRANSCRIPTOME \
	--fastqs=$BASEDIR/mkfastq/$FLOWCELL/outs/fastq_path \
	--sample=$SAMPLE \
  
  ```
  
