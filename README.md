# SeqTech2019 
## CSHL 2019 Advanced Sequencing Technologies and Applications Course Materials

These are the files produced by a NextSeq500 sequencing run:

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


```bash
$ head chromium_whitelist_3M_2018.txt 
AAACCCAAGAAACACT
AAACCCAAGAAACCAT
AAACCCAAGAAACCCA
AAACCCAAGAAACCCG
AAACCCAAGAAACCTG
AAACCCAAGAAACGAA
AAACCCAAGAAACGTC
AAACCCAAGAAACTAC
AAACCCAAGAAACTCA
AAACCCAAGAAACTGC
```

-------



**...R2... = Illumina Read 2.  This contains the gene body read (or Cell Hash / CITE-seq tags):**
```bash
@NB551387:259:H25HLBGXC:1:11101:7081:1841 2:N:0:GTTGCAGC
CTGCAAACCATCTCCTGTGCAGGGTCCTGCTGGCACCATGGTCTCACAGCCACCCG
+
A//AAEEE///EA<////</E/</EEAA/AA//EEE<<A<//EAA/6A/<EEEA</
```

-------









```bash
#!/bin/sh

SAMPLE=
TRANSCRIPTOME=/path/to/transcriptome/folder/


cellranger count \
  --id=$SAMPLE
  --jobmode=local
  --localcores=12
  --transcriptome=$TRANSCRIPTOME \
	--fastqs=$BASEDIR/mkfastq/$FLOWCELL/outs/fastq_path \
	--sample=$SAMPLE \
  
  ```
  
