# SeqTech2019 
## CSHL 2019 Advanced Sequencing Technologies and Applications Course Materials

This tutorial will be a guide through the first few steps of primary data analysis:
1. FASTQ generation with `cellranger mkfastq`
2. Mapping and count matrix generation with `cellranger count`
3. Combining two samples into a shared, normalized matrix with `cellranger aggr`
-------

### Generating a properly formatted set of FASTQs with cellranger mkfastq
**Congratulations!  You didn't screw up an experiment.  Now you might have data.**

You submitted your 10X libraries to your sequencing core and the run completed successfully.  If your core is nice enough to provide an Illumina quality score plot as part of your data delivery, it might look something like this:

![QC images](https://github.com/jpreall/SeqTech2019/blob/master/images/QC_scores.png "Example QC data from NextSeq")

Don't let those ugly spikes at the end of R1 and going on through the beginning of R2 worry you.  This is very typical, and comes from some common but tolerable artifacts of the library prep and sequencing.  If you look closely (or use a tool like FASTQC), you'll even be able to determine the sequence of an abundant "contaminating" signature at the start of Read2:

```bash
AAGCAGTGGTATCAACGCAGAGTACATGGG
```
As it turns out, this is the sequence of the 10X Template Switch Oligo (TSO)

### Mapping your reads with Cellranger mkfastq
You probably won't have to do this part yourself, but you might have to instruct your NGS core on how to generate properly formatted FASTQ files that will plug nicely into the subsequent `count` pipeline.  

```bash
cellranger mkfastq \
	--localcores=12 \
	--run=/path/to/basecalls/ \
	--samplesheet=/path/to/SampleSheet.csv \
```
**What the heck is a SampleSheet?**
A sample sheet tells the FASTQ generation pipeline how to break the reads out into separate folders based on their 8bp sample indices, and then how to name these files and organize them into folders.  Your Illumina sequencer will generate a SampleSheet.csv as part of the data generation process, but you may need to modify it in order to build a FASTQ file for 10X Genomics workflows:

**Example SampleSheet.csv**

```bash
[Header],,,,,,,
IEMFileVersion,4,,,,,,
Date,1/24/18,,,,,,
Workflow,GenerateFASTQ,,,,,,
Application,NextSeq FASTQ Only,,,,,,
Assay,TruSeq HT,,,,,,
Description,,,,,,,
Chemistry,Amplicon,,,,,,
,,,,,,,
[Reads],,,,,,,
26,,,,,,,
56,,,,,,,
[Settings],,,,,,,
,,,,,,,
[Data],,,:,,,,
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
300180,SeqCourse2018-10XGEX-LPLard,,,SI-GA-A1,SI-GA-A1,SeqTech2018,
300183,SeqCourse2018-10XGEX-LPcontrol,,,SI-GA-A2,SI-GA-A2,SeqTech2018,
```


10X Genomics uses a clever trick to make barcoding simple.  Each sample barcode is actually a carefully chosen pool of four unique 8bp indices.  This is why you won't see anything that looks like `AAGTCTGA` in these sample sheets, but rather something that looks like `SI-GA-C7`

If you use `cellranger mkfastq` to generate your FASTQs, it will translate `SI-GA-A2` into a set of four 8bp barcodes based on this [table provided by 10X Genomics.](https://s3-us-west-2.amazonaws.com/10x.files/supp/cell-exp/chromium-shared-sample-indexes-plate.csv).   For example:

```bash
SI-GA-A1,GGTTTACT,CTAAACGG,TCGGCGTC,AACCGTAA
SI-GA-A2,TTTCATGA,ACGTCCCT,CGCATGTG,GAAGGAAC
SI-GA-A3,CAGTACTG,AGTAGTCT,GCAGTAGA,TTCCCGAC
...etc
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
### What are all those files?

#### ...I1... = 8bp Sample index read:
```bash
@NB551387:259:H25HLBGXC:1:11101:7081:1841 1:N:0:GTTGCAGC
GTTGCAGC
+
AAAAAEEE

```
These files contain the Sample Index read information.  Note that the Cellranger mkfastq pipeline also writes the error-corrected sample index into the name line of the corresponding R1 and R2 reads.  See below.  

-------

#### ...R1... = Illumina Read 1.  This contains the cell barcode and UMI:
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
$ head -n 4 /path/to/CellRanger/cellranger-3.1.0/cellranger-cs/3.1.0/lib/python/cellranger/barcodes/chromium_whitelist_3M_2018.txt 
AAACCCAAGAAACACT
AAACCCAAGAAACCAT
AAACCCAAGAAACCCA
AAACCCAAGAAACCCG
```

-------



#### ...R2... = Illumina Read 2.  This contains the gene body read (or Cell Hash / CITE-seq tags):

Our current sequencing method of choice is the NextSeq500 HighOutput SE75 flow cell.  These kits allow us to do 56 base pair reads for the gene body.  This is because the 75 cycle kit actually ships enough extra reagent to sequence the Illumina indices, plus a little more to spare.  In reality, we can squeeze out 92bp from a SE75 kit. 

**56 (gene body) + 28 (cell barcode / UMI) + 8 (Sample Index) = 92 bases**
 
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

For each sample in your experiment, you'll need to run `cellranger count`.  

`cellranger count` will map the reads to the reference genome that you specify and count digital gene expression according to the transcriptome model that was used during building of the reference.  In most cases, we use the pre-built references for the [human and mouse genomes provided by 10X:](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)

`cellranger count` is widely used to generate the primary count table, and has become a de facto standard despite there being a few alternatives.  It requires a linux machine with 8GB memory per CPU (to hold the human reference genome in memory), and is quite astoundingly slow and prone to crashing.  Mapping and counting a single sample typically takes 4-8 hours on a single multi-core node with 12-16 CPU cores and 128GB memory. 

If you have a HPC cluster for distributed computing, good for you.  You'll probably get to know your sysadmin really well as you repeatedly crash jobs, hog memory, and wreak havoc on your HPC cluster.  

Unfortunately, Cellranger is still the best option for generating a robust, industry-standard count matrix plus easily shareable results files and a handsome QC summary HTML page.  If you are a power user and simply want an accurate count matrix as quickly as possible, consider switching to [STAR](https://github.com/alexdobin/STAR).  Cellranger actually relies on an older verion of STAR for the basic transcript alignment, but the newest versions of STAR supports generating count matrices from single cell libraries such as 10X Genomics, or any custom barcoding format of your own design.  You won't get the nice shareable Loupe file, but you will get your data after lunch instead of tomorrow morning.  

Here is how to run Cellranger in local mode, if you have only a single workstation computer or if you want to restrict the run to a single node on your cluster.  I find this is slow but reliable:

```bash
#!/bin/sh

SAMPLE=SeqCourse2018-10XGEX-LPLard
TRANSCRIPTOME=/path/to/transcriptome/folder/


cellranger count \
  --id=$SAMPLE \
  --jobmode=local \
  --localcores=12 \
  --transcriptome=$TRANSCRIPTOME \
	--fastqs=/path/to/folder/containing/your/fastqs/ \
	--sample=$SAMPLE \
  
  ```
In the FASTQ example above, the sample ID specified was `SeqCourse2018-10XGEX-LPLard`.  `cellranger count` will search through your FASTQ folders to find files whose names match this sample ID, and automatically recognize Read1 vs Read2 based on the file name.  If you ran multiple flowcells and wish to combine the data for additional depth, provide two paths to each fastq folder,separated by a comma with no spaces:
  
```bash
  	--fastqs=/path/to/fastq/folder1/,/path/to/fastq/folder2/
```
Make sure you used the same sample ID when preparing the fastq file from both flowcells.

### Congratulations, you successfully ran cellranger count!
How many tries did it take you?

Your results will now be stored in a folder called `Sample_ID/outs`.  There will also be a bunch of other files containing diagnostic information about the run that you can dig through if you are a masochist.  

In your `outs/` folder, you should see these files/folders:

```console
analysis/
filtered_feature_bc_matrix/
raw_feature_bc_matrix/
cloupe.cloupe
filtered_feature_bc_matrix.h5
metrics_summary.csv
molecule_info.h5
possorted_genome_bam.bam
possorted_genome_bam.bam.bai
raw_feature_bc_matrix.h5
web_summary.html
```
The first thing you should look at is the `web_summary.html`:
![Web Summary](https://github.com/jpreall/SeqTech2019/blob/master/images/CR_web_summary.png "Web Summary Preview")
   
Who am I kidding, the first thing you did was download and view the pretty Loupe file
