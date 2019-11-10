# SeqTech2019 
## CSHL 2019 Advanced Sequencing Technologies and Applications Course Materials
-------


### [Link to pre-baked data](http://34.239.1.158/workspace/Preall/SeqTech2018.tar.gz) (415MB .tar.gz file)

-------


This tutorial will be a guide through the first few steps of primary data analysis:
1. FASTQ generation with `cellranger mkfastq`
2. Mapping and count matrix generation with `cellranger count`
3. Combining two samples into a shared, normalized matrix with `cellranger aggr`

-------
*Example experiment:*

5' Gene Expression from murine intestinal lamina propria, gathered by 2018 SeqTech students:

**Lamina propria**
<img src="https://github.com/jpreall/SeqTech2019/blob/master/images/Oral_mucosa.png" width="300">

From [wikipedia:](https://en.wikipedia.org/wiki/Lamina_propria)
> The connective tissue of the lamina propria is loose and rich in cells. The cells of the lamina propria are variable and can include fibroblasts, lymphocytes, plasma cells, macrophages, eosinophilic leukocytes, and mast cells.[2] It provides support and nutrition to the epithelium, as well as the means to bind to the underlying tissue. Irregularities in the connective tissue surface, such as papillae found in the tongue, increase the area of contact of the lamina propria and the epithelium.[3]

We have intestinal lamina propria harvested from two groups of mice:
 * Normal chow
 * High-fat lard-based diet
 
Gross.  Let's dive in:

### Illumina sequencing output
*(This is taken care of for you this year.  But here is some useful information about this step anyway, in case it's ever your responsibility to do the FASTQ generation step.):*

**Congratulations!  You didn't screw up an experiment.  Now you might have data.**

You submitted your 10X libraries to your sequencing core and the run completed successfully.  If your core is nice enough to provide an Illumina quality score plot as part of your data delivery, it might look something like this:

![QC images](https://github.com/jpreall/SeqTech2019/blob/master/images/QC_scores.png "Example QC data from NextSeq")

Don't let those ugly spikes at the end of R1 and going on through the beginning of R2 worry you.  This is very typical, and comes from some common but tolerable artifacts of the library prep and sequencing.  If you look closely (or use a tool like [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)), you'll even be able to determine the sequence of an abundant "contaminating" signature at the start of Read2:

```bash
AAGCAGTGGTATCAACGCAGAGTACATGGG
```
As it turns out, this is the sequence of the 10X Template Switch Oligo (TSO).  Lots of reads contain with the TSO, due to artifacts of the library chemistry.  [Don't Panic.](https://en.wikipedia.org/wiki/Don%27t_Panic_(The_Hitchhiker%27s_Guide_to_the_Galaxy))

## Creating 10X-compatible FASTQ files with `cellranger mkfastq`
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

## Primary data analysis with `cellranger count`

For each sample in your experiment, you'll need to run `cellranger count`.  [(detailed map of the pipeline)](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/map/cr-counter)

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

For instance, this is what the BAM (possorted_genome_bam.bam) file looks like:
```bash
$samtools view possorted_genome_bam.bam 19:5795690-5802672 | head -n 1
NB551387:106:HFL3VBGX9:3:23601:23610:10303	0	19	5795691	255	56M	*	0	0	GAGAGATATTTAGTTTTTATTTCATAAAATCAAAGTATTCAATAAATAGTAGGAGG	AAAAAEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE6AEEAEEEEEEEE	NH:i:1	HI:i:1	AS:i:55	nM:i:0	TX:Z:ENSMUST00000172812,+6926,56M;ENSMUST00000173314,+500,56M	GX:Z:ENSMUSG00000092341	GN:Z:Malat1	RE:A:E	BC:Z:TTCCCGAC	QT:Z:AAAAAEEE	CR:Z:GTCTTCGCATCCTTGC	CY:Z:AAAAAEEEEEEEEEEE	CB:Z:GTCTTCGCATCCTTGC-1	UR:Z:GAACATCTTA	UY:Z:EEEEEEEEEE	UB:Z:GAACATCTTA	RG:Z:SeqCourse2018-10XGEX-LPLard-G3:MissingLibrary:1:HFL3VBGX9:3
```
*Note: The pre-baked data that I linked above doesn't include the BAM file, since it's about 13GB in size.*

If you want to learn about what all these columns and tags mean, check out [this guide](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam)

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
![Loupe snapshot](https://github.com/jpreall/SeqTech2019/blob/master/images/SeqTech_Loupe_Example.png "Your awesome Loupe file")

That's ok, we all do it.  But seriously, go back to the [Web Summary](https://github.com/jpreall/SeqTech2019/blob/master/files/web_summary.html).  We're going to talk over what all those values mean in class. 

*Instrumental Break*

#### Sequencing Saturation
The single most useful piece of information stored in this summary is the estimate of sequencing saturation.  This will tell you how deeply you have sequenced these libraries, and whether it would be worth your time and money to add additional lanes of sequencing to identify new transcripts improve count numbers for differential expression.  

Total saturation is listed on the summary page, with a more thorough view in the `analysis` tab:

<img src="https://github.com/jpreall/SeqTech2019/blob/master/images/SeqTech_Saturation.png" width="500">

Putting it all together, we can see these libraries are actually pretty lousy:
 * Median UMIs per cell: 1,218
 * Median genes per cell: 601
 * Sequencing saturation: 72.4%
 
This is kind of typical for 5' Gene Expression libraries on lymphocyte-rich samples.  T- and B-cells are transcriptionally quiescent, and always contain far fewer unique transcripts than other cell types.  

#### The Data Matrix
Two other files that will be of extreme value to you are the actual data matrices.  Cellranger packages what would otherwise be an enormous data file into a clean, compressed hierarchical data (HDF5) file format.  It creates two versions: one that has been filtered of "empty" cells based on its [filtering algorithm](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/algorithms/overview), and the raw matrix containing all 3M+ barcodes and their associated gene counts, regardless if they are likely to contain valid cells or not.

 * filtered_feature_bc_matrix.h5
 * raw_feature_bc_matrix.h5

The filtered and raw matrices are both also stored under a separate matrix market exchange (.mtx) file format along with separate .csv files listing the cell barcodes and gene names, which can be stiched together into a unified data matrix.  Look for these folders under `filtered_feature_bc_matrix` and `raw_feature_bc_matrix`, respectively.  You can open any of these files with [Seurat](https://satijalab.org/seurat/), or [Scanpy](https://scanpy.readthedocs.io/en/latest/), or potentially other 3rd party analysis packages.   


## Combining samples with Cellranger aggr

Let's combine both the control and the lard diet samples into a unified data matrix.  
Be careful not to accidentally bait a computation scientist into discussing the relative merits of the many different strategies for aggregating multiple data sets.  You will have to gnaw your foot off before they finally get to the punchline: 
*there is no single best way to jointly analyze multiple datasets*

10X Genomics has decided to sidestep the issue by providing a simple data aggregation pipeline that takes a conservative approach as a first step, and leaving the more sophisticated steps in your capable hands.  `cellranger aggr` combines two or more datasets by randomly downsampling (discarding) reads from the richer datasets, until all samples have approximately the same median number of reads per cell.  This helps mitigate one of the simplest and easy to fix batch-effects caused by sequencing depth, but will not correct for the zillions of other variables that injected unintended variation to your samples. 


#### Create an aggr.csv file:

First, tell cellranger which samples to aggregate by creating an aggr.csv file formatted thusly:

```bash
library_id,molecule_h5
LP_Lard,/path/to/LP_Lard/outs/molecule_info.h5
LP_Control,/path/to/LP_Control/outs/molecule_info.h5
```
`cellranger aggregate` uses the `molecule_info.h5` file as the primary data source to do its downsampling.  This file contains rich data about each unique cDNA detected, including the number of duplicated or redundant reads mapping to a common UMI.  It is cleaner to downsample sequencing data based on this data rather than a simplified count matrix, which has discarded any information about the library complexity, PCR duplications, etc.  Cellranger uses this richer data source, but other tools seem to work with the final count matrix just fine.  Again, don't ask a bioinformation about it if you have children to feed some time today.

#### Run cellranger aggr:

```bash
cellranger aggr --id=SeqTech2018_LP_combined \
	--localcores=12 \
	--csv=/path/to/aggr.csv \
	--normalize=mapped
```
`cellranger aggr` is significantly less memory and cpu intensive than `cellranger count`.  If you are aggregating only a few samples, this should take less than an hour.  

Once it's done, you can view the web summary to see what was done to normalize the libraries:

<img src="https://github.com/jpreall/SeqTech2019/blob/master/images/aggr_web_summary.png" width="800">

Let's take a look at that aggr Loupe file.  Each sample is now stored as a separate category under "LibraryID":

<img src="https://github.com/jpreall/SeqTech2019/blob/master/images/aggr_tsne.png" width="500">
