# Chapter 3 Juvenile CADO RNASeq Analysis
**Author: Megan Guidry**

All analysis was done on our lab shared server, KITT, made by [J. Puritz](https://github.com/jpuritz) in /home/mguidry/3_Juvenile-CADO/

Programs Installed/Needed for this Project:  
- HISAT2 
- StringTie 
- gffcompare 
- fastp, fastQC, multiqc
- samtools

File Naming and Information:

location of reference genome on KITT: `/RAID_STORAGE2/Shared_Data/Oyster_Genome/masked/masked.cvir.genome.fasta`

- `masked.cvir.genome.fasta`: Eastern Oyster genome
- `ref_C_virginica-3.0_top_level.gff3`: annotation file for the Eastern Oyster, I also got this from Erin, it has the matching header line convention to work with the above genome
- `prepDE.py`: python script for converting read count information from StringTie into a matrix fo DESeq2, full code [here](http://ccb.jhu.edu/software/stringtie/dl/prepDE.py)
- ``: this is an example of the naming convention of the samples


General notes: 

[**Easy Guide to tmux**](https://hamvocke.com/blog/a-quick-and-easy-guide-to-tmux/)
```
#create a new session and run desired command
tmux new -s my_session
your_command

#to detach from the session and have things run in the background
Ctrl + b, then d 
##or type
tmux detach

#list tmux sessions
tmux ls

#reattach to session
tmux attach -t my_session

#kill session (if needed)
tmux kill-session -t my_session
```


-------
## Steps:

1. Set up 
2. Raw read QC
3. Trimming with fastp
4. Trimmed read QC
5. Mapping/Aligning reads to genome (HISAT2)
    - considered HISAT2 and STAR (both common RNAseq alignment tools)
    - leaning toward going with HISAT2 because it is much faster and has a decent % of reads aligned and transcriptome coverage (when compared to other aligners)
    - Musich et al 2021: https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2021.657240/full
6. Assemble mapped reads into transcripts (StringTie)
7. Prep files for differential expression analysis
8. DESeq2
9. Gene ontology
10. PCA
11. heatmap
12. WGCNA - maybe?

________


## Before you begin...
Create conda environment & install packages
```{bash}
conda create -n Ch3 hisat2 star stringtie gffcompare fastp fastQC multiqc samtools
conda activate Ch3
```
`conda list` will allow you to look at all of the packages you have installed in your environment.


## 1. Set up
Reads were already de-multiplexed and assigned to each individual sample by Novogene. To get started, we'll first link the data to the current working directory. 
Should have 72 samples each with a F and R read (144 *.fq.gz* files).

Create symbolic link to fastq files in RAIDSTORAGE2
```
cd home/mguidry/3_Juvenile-CADO
mkdir raw_reads
cd raw_reads

#symlink reads from individual sample folders into raw_reads directory
for i in /RAID_STORAGE2/Raw_Data/ASMFC_RNA/usftp21.novogene.com/01.RawData/*; do
    ln -s "$i"/*.fq.gz .
done

#check number of files in directory 
ll *fq.gz | wc -l
```

### Check out read counts 

**Purpose:** Look at read counts across samples prior to trimming.

Look at the read counts for each file using a code from [this website](http://www.sixthresearcher.com/list-of-helpful-linux-commands-to-process-fastq-files-from-ngs-experiments/) made into a for-loop that went through all the files. It outputs the filename and the number of reads in that file. Takes 1-2mins/sample.

```
for fq in *.fq.gz
do
    echo "File: $fq"
    reads=$(zcat "$fq" | wc -l)
    echo "Reads: $((reads / 4))"
done > number_of_raw_reads.txt
```

## 2. Raw read quality check
**Purpose:** Establish a baseline of what the read quality looks like to help determine trimming parameters in the next step.

Checked out the quality of the reads in each .fq.gz file using MultiQC & FastQC. We'll go back and run the same on the trimmed data in the next couple of steps too. But first, we want to get an idea of what the data look like before trimming. *Took roughly 7-8hrs.*

```
cd /home/mguidry/3_Juvenile-CADO
mkdir fastqc_raw
fastqc -o /home/mguidry/3_Juvenile-CADO/fastqc_raw/ raw_reads/*fq.gz   #took about 8 hrs for these data 
cd raw_reads 
ll *.fastqc.zip | wc -l #check to make sure we have 144 files (2 per sample - fwd & rev)
cd ../fastqc
multiqc . #output a multiqc report
```

Tranferred `multiqc_report.html` file to computer and opened in web browser.

### Raw reads QC report

#### Raw FastQC: Sequence Counts 
Notes:
* lots of duplicate reads across all samples 
* 3.5-8M unique reads across samples

![sequence-counts](images/raw_qc/fastqc_sequence_counts_plot-3.png)

#### Raw FastQC: Raw sequence quality 
Mean quality value across each base position in the read. 

Notes:
* looks good overall 
* everything is between 35-40 across the board 
* J1_14_SC3_1_2 is the one read that is a little lower - not sure why, but it still has great quality scores across the read

![per-base-seq-quality](images/raw_qc/fastqc_per_base_sequence_quality_plot-3.png)

#### Raw FastQC: Per Base Sequence Content 
The proportion of each base position for which each of the four normal DNA bases has been called.

Notes:
* little messy 1-10bp, but this is to be expected bc of the low diversity in the adapter sequence making it challenging to resolve colors
* fairly straight lines from 10bp on 
* plot below representative of most samples 

![per-base-seq-content](images/raw_qc/per-base-seq-content.png)

### NOTES on Raw QC report: 
- lots of read duplication - typical for oysters bc of repeats throughout the genome? and they are highly related? and expressing similar suites of genes?
- avg. GC content isn't really normally distrubuted, is this something to investigate?
- adapter content plot (adapter contamination?) - one sample looks funky (J1_14SC3_1_2)
- does anything in particular catch your attention that I'm missing?
- proceeding with trimming - game plan 

## 3. Trimming raw reads with [`fastp`](https://github.com/OpenGene/fastp)
**Purpose:** Cleaning up the raw reads to improve quality of sequences used in future analyses. 

**Features of fastp:**
1. filter out "bad" reads (too low quality, too short, or too many missing basecalls) - doesn't apply to these data?
2. cut low quality bases in each read 5'-3' determined by mean quality in a sliding window
3. cut adapters - automatically detected by fastp and enabled by default
4. trim reads at front and tail - if desired
5. correct mismatched bp in overlapped regions of paired end reads - if one base has quality and is paired with another over very low quality -- again, dont think this applies here bc high quality data
6. trim polyX tailing 
7. report JSON format results for futher interpreting
8. visualize QC and filtering results to html 

`fastp` flags:\
`-i` = input file, forward\
`-I` = input file, reverse\
`-o` = output file, forward\
`-O` = output file, reverse\
`-f` or `-t` = front or tail trimming setting, forward\
`-F` or `-T` = front or tail trimming setting, reverse (if not specified, reverse settings will be default set to whatever the forward parameters are)\
`-u` = how many percent of bases are allowed to be unqualified (integer value 0-100); default is 40%


**Running `fastp`**
```{bash}
mkdir trimmed_reads

#loop through all files
for fq in *.1.fq.gz
do
    # get sample name (everything before ".1.fq.gz")
    sample=$(basename "$fq" .1.fq.gz)

    fastp -i ${sample}.1.fq.gz -I ${sample}.2.fq.gz -o trimmed_reads/${sample}.F.trim.fq.gz -O trimmed_reads/${sample}.R.trim.fq.gz -f 9 -j trimmed_reads/${sample}.json -h trimmed_reads/${sample}.html
done

#one sample
#fastp -i sampleA.1.fq.gz -I sampleA.2.fq.gz -o trimmed_reads/sampleA.F.trim.fq.gz -O trimmed_reads/sampleA.R.trim.fq.gz  -f 9 -j trimmed_reads/sampleA.json -h trimmed_reads/sampleA.html
```

Notes to self:   
-f (front trim) - should I focus on the quality score or the base content to remove that noise?   
Adapters are trimmed by default   
basecalls with Q<15 are removed and 40% of unqualified bases are allowed by default   

**Anything else to consider when trimming RNAseq data??**