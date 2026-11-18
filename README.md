# A Tutorial on Processing and Assembling 2bRAD data with STACKS
 
## Before starting
 
Read the following:
 
- [STACKS manual](https://catchenlab.life.illinois.edu/stacks/manual/)
- [Tutorial on what the major assembly parameters do](https://catchenlab.life.illinois.edu/stacks/param_tut.php)
- [Paris et al. (2017)](https://doi.org/10.1111/2041-210X.12775)
## Introduction
 
Assembly of RAD-seq data in Stacks requires multiple steps. Some are optional and others are required for adequate assembly. Whether you do optional steps will depend on what types of RAD-seq your data were generated with. I highly recommend considering what types of data will allow you to answer your questions when deciding on RAD-seq lab protocols.
 
**This tutorial is for 2bRAD data.**
 
This tutorial assumes a basic knowledge of working at the command line. If you need to submit jobs to a queue system, then you should be familiar with how to do that on your cluster. STACKS should be installed on the computer/cluster you will be using. This tutorial also assumes you have read the STACKS manual. So, you should know about the different STACKS programs/modules/scripts. You should also understand file formats for standard files like the popmap.

## Explanation of dependencies
Scripts used by this pipeline, specifically 2bRAD_stacks_prep.sh are in the helper-scripts folder.

You must install STACKS (https://catchenlab.life.illinois.edu/stacks/) and have executables in your PATH.

SeqKit (https://bioinf.shenwei.me/seqkit/) is required and should also be in your PATH. Can be installed via Conda 

```bash
conda install -c bioconda seqkit
```
 
## Steps
 
### 1) Demultiplex RAD-seq data (REQUIRED, but potentially done by sequencing facility or collaborator)
 
1. STACKS module `process_radtags` will demultiplex your data.
2. If a sequencing facility demultiplexes your data, you will likely still want to do some pre-processing to remove bad-quality reads and possibly adaptors. However, adequate study design (i.e., size selection) will ensure adaptors are not sequenced unless you are using 2bRAD.
3. I recommend not having the sequencing facility demultiplex your data for RAD-seq projects.
4. The following code is for demultiplexing 2bRAD data. For this tutorial, the file names are simple, but the file names from the sequencer are likely to be more complex.
   1. First make a copy of your sequencing file. This is "dummy" data because `process_radtags` assumes that dual-indexed data is paired end.
```bash
cp reads.R1.fastq.gz reads.dummy.fastq.gz
```
 
   2. Run `process_radtags`:
```bash
mkdir demultiplexed-reads
process_radtags -1 reads.R1.fastq.gz -2 reads.dummy.fastq.gz -b barcodes.txt -o ./demultiplexed-reads -c -q -r --index-index --retain-header --disable-rad-check --threads 20
```
 
   You will want to adjust output directories, file names, and `--threads` as appropriate for your data and machine.
 
### 2) Quality filtering
 
1. For 2bRAD data, scripts exist to trim reads that lack the 2bRAD cut sites.
2. For 2bRAD data, reads should also be trimmed to the RAD-locus to remove adapters.
3. `process_radtags` will **NOT** recognize 2bRAD cut sites because they are in the middle of the RAD-locus.
4. Run `2bRAD_stacks_prep.sh` in the folder with your demultiplexed reads. Adjust the script with the number of cores you have available and want to use.

### 3) Initial assembly to identify best-possible -m, -M and -n parameters (REQUIRED)
 
1. For moderately deep sequenced samples, I recommend `-m 3`, which means a minimum read depth of 3X.
2. Select a subsample of specimens. These should be individuals with an average number of reads and 2-3 individuals per population/site. The number of reads per individual can be found in the `process_radtags` log.
3. Run `denovo_map.pl` with a range of `-M` and `-n` values. The script `parameter-test_2bRAD.sh` is a SLURM batch script for doing this.
4. If you ran `denovo_map.pl` with `-X populations:"-R 0.80"` then you do not need to do the next step.
   1. Otherwise, run `populations` on each assembly with `-R 0.80` as the filter.
5. Run `paris_optimize.py` to determine best assembly parameters.
   1. Alternatively, you can do this by hand: examine `populations.log.distribs` for each metric and determine best assembly parameters following Paris et al. (2017). This is somewhat subjective.
### 4) Full assembly using chosen -m, -M, and -n parameters from above
 
### 5) Run populations with -R 0.80
 
1. Choose minimum minor allele frequency and missing data cutoffs.
2. Identify "bad apples" or individuals with high proportion of missing data. The definition of "high" is subjective and up to you, but some analyses are more sensitive to missing data than others.

### 6) Run populations for final times with a popmap that has bad apples removed
 
1. I recommend outputting a single-SNP-per-locus dataset and multiple-SNP-per-locus dataset. Requires running `populations` twice.
2. Filters should use `-p #ofCollectionSites -r #`. `-r` can vary, but 0.5-0.8 seems to usually work well.
3. Choose other filters as well like `--min-mac`, `--max-obs-het`.
4. Choose file outputs.
5. Much quicker than full assembly.

### 7) Perform population genomics analyses!
 
