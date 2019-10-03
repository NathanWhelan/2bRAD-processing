#This pipeline is for the processing of zipped illumina single end 75 bp files for AlfI 2b-RAD data.
#This pipeline was put together by Matt Galaska with the currently available skynet programs in mind. It has been modified by N. Whelan in Sept 2019 with some different scripts

###IT IS YOUR RESPONSIBILITY TO MAKE SURE SCRIPTS ARE IN YOUR PATH OR WORKING DIRECTORY
###This script requires BioPerl.
###This script requires GNU Parallel

####The AlfIExtract_mod.pl is a modified version by N. Whelan to recognize that first and base pairs were trimmed

#First, I suggest that you move all zipped files into a directory called /raw or other appropriate designation.
#This pipeline will unzip all files, truncate the files down to the AlfI contig size (minus a bp on either side where N's are more abundant).
#It will then discard any poor quality sequences and finally pull out only sequences containing the AlfI restriction site.
#Intermediate files are cleaned throughout the pipeline to keep directory size down as these are large files.

############################
##USES 10 cores by Default##
CORES=8  ##Can be changed##
############################

mkdir raw_backup
cp *.gz raw_backup/
rm Undetermined*.gz ##Removes file with reads that were not demultuiplexed correctly
for FILENAME in *.gz
        do
        gunzip $FILENAME
done


ls *.fastq | parallel -j $CORES 'TruncateFastq.pl -i {} -s 2 -e 36 -o {}.out' ##Elim Meyer Scripts  ##Modified to not trim last base as last base of Alf1 site as it is not the last sequenced base pair with 75b
rm -rf *.fastq
rename .fastq.out .fastq *.fastq.out
ls *.fastq |parallel -j $CORES 'perl /usr/local/genome/bin/QualFilterFastq.pl -i {} -m 20 -x 5 -o {}.qual' ##Fix directory if necessary
rm -rf *.fastq
rename .fastq.qual .fastq *.fastq.qual
ls *.fastq | parallel -j $CORES 'fastq_to_fasta {} > {}.fasta' #SHRiMP_2_2_3
rm -rf *.fastq
rename .fastq.fasta .fasta *.fastq.fasta


ls *.fasta | parallel -j $CORES 'AlfIExtract_mod.pl {} {}.alf.fasta'
rename fasta.alf.fasta .fasta *fasta.alf.fasta
