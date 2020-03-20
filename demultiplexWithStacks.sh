###Basic Script for demultiplexing data with process_radtags with dual indexed, single-end Illumina data The script was desinged for 2bRAD data, but could be adapted for other RAD types.

#The text file should have your taxon names and barcodes. See the process_radtags documentation for more details.

#Process_radtags was not designed to demultiplex dual-barcoded single-end data, but there is a workaround. Create a dummy paired reads file by make a 
# copy of you multiplexed sequencing file. Name it something like RealSequenceDataFileName.dummy.fastq.gz
# Script assumes process_radtags is in your path.
# This script will output your data into an output folder with the name of your sequence data with ".output" appended


###TO run: demultiplexWithStacks.sh <sequenceFile.fastq.gz> <dummySequenceFile.gz> <barcodeFile.txt> <number of threads>
##EXAMPLE: demultiplexWithStacks.sh L_compacta.fastq.gz L_compacta.dummy.fastq.gz L_compactaBarcodes.txt 30



mkdir $1.output
process_radtags -1 $1 -2 $2 -o $1.output -b $3 -c -r -t $4 -D --index_index --disable_rad_check
