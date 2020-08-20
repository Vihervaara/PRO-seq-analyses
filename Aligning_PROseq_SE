#====================
#====================

# shell script for aligning PROseq data
# generates unnormalized and coverage normalized bedgraph files that report the 3'-end of each read (active site of transcription) as well as the 5'-end of each read that enrich at the TSS.
# This script uses combined genome of hg19 and dm3. This allows the usage of extern normalization control from Drosophila cells, spiked into the samples derived from human cells before the run-on reaction. To discern the chromosomes of human and fly genomes in the combined hg19-dm3 genome, please add '_dm3' postfix to the dm3 chromosome names.

# Required:
					      
#1. appropriate reference genome, generated with: bowtie2-build [options]* <reference_in> <bt2_base> where reference_in is a comma-separated list of FASTA files containing the reference sequences to be aligned. bt2_base is the base name for the indexed files to write.
# Here, the reference_in for bowtie2-build is the combined hg19-dm3 genome with renamed dm3 chromosomes.

#2. fastx_tools	 - can be downloaded from http://hannonlab.cshl.edu/fastx_toolkit/download.html
#3. samtools - can be downloaded from http://www.htslib.org/download/
#4. bowtie2 - can be downlaoded from http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
#5. bedGraphTobBgWig script - can be downloaded from http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/
#6. bedtools - can be downloaded from https://bedtools.readthedocs.io/en/latest/content/installation.html
					      				              
#=============================================================================================
 
# Usage

#1. Place fastq files to be analyzed in fastq folder within the working directory. The output files will be written in this working directory.
#2. Make sure the Paths are correct and if not, edit the Paths.
#3. Edit the List; Put the names of the fastq files (without the .fastq ext) to be analyzed in the List within "", seperated by a space.  
#4. In order to save the terminal output in a log file, type <script log_align_PROseq_xxx.txt>. this will create a log file named log_align_PROseq_xxx.txt in the working directory which will save terminal output until it is terminated by pressing <control D> in the command line.
#5. When ready, run the shell script by typing: sh align_PROseq_hg198_dm3_shell.sh

#=============================================================================================

#Paths

genomehg19dm3="/yourPATH/hg19_dm3"
## bowtie2 built hg19 and dm3 combined genome

chromSizeshg19="/yourPATH/hg19chromSizes.txt"
## two columns: chromosome name and sizes for all chromosomes

chromSizesdm3="/yourPATH/dm3chromSizes.txt"
## two columns: chromosome name and sizes for all chromosomes
## remember to rename the chromosome names as in the combined hg19-dm3 genome, e.g. with the '_dm3' extension.

bgTObigWig="/yourPATH/mkbigWig/bedGraphToBigWig"
## converts bedgraph file to biwig file.

#=============================================================================================

## make a list of your file names, fithout the .fastq extension

List="unC_Scr_HS0 unC_Scr_HS30 unC_Scr_HS60 unC_sh1_HS0 unC_sh1_HS30 unC_sh1_HS60 preC_Scr_HS0 preC_Scr_HS30_bs preC_Scr_HS60 preC_sh1_HS0 preC_sh1_HS30 preC_sh1_HS60"


for x in ${List}

do	
	
	echo converting fastq to fasta file from ${x} :
	fastq_to_fasta -n -v -i ./fastq/${x}.fastq -o ${x}.fasta -Q33
	## converting fastq to fasta: makes file size smaller to handle later and outputs the number of reads. This is also important because the fastx_collapser only outputs fasta file not fastq. -v reports number of sequences and -n keeps sequences with unknown (N) nucleotides, -Q33 for sanger coded ascii

	echo removing adapter sequence from ${x} :
	fastx_clipper -v -n -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC -l 12 -i ${x}.fasta -o ${x}_clipped.fasta -Q33
	## removes the above given adapter sequence (which is the immediate sequence after the insert all the up to the barcode) all of it or part of it, from the sequencing reads. Removing adapter as the first step also helps in removing PCR duplicates that may have a difference in nucleotide after the adapter sequence, which have low accuracy. -v is for verbose and -n keeps seqeuences with N.
				
	echo removing PCR duplicates and collapsing reads from ${x} :
	fastx_collapser -v -i ${x}_clipped.fasta -o ${x}_collapsed.fasta -Q33
	## fastx_collapser only outputs fasta file
	
### Traditional proseq 3' adapter does not include a molecular barcode.
### One used here does, containing a sample-specific 6 nt barcode, and a constant 'C' at the ligation site. If your sample contains the barcode, execute the following command:
	#echo trimming the 7bp molecular barcodes on left from ${x} :
	#fastx_trimmer -f 8 -v -i ${x}_collapsed.fasta -o ${x}_q.fasta  -Q33
    ## trims the 7bp molecular barcode and leaves the first base as the 8th base (-f 8), -v is for verbose
			
            
	echo making reverse complement of the clipped sequences of ${x}:
	fastx_reverse_complement -i ${x}_collapsed.fasta -o ${x}_RC.fasta -Q33
    #or, if the barcode was removed:
    #fastx_reverse_complement -i ${x}_q.fasta -o ${x}_RC.fasta -Q33
	## Takes the _q.fasta file and flips the read because the illumnia sequences the DNA from 5'end so the adapters in PRO-seq is reversed. this -Q33 parameter instructs to use Illumina encoded quality scores, not Sanger encoding
	
	
	echo aligning reads from ${x} to the combined genome of hg19 and dm3 with bowtie2 :
	bowtie2 -f --end-to-end -p 7 -x ${genomehg19dm3} -U ${x}_RC.fasta | awk '$2 != 4 {print}' | sed '/XS:/d' | samtools view -S -b '-' > ${x}.bam
	## aligns reads to hg19-dm3 using bowtie2. --end-to-end required the whole read to align. I.e. no soft clipping allowed to maintain the nucleotide-resolution of PRO-seq. -p 7 uses 7 cores, -x ${genomehg19dm3} provides the bowtie2 index of genome, and -U ${x}_clipped_RC.fasta provides the input file.
	## the output sam file is piped to awk which removes the read that has no reported alignments (flag 4 in 2nd column of sam file), and the sed '/XS:/d' command gets rid of multimapped reads.
	## Then, samtools view command converts sam file to bam file, -S is for input is SAM, -b is for output BAM.

	
	samtools sort ${x}.bam ${x}_sorted
	## no .bam extension in sorted file because the samtools sort adds .bam on it.
	
	samtools index ${x}_sorted.bam
	## these set of commands create sorted and indexed bam file to be used in IGV.
	
	echo converting bam to bed file and sorting it:
	bamToBed -i ${x}.bam | sort -k1,1 -k2,2n > ${x}_sorted_bed.tmp
	## the bam file is converted to bed by bamToBed, and the bed file is piped to be sorted by chromosome and then start position.
	
	hg19_chr=$(awk '{print $1}' ${chromSizeshg19})
	## making a variable called hg19_chr by taking the first column of the chromSizes files
	
	dm3_chr=$(awk '{print $1}' ${chromSizesdm3})
	## this chromSizesdm3 is modified so that the chr names are not similar to hg19. Add e.g. '_dm3' to the chrX, chr4, and chrM  of dm3 chromosome names

	### IMPORTANT Note:  The steps below write repeatedly to the separated bed files.  Therefore, if you rerun this script without deleting or moving these files, they will be added to instead of overwritten. #######

	for y in ${hg19_chr}
	do
		echo ${y}
		awk '$1=="'${y}'" {print $0}' ${x}_sorted_bed.tmp >> ${x}_hg19.tmp
		## if the first column, which is chr name, matches to the entry in the list, then it writes the entire line {print $0} to ${x}_hg19.tmp. >> appends output to existing file.
	done

	for y in ${dm3_chr}
	do
		echo ${y}
		awk '$1=="'${y}'" {print $0}' ${x}_sorted_bed.tmp >> ${x}_dm3.tmp
	done 

	sort -k 1,1 -k2,2n ${x}_hg19.tmp > ${x}_sorted.bed
	## sorting the ${x}_hg19.tmp file
	
	n=$(grep -c ^ ${x}_dm3.tmp)
	## n is calculated by summing the reads, which is equivalent to the number of lines in the dm3 bed file. Since we only retained uniquely mapping reads (see the bowtie2 line), only reads that uniquely mapped to dm3 are counted here. Also, if a read mapped to both hg19 and dm3, it was discarded by the same sed '/XS:/d' selection.
	
	echo Number of reads mapping to spikeIn dm3 genome is: $n
    # prints the number of reads that uniquely map to dm3 (spike-in) genome
	
	echo generating non-normalized bedgraph files of ${x}:
	awk '$6 == "+"' ${x}_sorted.bed | genomeCoverageBed -i stdin -3 -bg -g ${chromSizeshg19} > ${x}_pl.bedgraph
	awk '$6 == "-"' ${x}_sorted.bed | genomeCoverageBed -i stdin -3 -bg -g ${chromSizeshg19} > ${x}_m.bedgraph
	awk '{$4=$4*-1; print}' ${x}_m.bedgraph > ${x}_mn.bedgraph
	## the first line takes all the sequence in the sorted bed file which has "+" in the 6th column and makes bedgraph file by collapsing the reads to 3'end
	## the second line takes all the sequence in the sorted bed file which has "-" in the 6th column and makes bedgraph file by collapsing the reads to 3'end ( remember, for reads in "-" strand, the 5' end is the 3rd column and not 2nd column as for reads in "+"strand)
	## the third line converts the value of 4th column, which is the number of reads mapped at the position specified by first three columns, to negative value in order to specify that they are in the "-" strand.
	
	echo making bigwig from non-normalized bedgraph files of ${x}:
	${bgTObigWig} ${x}_pl.bedgraph ${chromSizeshg19} ${x}_pl.bigWig
	${bgTObigWig} ${x}_mn.bedgraph ${chromSizeshg19} ${x}_mn.bigWig
	## makes bigwig file from the bedgraph files



### These lines generate coverage normalized bedgraph files for the first look of the data.
### Please note that for a proper normalization of PRO-seq data, we utilize the Spike-in count or sample-intrinsic normalisation count (e.g. ends of long genes for <1h heat shock time points, see Mahat et al., 2016, Mol Cell or Vihervaara et al., 2017, Nat Commun).
### In the spike-in normalization, the n (here spike-in count) of each sample is compared to an 'anchor' dataset, generally chosen to be the control condition, here the non-stressed cells. By that way, one can retain the coverage-normalized file of the control, and adjust the levels of each sample to the control. The result is essentially, adjusted RPM normalization.

	#############################
	
	echo generating PRMnormalized bedgraphs files of ${x}:
	
	c=$(grep -c ^ ${x}_sorted.bed)
    ### c is the total number of uniquely mappig reads to the hg19
    
	echo $c | awk '{ c="'$c'"; printf "%s\t%s\t%s\t%s\n", $1, $2, $3, ($4*1000000)/c}' ${x}_pl.bedgraph  > ${x}_RPMnorm_pl.bedgraph
	echo $c | awk '{ c="'$c'"; printf "%s\t%s\t%s\t%s\n", $1, $2, $3, ($4*1000000)/c}' ${x}_mn.bedgraph  > ${x}_RPMnorm_mn.bedgraph
	## generates the coverage normalized bedgraph files by dividing the fourth column of bedgraph file with total number of mapped reads and multiplying the resulting number with 1M (RPM - Reads Per Million)
	
	echo making binary files from RPMnormalized bedgraph files of ${x}:
	${bgTObi} -i ${x}_RPMnorm_pl.bedgraph -o ${x}_RPMnorm_pl.bi -c ${chromInfohg19}
	${bgTObi} -i ${x}_RPMnorm_mn.bedgraph -o ${x}_RPMnorm_mn.bi -c ${chromInfohg19}
	## generates binary file from normalized bedgraph files for use in Hojoong's kbro

	#############################
	
	

	echo finished aligning and generating files for ${x}
	
done
	
rm *.tmp

rm *_m.bedgraph
gzip ./fastq/*.fastq
## gzip fastq files in the fastq folder

rm *.fasta

## terminate the script command by typing <control D> in the terminal
