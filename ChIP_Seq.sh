### Reference Genome

nohup wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz &
tar -zxvf chromFa.tar.gz
cat *.fa > hg19.fa
bowtie2-build hg19.fa hg19


##################################################################################
#
# ---------------------------------- GSE19013 ---------------------------------- #
#
##################################################################################

#@ Step1: Get the data
nohup wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR035/SRR035186/SRR035186.fastq.gz &
nohup wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR035/SRR035187/SRR035187.fastq.gz &


#@ Step2: First FastQC & Fail
fastqc -o ./ -t 6 ../SRR035186.fastq.gz
fastqc -o ./ -t 6 ../SRR035187.fastq.gz


#@ Step3: First Trim
nohup java -jar /opt/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 10 -trimlog SRR035186_trim.log -summary SRR035186_summary.log SRR035186.fastq.gz SRR035186_trimed.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 &
nohup java -jar /opt/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 10 -trimlog SRR035187_trim.log -summary SRR035187_summary.log SRR035187.fastq.gz SRR035187_trimed.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 &


#@ Step4: Second FastQC & Pass
fastqc -o ./ -t 6 ../SRR035186_trimed.fastq
fastqc -o ./ -t 6 ../SRR035187_trimed.fastq


#@ Step5: Bowtie2 Alignment
nohup bowtie2 -p 10 -x hg19 -U SRR035186_trimed.fastq -S SRR035186_trimed.sam &
nohup bowtie2 -p 10 -x hg19 -U SRR035187_trimed.fastq -S SRR035187_trimed.sam &


#@ Step6: MACS2 Call Peak
nohup macs2 callpeak -t SRR035187_trimed.sam -c SRR035186_trimed.sam -f SAM -g hs --keep-dup 1 --outdir macs2_callpeak -n macs2 -B &>macs2_callpeak/MACS2.out &


#@ Step7: Get sequence of the peaks for motif searching
nohup bedtools getfasta -fi ../../hg19.fa -bed macs2_peaks.narrowPeak -fo macs2_peaks.fa &
nohup grep -c '^>' macs2_peaks.fa


### @ StepX: .sam file -> .bw file
#@ **PS** Commands after getting sam files:

# 1.	sam files --> bam files
samtools view -bS SRR015349.sam > SRR015349.bam
samtools view -bS SRR015350.sam > SRR015350.bam

# 2.	bam sort
samtools sort SRR015349.bam SRR015349_sort
samtools sort SRR015350.bam SRR015350_sort

# 3.	Build Index
samtools index SRR015349_sort.bam
samtools index SRR015350_sort.bam

# 4.	bamCoverage
bamCoverage -b SRR015349_sort.bam -o SRR015349_sort.bw
bamCoverage -b SRR015350_sort.bam -o SRR015350_sort.bw






##################################################################################
#
# ---------------------------------- GSE14664 ---------------------------------- #
#
##################################################################################


#@ Step1: Get the data
nohup wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR015/SRR015349/SRR015349.fastq.gz &
nohup wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR015/SRR015350/SRR015350.fastq.gz &


#@ Step2: First FastQC & Pass
nohup fastqc -o ./fastqc_out -t 6 ./SRR015349.fastq.gz &
nohup fastqc -o ./fastqc_out -t 6 ./SRR015350.fastq.gz &


#@ Step3: Bowtie2 Alignment
nohup bowtie2 -p 10 -x hg19 -U SRR015349.fastq.gz -S SRR015349.sam &
nohup bowtie2 -p 10 -x hg19 -U SRR015350.fastq.gz -S SRR015350.sam &


#@ Step4: MACS2 Call Peak
nohup macs2 callpeak -t SRR015350.sam -c SRR015349.sam -f SAM -g hs --keep-dup 1 --outdir macs2_callpeak -n macs2 -B &>macs2_callpeak/MACS2.out &


#@ Step5: Get sequence of the peaks for motif searching
nohup bedtools getfasta -fi ../../hg19.fa -bed macs2_peaks.narrowPeak -fo macs2_peaks.fa &
nohup grep -c '^>' macs2_peaks.fa


### @ StepX: .sam file -> .bw file
#@ **PS** Commands after getting sam files:

# 1.	sam files --> bam files
samtools view -bS SRR015349.sam > SRR015349.bam
samtools view -bS SRR015350.sam > SRR015350.bam

# 2.	bam sort
samtools sort SRR015349.bam SRR015349_sort
samtools sort SRR015350.bam SRR015350_sort

# 3.	Build Index
samtools index SRR015349_sort.bam
samtools index SRR015350_sort.bam

# 4.	bamCoverage
bamCoverage -b SRR015349_sort.bam -o SRR015349_sort.bw
bamCoverage -b SRR015350_sort.bam -o SRR015350_sort.bw

