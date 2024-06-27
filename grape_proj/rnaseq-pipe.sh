
for i in 095 096 097 098 099 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118
do
	mkdir ./ranseq_results-ncbi/result-SRR9701${i}
	cd /home/liuk/pcc-project/01RNA-seq/ranseq_results-ncbi/result-SRR9701${i} 
	fastq-dump SRR9701${i}
	fastqc -f fastq SRR9701${i}.fastq
	java -jar /home/liuk/software/Trimmomatic-0.38/trimmomatic-0.38.jar SE -phred33 -threads 4 SRR9701${i}.fastq SRR9701${i}.trim ILLUMINACLIP:/home/liuk/software/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:80
	fastqc -f fastq SRR9701${i}.trim
	hisat2 -p 4 --dta -x /home/liuk/grape-data/RNA_practice/data/ref/Vitis_vinifera/hisat2_index/GCF_000003745.3_12X_genomic.fna_gi -U SRR9701${i}.trim -S SRR9701${i}.sam
	python /home/liuk/grape-data/shiulab/1.1RNA-seq_data_processing/02_keep_reads_with_quality_60_and_unique_mapping_SE.py SRR9701${i}.sam
	samtools sort SRR9701${i}_quality_60_unique.sam -o SRR9701${i}_sorted.bam
	stringtie -p 10 -G /home/liuk/grape-data/RNA_practice/data/ref/Vitis_vinifera/GCF_000003745.3_12X_genomic.gff -o /home/liuk/pcc-project/01RNA-seq/ranseq_results-ncbi/stringtie-b-result/SRR9701${i}.gtf -l SRR9701${i} -e SRR9701${i}_sorted.bam #

done

