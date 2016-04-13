#   This file is part of the analysis for the publication 'CRISPR-directed mitotic recombination 
#   enables genetic mapping without crosses'. Code written by Joshua S Bloom.
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

# for downloading hiseq data
rsync --recursive --times --verbose --stats --progress --itemize-changes rsync://SxaQSEQsXbp048L1@pan.pellegrini.mcdb.ucla.edu/SxaQSEQsXbp048L1 L1/
rsync --recursive --times --verbose --stats --progress --itemize-changes rsync://SxaQSEQsXbp048L2@pan.pellegrini.mcdb.ucla.edu/SxaQSEQsXbp048L2 L2/

# 2 rapid runs 
cat L1_read1.fastq   L2_read1.fastq  | pigz  >  05_read1.fastq.gz
cat L1_read2.fastq   L2_read2.fastq  | pigz  >  05_read2.fastq.gz
cat L1_index1.fastq  L2_index1.fastq | pigz  >  05_index1.fastq.gz
cat L1_index2.fastq  L2_index2.fastq | pigz  >  05_index2.fastq.gz

#zcat 02_read1.fastq.gz | awk 'NR%4==1' - | head

# assuming n700 is index1 ... double check this
# N7* is index1
# N5* is index2
# to demultiplex
fastq-multx -B barcode_n700.txt 05_index1.fastq.gz 05_read1.fastq.gz 05_read2.fastq.gz 05_index2.fastq.gz \
    -o n/a read1.%.fq.gz  read2.%.fq.gz  index.%.fq.gz 

fastq-multx -m 1 -B N701_96_barcodes.txt index.N701.fq.gz read1.N701.fq.gz read2.N701.fq.gz \
    -o n/a fastq/%_R1.fq.gz fastq/%_R2.fq.gz &
fastq-multx -m 1 -B N702_96_barcodes.txt index.N702.fq.gz read1.N702.fq.gz read2.N702.fq.gz \
    -o n/a fastq/%_R1.fq.gz fastq/%_R2.fq.gz &
fastq-multx -m 1 -B N703_96_barcodes.txt index.N703.fq.gz read1.N703.fq.gz read2.N703.fq.gz \
    -o n/a fastq/%_R1.fq.gz fastq/%_R2.fq.gz &
fastq-multx -m 1 -B N704_96_barcodes.txt index.N704.fq.gz read1.N704.fq.gz read2.N704.fq.gz \
    -o n/a fastq/%_R1.fq.gz fastq/%_R2.fq.gz &


genome='/data/CRISPR_LOH/Reference/ucsc_sacCer3/sacCer3_CRISPR.fasta'
root_dir='/data/CRISPR_LOH/05/'
fastq_dir="${root_dir}Rfastq/"
bam_dir="${root_dir}bam/"
adapter_file='/data/CRISPR_LOH/Nextera_primers_PE.fa'
snames=(`cat "${root_dir}MS_chrVII_fm.list" `)

cputhreads=11
mkdir $bam_dir

for samp in ${snames[@]}
do
        trimmomatic PE -threads $cputhreads \
        ${fastq_dir}/${samp}_R1.fq.gz   ${fastq_dir}/${samp}_R2.fq.gz \
        ${fastq_dir}/${samp}_t_R1.fastq.gz ${fastq_dir}/${samp}_R1_up.fastq.gz \
        ${fastq_dir}/${samp}_t_R2.fastq.gz ${fastq_dir}/${samp}_R2_up.fastq.gz \
        ILLUMINACLIP:${adapter_file}:2:20:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:0:0 MINLEN:30 TOPHRED33
done

#for samp in ${snames[@]}
#do
#        trimmomatic PE -threads $cputhreads \
#        ${fastq_dir}/${samp}_R1.fq.gz   ${fastq_dir}/${samp}_R2.fq.gz \
#        ${fastq_dir}/${samp}_t_R1.fastq.gz ${fastq_dir}/${samp}_R1_up.fastq.gz \
#        ${fastq_dir}/${samp}_t_R2.fastq.gz ${fastq_dir}/${samp}_R2_up.fastq.gz \
#        ILLUMINACLIP:${adapter_file}:2:20:10 LEADING:0 TRAILING:0 SLIDINGWINDOW:0:0 MINLEN:30 TOPHRED33
#done


for t in  ${snames[@]}
    do 
    echo $t
    # note, this is for PE reads ... adjust for SE reads
    bwa mem -t $cputhreads ${genome} "${fastq_dir}${t}_t_R1.fastq.gz" "${fastq_dir}${t}_t_R2.fastq.gz" | \
    samtools view -hSuq 1 - | \
    samtools sort - ${bam_dir}${t}
done
#-----------------------------------------------------------------------------------------------------------------------------------

# Create a merged bam for all of your samples -------------------------------------------------------------
cd $root_dir
# create a header
sed MS_chrVII_fm.list -e 's/\(.*\)/@RG\tID:\1\tSM:\1/g' | cat ref_contigs_header.sam - > header.sam
# merge bam files, attach @RG flad to each read and attach new header
samtools merge -r -h header.sam MS_chrVII_fmR.bam $bam_dir/*.bam
# index bam file (necessary for variant calling)
samtools index MS_chrVII_fmR.bam

# 04/MS_chrVII_fm.bam comes from a miseq run
picard MergeSamFiles I=/data/CRISPR_LOH/04/MS_chrVII_fm.bam    \
                     I=/data/CRISPR_LOH/05/MS_chrVII_fmR.bam   \
                     O=/data/CRISPR_LOH/05/MS_chrVII_fmM.bam   \
                      CREATE_INDEX=TRUE ASSUME_SORTED=TRUE USE_THREADING=TRUE MAX_RECORDS_IN_RAM=1000000

#gatk -R  /data/CRISPR_LOH/Reference/ucsc_sacCer3/sacCer3_CRISPR.fasta \
#     -T UnifiedGenotyper \
#     -I MS_chrVII_fmM2.bam \
#     -o MS_chrVII_fmM2.vcf \
#     -stand_call_conf 50.0 \
#     -stand_emit_conf 10.0 \
#     -dcov 5000 \
#     -nt $cputhreads \
#     --sample_ploidy 2 \
#     -comp /data/CRISPR_LOH/Reference/RM.var.flt.MoreIndels.vcf 

#vcftools --vcf MS_chrVII_fmM2.vcf --out MS_chrVII_fmM2 --min-alleles 2 --max-alleles 2 --minQ 900 --minGQ 30  --keep-INFO-all --recode

#no merging of bam files
gatk -R  /data/CRISPR_LOH/Reference/ucsc_sacCer3/sacCer3_CRISPR.fasta \
     -T UnifiedGenotyper \
     -I MS_chrVII_fmM.bam \
     -o MS_chrVII_fmM.vcf \
     -stand_call_conf 50.0 \
     -stand_emit_conf 10.0 \
     -dcov 5000 \
     -nt 12 \
     --sample_ploidy 2 \
     --genotyping_mode GENOTYPE_GIVEN_ALLELES \
     -alleles /data/CRISPR_LOH/Reference/RM.var.flt.MoreIndels.vcf \
     -L     /data/CRISPR_LOH/Reference/RM.var.flt.MoreIndels.vcf \
     -comp /data/CRISPR_LOH/Reference/RM.var.flt.MoreIndels.vcf

# Added manual sanger sequencing instead of alignment optimization
# careful variant calling with intervals list
#gatk -R /data/CRISPR_LOH/Reference/ucsc_sacCer3/sacCer3_CRISPR.fasta \
#     -T HaplotypeCaller \
#    -I MS_chrVII_fmM.bam \
#    -o MS_chrVII_fmM_cut_site_denovo_variants.vcf \
#    -L /data/CRISPR_LOH/05/chrVII.cut.interval.list
