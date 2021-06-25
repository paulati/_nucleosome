# To install the required Bioconductor packages, run:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("genomation", "Gviz", "GenomicFeatures", "Rsamtools", "rtracklayer"))
# 
# # To install other packages used in this notebook, run:
# install.packages("ggplot2")
# install.packages("pryr")
# install.packages("RMariaDB")
# install.packages("caTools")
# install.packages("colorRamps")

library(genomation)
library(Gviz)
library(GenomicFeatures)
library(Rsamtools)
library(rtracklayer)
library(ggplot2)
library(pryr)
library(caTools)
library(colorRamps)

plot_reads_per_bp <- function(bam_file_path, plot_path) {
  
  bam_file <- Rsamtools::BamFile(bam_file_path)
  
  #scanBamWhat()
  
  # Import selected fields from the bam file
  fields_to_load <- c("qname", "rname", "strand", "pos", "isize")
  param <- ScanBamParam(what = fields_to_load)
  aln_list <- scanBam(bam_file, param = param)
  
  aln <- aln_list[[1]]
  
  # Get the total number of reads
  #length(aln$qname)
  # 93174202 reads que corresponden a 93174202 / 2 paired
  
  # Construct GRanges from the paired-end reads
  read_on_Watson_strand <- (aln$strand == '+') & (aln$isize > 0)
  reads <- GRanges(seqnames = Rle(aln$rname[read_on_Watson_strand]),
                   ranges = IRanges(start = aln$pos[read_on_Watson_strand],
                                    width = aln$isize[read_on_Watson_strand]),
                   seqinfo = seqinfo(bam_file))
  
  # Check the size of memory that are occupied by aln and reads objects
  #pryr::object_size(aln)
  #pryr::object_size(reads)
  
  rm(aln)
  
  nums <- c(1: 41)
  valid_seq_levels <- c(sapply(nums, function(x) paste0('TcChr', x, '-S')),
                        sapply(nums, function(x) paste0('TcChr', x, '-P')))
                        
  aln_seq_names <- levels(seqnames(reads))
  
  to_remove_mask <- sapply(aln_seq_names, function(x) ! is.element(x, valid_seq_levels)) 
                             
  to_remove <- levels(seqnames(reads))[to_remove_mask]
  
  #filter
  reads <- dropSeqlevels(reads, to_remove, pruning.mode="coarse")
  #seqlevels(reads)
  
  
  # Check the number of paired-end reads that were mapped to each chromosome
  no_of_reads <- table(seqnames(reads))
  #no_of_reads
  
  # Get the chromosome lengths
  chr_lengths <- seqlengths(seqinfo(bam_file))[valid_seq_levels]
  #chr_lengths
  
  # Compute the sequencing depth (number of reads / chromosome)
  no_of_reads_per_bp <- no_of_reads / chr_lengths
  #no_of_reads_per_bp
  
  png(plot_path, width = 350, height = 350)
  plot(no_of_reads_per_bp)  
  dev.off()
  
  return(no_of_reads_per_bp)
}

# CL1



out_base_path <- '/home/paula/2020/mnase_octubre2020/plots/CL1/filtered'

# clean
bam_file_path <- '/u03/paula/2020/amazon_jocampo/s3/MNase/alignments/clean_reads/bowtie2/TriTrypDB-46_TcruziCLBrenerEsmeraldo-Non-Esmeraldo-like_Genome/paired_end/CLJ_1_70U_comb.bam'
plot_name <- 'bowtie2_clean_es_noEs_reads_per_bp'
plot_path <- file.path(out_base_path, paste0(plot_name, '.png'))
no_of_reads_per_bp <- plot_reads_per_bp(bam_file_path, plot_path)
no_of_reads_per_bp_df <- data.frame(no_of_reads_per_bp)
colnames(no_of_reads_per_bp_df) <- c('chr', plot_name)

bam_file_path <- '/u03/paula/2020/amazon_jocampo/s3/MNase/alignments/clean_reads/bowtie2/TriTrypDB-46_TcruziCLBrenerAll_Genome/paired_end/CLJ_1_70U_comb.bam'
plot_name <- 'bowtie2_clean_all_reads_per_bp'
plot_path <- file.path(out_base_path, paste0(plot_name, '.png'))
no_of_reads_per_bp <- plot_reads_per_bp(bam_file_path, plot_path)
tmp <- data.frame(no_of_reads_per_bp)
colnames(tmp) <- c('chr', plot_name)
no_of_reads_per_bp_df <- merge(no_of_reads_per_bp_df, tmp, by='chr')


bam_file_path <- '/u03/paula/2020/amazon_jocampo/s3/MNase/alignments/clean_reads/hisat2/TriTrypDB-46_TcruziCLBrenerEsmeraldo-Non-Esmeraldo-like_Genome/paired_end/CLJ_1_70U_comb.bam'
plot_name <- 'hisat2_clean_es_noEs_reads_per_bp'
plot_path <- file.path(out_base_path, paste0(plot_name, '.png'))
no_of_reads_per_bp <- plot_reads_per_bp(bam_file_path, plot_path)
tmp <- data.frame(no_of_reads_per_bp)
colnames(tmp) <- c('chr', plot_name)
no_of_reads_per_bp_df <- merge(no_of_reads_per_bp_df, tmp, by='chr')

bam_file_path <- '/u03/paula/2020/amazon_jocampo/s3/MNase/alignments/clean_reads/hisat2/TriTrypDB-46_TcruziCLBrenerAll_Genome/paired_end/CLJ_1_70U_comb.bam'
plot_name <- 'hisat2_clean_all_reads_per_bp'
plot_path <- file.path(out_base_path, paste0(plot_name, '.png'))
no_of_reads_per_bp <- plot_reads_per_bp(bam_file_path, plot_path)
tmp <- data.frame(no_of_reads_per_bp)
colnames(tmp) <- c('chr', plot_name)
no_of_reads_per_bp_df <- merge(no_of_reads_per_bp_df, tmp, by='chr')

# raw
bam_file_path <- '/u03/paula/2020/amazon_jocampo/s3/MNase/alignments/raw_reads/bowtie2/TriTrypDB-46_TcruziCLBrenerEsmeraldo-Non-Esmeraldo-like_Genome/paired_end/CLJ_1_70U_comb.bam'
plot_name <- 'bowtie2_raw_es_noEs_reads_per_bp'
plot_path <- file.path(out_base_path, paste0(plot_name, '.png'))
no_of_reads_per_bp <-plot_reads_per_bp(bam_file_path, plot_path)
tmp <- data.frame(no_of_reads_per_bp)
colnames(tmp) <- c('chr', plot_name)
no_of_reads_per_bp_df <- merge(no_of_reads_per_bp_df, tmp, by='chr')

bam_file_path <- '/u03/paula/2020/amazon_jocampo/s3/MNase/alignments/raw_reads/bowtie2/TriTrypDB-46_TcruziCLBrenerAll_Genome/paired_end/CLJ_1_70U_comb.bam'
plot_name <- 'bowtie2_raw_all_reads_per_bp'
plot_path <- file.path(out_base_path, paste0(plot_name, '.png'))
no_of_reads_per_bp <- plot_reads_per_bp(bam_file_path, plot_path)
tmp <- data.frame(no_of_reads_per_bp)
colnames(tmp) <- c('chr', plot_name)
no_of_reads_per_bp_df <- merge(no_of_reads_per_bp_df, tmp, by='chr')

bam_file_path <- '/u03/paula/2020/amazon_jocampo/s3/MNase/alignments/raw_reads/hisat2/TriTrypDB-46_TcruziCLBrenerEsmeraldo-Non-Esmeraldo-like_Genome/paired_end/CLJ_1_70U_comb.bam'
plot_name <- 'hisat2_raw_es_noEs_reads_per_bp'
plot_path <- file.path(out_base_path, paste0(plot_name, '.png'))
no_of_reads_per_bp <- plot_reads_per_bp(bam_file_path, plot_path)
tmp <- data.frame(no_of_reads_per_bp)
colnames(tmp) <- c('chr', plot_name)
no_of_reads_per_bp_df <- merge(no_of_reads_per_bp_df, tmp, by='chr')

bam_file_path <- '/u03/paula/2020/amazon_jocampo/s3/MNase/alignments/raw_reads/hisat2/TriTrypDB-46_TcruziCLBrenerAll_Genome/paired_end/CLJ_1_70U_comb.bam'
plot_name <- 'hisat2_raw_all_reads_per_bp'
plot_path <- file.path(out_base_path, 'hisat2_raw_all_reads_per_bp.png')
no_of_reads_per_bp <- plot_reads_per_bp(bam_file_path, plot_path)
tmp <- data.frame(no_of_reads_per_bp)
colnames(tmp) <- c('chr', plot_name)
no_of_reads_per_bp_df <- merge(no_of_reads_per_bp_df, tmp, by='chr')

no_of_reads_per_bp_file_path <- file.path(out_base_path, 'no_of_reads_per_bp.csv')
write.table(no_of_reads_per_bp_df, no_of_reads_per_bp_file_path, quote = FALSE,  
            row.names = FALSE, col.names = TRUE)



# ---------------------------------------------------

# pairs(~no_of_reads_per_bp_df$bowtie2_clean_es_noEs_reads_per_bp+
#         no_of_reads_per_bp_df$bowtie2_clean_all_reads_per_bp+
#         no_of_reads_per_bp_df$hisat2_clean_es_noEs_reads_per_bp+
#         no_of_reads_per_bp_df$hisat2_clean_all_reads_per_bp+
#         no_of_reads_per_bp_df$bowtie2_raw_es_noEs_reads_per_bp+
#         no_of_reads_per_bp_df$bowtie2_raw_all_reads_per_bp+
#         no_of_reads_per_bp_df$hisat2_raw_es_noEs_reads_per_bp+
#         no_of_reads_per_bp_df$hisat2_raw_all_reads_per_bp,
#       data = no_of_reads_per_bp_df,
#       main = "bowtie2 - hisat2")

pairs(~bowtie2_clean_es_noEs_reads_per_bp+
        bowtie2_clean_all_reads_per_bp+
        hisat2_clean_es_noEs_reads_per_bp+
        hisat2_clean_all_reads_per_bp+
        bowtie2_raw_es_noEs_reads_per_bp+
        bowtie2_raw_all_reads_per_bp+
        hisat2_raw_es_noEs_reads_per_bp+
        hisat2_raw_all_reads_per_bp,
      data = no_of_reads_per_bp_df,
      main = "bowtie2 - hisat2")


# install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)

my_data <- no_of_reads_per_bp_df[,-1]
chart.Correlation(my_data, histogram=TRUE, pch=19)






