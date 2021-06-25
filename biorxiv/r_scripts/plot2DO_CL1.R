get_command_line_args_base <- function(data_file_path) {
  
  option_file <- paste0("--file=", data_file_path)
  option_type <- "--type=occ"
  option_minLength <- "--minLength=90"
  option_maxLength <- "--maxLength=180"
  option_genome <- "--genome=TcruziCLBrenerEsmeraldo"
  option_align <- "--align=fivePrime"
  command_line_args_base = c(option_file, option_type, option_minLength, option_maxLength, option_genome, option_align)
  
  return(command_line_args_base)
}

base_path <- "/home/paula/2020/mnase_octubre2020/plot2DO"

setwd(base_path)

suppressPackageStartupMessages({
  library(yaml)
})

# Load default paths from config.yaml
config <- yaml.load_file("config/config.yaml")

sourceBasePath <- config$application$paths$source
if(is.null(sourceBasePath)) {
  sourceBasePath <- file.path(getwd(), "source")
} 

readsBasePath <- config$application$paths$reads
if(is.null(readsBasePath)) {
  readsBasePath <- getwd()
}

outputBasePath <- config$application$paths$output
if(is.null(outputBasePath)) {
  outputBasePath <- file.path(getwd(), "output/50bp")
}

annotationsBasePath <- config$application$paths$annotations
if(is.null(annotationsBasePath)) {
  annotationsBasePath <- getwd()
}

configBasePath <- config$application$paths$config
if(is.null(configBasePath)) {
  configBasePath <- file.path(getwd(), "config")
}

main <- file.path(sourceBasePath, "plot2DO_main.R")  

source(main)

data_base_path <- '/u03/paula/2020/amazon_jocampo/s3/MNase/alignments'


#### bowtie all clean
#data_file_path <- file.path(data_base_path, 'clean_reads/bowtie2/TriTrypDB-46_TcruziCLBrenerAll_Genome/paired_end/CLJ_1_70U_comb.bam')
#outputBasePath <- file.path(outputBasePath, 'clean_reads-bowtie2-TriTrypDB-46_TcruziCLBrenerAll_Genome-paired_end')

data_file_path <- file.path(data_base_path, '50bp/clean_reads/bowtie2/TriTrypDB-46_TcruziCLBrenerAll_Genome/paired_end/CLJ_1_70U_comb.bam')
outputBasePath <- file.path(outputBasePath, '50bp/clean_reads-bowtie2-TriTrypDB-46_TcruziCLBrenerAll_Genome-paired_end')


###########################

#### bowtie all raw
#data_file_path <- file.path(data_base_path, 'raw_reads/bowtie2/TriTrypDB-46_TcruziCLBrenerAll_Genome/paired_end/CLJ_1_70U_comb.bam')
#outputBasePath <- file.path(outputBasePath, 'raw_reads-bowtie2-TriTrypDB-46_TcruziCLBrenerAll_Genome-paired_end')

###########################

#### bowtie es_noEs clean
#data_file_path <- file.path(data_base_path, 'clean_reads/bowtie2/TriTrypDB-46_TcruziCLBrenerEsmeraldo-Non-Esmeraldo-like_Genome/paired_end/CLJ_1_70U_comb.bam')
#outputBasePath <- file.path(outputBasePath, 'clean_reads-bowtie2-TriTrypDB-46_TcruziCLBrenerEsmeraldo-Non-Esmeraldo-like_Genome-paired_end')

#data_file_path <- file.path(data_base_path, '50bp/clean_reads/bowtie2/TriTrypDB-46_TcruziCLBrenerEsmeraldo-Non-Esmeraldo-like_Genome/paired_end/CLJ_1_70U_comb.bam')
#outputBasePath <- file.path(outputBasePath, '50bp/clean_reads-bowtie2-TriTrypDB-46_TcruziCLBrenerEsmeraldo-Non-Esmeraldo-like_Genome-paired_end')

###########################

#### bowtie es_noEs raw
#data_file_path <- file.path(data_base_path, 'raw_reads/bowtie2/TriTrypDB-46_TcruziCLBrenerEsmeraldo-Non-Esmeraldo-like_Genome/paired_end/CLJ_1_70U_comb.bam')
#outputBasePath <- file.path(outputBasePath, 'raw_reads-bowtie2-TriTrypDB-46_TcruziCLBrenerEsmeraldo-Non-Esmeraldo-like_Genome-paired_end')
###########################

#### hisat all clean
#data_file_path <- file.path(data_base_path, 'clean_reads/hisat2/TriTrypDB-46_TcruziCLBrenerAll_Genome/paired_end/CLJ_1_70U_comb.bam')
#outputBasePath <- file.path(outputBasePath, 'clean_reads-hisat2-TriTrypDB-46_TcruziCLBrenerAll_Genome-paired_end')
###########################

#### hisat all raw
#data_file_path <- file.path(data_base_path, 'raw_reads/hisat2/TriTrypDB-46_TcruziCLBrenerAll_Genome/paired_end/CLJ_1_70U_comb.bam')
#outputBasePath <- file.path(outputBasePath, 'raw_reads-hisat2-TriTrypDB-46_TcruziCLBrenerAll_Genome-paired_end')
###########################

#### hisat es noEs clean
#data_file_path <- file.path(data_base_path, 'clean_reads/hisat2/TriTrypDB-46_TcruziCLBrenerEsmeraldo-Non-Esmeraldo-like_Genome/paired_end/CLJ_1_70U_comb.bam')
#outputBasePath <- file.path(outputBasePath, 'clean_reads-hisat2-TriTrypDB-46_TcruziCLBrenerEsmeraldo-Non-Esmeraldo-like_Genome')
###########################

#### hisat es noEs raw
#data_file_path <- file.path(data_base_path, 'raw_reads/hisat2/TriTrypDB-46_TcruziCLBrenerEsmeraldo-Non-Esmeraldo-like_Genome/paired_end/CLJ_1_70U_comb.bam')
#outputBasePath <- file.path(outputBasePath, 'raw_reads-hisat2-TriTrypDB-46_TcruziCLBrenerEsmeraldo-Non-Esmeraldo-like_Genome')
###########################



command_line_args_base <- get_command_line_args_base(data_file_path)
sites_file_name <- 'annotations/Corrida_UTRme-Epi_Li_46-5-best-score.bed'
sites_arg <-  paste0("--sites=", sites_file_name, "")

command_line_args <- c(command_line_args_base, sites_arg)

Main(command_line_args)



