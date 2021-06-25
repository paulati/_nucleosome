#!/bin/bash

#genome='TriTrypDB-46_TcruziCLBrenerEsmeraldo-like_Genome' 
#genome='TriTrypDB-46_TcruziCLBrenerNon-Esmeraldo-like_Genome' 
#genome='TriTrypDB-46_TcruziCLBrener_Genome'
genome='TriTrypDB-46_TcruziCLBrenerAll_Genome'
tool='hisat2' #options: hisat2 bowtie2
proc_num=8
index_base_path='/u03/paula/2020/amazon_jocampo/s3/genome/index'
reads_base_path='/u03/paula/2020/amazon_jocampo/s3/MNase/reads/clean'
alignments_base_path='/u03/paula/2020/amazon_jocampo/s3/MNase/alignments'


#TODO: algo que checkee estructura de directorios y cree carpetas

index_param=$index_base_path'/'$tool'/'$genome'/'$genome

reads_base_name='CLJ_1_70U_comb'

paired1_param=$reads_base_path'/'$reads_base_name'_R1.fastq.gz'
paired2_param=$reads_base_path'/'$reads_base_name'_R2.fastq.gz'

paired_suffix=''
paired_discordant_suffix='_discordant'
paired_discordant_mixed_suffix='_discordant_mixed'

alignment_paired_param=$alignments_base_path'/'$tool'/'$genome'/paired_end/'$reads_base_name'.sam'
alignment_paired_discordant_param=$alignments_base_path'/'$tool'/'$genome'/paired_end/'$reads_base_name$paired_discordant_suffix'.sam'
alignment_paired_discordant_mixed_param=$alignments_base_path'/'$tool'/'$genome'/paired_end/'$reads_base_name$paired_discordant_mixed_suffix'.sam'

alignment_single1_param=$alignments_base_path'/'$tool'/'$genome'/single_read/'$reads_base_name'_R1.sam'
alignment_single2_param=$alignments_base_path'/'$tool'/'$genome'/single_read/'$reads_base_name'_R2.sam'

summary_param=$alignments_base_path'/'$tool'/'$genome'/paired_end/'$reads_base_name'_summary.txt'
summary_discordant_param=$alignments_base_path'/'$tool'/'$genome'/paired_end/'$reads_base_name$paired_discordant_suffix'_summary.txt'
summary_discordant_mixed_param=$alignments_base_path'/'$tool'/'$genome'/paired_end/'$reads_base_name$paired_discordant_mixed_suffix'_summary.txt'

summary_single1_param=$alignments_base_path'/'$tool'/'$genome'/single_read/'$reads_base_name'_R1_summary.txt'
summary_single2_param=$alignments_base_path'/'$tool'/'$genome'/single_read/'$reads_base_name'_R2_summary.txt'

log_param=$alignments_base_path'/'$tool'/'$genome'/paired_end/'$reads_base_name'.log'
log_discordant_param=$alignments_base_path'/'$tool'/'$genome'/paired_end/'$reads_base_name$paired_discordant_suffix'.log'
log_discordant_mixed_param=$alignments_base_path'/'$tool'/'$genome'/paired_end/'$reads_base_name$paired_discordant_mixed_suffix'.log'

log_single1_param=$alignments_base_path'/'$tool'/'$genome'/single_read/'$reads_base_name'_R1.log'
log_single2_param=$alignments_base_path'/'$tool'/'$genome'/single_read/'$reads_base_name'_R2.log'

command_paired=$tool' --max-seeds 20 --no-spliced-alignment --score-min L,-4,-0.5 --phred33 -p '$proc_num' -X 1000 --no-mixed --no-discordant --no-unal -q ''-x '$index_param' -1 '$paired1_param' -2 '$paired2_param' -S '$alignment_paired_param' --summary-file '$summary_param' > '$log_param

command_paired_discordant=$tool' --max-seeds 20 --no-spliced-alignment --score-min L,-0.6,-0.6 -p '$proc_num' -X 1000 --no-mixed --no-unal -q ''-x '$index_param' -1 '$paired1_param' -2 '$paired2_param' -S '$alignment_paired_discordant_param' --summary-file '$summary_discordant_param' > '$log_discordant_param

command_paired_discordant_mixed=$tool' --max-seeds 20 --no-spliced-alignment --score-min L,-0.6,-0.6 -p '$proc_num' -X 1000 --no-unal -q ''-x '$index_param' -1 '$paired1_param' -2 '$paired2_param' -S '$alignment_paired_discordant_mixed_param' --summary-file '$summary_discordant_mixed_param' > '$log_discordant_mixed_param

command_single1=$tool' --max-seeds 20 --no-spliced-alignment --score-min L,-0.6,-0.6 -p '$proc_num' --no-unal -q ''-x '$index_param' -U '$paired1_param' -S '$alignment_single1_param' --summary-file '$summary_single1_param' > '$log_single1_param

command_single2=$tool' --max-seeds 20 --no-spliced-alignment --score-min L,-0.6,-0.6 -p '$proc_num' --no-unal -q ''-x '$index_param' -U '$paired2_param' -S '$alignment_single2_param' --summary-file '$summary_single2_param' > '$log_single2_param

#command='ls -l'

#echo $command_paired

eval $command_paired
#eval $command_paired_discordant
#eval $command_paired_discordant_mixed
#eval $command_single1
#eval $command_single2
#echo $command_single1





