base_path <- '/u03/paula/2020/amazon_jocampo/s3/MNase/UTRme_2020-02-27_14-31-50'

file_name <- 'Corrida_UTRme-Epi_Li_46-5-best-score.gff'
file_path <- file.path(base_path, file_name)

out_file_name <- 'Corrida_UTRme-Epi_Li_46-5-best-score.bed'
out_file_path <- file.path(base_path, out_file_name)

data <- read.csv(file_path, sep='\t', comment.char = '#',
                 header = FALSE, stringsAsFactors = FALSE)

five_prime_UTR_mask <- data$V3 == 'five_prime_UTR'

columns <- c('V1', 'V4', 'V5', 'V7')

data_five_prime_UTR <- data[five_prime_UTR_mask, columns]
# hay que dar vuelta los del strand - ?

write.table(data_five_prime_UTR, out_file_path, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')




