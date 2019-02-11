rm(list=ls(all=TRUE))
graphics.off()

library(ggplot2)
library(grid)
library(gridExtra)

# ### Non Unique
# input_folder <- 'delta_2018_01_03_16_16/output_data/'
#
# depth_folder <- 'depth_2018_01_09_12_06/output_data/'
#
# dar_cab_diffs <- read.table(paste0('parents_2018_01_08_12_25/output_data',
#     '/darmor_cabriolet.csv'), sep='\t', stringsAsFactors=FALSE, header=TRUE)
#
# dar_dar_diffs <- read.table(paste0('parents_2018_01_08_12_25/output_data',
#     '/published_darmor_darmor.csv'), sep='\t', stringsAsFactors=FALSE,
#     header=TRUE)

### Unique
input_folder <- 'delta_unique_2018_04_23_13_47/output_data/'

depth_folder <- 'depth_unique_2018_04_24_11_42/output_data/'

dar_cab_diffs <- read.table(paste0('parents_unique_2018_04_23_18_31/output_data',
    '/darmor_cabriolet.csv'), sep='\t', stringsAsFactors=FALSE, header=TRUE)

dar_dar_diffs <- read.table(paste0('parents_unique_2018_04_23_18_31/output_data',
    '/published_darmor_darmor.csv'), sep='\t', stringsAsFactors=FALSE,
    header=TRUE)



samples <- list.files(input_folder)

results_table_list <- sapply(samples, function(file_name) {
    read.table(paste0(input_folder, file_name), sep='\t',
    stringsAsFactors=FALSE, header=TRUE)}, simplify=FALSE, USE.NAMES=TRUE)

depth_table_list <- sapply(samples, function(file_name) {
    read.table(paste0(depth_folder, file_name), sep='\t',
    stringsAsFactors=FALSE, header=TRUE)}, simplify=FALSE, USE.NAMES=TRUE)

chromosomes <- unique(results_table_list[[1]]$chromosome)


threshold_table <- c()
for (file_name in samples) {
    results <- results_table_list[[file_name]]
    threshold <- quantile(abs(results$delta_snp_index), 0.99, names=FALSE)
    threshold_table <- rbind(threshold_table,
        data.frame(file_name=rep(file_name, 3),
        value=c(mean(results$delta_snp_index), -threshold, threshold),
        type=c('mean', 'lo', 'hi')))
    filtered_results <- results[results$delta_snp_index <=
            threshold_table$value[threshold_table$file_name==file_name &
            threshold_table$type=='lo'] |
            threshold_table$value[threshold_table$file_name==file_name &
            threshold_table$type=='hi'] <= results$delta_snp_index,]
    gct_output <- data.frame(Name=paste(filtered_results$chromosome,
        filtered_results$position, filtered_results$reference_parent_allele,
        filtered_results$other_parent_allele, sep='_'),
        Description=paste0('|@', filtered_results$chromosome, ':',
        filtered_results$position, '-', filtered_results$position + 1, '|'))
    gct_output[file_name] <- filtered_results$delta_snp_index

    filepath <- paste0('snps_of_interest_unique/', gsub('.csv', '.gct', file_name))

#     file_connection <- file(filepath)
#     writeLines(c('#1.2', paste(nrow(gct_output), 1, collapse='\t')),
#         file_connection)
#     close(file_connection)
#     write.table(gct_output, file=filepath, quote=FALSE,
#         sep='\t', row.names=FALSE, append=TRUE)
}

snp_windows <- c(100, 1000)

for (chromosome in chromosomes) {
    file_path <- paste0('output_plots_unique/', chromosome, '.pdf')
    if (!file.exists(file_path)) {
    pdf(file_path, height=11, width=8)

    print(chromosome)
    combined_table <- c()

    combined_depth_table <- c()

    for (file_name in samples) {
        results <- results_table_list[[file_name]]
        results$threshold <- c(results$delta_snp_index <=
            threshold_table$value[threshold_table$file_name==file_name &
            threshold_table$type=='lo'] |
            threshold_table$value[threshold_table$file_name==file_name &
            threshold_table$type=='hi'] <= results$delta_snp_index)
        combined_table <- rbind(combined_table, cbind(results[results$chromosome==chromosome,],
            file_name=rep(file_name, nrow(results[results$chromosome==chromosome,]))))

        depth_table <- depth_table_list[[file_name]]

        valid_positions <- as.numeric(unlist(apply(
            depth_table[depth_table$chromosome==chromosome,], 1,
            function(vec) as.numeric(vec['start']):as.numeric(vec['end']))))

        combined_depth_table <- rbind(combined_depth_table,
            data.frame(chromosome=rep(chromosome, length(valid_positions)),
            position=valid_positions,
            file_name=rep(file_name, length(valid_positions))))
    }

    combined_sliding_window_table <- c()

    for (snp_window in snp_windows) {
        for (file_name in samples) {
            combined_table_subset <- combined_table[combined_table$file_name==file_name,]
            if (nrow(combined_table_subset) > snp_window) {
                combined_sliding_window_table <- rbind(combined_sliding_window_table,
                    do.call('rbind', lapply(1:(nrow(combined_table_subset) - snp_window),
                    function(row_idx) {return(data.frame(chromosome=chromosome,
                        position=mean(combined_table_subset[row_idx:(row_idx + snp_window),
                        ]$position),
                        average_delta_snp_index=mean(combined_table_subset[
                        row_idx:(row_idx + snp_window),]$delta_snp_index),
                        file_name=file_name,
                        window_size=snp_window))})))
            }
        }
    }

    combined_sliding_window_table$window_size <- as.factor(
        combined_sliding_window_table$window_size)
    print(paste(chromosome, 'first plot start'))
    p <- ggplot()
#     print(combined_table)
    p <- p + geom_point(aes(x=position, y=delta_snp_index, colour=threshold),
        size=0.5, data=combined_table)
#     print(combined_sliding_window_table)
    if (length(combined_sliding_window_table$window_size) > 0) {
        p <- p + geom_line(aes(x=position, y=average_delta_snp_index,
            colour=window_size), data=combined_sliding_window_table)
    }
    p <- p + facet_grid(file_name~chromosome)
    p <- p + scale_y_continuous(limits=c(-1, 1))
#     print(threshold_table)
    p <- p + geom_hline(aes(yintercept=value), data=threshold_table)

    bin_width <- diff(range(combined_table$position)) / 100

    print(paste(chromosome, 'second plot start'))
    dar_cab_subset <- dar_cab_diffs[dar_cab_diffs$chromosome==chromosome,]
    pd <- ggplot(dar_cab_subset)
    pd <- pd + geom_histogram(aes(x=position, fill=type, colour=type),
        binwidth=bin_width, alpha=0.1)
    pd <- pd + facet_grid(type~chromosome, scale='free_y')
    pd <- pd + ggtitle(paste0(
        'JIC_Darmor - Cabriolet Differences : Bin width - ',
        round(bin_width, 1) / 1000, ' kb'))

    print(paste(chromosome, 'third plot start'))
    dar_dar_subset <- dar_dar_diffs[dar_dar_diffs$chromosome==chromosome,]
    pdd <- ggplot(dar_dar_subset)
    pdd <- pdd + geom_histogram(aes(x=position, fill=type, colour=type),
        binwidth=bin_width, alpha=0.1)
    pdd <- pdd + facet_grid(type~chromosome, scale='free_y')
    pdd <- pdd + ggtitle(paste0(
        'Darmor - JIC_Darmor Differences : Bin width - ',
        round(bin_width, 1) / 1000, ' kb'))

    print(paste(chromosome, 'fourth plot start'))
    d <- ggplot(combined_depth_table)
    d <- d + geom_histogram(aes_string(x='position',
        y=paste0('(..count../', bin_width, ')*100')), binwidth=bin_width)
    d <- d + ylab('Percentage of bases in bin with adequate depth')
    d <- d + facet_grid(file_name~chromosome)

    print(paste(chromosome, 'combining plots start'))
    g1 <- ggplotGrob(p)
    g2 <- ggplotGrob(pd)
    g3 <- ggplotGrob(pdd)
    g4 <- ggplotGrob(d)
    colnames(g1) <- paste0(seq_len(ncol(g1)))
    colnames(g2) <- paste0(seq_len(ncol(g2)))
    colnames(g3) <- paste0(seq_len(ncol(g3)))
    colnames(g4) <- paste0(seq_len(ncol(g4)))

    print(paste(chromosome, 'drawing start'))
    grid.draw(p)
    grid.newpage()
    grid.draw(combine(g1, g2, g3, g4, along=2))

    dev.off()
    }
}

