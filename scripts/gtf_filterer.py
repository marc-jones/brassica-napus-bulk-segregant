#!/bin/python

import re
import os
import csv

gtf_filepath = ('/home/jonesd/morris_local' +
    '/2016_08_02_final_sequencing_results/out_2016_08_15_19_05' +
    '/merged_asm/merged.gtf')

annotation_filepath = ('/home/jonesd/morris_local' +
    '/2016_08_02_final_sequencing_results/out_2016_08_15_19_05' +
    '/merged_genes_annotation.tsv')

flowering_filepaths = ['/home/jonesd/Dropbox/JIC/2016' +
    '/2016_06_02_differential_expression_analysis/flor_id' +
    '/flowering_time_2016_08_19.csv',
    '/home/jonesd/Dropbox/JIC/2016' +
    '/2016_06_02_differential_expression_analysis/flor_id' +
    '/flower_development_2016_08_19.csv']

def get_annotation_dict(arabidopsis_gene_list=[]):
    # Only gets the top hit for each gene
    annotation = open(annotation_filepath)
    annotation.readline()
    annotation_dict = {}
    seen_genes = []
    for line in annotation:
        line = line.strip().split('\t')
        if not line[0] in seen_genes:
            seen_genes.append(line[0])
            if (line[8][0:9] in arabidopsis_gene_list or
                arabidopsis_gene_list == []):
                if len(line) == 10:
                    symbols = line[9]
                elif len(line) == 9:
                    symbols == ''
                else:
                    sys.exit('Unexpected annotation line length encountered')
                annotation_dict.setdefault(line[0], [line[8][0:9], symbols])
    return(annotation_dict)

def get_flowering_dict():
    output_dict = {}
    for file_path in flowering_filepaths:
        with open(file_path, 'rb') as csv_file:
            csv_reader = csv.reader(csv_file)
            headers = csv_reader.next()
            for row in csv_reader:
                output_dict.setdefault(row[6], row[0:2])
    return(output_dict)

if __name__ == '__main__':
    annotated_gtf_filepath = 'gtf_files/annotated_merged_asm.gtf'
    flowering_gtf_filepath = 'gtf_files/flowering_merged_asm.gtf'

    annotated_gtf_file_hook = open(annotated_gtf_filepath, 'w')
    flowering_gtf_file_hook = open(flowering_gtf_filepath, 'w')

    bra_to_ara_dict = get_annotation_dict()

    flowering_dict = get_flowering_dict()
    bra_to_ara_flowering_dict = get_annotation_dict(flowering_dict.keys())

    with open(gtf_filepath) as gtf:
        for line in gtf:
            row = line.strip().split('\t')
            brassica_gene = re.findall('XLOC_[0-9]+', row[8])[0]

            arabidopsis_gene_agi = ''
            arabidopsis_gene_symbols = ''

            if brassica_gene in bra_to_ara_dict.keys():
                arabidopsis_gene_agi = bra_to_ara_dict[brassica_gene][0]
                arabidopsis_gene_symbols = (
                    bra_to_ara_dict[brassica_gene][1].replace(' ', '_'))

            annotated_gtf_file_hook.write(line.replace('\n',
                ' AGI "' + arabidopsis_gene_agi + '"; symbols "' +
                arabidopsis_gene_symbols + '";\n'))

            if brassica_gene in bra_to_ara_flowering_dict.keys():
                flowering_gtf_file_hook.write(line.replace('\n',
                    ' AGI "' + arabidopsis_gene_agi + '"; symbols "' +
                    arabidopsis_gene_symbols + '";\n'))

    annotated_gtf_file_hook.close()
    flowering_gtf_file_hook.close()
