import sys
import re
from pileup import *

# inspiration http://blog.nextgenetics.net/?e=56

#filepath = 'mpileup_2017_12_20_19_08/output_data/combined.pileup'
filepath = 'mpileup_unique_2018_04_20_15_02/output_data/combined.pileup'

def determine_type_of_change(ref_allele_str, allele_str):
    if ref_allele_str in ['A', 'T', 'C', 'G']:
        if len(allele_str) > 1:
            return('insertion')
        elif allele_str == '*':
            return('deletion')
        elif allele_str in ['A', 'T', 'C', 'G']:
            return('snp')
        else:
            sys.exit('Unexpected change observed: ' + ref_allele_str + ':' +
                allele_str)
    elif ref_allele_str == '*':
        if len(allele_str) > 1:
            return('insertion')
        elif allele_str in ['A', 'T', 'C', 'G']:
            return('insertion')
        else:
            sys.exit('Unexpected change observed: ' + ref_allele_str + ':' +
                allele_str)
    elif ref_allele_str == 'N':
        if allele_str in ['A', 'T', 'C', 'G']:
            return('replacement')
        else:
            sys.exit('Unexpected change observed: ' + ref_allele_str + ':' +
                allele_str)
    elif not len(ref_allele_str) == len(allele_str):
        if len(ref_allele_str) < len(allele_str):
            return('insertion')
        elif len(allele_str) < len(ref_allele_str):
            return('deletion')
        else:
            sys.exit('Unexpected change observed: ' + ref_allele_str + ':' +
                allele_str)
    elif len(ref_allele_str) == len(allele_str) and len(allele_str) > 1:
        number_of_diffs = sum([int(ref_allele_str[idx] == allele_str[idx])
            for idx in range(len(allele_str))])
        if number_of_diffs == 1:
            return('snp')
        else:
            sys.exit('Unexpected change observed: ' + ref_allele_str + ':' +
                allele_str)
    else:
        sys.exit('Unexpected reference allele observed (' + ref_allele_str +
            ') with test allele (' + allele_str + ')')

if __name__ == '__main__':
    min_depth = 20
    min_snp_percentage = 0.95
    min_snp_idx = 0.3
    ref_par_num = int(sys.argv[3])
    oth_par_num = int(sys.argv[4])
    ref_ref_output_file = sys.argv[1]
    ref_oth_output_file = sys.argv[2]
    ref_ref_output_file_hook = open(ref_ref_output_file, 'w')
    ref_ref_output_file_hook.write('\t'.join(['chromosome',
        'position', 'published_reference', 'reference_parent_allele',
        'reference_parent_depth', 'type']) + '\n')
    ref_oth_output_file_hook = open(ref_oth_output_file, 'w')
    ref_oth_output_file_hook.write('\t'.join(['chromosome',
        'position', 'reference_parent_allele',
        'other_parent_allele', 'reference_parent_depth',
        'other_parent_depth', 'type']) + '\n')
    with open(filepath) as pileup:
        for line in pileup:
            line_list = line.strip().split()
            (chromosome, position, reference) = line_list[0:3]
            (ref_par_depth, ref_par_reads, ref_par_quals) = (
                extract_relevant_columns(line_list, ref_par_num))
            (oth_par_depth, oth_par_reads, oth_par_quals) = (
                extract_relevant_columns(line_list, oth_par_num))

            # check what the reference depth is
            if int(ref_par_depth) >= min_depth:

                ref_par_reads_list = parse_read_column(ref_par_reads,
                    reference, int(ref_par_depth))

                # check if the reference parent allele is above threshold
                if (return_major_allele_freq(ref_par_reads_list) >=
                    min_snp_percentage):

                    ref_allele = return_major_allele(ref_par_reads_list)

                    # if the reference parent allele is different to the
                    # reference
                    if not ref_allele == reference:

                        ref_ref_output_file_hook.write('\t'.join([str(val)
                            for val in [chromosome, position, reference,
                            ref_allele, ref_par_depth,
                            determine_type_of_change(reference, ref_allele)]])
                            + '\n')


            # check for min read depth
            if (int(ref_par_depth) >= min_depth and
                int(oth_par_depth) >= min_depth):

                ref_par_reads_list = parse_read_column(ref_par_reads,
                    reference, int(ref_par_depth))
                oth_par_reads_list = parse_read_column(oth_par_reads,
                    reference, int(oth_par_depth))

                # check that parental alleles are above threshold
                if (return_major_allele_freq(ref_par_reads_list) >=
                    min_snp_percentage and return_major_allele_freq(
                    oth_par_reads_list) >= min_snp_percentage):

                    ref_allele = return_major_allele(ref_par_reads_list)
                    oth_allele = return_major_allele(oth_par_reads_list)

                    # check if ref and oth are different
                    if not ref_allele == oth_allele:

                        ref_oth_output_file_hook.write('\t'.join([str(val)
                            for val in [chromosome,
                            position, ref_allele, oth_allele,
                            ref_par_depth, oth_par_depth,
                            determine_type_of_change(ref_allele, oth_allele)]])
                            + '\n')

    ref_ref_output_file_hook.close()
    ref_oth_output_file_hook.close()
