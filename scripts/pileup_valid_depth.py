import sys
import re

# inspiration http://blog.nextgenetics.net/?e=56

filepath = 'mpileup_2017_12_20_19_08/output_data/combined.pileup'

def parse_read_column(read_str, reference, depth):
    read_str = read_str.upper().replace(',', '.')
    output_list = []
    current_out_base = len(output_list)
    idx = 0
    while idx < len(read_str):
        read = read_str[idx]
        if read == '^': # if it's the start of a read, skip it
            idx += 2
        elif read == '$': # if it's the end of a read, skip it
            idx += 1
        elif read == '.':
            output_list.append(reference)
            current_out_base = len(output_list)
            idx += 1
        elif read == '*':
            output_list.append('*')
            current_out_base = len(output_list)
            idx += 1
        elif read == '-': # skip deletions, they're covered by asterisk
            del_len = int(re.match('[0-9]+', read_str[(idx + 1):]).group(0))
            idx += 1 + len(str(del_len)) + del_len
        elif read == '+': # add insertions to the last base
            ins_len = int(re.match('[0-9]+', read_str[(idx + 1):]).group(0))
            insert = read_str[(idx + 2):(idx + 2 + ins_len)]
            output_list[current_out_base - 1] += insert
            idx += 1 + len(str(ins_len)) + ins_len
        elif read in ['A', 'T', 'C', 'G']:
            output_list.append(read)
            current_out_base = len(output_list)
            idx += 1
        else:
            sys.exit('Unexpected character in read: ' + read)
    if not (len(output_list) == int(depth) or int(depth) == 0):
        sys.exit('Depth did not match: ' + chromosome + ' ' + depth)
    return(output_list)

def extract_relevant_columns(line_list, num):
    return(line_list[(3 * num + 3):(3 * num + 6)])

def return_major_allele_freq(reads_list):
    max_count = 0.0
    for val in list(set(reads_list)):
        max_count = max(max_count, float(reads_list.count(val)))
    return(max_count / len(reads_list))

def return_major_allele(reads_list):
    # caution, only returns the first allele with that frequency
    max_count = 0.0
    count_dict = {}
    for val in list(set(reads_list)):
        max_count = max(max_count, float(reads_list.count(val)))
        count_dict[val] = float(reads_list.count(val))
    return(count_dict.keys()[count_dict.values().index(max_count)])

if __name__ == '__main__':
    min_depth = 20
    min_snp_percentage = 0.95
    min_snp_idx = 0.3
    ref_par_num = int(sys.argv[2])
    oth_par_num = int(sys.argv[3])
    output_file = sys.argv[1]
    bulk1_num = int(sys.argv[4])
    bulk2_num = int(sys.argv[5])
    output_file_hook = open(output_file, 'w')
    output_file_hook.write('\t'.join(['chromosome', 'region_no',
        'start', 'end']) + '\n')
    with open(filepath) as pileup:
        current_chr = ''
        current_start = 0
        current_end = 0
        current_region_no = 0
        for line in pileup:
            line_list = line.strip().split()
            (chromosome, position, reference) = line_list[0:3]

            if not current_chr == chromosome:
                if not current_chr == '':
                    output_file_hook.write('\t'.join([str(val) for val in
                        [current_chr, current_region_no, current_start,
                        current_end]]) + '\n')
                current_chr = chromosome
                current_start = int(position)
                current_end = int(position)
                current_region_no = 1

            (ref_par_depth, ref_par_reads, ref_par_quals) = (
                extract_relevant_columns(line_list, ref_par_num))
            (oth_par_depth, oth_par_reads, oth_par_quals) = (
                extract_relevant_columns(line_list, oth_par_num))
            (bulk1_depth, bulk1_reads, bulk1_quals) = (
                extract_relevant_columns(line_list, bulk1_num))
            (bulk2_depth, bulk2_reads, bulk2_quals) = (
                extract_relevant_columns(line_list, bulk2_num))

            # check for min read depth
            if (int(ref_par_depth) >= min_depth and
                int(oth_par_depth) >= min_depth and
                int(bulk1_depth) >= min_depth and
                int(bulk2_depth) >= min_depth):

                if int(position) == (current_end + 1):
                    current_end = int(position)
                else:
                    output_file_hook.write('\t'.join([str(val) for val in
                        [current_chr, current_region_no, current_start,
                        current_end]]) + '\n')
                    current_start = int(position)
                    current_end = int(position)
                    current_region_no += 1
    output_file_hook.close()
