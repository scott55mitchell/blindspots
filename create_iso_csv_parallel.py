#!/usr/bin/env python2.7
import argparse
import csv
import os
from multiprocessing import Pool

parser = argparse.ArgumentParser(description="This script outputs two files. The first is a file of file names that "
                                             "will be used as input for make-snp-fasta.py named 'vcfs.fofn'. The "
                                             "second is a csv file named position_isolates.csv for the "
                                             "parse_newick.py script. The csv is in the format of position,isolate1,"
                                             "isolate2,isolate3....etc. This script has several requirements; the "
                                             "following files must all be in the same directory. For each genome in "
                                             "the dataset, three files are needed. First, a vcf file is needed to call "
                                             "variants, which is needed to create the phylogeny This file must have "
                                             "the file extension of '.vcf.gz'. The second is a file with a single "
                                             "value of the average coverage for all positions in the genome, which "
                                             "will be used to exclude genomes with potentially unreliable assemblies. "
                                             "This file must have the file extension of .average_coverage.txt. The "
                                             "third and last file is a space separated text file with each low "
                                             "coverage position (determined by your criteria) in the first column of "
                                             "each line, and the coverage for that position in the second column of "
                                             "each line. This file must have the file extension of .LC_positions.txt ")
parser.add_argument("-i", "--input", dest="variant_directory", help="Enter directory that contains all of the required "
                                                                    "files described in the main description",
                    required=True)
parser.add_argument("--nproc", dest="number_cores", help="Enter the number of cores desired for processing, default "
                                                         "is 1.", required=False, default=1)
parser.add_argument("-c", dest='coverage_cutoff', help='Enter the desired cutoff for average coverage of the genome. '
                                                       'Genomes below this cutoff will not be included.', default=37)
parser.add_argument("-o", dest='output_file', help='Enter the desired output filename.')
parser.add_argument("-v", dest='vcf_fofn', help="Enter the desired name for the file of vcf file names.",
                    default='vcf_fofn')
args = parser.parse_args()
vcf_fofn = args.vcf_fofn
coverage_threshold = args.coverage_cutoff
output_file = args.output_file
variant_dir = args.variant_directory
nproc = args.number_cores


def check_avg_coverage(coverage_file_path):
    for line in file(coverage_file_path):
        line = line.strip()
        if float(line) > int(coverage_threshold):
            return True
        else:
            return False


def get_fofn(blindspot_dir):
    vcf_file_list = []
    position_file_list = []
    coverage_list = []
    for position_file in os.listdir(blindspot_dir):
        if position_file.endswith('.LC_positions.txt'):
            isolate = position_file.split('.')[0]
            if os.path.exists(blindspot_dir + isolate + '.vcf.gz'):
                if check_avg_coverage(blindspot_dir + isolate + '.average_coverage.txt') is True:
                    coverage_list.append('okay')
                    position_file_list.append(blindspot_dir + position_file)
                    vcf_file_list.append(blindspot_dir + isolate + '.vcf.gz')
                else:
                    coverage_list.append('low')
    print 'Your average genome coverage threshold excluded ' , coverage_list.count('low') , 'genomes. Your dataset is ' \
                                                            , coverage_list.count('okay') , 'genomes.'

    return position_file_list, vcf_file_list


def write_fofn_files():
    pos_vcf_lists = get_fofn(variant_dir)
    position_file_list = pos_vcf_lists[0]
    vcf_file_list = pos_vcf_lists[1]
    with open('lc_pos_file.fofn', 'w') as file_handle1:
        for file_path1 in position_file_list:
            file_handle1.write(file_path1+'\n')
    with open(vcf_fofn, 'w') as file_handle2:
        for file_path2 in vcf_file_list:
            file_handle2.write(file_path2+'\n')

write_fofn_files()

def get_all_positions(variant_fofn):
    all_position_list = []
    var_fofn = file(variant_fofn)
    for variant_file in var_fofn:
        variant_file = variant_file.strip()
        for line in file(variant_file):
            position = line.strip().split('\t')[0]
            all_position_list.append(position)
    unique_position_list = set(sorted(all_position_list))
    return unique_position_list


def get_blindspot_dict(variant_fofn):
    isolate_blindspots_dict = {}
    var_fofn = file(variant_fofn)
    for variant_file in var_fofn:
        isolate = variant_file.split('/')[6].split('.')[0]
        position_list = []
        variant_file = variant_file.strip()
        for line in file(variant_file):
            position = line.strip().split('\t')[0]
            position_list.append(position)
        isolate_blindspots_dict[isolate] = position_list
    return isolate_blindspots_dict


iso_bs_dict = get_blindspot_dict('lc_pos_file.fofn')


def get_output_dict(position):
    output_dict = {}
    isolate_list = []
    for iso, pos_list in iso_bs_dict.items():
        if position in pos_list:
            isolate_list.append(iso)
    output_dict[position] = isolate_list
    return output_dict


def write_files(output_dictionary_list):
    with open(output_file, 'w') as csvfile:
        writer = csv.writer(csvfile)
        for output_dictionary in sorted(output_dictionary_list):
            for position, isolates in sorted(output_dictionary.items()):
                writer.writerow([position, isolates])


def main():
    load_pool = Pool(processes=int(nproc))
    results = load_pool.map(get_output_dict, get_all_positions('lc_pos_file.fofn'))
    load_pool.close()
    load_pool.join()
    write_files(results)
    with open('position_isolates_cleaned.csv', 'wb') as file_guy:
        for line in file('position_isolates.csv'):
            line = line.replace('[', '').replace(']', '').replace('(', '').replace("'", '').\
                replace(')', '').replace('"', '').replace(' ', '')
            file_guy.writelines(line)
    os.remove('lc_pos_file.fofn')
    os.rename('position_isolates_cleaned.csv', output_file)


if __name__ == '__main__':
    main()
