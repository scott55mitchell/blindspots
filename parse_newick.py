#!/usr/bin/env python2.7
import argparse
from ete3 import Tree
from multiprocessing import Pool
import os.path
import os

parser = argparse.ArgumentParser(description="Outputs a csv file named cladistics.csv, with the low coverage position "
                                             "of interest in the first column of each line, then the number of "
                                             "isolates with that variant, and if the isolates make up a "
                                             "monophyletic group on the tree provided.")
parser.add_argument("--nproc", dest="number_cores", help="Enter the number of cores desired for processing, default "
                                                         "is 1.", required=False, default=1)
parser.add_argument("-t", "--tree", dest="tree_file", help="Input tree file in newick format", required=True)
parser.add_argument("-i", "--csv", dest="pos_csv", help="Input csv file output by create_iso_csv_parallel.py",
                    required=True)
parser.add_argument("-o", "--out", dest="output_file", help="Input desired output filename", required=True)
args = parser.parse_args()
ete3_tree = Tree(args.tree_file)
pos_csv = args.pos_csv
output_file = args.output_file
nproc = args.number_cores
results_list = []
chunk_list = []


lines_per_file = 10000
chunkfile = None
if not os.path.isdir('chunked_input_files/'):
    os.mkdir('chunked_input_files/')
with open(pos_csv) as bigfile:
    for lineno, line in enumerate(bigfile):
        if lineno % lines_per_file == 0:
            if chunkfile:
                chunkfile.close()
            small_filename = 'chunked_input_files/chunked_file_{}.csv'.format(lineno + lines_per_file)
            chunkfile = open(small_filename, "w")
        chunkfile.write(line)
    if chunkfile:
        chunkfile.close()


for chunk_files in os.listdir('chunked_input_files/'):
    if chunk_files.endswith('.csv'):
        chunk_list.append('chunked_input_files/' + chunk_files)


def get_variant_dict(variant_file):
    variant_dict = {}
    for lines in file(variant_file):
        lines = lines.strip()
        lines = lines.replace('[', '').replace(']', '').replace("'", '').replace('"', '')
        variant = lines.strip().split(',')[0]
        isolates = lines.strip().split(',')[1:]
        variant_dict[variant] = isolates
    return variant_dict


def get_ete3_monophyletic_status(input_list):
    results_dict = {}
    for variant, clades in get_variant_dict(input_list).items():
        isolates_monophyly_list = []
        number_isolates_with_variant = len(clades)
        results = ete3_tree.check_monophyly(values=clades, target_attr='name', ignore_missing=True)
        monophyletic_status = results[1]
        isolates_for_monophyly = results[2]
        for isolate in isolates_for_monophyly:
            isolate = str(isolate).strip().split('--')[1]
            if isolate not in isolates_monophyly_list:
                isolates_monophyly_list.append(isolate)
        if len(isolates_monophyly_list) == 0:
            isolates_monophyly_list = 'None'
        if len(isolates_monophyly_list) > 30:
            isolates_monophyly_list = 'More than Thirty Isolates Needed for Monophyly'
        if number_isolates_with_variant == 1:
            isolates_monophyly_list = 'Single Isolate With Blindspot'
        results_dict[variant] = monophyletic_status, number_isolates_with_variant, isolates_monophyly_list, \
                            len(isolates_monophyly_list)
    return results_dict


def main():
    load_pool = Pool(processes=int(nproc))
    results = load_pool.map(get_ete3_monophyletic_status, chunk_list)
    load_pool.close()
    load_pool.join()
    with open(output_file, 'w') as f:
        for output_dictionary in sorted(results):
            for position, isolates in output_dictionary.items():
                f.writelines([str(position), str(isolates)])
                f.writelines('\n')
    with open('cladistics_cleaned.csv', 'w') as g:
        for file_line in file('cladistics.csv'):
            file_line = file_line.replace('[', '').replace(']', '').replace('(', ',').replace(')', '').\
                replace("'", '').replace(' ', '')
            g.writelines(file_line)
    os.rename('cladistics_cleaned.csv', output_file)


if __name__ == '__main__':
    main()
