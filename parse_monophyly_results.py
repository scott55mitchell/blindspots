#!/usr/bin/env python2.7
import argparse
import os
from ete3 import Tree
import shutil


parser = argparse.ArgumentParser(description="")
parser.add_argument("-c", "--cutoff", dest="genome_cutoff", help="Input the threshold for number of genomes"
                                                                         "that must have the low coverage position. "
                                                                         "Low coverage positions that are in less than "
                                                                         "the threshold number of genomes will be "
                                                                         "excluded.", required=True)
parser.add_argument("-i", "-input", dest="pos_iso", help="Input the csv file output by create_iso_csv_parallel.py")
parser.add_argument("-t", "--tree", dest="tree_file", help="Input tree file in newick format", required=True)
args = parser.parse_args()
pos_iso = args.pos_iso
t = Tree(args.tree_file, format=1)
threshold = int(args.genome_cutoff)
polyphyly_list = []
paraphyly_list = []
monophyly_list = []


with open('monophyletic_positions.csv', 'wb') as f:
    with open('paraphyletic_positions.csv', 'wb') as g:
        with open('polyphyletic_positions.csv', 'wb') as h:
            for line in file('cladistics.csv'):
                line = line.strip()
                position = line.split(',')[0]
                isolates = line.split(',')[2]
                if int(isolates) >= threshold:
                    if 'monophyletic' in line:
                        monophyly_list.append(position)
                        line = line.replace(',None,4', '')
                        if line.endswith('Monophyly, 43'):
                            line = line.replace(',43', '')
                        f.writelines(line)
                        f.writelines('\n')
                    if 'paraphyletic' in line:
                        paraphyly_list.append(line.split(',')[0])
                        g.writelines(line)
                        g.writelines('\n')
                    if 'polyphyletic' in line:
                        polyphyly_list.append(line.split(',')[0])
                        h.writelines(line)
                        h.writelines('\n')


with open('polyphyletic_position_isolates.csv', 'wb') as f:
    with open('paraphyletic_position_isolates.csv', 'wb') as g:
        with open('monophyletic_position_isolates.csv', 'wb') as h:
            for line in file(pos_iso):
                position = line.strip().split(',')[0]
                isolates = line.strip().split(',')[1:]
                if len(isolates) >= int(threshold):
                    if position in polyphyly_list:
                        f.writelines(line)
                    if position in paraphyly_list:
                        g.writelines(line)
                    if position in monophyly_list:
                        h.writelines(line)

os.rename('polyphyletic_position_isolates.csv', 'polyphyletic_positions.csv')
os.rename('paraphyletic_position_isolates.csv', 'paraphyletic_positions.csv')
os.rename('monophyletic_position_isolates.csv', 'monophyletic_positions.csv')


def search_by_size(node, size):
    node_descendants_dict = {}
    count = 0
    for node in node.traverse():
        if len(node) > size:
            descendant_list = []
            descendants = node.get_descendants()
            for descendant in descendants:
                if descendant.name.startswith('SRR'):
                    descendant_list.append(descendant.name)
            node_descendants_dict[count] = descendant_list
            count += 1
    return node_descendants_dict

node_descendants_dict = search_by_size(t, size=4)

# Creates the platform and library information for each position
unique_instrument_list = []
count_instrument_list = []
unique_prep_list = []
count_prep_list = []
library_prep_dict = {}
instrument_dict = {}
for line in file('/home/smitchel/blindspots/SRA_df.csv'):
    line = line.strip().replace('"','').split(',')
    isolate = line[1]
    library_prep = line[4]
    instrument = line[16]
    if 'Only' in library_prep:
        library_prep = 'Unknown'
    if library_prep != 'Library_prep':
        library_prep_dict[isolate] = library_prep
    if instrument == '0':
        instrument = line[17]
    if instrument != 'Instrument':
        instrument_dict[isolate] = instrument
    if instrument not in unique_instrument_list:
        unique_instrument_list.append(instrument)
    if library_prep not in unique_prep_list:
        unique_prep_list.append(library_prep)
    count_prep_list.append(library_prep)
    count_instrument_list.append(instrument)

position_instrument_dict = {}
unique_instrument_list = []
for line in file('position_isolates.csv'):
    instrument_list = []
    line = line.strip().split(',')
    position = line[0]
    isolates = line[1:]
    for isolate in isolates:
        instrument = instrument_dict[isolate]
        if instrument not in instrument_list:
            instrument_list.append(instrument)
        if instrument not in unique_instrument_list:
            unique_instrument_list.append(instrument)
    position_instrument_dict[position] = instrument_list

position_prep_dict = {}
unique_prep_list = []
for line in file('position_isolates.csv'):
    prep_list = []
    line = line.strip().split(',')
    position = line[0]
    isolates = line[1:]
    for isolate in isolates:
        prep = library_prep_dict[isolate]
        if prep not in prep_list:
            prep_list.append(prep)
        if prep not in unique_prep_list:
            unique_prep_list.append(prep)
    position_prep_dict[position] = prep_list

with open('blindspots_library_prep.csv','wb') as f:
    for position, library_preps in position_prep_dict.items():
        f.writelines(position)
        f.writelines(',')
        for library_prep in library_preps:
            f.writelines(library_prep)
            f.writelines(',')
        f.writelines('\n')

with open('blindspots_instrument.csv','wb') as f:
    for position, instruments in position_instrument_dict.items():
        f.writelines(position)
        f.writelines(',')
        for instrument in instruments:
            f.writelines(instrument)
            f.writelines(',')
        f.writelines('\n')

# Creates library/prep combination info
mono_list = []
for line in file('monophyletic_positions.csv'):
    line = line.strip()
    position = line.split(',')[0]
    number_of_isolates = line.split(',')[1:]
    if len(number_of_isolates) > threshold:
        mono_list.append(position)

platform_list = ['Illumina HiSeq 2000','Illumina HiSeq 2500','NextSeq 500','Illumina MiSeq','Illumina MiniSeq']
platform_dict = {}
position_isolate_list_dict = {}
position_list = []

for line in file('position_isolates.csv'):
    position = line.strip().split(',')[0]
    isolates = line.strip().split(',')[1:]
    position_list.append(position)
    position_isolate_list_dict[position] = isolates

pos_prep_dict = {}
for line in file('blindspots_library_prep.csv'):
    line = line.strip().replace('[','').replace(']','').replace("'",'').split(',')[0:-1]
    file_position = line[0]
    prep = line[1:]
    pos_prep_dict[file_position] = prep


pos_platform_dict = {}
for line in file('blindspots_instrument.csv'):
    line = line.strip().replace('[','').replace(']','').replace("'",'').split(',')[0:-1]
    file_position = line[0]
    platform = line[1:]
    pos_platform_dict[file_position] = platform


for platform in platform_list:
    platform_position_list = []
    for line in file('blindspots_instrument.csv'):
        line = line.strip().replace('[','').replace(']','').replace("'",'')
        line = line.split(',')[0:-1]
        position = line[0]
        instruments = line[1:]
        if platform in instruments:
            platform_position_list.append(position)
    platform_dict[platform] = platform_position_list

platform_prep_dict = {}
for platform in platform_list:
    prep_list = []
    for line in file('/home/smitchel/blindspots/SRA_df.csv'):
        line = line.strip().replace('"','').split(',')
        library_prep = line[4]
        instrument = line[16]
        if instrument == '0':
            instrument = line[17]
        if library_prep != 'Library_prep':
            if 'Only' not in library_prep:
                if 'house' not in library_prep:
                    if library_prep not in prep_list:
                        prep_list.append(library_prep)
    platform_prep_dict[platform] = prep_list


isolate_platform_prep_dict = {}
platform_prep_dict2 = {}
for platform, prep_list in platform_prep_dict.items():
    for prep in prep_list:
        platform_prep_isolate_list = []
        for line in file('/home/smitchel/blindspots/SRA_df.csv'):
            line = line.strip().replace('"', '').split(',')
            isolate = line[1]
            library_prep = line[4]
            instrument = line[16]
            if instrument == '0':
                instrument = line[17]
            isolate_platform_prep_dict[isolate] = [platform, library_prep]
            if platform == instrument:
                if prep == library_prep:
                    platform_prep_isolate_list.append(isolate)
        if len(platform_prep_isolate_list) > 19:
            platform_prep_dict2[(platform,prep)] = platform_prep_isolate_list

if os.path.exists('library_combos/'):
    shutil.rmtree('library_combos/')
os.mkdir('library_combos/')
for platform_prep, isolates_list in platform_prep_dict2.items():
    output_dict = {}
    platform = platform_prep[0]
    prep = platform_prep[1]
    for position in position_list:
        mono_poly_list = []
        isos = list(set(isolates_list)& set(position_isolate_list_dict[position]))
        for node_number, descendants in node_descendants_dict.items():
            if set(descendants).issubset(set(isos)):
                for descendant in descendants:
                    if descendant not in mono_poly_list:
                        mono_poly_list.append(descendant)
        number_mono_poly_isos = len(mono_poly_list)
        num_isos = len(isos) - number_mono_poly_isos
        if platform in pos_platform_dict[position]:
            if prep in pos_prep_dict[position]:
                print platform_prep, position, num_isos, len(isolates_list), number_mono_poly_isos, (len(isolates_list) - number_mono_poly_isos)
                if num_isos > 0:
                    if position not in mono_list:
                        output_dict[position] = [num_isos,(len(isolates_list)-number_mono_poly_isos)]

    with open('library_combos/' + platform + '_' + prep + '.csv','wb') as f:
        f.writelines('#position,number_of_isolates_with_low_coverage_position,number_of_isolates_with_library_prep_combo')
        f.writelines('\n')
        for position, number_isos_count in sorted(output_dict.iteritems(),key=lambda x:int(x[0])):
            pos_count = number_isos_count[0]
            plat_lib_count = number_isos_count[1]
            f.writelines([str(position),',',str(pos_count),',',str(plat_lib_count)])
            f.writelines('\n')



