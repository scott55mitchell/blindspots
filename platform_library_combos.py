#!/usr/bin/env python2.7
from ete3 import Tree
from multiprocessing import Pool
import os.path


chunk_list = []
lines_per_file = 10000
chunkfile = None
with open('/home/smitchel/blindspots/blindspots_isolates_cleaned.csv') as bigfile:
    for lineno, line in enumerate(bigfile):
        if lineno % lines_per_file == 0:
            if chunkfile:
                chunkfile.close()
            small_filename = '/home/smitchel/blindspots/polyphyletic/chunked_input_files/chunked_file_{}.csv'.format\
                (lineno + lines_per_file)
            chunkfile = open(small_filename, "w")
        chunkfile.write(line)
    if chunkfile:
        chunkfile.close()

for chunk_files in os.listdir('/home/smitchel/blindspots/polyphyletic/chunked_input_files/'):
    if chunk_files.endswith('.csv'):
        chunk_list.append('/home/smitchel/blindspots/polyphyletic/chunked_input_files/' + chunk_files)

t = Tree('/home/smitchel/blindspots/phylogeny/blindspots.bipartitions', format=1)


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


mono_list = []
for line in file('/home/smitchel/blindspots/monophyletic/monophyletic_blindspots.csv'):
    line = line.strip()
    position = line.split(',')[0]
    number_of_isolates = line.split(',')[2]
    if int(number_of_isolates) > 4:
        mono_list.append(position)

def get_position_dicts():
    mono_position_dict = {}
    para_position_dict = {}
    for line in file('/home/smitchel/blindspots/monophyletic/monophyletic_blindspots.csv'):
        line = line.strip()
        position = line.split(',')[0]
        number_isolates = line.split(',')[2]
        if int(number_isolates) > 4:
            mono_position_dict[position] = number_isolates
    for line in file('/home/smitchel/blindspots/paraphyletic/paraphyletic_blindspots.csv'):
        line = line.strip()
        position = line.split(',')[0]
        number_isolates = line.split(',')[2]
        para_position_dict[position] = number_isolates
    return mono_position_dict, para_position_dict

mono_position_dict = get_position_dicts()[0]
platform_list = ['Illumina HiSeq 2000','Illumina HiSeq 2500','NextSeq 500','Illumina MiSeq','Illumina MiniSeq']
platform_dict = {}
position_isolate_list_dict = {}
position_list = []



for line in file('/home/smitchel/blindspots/blindspots_isolates_cleaned.csv'):
    position = line.strip().split(',')[0]
    isolates = line.strip().split(',')[1:]
    position_list.append(position)
    position_isolate_list_dict[position] = isolates




pos_prep_dict = {}
for line in file('/home/smitchel/blindspots/blindspots_library_prep.csv'):
    line = line.strip().replace('[','').replace(']','').replace("'",'').split(',')[0:-1]
    file_position = line[0]
    prep = line[1:]
    pos_prep_dict[file_position] = prep


pos_platform_dict = {}
for line in file('/home/smitchel/blindspots/blindspots_instrument.csv'):
    line = line.strip().replace('[','').replace(']','').replace("'",'').split(',')[0:-1]
    file_position = line[0]
    platform = line[1:]
    pos_platform_dict[file_position] = platform


for platform in platform_list:
    platform_position_list = []
    for line in file('/home/smitchel/blindspots/blindspots_instrument.csv'):
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


def get_max_mono(input_list):
    position_iso_list_dict = {}
    for line in file(input_list):
        line = line.strip().split(',')
        position = line[0]
        print position
        isolates = line[1:]
        isolate_list = []
        if len(isolates) > 4:
            for node_number, descendants in node_descendants_dict.items():
                if set(descendants).issubset(set(isolates)):
                    for isolate in descendants:
                        if isolate not in isolate_list:
                            isolate_list.append(isolate)
        if len(isolate_list) > 4:
            position_iso_list_dict[position] = isolate_list
    return position_iso_list_dict


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

    with open('/home/smitchel/blindspots/library_combos/' + platform + '_' + prep + '.csv','wb') as f:
        f.writelines('#position,number_of_isolates_with_low_coverage_position,number_of_isolates_with_library_prep_combo')
        f.writelines('\n')
        for position, number_isos_count in sorted(output_dict.iteritems(),key=lambda x:int(x[0])):
            pos_count = number_isos_count[0]
            plat_lib_count = number_isos_count[1]
            f.writelines([str(position),',',str(pos_count),',',str(plat_lib_count)])
            f.writelines('\n')





