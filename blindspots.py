#!/usr/bin/env python2.7
from Bio.SeqRecord import SeqRecord
import os
from Bio import SeqIO
import os.path
import csv
import shutil
import argparse
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Blast.Applications import NcbiblastnCommandline



parser = argparse.ArgumentParser(description="")
parser.add_argument("-c", "--coverage", dest="average_coverage", required=False, default=37,
                    help="Input threshold for a genome's average coverage? Any genome with an average coverage "
                         "below this threshold will not be processed. Default is 37")
parser.add_argument("-i", "--input", dest="input_dir", help="Input the directory with all of your .LC_position.txt, "
                                                            ".vcf.gz and average-coverage.txt files? For more "
                                                            "information please see the instructions.txt file",
                    required=True)
parser.add_argument("--nproc", dest="number_cores", help="Enter the number of cores desired for processing, default"
                   "is 1.", required=False, default=1)
parser.add_argument("--GTRGAMMA", dest='raxml_model', help='Input --GTRGAMMA if you would like the ML search performed'
                                                           'with GTRGAMMA instead of GTRCAT. GTRCAT is much faster if '
                                                           'the dataset includes many genomes', required=False,
                    default='GTRCAT')
parser.add_argument('-t', '--tree', dest='tree', help='Input the newick file output by RAxML', required=True)

args = parser.parse_args()
input_dir = args.input_dir
avg_coverage = args.average_coverage
nproc = args.number_cores
tree = args.tree
raxml_model = args.raxml_model

if os.path.exists('blindspots/'):
    shutil.rmtree('blindspots/')
os.mkdir('blindspots/')

# Takes directory with .LC_positions files, .vcf.gz files, and .average_coverage files and creates
# 'blindspots_isolates.fofn' and 'blindspots_isolates.csv'
run_create_iso_csv_parallel = './src/create_iso_csv_parallel.py -i ' +  input_dir + ' -o position_isolates.csv --nproc ' + \
    nproc + ' -c ' + avg_coverage + ' -v isolate_vcfs.fofn'
os.system(run_create_iso_csv_parallel)


# Creates multi-fasta file to be used by RAxML
run_make_snp_fasta = './src/make-snp-fasta.py ' + ' -i isolate_vcfs.fofn -o snps.fasta --outgroups ' \
                                              '--og_type ' + "'bovis' 'canetti' " + '--nproc ' + str(nproc)
os.system(run_make_snp_fasta)


# Generates newick file to be used for phylogenetic analysis.
run_raxml = 'raxmlHPC-PTHREADS -T ' + str(nproc)  + ' -s snps.fasta -n snps.out -m ' \
            + raxml_model + ' -f a -x 12345 -N 100 -p 12345 -o mbovis,mcanetti'
os.system(run_raxml)
os.mkdir('blindspots/phylogeny/')
os.rename('RAxML_bipartitionsBranchLabels.snps.out','blindspots/phylogeny/RAxML_bipartitionsBranchLabels.snps.out')
os.rename('RAxML_bipartitions.snps.out','blindspots/phylogeny/snps.newick')
os.rename('RAxML_bootstrap.snps.out','blindspots/phylogeny/RAxML_bootstrap.snps.out')
os.rename('RAxML_info.snps.out','blindspots/phylogeny/RAxML_info.snps.out')
os.rename('snps.fasta','blindspots/phylogeny/snps.fasta')
os.rename('isolate_vcfs.fofn','blindspots/phylogeny/isolate_vcfs.fofn')


# Performs phylogenetic analysis on blindspots to define possible monophyletic artifacts.
run_parse_newick = './src/parse_newick.py --nproc ' + nproc + ' -t blindspots/phylogeny/snps.newick ' \
                                                              '-i position_isolates.csv -o cladistics.csv'
os.system(run_parse_newick)
shutil.rmtree('chunked_input_files/')

# Separates the cladistics.csv file into mono, para and polyphyletic files and lists the corresponding genomes,
# as well as creates files for each platform/library prep combination that include low coverage positions, isolates
# with that low coverage position, and the total number of isolates with that library prep/platform comination
run_parse_monophyly_results = './src/parse_monophyly_results.py -c ' + avg_coverage + \
                              ' -i position_isolates.csv -t blindspots/phylogeny/snps.newick'
os.system(run_parse_monophyly_results)

# Creates similar files to the 'run_parse_monophyly_results' command above, but this uses a subset of isolates for the
# combinations above.
run_platform_library_combos = './src/platform_library_combos.py -t blindspots/phylogeny/snps.newick'


