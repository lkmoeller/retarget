#######################################################################################################################
#
# Search sequence window for loci suitable for Recursive Editing
#
# Author: Lukas Moeller - 11/2021
#
#######################################################################################################################



import os
import argparse
from recursive_editing.hdr_optimization import run_retarget
from recursive_editing.utils import BASE_PATH
from recursive_editing import config



# usage: python sequence_search/sequence_search.py -q [add your sequence]
# parse command line arguments
parser = argparse.ArgumentParser(description='selecting best locus for editing')
parser.add_argument('-n', '--name', type=str, help='path to output folder', default='REtarget')
parser.add_argument('-q', '--seq', required=True, type=str, help='sequence to search', default='')
parser.add_argument('-s', '--start', type=int, help='start index', default=-1)
parser.add_argument('-e', '--end', type=int, help='end index', default=-1)
parser.add_argument('-f', '--database', required=False, type=str, help='Name of FleshFry database', default='chr22_cas9ngg_database')
parser.add_argument('-t', '--template', required=False, type=str, help='Sequence of HDR template', default='atccacagcagacccaccctcaccttcagttttattgttttgctccaaacGATataactctgctgcttccactgctctggggctggtaaaaatgagtcccccg')
parser.add_argument('-b', '--blacklist', required=False, type=list, help='gRNA Blacklist: List of gRNAs that should not be used by REtarget', default=[])
args = parser.parse_args()

if str(args.database) != '':
    config.FF_DATABASE_NAME = str(args.database)
if str(args.template) != '':
    config.HDR_TEMPLATE = str(args.template).upper()
if len(list(args.blacklist)) > 0:
    config.BLACKLIST = [str(g).upper() for g in list(args.blacklist)]

# load sequence data
gen_dict = args.seq

# set start and end positions for search
if args.start == -1:
    pos_start = 0
else:
    pos_start = args.start
if args.end == -1:
    pos_end = len(args.seq)
else:
    pos_end = args.end

# set output path
out_path = os.path.join(BASE_PATH, 'data/output_files', str(args.name))
os.makedirs(out_path, exist_ok=True)

# set parameter for search
cutsite = 50
pam_fwd, pam_rev = 'GG', 'CC'

# conduct search
for i in range(pos_start + 50, pos_end - 50):
    
    pam = str(gen_dict[i:i + 2])

    try:
        # print(i/(pos_end-pos_start)*100)
        if pam == pam_fwd:
            seq = str(gen_dict[i - 4 - 50 : i - 4 + 50]).upper()
            seq_type = 'fwd'
            edit_df, guide_df, res_df, layer_num = run_retarget(seq, cutsite, seq_type, store_res=True, store_name=str(args.name), res_id='id_' + str(i - 4) + '_f')
        elif pam == pam_rev:
            seq = str(gen_dict[i + 6 - 50 : i + 6 + 50]).upper()
            seq_type = 'rev'
            edit_df, guide_df, res_df, layer_num = run_retarget(seq, cutsite, seq_type, store_res=True, store_name=str(args.name), res_id='id_' + str(i + 6) + '_r')
    
    except Exception as e:
        pass
