#######################################################################################################################
#
# Additional functions for Recursive Editing
#
#######################################################################################################################



import os
import pandas as pd
import numpy as np
from Bio.Seq import Seq



BASE_PATH = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
FLASHFRY_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))), 'flash_fry')



def append_edit_df(df, cutting_score, seq_type, guide_id, level_num, cutsite):
    """
    append additional columns to edit data frame
    """
    df['REtarget_Score_Contribution'] = df['Predicted frequency'].values / 100 * cutting_score
    df['Seq_Type'] = seq_type
    df.loc[df['Category'] == 'ins', 'Genotype position'] = 2
    two_bp_ins_index = [index for index in range(len(df['Inserted Bases'].values)) if len(str(df['Inserted Bases'].values[index])) == 2]
    df.iloc[two_bp_ins_index, 1] = -2
    df['Cutsite'] = cutsite - (df['Length'].values - df['Genotype position'].values)
    df['Edit_ID'] = np.array([f'{guide_id}_{k}' for k in range(1, len(df) + 1)])
    df['Guide_ID'] = guide_id
    df['Level'] = level_num
    return df


def append_guide_df(df, guide_id, guide, cutting_score, dscore, seq, seq_type, cutsite, target_origin, edit_df, old_guide_df, task_c):
    """
    append new row to guide data frame
    """
    from recursive_editing import config

    edit_df = edit_df[edit_df['Guide_ID'] == guide_id]
    editing_score = np.sum(edit_df['Predicted frequency'].tolist())*0.01
    REtarget_score = np.sum(edit_df['REtarget_Score_Contribution'].tolist())
    level = int(len(guide_id.split('_')))
    if guide in old_guide_df['Guide_RNA'].values:
        multi_occurance = 1
        mga_REtarget_score = REtarget_score + np.sum(old_guide_df.loc[old_guide_df['Guide_RNA'] == guide, 'REtarget_Score'].values)
    else:
        multi_occurance = 0
        mga_REtarget_score = REtarget_score
    if config.CHECK_OFFT_GENOME is True and config.CHECK_OFFT_SCORES is True:
        flash_fry_results = pd.read_csv(os.path.join(FLASHFRY_PATH, 'level_'+str(task_c)+'.output.scored'), delimiter='\t', index_col=None)
        guide_scores = flash_fry_results[flash_fry_results['contig'] == guide]
        d_cfd_otmax = guide_scores['DoenchCFD_maxOT'].values[0]
        d_cfd_spec = guide_scores['DoenchCFD_specificityscore'].values[0]
        hsu_score = guide_scores['Hsu2013'].values[0]
        num_ot = guide_scores['0-1-2-3-4_mismatch'].values[0]
        closest_ot = guide_scores['basesDiffToClosestHit'].values[0]
        num_closest_ot = guide_scores['closestHitCount'].values[0]
        ff_guide_orientation = guide_scores['orientation'].values[0]
        ff_pam = guide_scores['target'].values[0][-3:]
    elif config.CHECK_OFFT_GENOME is True:
        flash_fry_results = pd.read_csv(os.path.join(FLASHFRY_PATH, 'level_'+str(task_c)+'.output.scored'), delimiter='\t', index_col=None)
        guide_scores = flash_fry_results[flash_fry_results['contig'] == guide]
        num_ot = guide_scores['0-1-2-3-4_mismatch'].values[0]
        closest_ot = guide_scores['basesDiffToClosestHit'].values[0]
        num_closest_ot = guide_scores['closestHitCount'].values[0]
        ff_guide_orientation = guide_scores['orientation'].values[0]
        ff_pam = guide_scores['target'].values[0][-3:]
        d_cfd_otmax, d_cfd_spec, hsu_score = -1, -1, -1
    elif config.CHECK_OFFT_SCORES is True:
        flash_fry_results = pd.read_csv(os.path.join(FLASHFRY_PATH, 'level_'+str(task_c)+'.output.scored'), delimiter='\t', index_col=None)
        guide_scores = flash_fry_results[flash_fry_results['contig'] == guide]
        d_cfd_otmax = guide_scores['DoenchCFD_maxOT'].values[0]
        d_cfd_spec = guide_scores['DoenchCFD_specificityscore'].values[0]
        hsu_score = guide_scores['Hsu2013'].values[0]
        ff_guide_orientation = guide_scores['orientation'].values[0]
        ff_pam = guide_scores['target'].values[0][-3:]
        num_ot, closest_ot, num_closest_ot = '-1,-1,-1,-1,-1', -1, -1
    else:
        d_cfd_otmax, d_cfd_spec, hsu_score, num_ot, closest_ot, num_closest_ot, ff_guide_orientation, ff_pam = -1, -1, -1, '-1,-1,-1,-1,-1', -1, -1, seq_type, 'NaN'
    new_row = pd.DataFrame(list(zip([guide_id], [guide], [level], [REtarget_score], [mga_REtarget_score], [cutting_score], [editing_score], [dscore], [d_cfd_otmax], [d_cfd_spec],\
                                    [hsu_score], [num_ot], [closest_ot], [num_closest_ot], [ff_guide_orientation], [ff_pam], [seq], [seq_type], [cutsite], [target_origin],\
                                    [multi_occurance])),
                                columns=['Guide_ID', 'Guide_RNA', 'Level', 'REtarget_Score', 'MGA_REtarget_Score', 'Cutting_Score', 'Editing_Score', 'Doench2014_Score',\
                                        'DoenchCFD_MaxOT', 'DoenchCFD_SpecScore', 'Hsu2013_Score', 'Number_Off_Targets', 'Mismatches_Closest_OT', 'Number_Closest_OT',\
                                        'Guide_Orientation', 'PAM', 'Target', 'Target_Type', 'Cut_Position', 'Target_Origin', 'Multi'])
    df = df.append(new_row)
    return df


def append_res_df(df, res_id, level_num, guide_df):
    """
    append new row to result data frame
    """
    from recursive_editing import config

    guide_num = len(guide_df['Guide_ID'])
    guide_rnas = '\n'.join(guide_df['Guide_RNA'].tolist())

    # calculate total REtarget score
    level_REtarget_score = np.zeros((level_num))
    for i in range(1, level_num + 1):
        tmp_df = guide_df[guide_df['Level'] == i]
        level_REtarget_score[i - 1] = (config.LEVEL_SCORE_FACTOR**(i - 1)) * np.sum(tmp_df['REtarget_Score'].tolist())
    tot_REtarget_score = np.sum(level_REtarget_score)

    if level_num >= 2:
        new_row = pd.DataFrame(list(zip([res_id], [level_num], [guide_num], [guide_rnas], [tot_REtarget_score], [level_REtarget_score[0]], [level_REtarget_score[1]])),
            columns=['Res_ID', 'Level_Number', 'Guide_Number', 'Guide_RNAs', 'Total_REtarget_Score', 'REtarget_Score_L1', 'REtarget_Score_L2'])
    else:
        new_row = pd.DataFrame(list(zip([res_id], [level_num], [guide_num], [guide_rnas], [tot_REtarget_score], [level_REtarget_score[0]], [level_REtarget_score[0]])),
            columns=['Res_ID', 'Level_Number', 'Guide_Number', 'Guide_RNAs', 'Total_REtarget_Score', 'REtarget_Score_L1', 'REtarget_Score_L2'])

    df = df.append(new_row)
    return df


def get_rc_seq(seq):
    """
    returns reverse complement of sequence seq as string
    """
    return str(Seq(str(seq)).reverse_complement())


def get_rc_cut(seq, cut):
    """
    returns new cutsite for reverse complement of seq
    """
    return len(seq) - cut


def get_guide(seq, cutsite, seq_type):
    """
    get gRNA from sequence, cutsite context
    """
    if seq_type == 'fwd':
        return seq[cutsite - 17 : cutsite + 3]
    else:
        return str(Seq(seq[cutsite - 3 : cutsite + 17]).reverse_complement())


def get_pam(seq, cutsite, seq_type):
    """
    get gRNA from sequence, cutsite context
    """
    if seq_type == 'fwd':
        return seq[cutsite + 3 : cutsite + 6]
    else:
        return str(Seq(seq[cutsite - 6 : cutsite - 3]).reverse_complement())
