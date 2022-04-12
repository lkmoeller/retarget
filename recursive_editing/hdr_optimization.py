#######################################################################################################################
#
# Recursive optimization to find gRNAs suitable for Recursive Editing
#
#######################################################################################################################



import pandas as pd
import numpy as np
import os, re, subprocess, logging
from prediction_tools import predict_lindel, predict_indelphi
from recursive_editing import config
from recursive_editing.doench_score import calcDoenchScore
from prediction_tools.predict_indelphi import init_indelphi
from recursive_editing.utils import (
    append_edit_df,
    append_guide_df,
    append_res_df,
    get_rc_seq,
    get_rc_cut,
    get_guide,
    get_pam,
    BASE_PATH,
    FLASHFRY_PATH
)



# set logging configuration
if config.VERBOSE == True:
    logging.basicConfig(format='%(asctime)s: %(message)s', level=logging.INFO)
else:
    logging.basicConfig(format='%(asctime)s: %(message)s', level=logging.WARNING)



def run_prediction(model, *args):
    """
    run indel prediction tools
    """
    if model == 0:
        # predicted edit outcomes with Lindel
        return predict_lindel.run_lindel(*args)
    elif model == 1:
        # predicted edit outcomes with inDelphi
        return predict_indelphi.run_indelphi(*args)


def check_guide(guide, guide_scores=None, wt_seq=None, check_offt_prox=True, l1=False):
    """
    function checks quality of guides and returnes 2014 Doench scores
    """
    
    # ensure guide is str in uppercase letters
    guide = str(guide).upper()

    # 1. check on-target scores
    # check activity of sgRNA with Doench 2014 score
    dscore = calcDoenchScore(guide)
    if dscore < config.D2014ON_MIN:
        return False, dscore, 'Doench 2014 score of guide {} ({}) is below threshold of {}.'.format(guide, dscore, config.D2014ON_MIN)

    # 2. check if guide cuts WT sequence in proximity of first cut
    if (config.CHECK_OFFT_PROX == True) and (check_offt_prox == True):
        if wt_seq is not None:
            search_str = re.compile(guide + '.GG')
            wt_seq = str(wt_seq).upper()
            if bool(search_str.search(wt_seq)) is True:
                return False, dscore, 'Guide {} cuts WT sequence in proximity of initial cutsite.'.format(guide)
        else:
            return False, dscore, 'WT sequence not specified. Checking of off-targets in proximity of initial cut site not possible.'
    
    # 3. check off-targets
    if config.CHECK_OFFT_GENOME == True:
        if guide_scores is not None:
            off_targets = guide_scores['0-1-2-3-4_mismatch'].values
            if len(off_targets) > 0:
                off_targets = np.array(off_targets[0].split(','), dtype=int)
                off_check = off_targets - np.array(config.MAX_OFF_TAR, dtype=int)
                if l1 == True:
                    off_check[0] -= 1
                if np.sum(np.where(off_check > 0, 1, 0)) > 0:
                    return False, dscore, 'Guide {} has more off-targets than allowed. Number of off-targets with 0,1,2,3,4 mismatches: {}. Allowed off-targets: {}.'.format(guide, off_targets, config.MAX_OFF_TAR)
            else:
                return False, dscore, 'Guide {} sorted out due to error in off-target calculation: Flash Fry overflow'.format(guide)
        else:
            return False, dscore, 'Guide {} sorted out due to error in off-target calculation: Flash Fry overflow'.format(guide)
    
    # 4. check off-target scores
    if config.CHECK_OFFT_SCORES == True:
        if guide_scores is not None:
            # Doench 2016 CFD score - highest off-target score
            score_value = guide_scores['DoenchCFD_maxOT'].values
            if len(score_value) > 0:
                if np.float(score_value[0]) > config.DCFD_MAXOT_MAX:
                    return False, dscore, 'Doench CFD MaxOT score of guide {} ({}) is above threshold of {}.'.format(guide, score_value, config.DCFD_MAXOT_MAX)
            else:
                return False, dscore, 'Guide {} sorted out due to error in off-target score calculation: Flash Fry overflow'.format(guide)
            # Doench 2016 CFD score - specificity score
            score_value = guide_scores['DoenchCFD_specificityscore'].values
            if len(score_value) > 0:
                if np.float(score_value[0]) < config.DCFD_SPEC_MIN:
                    return False, dscore, 'Doench CFD specificity score of guide {} ({}) is below threshold of {}.'.format(guide, score_value, config.DCFD_SPEC_MIN)
            else:
                return False, dscore, 'Guide {} sorted out due to error in off-target score calculation: Flash Fry overflow'.format(guide)
            # Hsu 2013 score
            score_value = guide_scores['Hsu2013'].values
            if len(score_value) > 0:
                if np.float(score_value[0]) < config.HSU2013_MIN:
                    return False, dscore, 'Hsu 2013 score of guide {} ({}) is below threshold of {}.'.format(guide, score_value, config.HSU2013_MIN)
            else:
                return False, dscore, 'Guide {} sorted out due to error in off-target score calculation: Flash Fry overflow'.format(guide)
        else:
            return False, dscore, 'Guide {} sorted out due to error in off-target score calculation: Flash Fry overflow'.format(guide)
    
    # 5. check GC content
    if config.CHECK_GC == True:
        GC_content = (guide.count('C') + guide.count('G')) / 20
        if (GC_content < config.GC_LOW):
            return False, dscore, 'GC-content of guide {} ({}) below threshold of {}.'.format(guide, GC_content, config.GC_LOW) 
        if (GC_content > config.GC_HIGH):
            return False, dscore, 'GC-content of guide {} ({}) above threshold of {}.'.format(guide, GC_content, config.GC_HIGH)
    
    # 6. check homopolymers
    if config.CHECK_POLY_N == True:
        homopolymer_num = guide.count('AAAA') + guide.count('TTTT') + guide.count('GGGG') + guide.count('CCCC')
        if homopolymer_num > 0:
           return False, dscore, 'Found homopolymer sequence in guide {}.'.format(guide)
    
    # 7. check uracil triplets
    if config.CHECK_TTT == True:
        uracil_triplet_num = guide.count('TTT')
        if uracil_triplet_num > 0:
            return False, dscore, 'Found TTT sequence in guide {}.'.format(guide)
    
    # 8. check if guide cuts HDR template
    search_str = re.compile(guide + '.GG')
    if config.CHECK_TEMPLATE == True:
        try:
            hdr_template = str(config.HDR_TEMPLATE).upper()
            if bool(search_str.search(hdr_template)) is True or bool(search_str.search(get_rc_seq(hdr_template))) is True:
                return False, dscore, 'Guide {} cuts HDR template.'.format(guide)
        except Exception as e:
            return False, dscore, 'Error while checking if gRNA cuts HDR template. Please specify valid HDR template (no whitespaces allowed)'

    # 9. check if guide blacklisted
    if config.CHECK_BLACKLIST == True:
        try:
            blacklist = [str(g).upper() for g in config.BLACKLIST]
            if len(np.intersect1d(guide, blacklist)) > 0:
                return False, dscore, 'Blacklisted guide {} was excluded from analysis.'.format(guide)
        except Exception as e:
            return False, dscore, 'Error while checking if gRNA blacklist. Please specify valid blacklist (format: ["gRNA1","gRNA2", ...])'
    
    # if all criteria met: return valid guide
    return True, dscore, 'Valid guide.'


def find_all(genotype, pam_seq, start, end):
    """
    find all occurrences of pam_seq in genotype & returns list with pam-positions if list(find_all()) used
    """

    while True:
        start = genotype.find(pam_seq, start, end)
        if start == -1:
            return
        yield start
        start += 1


def run_flash_fry(guide_list, database_name, task_c, use_taskset=True):
    """
    run FlashFry to analyze if gRNAs have off-targets (genome-wide off-target analysis)
    """
    # generate fasta file with all guides of level
    with open(os.path.join(FLASHFRY_PATH, 'level_'+str(task_c)+'.fasta'), 'w+') as f_out:
        for guide_pam in guide_list:
            f_out.write('>'+str(guide_pam[:-3])+'\n')
            f_out.write(str(guide_pam)+'\n')
    
    if (os.path.exists(os.path.join(FLASHFRY_PATH, str(database_name))) == True) and (os.path.exists('../flash_fry/FlashFry-assembly-1.12.jar') == True) and (os.path.exists(config.JDK_PATH) == True):
        logging.info('Running FlashFry for off-target analysis')
        if use_taskset == True:
            # run FlashFry
            subprocess.run(['taskset', '-c', str(task_c), config.JDK_PATH, '-Xmx4g', '-jar', '../flash_fry/FlashFry-assembly-1.12.jar', 'discover', '--database', os.path.join(FLASHFRY_PATH, str(database_name)), '--fasta', '../flash_fry/level_'+str(task_c)+'.fasta', '--output', '../flash_fry/level_'+str(task_c)+'.output', '--maximumOffTargets', '5000', '--maxMismatch', '4'], stdout=subprocess.DEVNULL)
            subprocess.run(['taskset', '-c', str(task_c), config.JDK_PATH, '-Xmx4g', '-jar', '../flash_fry/FlashFry-assembly-1.12.jar', 'score', '--input', '../flash_fry/level_'+str(task_c)+'.output', '--output', '../flash_fry/level_'+str(task_c)+'.output.scored', '--scoringMetrics', 'doench2016cfd,hsu2013,minot', '--database', os.path.join(FLASHFRY_PATH, str(database_name))], stdout=subprocess.DEVNULL)
        else:
            # run FlashFry
            subprocess.run([config.JDK_PATH, '-Xmx4g', '-jar', '../flash_fry/FlashFry-assembly-1.12.jar', 'discover', '--database', os.path.join(FLASHFRY_PATH, str(database_name)), '--fasta', '../flash_fry/level_'+str(task_c)+'.fasta', '--output', '../flash_fry/level_'+str(task_c)+'.output', '--maximumOffTargets', '5000', '--maxMismatch', '4'], stdout=subprocess.DEVNULL)
            subprocess.run([config.JDK_PATH, '-Xmx4g', '-jar', '../flash_fry/FlashFry-assembly-1.12.jar', 'score', '--input', '../flash_fry/level_'+str(task_c)+'.output', '--output', '../flash_fry/level_'+str(task_c)+'.output.scored', '--scoringMetrics', 'doench2016cfd,hsu2013,minot', '--database', os.path.join(FLASHFRY_PATH, str(database_name))], stdout=subprocess.DEVNULL)


def select_pam(guide_list, pam_list, flash_fry_results, wt_seq):
    """
    check quality of guides, sort out unsuitable guides, and rank guides according to on-target efficiency
    """

    pam_sorted = []
    for guide, j in zip(guide_list, pam_list):
        if ((config.CHECK_OFFT_GENOME == True) or (config.CHECK_OFFT_SCORES == True)) and (flash_fry_results is not None):
            if guide[:-3] in flash_fry_results['contig'].values:
                # get scores for individual guide
                guide_scores = flash_fry_results[flash_fry_results['contig'] == guide[:-3]]
            else:
                guide_scores = None
        else:
            guide_scores = None

        # check scores for individual guide
        keep_guide, dscore, guide_message = check_guide(guide[:-3], guide_scores=guide_scores, wt_seq=wt_seq)
        logging.info('Evaluating gRNA: ' + str(guide[:-3]) + '. PAM: ' + str(guide[-3:]) + '. ' + str(guide_message))
        if keep_guide == True:
            pam_sorted.append([dscore, j - 4, guide[:-3]])

    # sort according to Doench 2014 score
    pam_sorted = sorted(pam_sorted, key=lambda x: x[0], reverse=True)
    
    # return best-ranking guides
    return pam_sorted[:min(config.MAX_PAM_NUM, len(pam_sorted))]


def find_pam(seq_list, old_cutsite_list, guide_id_list, task_c):
    """
    find all PAMs for potential new guides & run FlashFry
    """
    
    # initialize guide list
    guide_list = []
    # initialize dictionaries to store results
    guide_dict_fwd, pam_dict_fwd, guide_dict_rev, pam_dict_rev = {}, {}, {}, {}

    for seq, old_cutsite, guide_id in zip(seq_list, old_cutsite_list, guide_id_list):
        # get reverse complement of seq
        rc_seq = get_rc_seq(seq)
        rc_old_cutsite = get_rc_cut(seq, old_cutsite)
        
        # get all possible PAM sequence positions
        pam_list_fwd = list(find_all(seq, 'GG', int(old_cutsite + config.MIN_OVERLAP), int(old_cutsite + config.MAX_OVERLAP)))
        pam_dict_fwd[guide_id] = pam_list_fwd
        pam_list_rev = list(find_all(rc_seq, 'GG', int(rc_old_cutsite + config.MIN_OVERLAP), int(rc_old_cutsite + config.MAX_OVERLAP)))
        pam_dict_rev[guide_id] = pam_list_rev
        
        # get potential new guides + PAM
        guide_list_fwd = [str(seq[j - 21:j - 1])+str(get_pam(seq, j - 4, 'fwd')) for j in pam_list_fwd]
        guide_dict_fwd[guide_id] = guide_list_fwd
        guide_list += guide_list_fwd # append list to list
        guide_list_rev = [str(rc_seq[j - 21:j - 1])+str(get_pam(rc_seq, j - 4, 'fwd')) for j in pam_list_rev] # 'fwd' correct here!
        guide_dict_rev[guide_id] = guide_list_rev
        guide_list += guide_list_rev # append list to list
    
    # make sure to list every guide-PAM combination only once
    guide_list = list(np.unique(np.array(guide_list)))

    if (config.CHECK_OFFT_GENOME == True) or (config.CHECK_OFFT_SCORES == True):
        # run flash fry for all potential guides in level
        run_flash_fry(guide_list, config.FF_DATABASE_NAME, task_c, use_taskset=config.USE_TASKSET_JKD)
        # get path to results
        ffpath = os.path.join(FLASHFRY_PATH, 'level_'+str(task_c)+'.output.scored')
        # check if flash fry run successful
        if os.path.exists(ffpath) == True:
            # get flash fry results
            flash_fry_results = pd.read_csv(ffpath, delimiter='\t', index_col=None)
        else:
            flash_fry_results = None
    else:
        flash_fry_results = None
    
    return guide_dict_fwd, pam_dict_fwd, guide_dict_rev, pam_dict_rev, flash_fry_results


def predict_single_guide(seq, cutting_score, seq_type, guide_id, level_num, pam_list, best_edit_score, best_edit_df, best_guide, best_guide_id, best_guide_dscore, best_seq_type, best_cutsite):
    """
    generate predictions for single guide
    """

    for dscore, cutsite, guide in pam_list:
        tmp_edit_df = run_prediction(config.MODEL, seq, cutsite)
        if not tmp_edit_df.empty:
            edit_score = np.sum(tmp_edit_df['Predicted frequency'].values)
            if edit_score > best_edit_score:
                best_edit_df = append_edit_df(tmp_edit_df, cutting_score, seq_type, guide_id, level_num, cutsite)
                best_edit_score, best_guide, best_guide_id, best_guide_dscore, best_seq_type, best_cutsite = edit_score, guide, guide_id, dscore, seq_type, cutsite

    return best_edit_score, best_edit_df, best_guide, best_guide_id, best_guide_dscore, best_seq_type, best_cutsite


def predict_level(seed_df, level_num, wt_seq, guide_df, task_c):
    """
    generate predictions for all sequences in level
    """
    
    # extract information from seed_df to generate predictions
    seq_list = seed_df['Genotype'].values
    old_cutsite_list = seed_df['Cutsite'].values
    cutting_score_list = seed_df['REtarget_Score_Contribution'].values * (config.LEVEL_WEIGHT_FACTOR**(level_num - 1))
    seq_type_list = seed_df['Seq_Type'].values
    guide_id_list = seed_df['Edit_ID'].values
    target_origin_list = seed_df['Guide_ID'].values

    # initialize empty data frame for edit outcomes
    new_edit_df = pd.DataFrame(columns=['Category', 'Genotype position', 'Inserted Bases', 'Length', 'Predicted frequency', 'Genotype',\
                                        'REtarget_Score_Contribution', 'Seq_Type', 'Cutsite', 'Edit_ID', 'Guide_ID', 'Level'])
    new_guide_df = pd.DataFrame(columns=['Guide_ID', 'Guide_RNA', 'Level', 'REtarget_Score', 'MGA_REtarget_Score', 'Cutting_Score', 'Editing_Score', 'Doench2014_Score',\
                                        'DoenchCFD_MaxOT', 'DoenchCFD_SpecScore', 'Hsu2013_Score', 'Number_Off_Targets', 'Mismatches_Closest_OT', 'Number_Closest_OT',\
                                        'Guide_Orientation', 'PAM', 'Target', 'Target_Type', 'Cut_Position', 'Target_Origin', 'Multi'])

    # find PAMs & run FlashFry to check off-targets
    guide_dict_fwd, pam_dict_fwd, guide_dict_rev, pam_dict_rev, flash_fry_results = find_pam(seq_list, old_cutsite_list, guide_id_list, task_c)

    # loop over new seeds
    for seq, cutting_score, seq_type, guide_id, target_origin in zip(seq_list, cutting_score_list, seq_type_list, guide_id_list, target_origin_list):
        
        # sort out guides that do not meet selection criteria (score thresholds)
        pam_list_fwd = select_pam(guide_dict_fwd[guide_id], pam_dict_fwd[guide_id], flash_fry_results, wt_seq)
        pam_list_rev = select_pam(guide_dict_rev[guide_id], pam_dict_rev[guide_id], flash_fry_results, wt_seq)
    
        # initialize edit score for single guide
        best_edit_score = 0
        best_edit_df, best_guide, best_guide_id, best_guide_dscore, best_seq_type, best_cutsite = pd.DataFrame(), None, None, None, None, None

        # run prediction for single guides
        if len(pam_list_fwd) > 0:
            best_edit_score, best_edit_df, best_guide, best_guide_id, best_guide_dscore, best_seq_type, best_cutsite = predict_single_guide(seq, cutting_score, seq_type, guide_id, level_num, pam_list_fwd, best_edit_score, best_edit_df, best_guide, best_guide_id, best_guide_dscore, best_seq_type, best_cutsite)
        if len(pam_list_rev) > 0:
            rc_seq = get_rc_seq(seq)
            if seq_type == 'fwd':
                rc_seq_type = 'rev'
            else:
                rc_seq_type = 'fwd'
            best_edit_score, best_edit_df, best_guide, best_guide_id, best_guide_dscore, best_seq_type, best_cutsite = predict_single_guide(rc_seq, cutting_score, rc_seq_type, guide_id, level_num, pam_list_rev, best_edit_score, best_edit_df, best_guide, best_guide_id, best_guide_dscore, best_seq_type, best_cutsite)
        
        # create output data frame for level
        if best_edit_score > 0:
            new_edit_df = new_edit_df.append(best_edit_df)
            if (best_seq_type == 'fwd' and seq_type == 'fwd') or (best_seq_type == 'rev' and seq_type == 'rev'):
                best_seq = seq
            else:
                best_seq = rc_seq
            new_guide_df = append_guide_df(new_guide_df, best_guide_id, best_guide, cutting_score, best_guide_dscore, best_seq, best_seq_type, best_cutsite, target_origin, best_edit_df, guide_df, task_c)

    return new_edit_df, new_guide_df


def check_stop_criteria(guide_df, level_num):
    """
    check if stopping criteria met
    """

    # df with guides of last level
    stop_df = guide_df[guide_df['Level'] == level_num]

    # test if overall optimum
    if stop_df.empty:
        logging.info('No gRNAs in current level satisfy selection criteria. Stopping recursive optimization.')
        return True
    
    # test if stopping criterion met
    if np.sum(stop_df['REtarget_Score'].values) < config.LEVEL_FREQ_CUT:
        logging.info('Level REtarget score too low: ' + str(np.sum(stop_df['REtarget_Score'].values)))
        return True
    if level_num == 1 and (np.sum(stop_df['REtarget_Score'].values) < config.L1_CUT):
        logging.info('Level 1 REtarget score too low: ' + str(np.sum(stop_df['REtarget_Score'].values)))
        return True
    if level_num == 2 and (np.sum(stop_df['REtarget_Score'].values) < config.L2_CUT):
        logging.info('Level 2 REtarget score too low: ' + str(np.sum(stop_df['REtarget_Score'].values)))
        return True

    # test if max number of levels reached
    if level_num >= config.MAX_LEVEL_NUM:
        logging.info('Maximal number of levels reached. Stopping recursive optimization.')
        return True

    # stopping criteria not yet met
    return False


def select_guides(new_edit_df, new_guide_df, new_guide_num, edit_df, guide_df, guide_num):
    """
    select guides from existing and new guides with highest REtarget scores, respectively MGA REtarget scores (MGA = multiple guide appearance: taking into consideration that
    same guide could appear multiple times if editing restores sequences previously present)
    """

    # (a): if number of existing and new guides smaller than maximal guide number: keep all new guides
    if guide_num + new_guide_num <= config.MAX_GUIDE_NUM:
        # update existing guide, edit df
        guide_df = guide_df.append(new_guide_df)
        edit_df = edit_df.append(new_edit_df)
        # update new guide, edit df
        new_guide_df, new_edit_df = pd.DataFrame(), pd.DataFrame()
        # update guide number
        guide_num += new_guide_num
    
    # (b): if (a) not met but number of existing guides smaller than maximal guide number: keep new guides with highest REtarget scores
    elif guide_num < config.MAX_GUIDE_NUM:
        highest_guide_df = new_guide_df.head(config.MAX_GUIDE_NUM - guide_num)
        highest_guide_ids = highest_guide_df['Guide_ID'].values
        guide_df = guide_df.append(highest_guide_df)
        edit_df = edit_df.append(new_edit_df[new_edit_df['Guide_ID'].isin(highest_guide_ids)])
        new_guide_df = new_guide_df[~new_guide_df['Guide_ID'].isin(highest_guide_ids)]
        new_edit_df = new_edit_df[~new_edit_df['Guide_ID'].isin(highest_guide_ids)]
        # special case: guide that already occurred is in new_guide_df
        if 1 in new_guide_df['Multi'].values:
            highest_guide_df = new_guide_df[new_guide_df['Multi'] == 1]
            highest_guide_ids = highest_guide_df['Guide_ID'].values
            guide_df = guide_df.append(highest_guide_df)
            edit_df = edit_df.append(new_edit_df[new_edit_df['Guide_ID'].isin(highest_guide_ids)])
            new_guide_df = new_guide_df[~new_guide_df['Guide_ID'].isin(highest_guide_ids)]
            new_edit_df = new_edit_df[~new_edit_df['Guide_ID'].isin(highest_guide_ids)]
        guide_num = config.MAX_GUIDE_NUM
    
    # if maximal guide number reached and new guides left: keep guides with highest REtarget scores from new and existing guides
    if guide_num == config.MAX_GUIDE_NUM and not new_guide_df.empty:
        # initialize dict with guide with lowest REtarget score: 0 - guide, 1 - score, 2 - if guide occurrs multiple times
        old_min_score = {0: '', 1: 0, 2: 0}
        # find guide index of guide with lowest REtarget score
        index_min_score = np.argmin(guide_df['MGA_REtarget_Score'].values)
        # store guide with lowest REtarget score in dict
        min_score = {0: guide_df['Guide_ID'].values[index_min_score], 1: guide_df['MGA_REtarget_Score'].values[index_min_score], 2: guide_df['Multi'].values[index_min_score]}
        # remove existing guide and add new to df as long as new guide in each iteration with minimal REtarget score
        while old_min_score[0] != min_score[0]:
            old_min_score = min_score
            new_guide_df = new_guide_df[new_guide_df['MGA_REtarget_Score'] > min_score[1]]
            if not new_guide_df.empty:
                # append new guide
                highest_guide_ids = new_guide_df['Guide_ID'].values[0]
                guide_df = guide_df.append(new_guide_df.head(1))
                edit_df = edit_df.append(new_edit_df[new_edit_df['Guide_ID'] == highest_guide_ids])
                # delete old guide
                guide_df = guide_df[guide_df['Guide_ID'] != min_score[0]]
                edit_df = edit_df[edit_df['Guide_ID'] != min_score[0]]
                # special deletion case: multiple occurrances of same guide
                if min_score[2] == 1:
                    guide_df = guide_df[guide_df['MGA_REtarget_Score'] != min_score[1]]
                    edit_df = edit_df[edit_df['Guide_ID'].isin(guide_df['Guide_ID'].values)]
                # delete added guide from new df
                new_guide_df = new_guide_df[new_guide_df['Guide_ID'] != highest_guide_ids]
                new_edit_df = new_edit_df[new_edit_df['Guide_ID'] != highest_guide_ids]
                # update min guide score
                index_min_score = np.argmin(guide_df['MGA_REtarget_Score'].values)
                min_score = {0: guide_df['Guide_ID'].values[index_min_score], 1: guide_df['MGA_REtarget_Score'].values[index_min_score], 2: guide_df['Multi'].values[index_min_score]}
    
    return edit_df, guide_df, guide_num


def final_analysis(edit_df, guide_df, level_num, wt_seq, task_c):
    """
    predict gRNAs with maximal on-target efficiency for last level after stopping criteria met (recalculation of last level)
    """

    # get seeds for new level
    seed_df = edit_df[edit_df['Level'] == level_num]
    
    # update level number
    level_num += 1

    # extract information from seed_df to generate predictions
    seq_list = seed_df['Genotype'].values
    old_cutsite_list = seed_df['Cutsite'].values
    cutting_score_list = seed_df['REtarget_Score_Contribution'].values * (config.LEVEL_WEIGHT_FACTOR**(level_num - 1))
    seq_type_list = seed_df['Seq_Type'].values
    guide_id_list = seed_df['Edit_ID'].values
    target_origin_list = seed_df['Guide_ID'].values
    
    # find PAMs & run FlashFry to check off-targets
    guide_dict_fwd, pam_dict_fwd, guide_dict_rev, pam_dict_rev, flash_fry_results = find_pam(seq_list, old_cutsite_list, guide_id_list, task_c)

    # loop over new seeds
    for seq, cutting_score, seq_type, guide_id, target_origin in zip(seq_list, cutting_score_list, seq_type_list, guide_id_list, target_origin_list):
        
        # sort out guides that do not meet selection criteria (score thresholds)
        pam_list_fwd = select_pam(guide_dict_fwd[guide_id], pam_dict_fwd[guide_id], flash_fry_results, wt_seq)
        pam_list_rev = select_pam(guide_dict_rev[guide_id], pam_dict_rev[guide_id], flash_fry_results, wt_seq)
    
        # initialize edit score for single guide
        best_edit_score = 0
        best_edit_df, best_guide, best_guide_id, best_guide_dscore, best_seq_type, best_cutsite = pd.DataFrame(), None, None, 0, None, None

        # run prediction for single guides
        if len(pam_list_fwd) > 0:
            best_edit_score, best_edit_df, best_guide, best_guide_id, best_guide_dscore, best_seq_type, best_cutsite = predict_single_guide(seq, cutting_score, seq_type, guide_id, level_num, [pam_list_fwd[0]], best_edit_score, best_edit_df, best_guide, best_guide_id, best_guide_dscore, best_seq_type, best_cutsite)
        if len(pam_list_rev) > 0 and best_guide_dscore < pam_list_rev[0][0]:
            rc_seq = get_rc_seq(seq)
            if seq_type == 'fwd':
                rc_seq_type = 'rev'
            else:
                rc_seq_type = 'fwd'
            best_edit_score, best_edit_df, best_guide, best_guide_id, best_guide_dscore, best_seq_type, best_cutsite = predict_single_guide(rc_seq, cutting_score, rc_seq_type, guide_id, level_num, [pam_list_rev[0]], best_edit_score, best_edit_df, best_guide, best_guide_id, best_guide_dscore, best_seq_type, best_cutsite)
        
        if not best_edit_df.empty:
            edit_df = edit_df.append(best_edit_df)
            if (best_seq_type == 'fwd' and seq_type == 'fwd') or (best_seq_type == 'rev' and seq_type == 'rev'):
                best_seq = seq
            else:
                best_seq = rc_seq
            guide_df = append_guide_df(guide_df, best_guide_id, best_guide, cutting_score, best_guide_dscore, best_seq, best_seq_type, best_cutsite, target_origin, best_edit_df, guide_df, task_c)

    return edit_df, guide_df, level_num


def recursive_hdr_optimization(edit_df, guide_df, wt_seq, level_num=1, guide_num=1, task_c=1):
    """
    keep optimal guides and sort out remaining
    """

    # 1. check if stopping criteria met
    if check_stop_criteria(guide_df, level_num) == True:
        # sort out last level that did not fulfill selection criteria
        guide_df = guide_df[guide_df['Level'] < level_num]
        edit_df = edit_df[edit_df['Level'] < level_num]

        # redo guide selection for last level: select guides according to guide scores, not potential for retargeting
        if config.FINAL_BY_D2014ON is True:
            logging.info('Recalculation of last level. Selecting gRNAs with maximal on-target efficiency.')
            config.SINGLE_FREQ_CUT = 0
            config.D2014ON_MIN = 0
            edit_df, guide_df, level_num = final_analysis(edit_df, guide_df, level_num - 1, wt_seq, task_c)

        # terminate recursive search
        return edit_df, guide_df, (level_num - 1) # level num only considers Recursive Editing layers for decision if results will be saved or not

    # 2. predictions for new level
    # get seeds for new level
    seed_df = edit_df[edit_df['Level'] == level_num]
    # update level number
    level_num += 1
    logging.info('Recursive optimization level ' + str(level_num))
    # predict guides for new level
    new_edit_df, new_guide_df = predict_level(seed_df, level_num, wt_seq, guide_df, task_c)

    # 3. selection of guides with highest REtarget score
    if not new_guide_df.empty:
        # sort new guides by REtarget scores
        new_guide_df = new_guide_df.sort_values(by='MGA_REtarget_Score', ascending=False)
        # get number of new guides
        new_guide_num = len(new_guide_df.loc[new_guide_df['Multi'] != 1, 'Guide_RNA'].values)
        # select guides
        edit_df, guide_df, guide_num = select_guides(new_edit_df, new_guide_df, new_guide_num, edit_df, guide_df, guide_num)
        
    # run function optimize_hdr recursively
    return recursive_hdr_optimization(edit_df, guide_df, wt_seq, level_num, guide_num, task_c)


def initial_analysis(seq, init_cutsite, seq_type, init_guide, dscore, task_c):
    """
    initial analysis step for given cutsite/guide
    """
    if seq_type == 'rev':
        seq = get_rc_seq(seq)
        init_cutsite = get_rc_cut(seq, init_cutsite)

    edit_df = run_prediction(config.MODEL, seq, init_cutsite)
    
    if edit_df.empty:
        return edit_df, None
    else:
        edit_df = append_edit_df(edit_df, 1, seq_type, '1', 1, init_cutsite)
        guide_df = pd.DataFrame(columns=['Guide_ID', 'Guide_RNA', 'Level', 'REtarget_Score', 'MGA_REtarget_Score', 'Cutting_Score', 'Editing_Score', 'Doench2014_Score',\
                                        'DoenchCFD_MaxOT', 'DoenchCFD_SpecScore', 'Hsu2013_Score', 'Number_Off_Targets', 'Mismatches_Closest_OT', 'Number_Closest_OT',\
                                        'Guide_Orientation', 'PAM', 'Target', 'Target_Type', 'Cut_Position', 'Target_Origin', 'Multi'])
        guide_df = append_guide_df(guide_df, '1', init_guide, 1, dscore, seq, seq_type, init_cutsite, 'wt', edit_df, guide_df, task_c)

        return edit_df, guide_df


def run_retarget(seq, cutsite, seq_type, res_df=None, res_id=1, start_guide=None, store_res=False, store_name='', celltype=config.CELLTYPE_MODEL, task_c=1):
    """
    run complete hdr optimization script for single locus
    """

    logging.info('==> Running REtarget for sample with ID ' + str(res_id))

    # get fwd and reverse WT sequence to check if guides cut WT
    wt_seq = seq + ' ' + get_rc_seq(seq)

    # initialize selected model
    if config.MODEL == 1:
        init_indelphi(celltype)
    model_dict = {0: 'Lindel', 1: 'inDelphi'}
    logging.info('Selected indel prediction tool: ' + str(model_dict[config.MODEL]))

    # get initial guide from seq, cutsite and check if identical to initial guide specified by user
    init_guide = get_guide(seq, cutsite, seq_type)
    if start_guide != None:
        if start_guide != init_guide:
            logging.info('Error: Specified initial gRNA {} does not match gRNA at specified cutsite {}. Please select correct initial gRNA.'.format(start_guide, init_guide))
            return None, None, None, 'Specified initial gRNA {} does not match gRNA at specified cutsite {}.'.format(start_guide, init_guide)
    
    # check if initial guide meets search criteria
    if config.CHECK_OFFT_GENOME is True or config.CHECK_OFFT_SCORES is True:
        # get pam
        init_pam = get_pam(seq, cutsite, seq_type)

        logging.info('Evaluating gRNA ' + str(init_guide) + '. PAM: ' + str(init_pam))

        # run flashfry for all potential guides in level
        run_flash_fry([str(init_guide)+str(init_pam)], config.FF_DATABASE_NAME, task_c, use_taskset=config.USE_TASKSET_JKD)
        # get flash fry results
        flash_fry_results = pd.read_csv(os.path.join(FLASHFRY_PATH, 'level_'+str(task_c)+'.output.scored'), delimiter='\t', index_col=None)
    else:
        logging.info('Evaluating gRNA ' + str(init_guide))
        flash_fry_results = None
    do_analysis, dscore, message = check_guide(init_guide, guide_scores=flash_fry_results, wt_seq=wt_seq, check_offt_prox=False, l1=True)
    if do_analysis is False:
        logging.info('Optimization stopped for sample with ID ' + str(res_id) + '. ' + str(message))
        return None, None, None, message

    # initial analysis for given first guide
    edit_df, guide_df = initial_analysis(str(seq), int(cutsite), seq_type, init_guide, dscore, task_c)
    # if editing df is empty: no editing outcome matches criteria for further search --> terminate analysis
    if edit_df.empty:
        logging.info('Optimization stopped for sample with ID ' + str(res_id) + '. No retargetable editing outcome resulting from initial gRNA editing. Consider adjusting search parameters.')
        return None, None, None, 'No editing outcome resulting from first gRNA suitable for Recursive Editing according to specified search parameters.'

    # run hdr optimization
    edit_df, guide_df, level_num = recursive_hdr_optimization(edit_df, guide_df, wt_seq, level_num=1, guide_num=1, task_c=task_c)

    # summarize results in res_df
    if res_df == None:
        res_df = pd.DataFrame(columns=['Res_ID', 'Level_Number', 'Guide_Number', 'Guide_RNAs', 'Total_REtarget_Score', 'REtarget_Score_L1', 'REtarget_Score_L2'])
    res_df = append_res_df(res_df, res_id, level_num, guide_df)
    if store_res == True and level_num >= config.L_MIN_SAVE:
        # create output directory
        store_path = os.path.join(BASE_PATH, 'data', 'output_files', store_name, str(level_num))
        os.makedirs(store_path, exist_ok=True)
        logging.info('ID ' + str(res_id) + ' suitable for Recursive Editing according to specified parameters.')
        logging.info('Saving results to path:' + str(store_path))

        # save results
        edit_df.to_csv(os.path.join(store_path, str(res_id) + '_edit_df.csv'), mode='w+', index=False)
        guide_df.to_csv(os.path.join(store_path, str(res_id) + '_guide_df.csv'), mode='w+', index=False)
        res_df.to_csv(os.path.join(os.path.dirname(store_path), 'res_df.csv'), mode='a+', index=False, header=None)
    else:
        logging.info('Number of resulting levels smaller than config.L_MIN_SAVE: results for {} not saved. Consider adjusting search parameters.'.format(str(res_id)))
        level_num = 'Number of resulting levels smaller than config.L_MIN_SAVE: results for {} not saved'.format(str(res_id))
    
    # delete temporary files if created
    if (config.CHECK_OFFT_GENOME == True) or (config.CHECK_OFFT_SCORES == True):
        try:
            # remove temporary files
            os.remove(os.path.join(FLASHFRY_PATH, 'level_'+str(task_c)+'.fasta'))
            os.remove(os.path.join(FLASHFRY_PATH, 'level_'+str(task_c)+'.output'))
            os.remove(os.path.join(FLASHFRY_PATH, 'level_'+str(task_c)+'.output.scored'))
        except Exception as e:
            pass

    return edit_df, guide_df, res_df, level_num



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='plot')
    parser.add_argument('-s', '--seq', required=False, type=str, help='Input sequence. Specify forward strand.', default='gaatttagtctcccagcaggagaaataagagtaaataaatgtccagtggaatccacagcagacccaccctcaccttcagttttattgttttgctccaaacataactctgctgcttccactgctctggggctggtaaaaatgagtcccccgtaatcttcaggatgagaaagctgcacaccaaaaagcaataaagacatttt')
    parser.add_argument('-c', '--cutsite', required=False, type=int, help='Initial cut site', default='100')
    parser.add_argument('-o', '--orientation', required=False, type=str, help='gRNA orientation: "fwd" if gRNA and seq have same orientation, "rev" if gRNA and seq have different orientation', default='rev')
    parser.add_argument('-n', '--name', required=False, type=str, help='Name of output folder', default='REtarget')
    parser.add_argument('-f', '--database', required=False, type=str, help='Name of FleshFry database', default='chr22_cas9ngg_database')
    parser.add_argument('-t', '--template', required=False, type=str, help='Sequence of HDR template', default='CGGGGGACTCATTTTTACCAGCCCCAGAGCAGTGGAAGCAGCAGAGTTATgatGTTTGGAGCAAAACAATAAAACTGAAGGTGAGGGTGGGTCTGCTGTGGAT')
    parser.add_argument('-b', '--blacklist', required=False, type=list, help='gRNA Blacklist: List of gRNAs that should not be used by REtarget', default=[])
    args = parser.parse_args()

    if str(args.database) != '':
        config.FF_DATABASE_NAME = str(args.database)
    if str(args.template) != '':
        config.HDR_TEMPLATE = str(args.template).upper()
    if len(list(args.blacklist)) > 0:
        config.BLACKLIST = [str(g).upper() for g in list(args.blacklist)]
    
    edit_df, guide_df, res_df, level_num = run_retarget(str(args.seq).upper(), int(args.cutsite), str(args.orientation), store_res=True, store_name=str(args.name))
