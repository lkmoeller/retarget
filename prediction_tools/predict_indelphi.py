#######################################################################################################################
#
# Generate and filter inDelphi predictions
#
# Author: Lukas Moeller - 11/2021
#
#######################################################################################################################
 


import warnings
import pandas as pd
import inDelphi
from recursive_editing import config



def warn(*args, **kwargs):
    """
    silence warnings that are caused by loading/running imported prediction models
    """
    pass


def init_indelphi(cell_type):
    # silence inDelphi warnings
    warnings.warn = warn
    # initialize inDelphi model
    inDelphi.init_model(celltype=cell_type)


def run_indelphi(seq, cutsite):
    """
    function to get predictions for edit outcomes by inDelphi
    """
    try:
        # generation of predictions
        tmp_df, stats = inDelphi.predict(str(seq).upper(), int(cutsite))
        tmp_df = inDelphi.add_genotype_column(tmp_df, stats)
        # sorting of predictions
        tmp_df = tmp_df[tmp_df.Genotype.notna()]
        tmp_df = tmp_df.sort_values(by='Predicted frequency', ascending=False)
        tmp_df = tmp_df.head(config.MAX_GUIDES)
        tmp_df = tmp_df[tmp_df['Predicted frequency'] > config.SINGLE_FREQ_CUT]
        return tmp_df
    except ValueError:
        return pd.DataFrame()
