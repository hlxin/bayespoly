# -*- coding: utf-8 -*-
"""
Test data
"""

# Imports

import pandas as pd
from edbo.feature_utils import build_experiment_index

# Build data sets from indices


def Al(YR3='dft', R1R2='dft', Temp='dft', Ratio='dft'):
    """
    Load Al data with different features.
    """

    # SMILES index

    index = pd.read_csv('../BO/data/Al/experiment_index_86.csv')

    # Electrophile features

    if YR3 == 'mordred':
        YR3_features = pd.read_csv('../BO/data/Al/YR3_mordred.csv')
    elif YR3 == 'dft':
        YR3_features = pd.read_csv('../BO/data/Al/YR3_dft.csv')
    elif YR3 == 'cm':
        YR3_features = pd.read_csv('../BO/data/Al/YR3_cm.csv')
    elif YR3 == 'ei':
        YR3_features = pd.read_csv('../BO/data/Al/YR3_ei.csv')

    # Nucleophile features

    if R1R2 == 'mordred':
        R1R2_features = pd.read_csv('../BO/data/Al/R1R2_mordred.csv')
    elif R1R2 == 'dft':
        R1R2_features = pd.read_csv('../BO/data/Al/R1R2_dft.csv')
    elif R1R2 == 'cm':
        R1R2_features = pd.read_csv('../BO/data/Al/R1R2_cm.csv')
    elif R1R2 == 'ei':
        R1R2_features = pd.read_csv('../BO/data/Al/R1R2_ei.csv')


    if Temp == 'mordred':
        Temp_features = pd.read_csv('../BO/data/Al/Temp_mordred.csv')
    elif Temp == 'dft':
        Temp_features = pd.read_csv('../BO/data/Al/Temp_dft.csv')
    elif Temp == 'cm':
        Temp_features = pd.read_csv('../BO/data/Al/Temp_cm.csv')
    elif Temp == 'ei':
        Temp_features = pd.read_csv('../BO/data/Al/Temp_ei.csv')

    # Base features

    if Ratio == 'mordred':
        Ratio_features = pd.read_csv('../BO/data/Al/Ratio_mordred.csv')
    elif Ratio == 'dft':
        Ratio_features = pd.read_csv('../BO/data/Al/Ratio_dft.csv')
    elif Ratio == 'cm':
        Ratio_features = pd.read_csv('../BO/data/Al/Ratio_cm.csv')
    elif Ratio == 'ei':
        Ratio_features = pd.read_csv('../BO/data/Al/Ratio_ei.csv')


    # Build the descriptor set

    index_list = [index['R1R2_code'],
                  index['YR3_code'],
                  index['Temp_code'],
                  index['Ratio_code']]

    lookup_table_list = [R1R2_features,
                         YR3_features,
                         Temp_features,
                         Ratio_features]

    lookup_list = ['R1R2_code',
                   'YR3_code',
                   'Temp_code',
                   'Ratio_code']

    experiment_index = build_experiment_index(index['entry'],
                                              index_list,
                                              lookup_table_list,
                                              lookup_list)

    experiment_index['yield'] = index['yield']

    return experiment_index
