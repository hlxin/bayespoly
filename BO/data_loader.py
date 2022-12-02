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

    index = pd.read_csv('data/Al/experiment_index.csv')

    # Electrophile features

    if YR3 == 'mordred':
        YR3_features = pd.read_csv('data/Al/YR3_mordred.csv')
    elif YR3 == 'dft':
        YR3_features = pd.read_csv('data/Al/YR3_dft.csv')
    elif YR3 == 'cm':
        YR3_features = pd.read_csv('data/Al/YR3_cm.csv')
    elif YR3 == 'ei':
        YR3_features = pd.read_csv('data/Al/YR3_ei.csv')

    # Nucleophile features

    if R1R2 == 'mordred':
        R1R2_features = pd.read_csv('data/Al/R1R2_mordred.csv')
    elif R1R2 == 'dft':
        R1R2_features = pd.read_csv('data/Al/R1R2_dft.csv')
    elif R1R2 == 'cm':
        R1R2_features = pd.read_csv('data/Al/R1R2_cm.csv')
    elif R1R2 == 'ei':
        R1R2_features = pd.read_csv('data/Al/R1R2_ei.csv')


    if Temp == 'mordred':
        Temp_features = pd.read_csv('data/Al/Temp_mordred.csv')
    elif Temp == 'dft':
        Temp_features = pd.read_csv('data/Al/Temp_dft.csv')
    elif Temp == 'cm':
        Temp_features = pd.read_csv('data/Al/Temp_cm.csv')
    elif Temp == 'ei':
        Temp_features = pd.read_csv('data/Al/Temp_ei.csv')

    # Base features

    if Ratio == 'mordred':
        Ratio_features = pd.read_csv('data/Al/Ratio_mordred.csv')
    elif Ratio == 'dft':
        Ratio_features = pd.read_csv('data/Al/Ratio_dft.csv')
    elif Ratio == 'cm':
        Ratio_features = pd.read_csv('data/Al/Ratio_cm.csv')
    elif Ratio == 'ei':
        Ratio_features = pd.read_csv('data/Al/Ratio_ei.csv')


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

def Al_ohe(Y='ohe', R3='ohe', R1R2='ohe', Temp='ohe', Ratio='ohe'):
    """
    Load Al data with different features.
    """

    # SMILES index

    index = pd.read_csv('data/Al/experiment_index_ohe.csv')

    # Electrophile features

    if Y == 'ohe':
        Y_features = pd.read_csv('data/Al/Y_ohe.csv')

    # Nucleophile features

    if R3 == 'ohe':
        R3_features = pd.read_csv('data/Al/R3_ohe.csv')

    # Base features

    if R1R2 == 'ohe':
        R1R2_features = pd.read_csv('data/Al/R1R2_ohe.csv')

    if Temp == 'ohe':
        Temp_features = pd.read_csv('data/Al/Temp_ohe.csv')

    # Base features

    if Ratio == 'ohe':
        Ratio_features = pd.read_csv('data/Al/Ratio_ohe.csv')


    # Build the descriptor set

    index_list = [index['Y_code'],
                  index['R3_code'],
                  index['R1R2_code'],
                  index['Temp_code'],
                  index['Ratio_code']]

    lookup_table_list = [Y_features,
                         R3_features,
                         R1R2_features,
                         Temp_features,
                         Ratio_features]

    lookup_list = ['Y_code',
                   'R3_code',
                   'R1R2_code',
                   'Temp_code',
                   'Ratio_code']

    experiment_index = build_experiment_index(index['entry'],
                                              index_list,
                                              lookup_table_list,
                                              lookup_list)

    experiment_index['yield'] = index['yield']

    return experiment_index
