import pandas as pd
from autoqchem.gaussian_log_extractor import *
import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
import statistics
import os

#generate DFT descriptors for rings as an example

#generate path list for 16 rings
HomeDir = os.getcwd()
base_path = HomeDir+'/log_files/ring/'
ring_number = np.arange(17)
ring_list = []
for rn in ring_number:
    if rn!=0:
        ring = base_path+str(rn)+'.log'
        ring_list.append(ring)

'''
#generate path list for 36 linkers
linker_number_1 = np.arange(13).tolist()
linker_number_1.remove(0)
linker_number_2 = np.arange(4).tolist()
linker_number_2.remove(0)
linker_list = []
for lk_2 in linker_number_2:
    for lk_1 in linker_number_1:
        linker = base_path+str(lk_1)+str(lk_2)+'.log'
        linker_list.append(linker)
'''

#generate DFT descriptors
features_list = []
for path in ring_list:
    feature = []
    descriptors = gaussian_log_extractor(path).get_descriptors()
    for des_key in descriptors['descriptors'].keys():
        feature.append(descriptors['descriptors'][des_key])
    for ades_key in ['VBur', 'Mulliken_charge', 'APT_charge']:
        full_list = descriptors['atom_descriptors'][ades_key]
        fl_list = []
        for f in full_list:
            fl_list.append(float(f))
        feature.append(statistics.mean(fl_list))
        feature.append(max(fl_list))
        feature.append(min(fl_list))
    for modes_key in descriptors['modes'].keys():
        full_list = descriptors['modes'][modes_key]
        fl_list = []
        for f in full_list:
            fl_list.append(float(f))
        feature.append(statistics.mean(fl_list))
        feature.append(max(fl_list))
        feature.append(min(fl_list))
    features_list.append(feature)

data = pd.DataFrame(features_list)
name_list = list(descriptors['descriptors'].keys())
name_list = name_list + ['VBur_m', 'VBur_max', 'VBur_min', 'Mulliken_charge_m', 'Mulliken_charge_max', 'Mulliken_charge_min',
                            'APT_charge_m', 'APT_charge_max', 'APT_charge_min'] + ['Frequencies_m', 'Frequencies_max', 'Frequencies_min',
                            'Red masses_m', 'Red masses_max', 'Red masses_min',
                            'Frc consts_m', 'Frc consts_max', 'Frc consts_min',
                            'IR Inten_m', 'IR Inten_max', 'IR Inten_min']
data.columns = name_list
data.to_csv("DFT_rings.csv")
