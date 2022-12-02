from data_loader import Al
from data_loader import Al_ohe
from edbo.utils import Data
from edbo.bro import BO
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from edbo.models import GP_Model
import os


#load descriptors and generate features for DFT, Mordred, CM, EI
reaction = Al(YR3='dft',
                R1R2='dft',
                Temp='dft',
                Ratio='dft')
'''
#load descriptors and generate features for OHE
reaction = Al_ohe(Y='ohe',
                    R3='ohe',
                    R1R2='ohe',
                    Temp='ohe',
                    Ratio='ohe')
'''

reaction = Data(reaction)
reaction.clean()
reaction.drop(['entry', 'YR3_code', 'R1R2_code','Temp_code','Ratio_code'])
#reaction.drop(['entry', 'Y_code', 'R3_code', 'R1R2_code','Temp_code','Ratio_code'])
reaction.standardize(scaler='minmax')
reaction.uncorrelated(threshold=0.95)

#the number of points desired to be proposed by the model
n_suggest = 40

#training dataset for the current iteration
n_train = 56

#Pm or Pr as the maximization target
Pmr = 'Pm'

bo = BO(exindex=reaction.data,
        model=GP_Model,
        domain=reaction.data.drop('yield', axis=1),
        results=reaction.data.iloc[list(range(n_train))],
        init_method='external',
        batch_size=n_suggest,
        acquisition_function='EI',
        fast_comp=True)

bo.run()

#save original dataset proposed by the model
orig_path = Pmr+'_'+str(n_suggest)+'_'+str(n_train)+'_orig.csv'
bo.acquisition_summary().to_csv(orig_path)

#organize the model-proposed dataset and clearly list all the important information including catalyst molecules, predicted mean & variance and synthetic scales
HomeDir = os.getcwd()
overall_pd = pd.read_csv(HomeDir + '/data/Al/experiment_index.csv')
overall_df = pd.DataFrame(overall_pd)
suggest_pd = pd.read_csv(orig_path)
suggest_dic = {}
suggest_dic['entry'] = np.array(suggest_pd.iloc[:,0])
suggest_dic['mean'] = np.array(suggest_pd.iloc[:,-2])
suggest_dic['variance'] = np.array(suggest_pd.iloc[:,-1])
suggest_df = pd.DataFrame(suggest_dic)
line_list = []
for entry_index in suggest_dic['entry'].tolist():
    line_list.append(overall_df.loc[overall_df['entry'] == entry_index])
line_df = pd.concat(line_list).reset_index()
cat_df = pd.concat([line_df,suggest_df],axis=1)

synthesis_pd = pd.read_csv(HomeDir+'/data/Al/synthesis_scale_576.csv')
cat_line_list = []
for cat in list(cat_df['catalyst']):
    cat_line_list.append(synthesis_pd.loc[synthesis_pd['catalyst'] == cat])
cat_line_df = pd.concat(cat_line_list).reset_index()
temp_df = pd.concat([cat_line_df,cat_df],axis=1)

final_dic = {'catalyst':list(temp_df['catalyst'].iloc[:,-1]),'syn_scale':np.array(temp_df['AB_synth'])+np.array(temp_df['YB_synth'])+2,'mean':np.array(temp_df['mean']),'variance':np.array(temp_df['variance'])}
df = pd.DataFrame(final_dic)
df.to_csv(Pmr+'_sugg_'+str(n_suggest)+'.csv')
