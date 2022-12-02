from data_loader import Al
from edbo.utils import Data
from edbo.bro import BO
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
from edbo.models import GP_Model


#load data
desc = 'dft'
reaction = Al(YR3=desc,
                R1R2=desc,
                Temp=desc,
                Ratio=desc)
'''
reaction = Al(Y='ohe',
                R3='ohe',
                R1R2='ohe',
                Temp='ohe',
                Ratio='ohe')
'''



reaction = Data(reaction)
reaction.clean()
reaction.standardize(scaler='minmax')
reaction.uncorrelated(threshold=0.95)
reaction.drop(['entry', 'YR3_code', 'R1R2_code'])
#reaction.drop(['entry', 'Y_code', 'R3_code', 'R1R2_code'])

reaction.data = reaction.data[:56]
reaction.data = reaction.data.sample(frac=1).reset_index(drop=True)

#get Pm or Pr values as the target vector
yd = np.array(reaction.data['yield'].tolist())[:56]

#5-fold cross validation to obtain the mean absolute errors for each descriptor
MAE_list = []
for i in range(5):
    bo = BO(exindex=reaction.data,
            model=GP_Model,
            domain=reaction.data.drop('yield', axis=1),
            results=reaction.data.iloc[np.array(list(set(np.arange(56).tolist()).difference(set(np.arange(11*i,11*(i+1)).tolist()))))],
            init_method='external',
            batch_size=0,
            acquisition_function='EI',
            fast_comp=True)

    bo.run()
    abs_diff = abs(bo.obj.scaler.unstandardize(bo.model.predict(bo.obj.domain))-reaction.data['yield'].values)
    test_list = np.arange(11*i,11*(i+1)).tolist()
    MAE = abs_diff[test_list].mean()
    MAE_list.append(MAE)

print(np.array(MAE_list).mean(),np.array(MAE_list).std())



#randomly select one fold and draw parity plot between model-predicted and observed Pm or Pr values
slice = random.sample(np.arange(56).tolist(), 45)
c_slice = list(set(np.arange(56).tolist()).difference(set(np.array(slice))))

bo = BO(exindex=reaction.data,
        model=GP_Model,
        domain=reaction.data.drop('yield', axis=1),
        results=reaction.data.iloc[np.array(slice)],
        init_method='external',
        batch_size=0,
        acquisition_function='EI',
        fast_comp=True)
bo.run()

y_train = yd[slice]
y_test = yd[c_slice]
predict_train = (bo.obj.scaler.unstandardize(bo.model.predict(bo.obj.domain))[slice],(np.sqrt(bo.model.variance(bo.obj.domain))[slice] * bo.obj.scaler.std)**2)
predict_test = (bo.obj.scaler.unstandardize(bo.model.predict(bo.obj.domain))[c_slice],(np.sqrt(bo.model.variance(bo.obj.domain))[c_slice] * bo.obj.scaler.std)**2)

x_line = np.linspace(0,1,100)
y_line = x_line
plt.xlim = ([0,1])
plt.ylim = ([0,1])
plt.plot(x_line, y_line, '--', color = 'blue')
plt.scatter(y_train, predict_train[0], s=50, label = 'training',color = 'green', alpha = 0.7)
plt.scatter(y_test, predict_test[0], s=50, label = 'test',color = 'purple', alpha = 0.7)
plt.errorbar(y_train, predict_train[0], yerr = predict_train[1]**0.5,fmt = 'o', elinewidth=2,capsize = 5,color = 'green', alpha = 0.7)
plt.errorbar(y_test, predict_test[0], yerr = predict_test[1]**0.5,fmt = 'o', elinewidth=2,capsize = 5,color = 'purple', alpha = 0.7)
plt.xlabel('Observed Pm')
plt.ylabel('Predicted Pm')
plt.legend()
plt.savefig('dft_regression_parity.png', format='png',transparent = True)
plt.show()
