from data_loader import Al
from edbo.utils import Data
from edbo.bro import BO
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
from edbo.models import GP_Model
from gpytorch.priors import GammaPrior
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score, mean_absolute_error
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern
import shap

#SHAP analysis

#load data

desc = 'dft'
reaction = Al(YR3=desc,
                   R1R2=desc,
                   Temp=desc,
                   Ratio=desc)
reaction = Data(reaction)
reaction.clean()
reaction.drop(['entry', 'YR3_code', 'R1R2_code', 'Temp_code', 'Ratio_code'])
reaction.standardize(scaler='minmax')
reaction.uncorrelated(threshold=0.95)

#prepare faeture matrix and target vector for GPR training
data_df = reaction.data
data_y = data_df[['yield']]
data_y = data_y.values
data_y = data_y.reshape(86,)
data_df.drop(columns=['yield'], inplace = True)
data_X = data_df
kernel = Matern(nu=2.5)
gpr = GaussianProcessRegressor(kernel=kernel, random_state=0).fit(data_X, data_y)

#SHAP
explainer = shap.Explainer(gpr.predict,data_X)
shap_values_gpr = explainer(data_X)
shap.summary_plot(shap_values_gpr, data_X.values,show=False)
plt.savefig('shap_summarize_final.png',bbox_inches='tight',dpi=300)

#multivariant linear regression

#load the data
feature_pd = pd.read_csv('DFT_descriptors_86.csv')

#select the variables for linear regression
feature_columns = ['homo_energy_x','VBur_max_x','homo_energy_y','electronegativity_y','VBur_min_y']
X = np.array(pd.DataFrame(feature_pd, columns = feature_columns))
y = np.array(pd.DataFrame(feature_pd, columns = ['Pm']))

#fit the linear regression model
reg = LinearRegression().fit(X, y)

#print out the model performance and all the coefficients
print('R2 score is: ', r2_score(y,reg.predict(X)))
print('MAE is ', mean_absolute_error(y,reg.predict(X)))

title_string = 'Pm = '
for i in range(len(feature_columns)):
    if reg.coef_[0][i] > 0:
        title_string = title_string + '+' + str(reg.coef_[0][i]) + '*' + feature_columns[i] + ' '
    if reg.coef_[0][i] < 0:
        title_string = title_string + str(reg.coef_[0][i]) + '*' + feature_columns[i] + ' '
title_string = title_string + '+' + str(reg.intercept_[0])
print(title_string)


#plot the linear regression parity plot
x_line = np.linspace(0,1,100)
y_line = x_line
plt.xlim = ([0,1])
plt.ylim = ([0,1])
plt.plot(x_line, y_line, '--', color = 'tomato',linewidth = 3)
plt.scatter(reg.predict(X), y, s=100,marker='^', edgecolor='black',linewidth = 2,label = 'training',color = 'darkseagreen', alpha = 1)

plt.xlabel(r'Predicted $P_m$')
plt.ylabel(r'Observed $P_m$')
ax=plt.gca()
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.spines['right'].set_linewidth(2)
ax.spines['top'].set_linewidth(2)
plt.savefig('linear_regression.png', format='png',transparent = True)
plt.show()
