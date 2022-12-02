from data_loader import Al
from edbo.utils import Data
from edbo.bro import BO
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
from edbo.models import GP_Model, RF_Model, Random
from edbo.plot_utils import average_convergence, plot_avg_convergence
import math
from random import choice


def simulate(data,                 # Data container (vide supra)
             acq_func,             # Acquisition function: 'EI', 'PI', 'UCB', 'TS', ...
             model=GP_Model,       # Surrogate model: GP_Model or RF_Model
             init_method='rand',   # Initialization method: 'rand', 'pam', 'kmeans'
             batch_size=5,         # Parallel acquisition: int
             iterations=9,         # Number of iterations to run each simulation for
             average_of=10,        # Number of simulations to average
             export_path=None,     # Export the simulation results to a CSV file
             plot=False):          # Plot the average and standard deviation averaged convergence
    """
    Simulation function which averages BO runs with given optimization paramters.
    """

    # Average N optimizations with different random initializations
    results = []
    for i in range(average_of):

        # Use random function if the acquisition function is random selection
        if acq_func == 'rand':
            bo = BO(exindex=data.data,
                    domain=data.data.drop('yield', axis=1),
                    init_method=init_method,
                    model=Random,
                    batch_size=batch_size,
                    acquisition_function=acq_func)
        # Otherwise use specified model
        else:
            bo = BO(exindex=data.data,
                    domain=data.data.drop('yield', axis=1),
                    model=model,
                    init_method=init_method,
                    batch_size=batch_size,
                    acquisition_function=acq_func,
                    fast_comp=True)

        # Simulate
        bo.init_seq.visualize = False
        bo.simulate(iterations=iterations, seed=i)

        # Append results to record
        results.append(bo.obj.results_input()['yield'].values)

    # Save the results to a CSV file
    results = pd.DataFrame(results)
    if export_path != None:
        results.to_csv(export_path)

    # Average performance
    index, mean, std = average_convergence(results, batch_size)

    # Plot
    if plot:
        plot_avg_convergence(results, batch_size)

    return results, mean, std


#load data
reaction = Al(YR3='dft',
                   R1R2='dft',
                   Temp='dft',
                   Ratio='dft')

reaction = Data(reaction)
reaction.clean()
reaction.standardize(scaler='minmax')
reaction.uncorrelated(threshold=0.95)
reaction.drop(['entry', 'YR3_code', 'R1R2_code', 'Temp_code', 'Ratio_code'])

reaction.data = reaction.data[:56]


#BO on 56 literature data points, aiming to efficiently identify Pm or Pr maximum

results, mean, std = simulate(reaction,           # Data from reaction
                              'EI',                 # Probability of improvement
                              model=GP_Model,       # Gaussian process model
                              init_method='rand',   # Random selection initialization
                              batch_size=3,        # Choose 3 experiments at a time
                              iterations=11,         # Run for 12 iterations including initialization
                              average_of=10,        # Run with 10 different random seeds
                              export_path=None,     # Don't save the results
                              plot=True)

#Organize generated dataset and generate data for optimization curve

doc_name = 'Pm_BO_56.csv'
results.to_csv(doc_name)
N_pd = pd.read_csv(doc_name)
df = pd.DataFrame(N_pd)
mx_dic = {}
for i in np.arange(10):
    mx_dic[str(i)]=[]
    mx_list = mx_dic[str(i)]
    N_list = df.values.tolist()[i][1:]
    mx = max(N_list[:3])
    mx_list.append(mx)
    for j in range(11):
        if max(N_list[3*(j+1):3*(j+2)]) > mx:
            mx = max(N_list[3*(j+1):3*(j+2)])
        mx_list.append(mx)

mean_list = []
std_list = []
for i in range(12):
    h_list = [mx_dic[str(j)][i] for j in range(10)]
    mean_list.append(np.mean(h_list))
    std_list.append(np.std(h_list))
Pm_BO_mean_list = mean_list
Pm_BO_std_list = std_list


#Generate optimization curve for random search on the same dataset

mx_dic = {}
for i in range(10):
    mx_dic[str(i)]=[]
    mx_list = mx_dic[str(i)]
    p_list = reaction.data['yield'].tolist()
    slice = random.sample(p_list,3)
    mx = max(slice)
    mx_list.append(mx)
    for j in range(11):
        for k in range(3):
            p_list.remove(slice[k])
        slice = random.sample(p_list,3)
        if max(slice) > mx:
            mx = max(slice)
        mx_list.append(mx)

mean_list = []
std_list = []
for i in range(12):
    h_list = [mx_dic[str(j)][i] for j in range(10)]
    mean_list.append(np.mean(h_list))
    std_list.append(np.std(h_list))
Pm_RS_mean_list = mean_list
Pm_RS_std_list = std_list


#plot optimization curves of both BO and RS
X_list = np.linspace(1,12,12)
plt.plot(X_list, Pm_BO_mean_list, 'o-', color="red", label="BO",linewidth=2, markersize=5)
plt.plot(X_list, Pm_RS_mean_list, 'o-', color="blue", label="RS",linewidth=2, markersize=5)
plt.errorbar(X_list, Pm_BO_mean_list, yerr = Pm_BO_std_list,fmt = 'o', elinewidth=2,capsize = 5,color="red")
plt.errorbar(X_list, Pm_RS_mean_list, yerr = Pm_RS_std_list,fmt = 'o', elinewidth=2,capsize = 5,color="blue")
plt.xlabel('Iteration')
plt.ylabel('Max Pr observed')
plt.xticks([2,4,6,8,10,12])
plt.legend(frameon=False, loc = 'upper right', bbox_to_anchor=(0.35,0.25))
plt.title('Pm')
plt.savefig('Pm_BO_vs_RS_56.png', format='png',transparent = True)
plt.show()
