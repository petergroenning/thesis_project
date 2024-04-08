import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


s1 = pd.read_csv('predictions/linear/sim2.csv')
s2 = pd.read_csv('predictions/nonlinear/sim2.csv')
true = pd.read_csv('predictions/linear/test_data.csv')

d = np.linspace(0, 15, 16)

def get_s(idx):
    return s1.iloc[idx,1:].values
def get_pred(idx):
    return s2.iloc[idx,1:].values
def get_predw(idx):
    return pw.iloc[idx,1:].values
def get_pred2w(idx):
    return pw2.iloc[idx,1:].values
def get_true(idx):
    return true.iloc[idx, 2:18].values
def get_time(idx):
    return true.X[idx]
def get_sim(idx):
    return sim.iloc[idx, 1:].values

def getControl(idx):
    return true.loc[idx,['Fbot','Fmid','Ftop']]

fig, ax = plt.subplots(2,1,figsize = (10, 6))
def animate(i):
    
    p = get_s(i)
    pd = get_pred(i)

    # pw = get_predw(i)
    # pw2 = get_pred2w(i)
    # s = get_sim(i)
    # s, sw = get_sd(i)p
    ax[0].clear()
    ax[0].set_xlim(-0.5, 16.5)
    ax[0].set_ylim(0, 90)
    
    ax[0].grid(True, linestyle='--', alpha=0.5)
    ax[0].set_xlabel('Depth (m)')
    ax[0].set_ylabel('Temperature (C)')
    
    ax[0].set_title(f'Water temperature at different depths {get_time(i)}')
    ax[0].plot(d, get_true(i), label='True', marker = 'o')
    ax[0].plot(d, p, label='Model (linear)', marker = 'o')
    ax[0].plot(d, pd, label='Prediction (non linear)', marker = 'o',alpha=0.5)

    
    # ax[0].plot(d, pw2, label='Prediction (14d)', marker = 'o')
    ax[0].legend()

    c = getControl(i)
    
    ax[1].clear()
    ax[1].set_ylim(-250,250)
    ax[1].bar(['Fbot','Fmid','Ftop'], c)
    ax[1].vlines(['Fbot','Fmid','Ftop'], -200, 200, linestyle='--', alpha=0.5)
    ax[1].hlines(0, -0.5, 2.5, linestyle='--', alpha=0.5)



# # print(simw)
# print(sim)
ani = FuncAnimation(fig, animate, frames=len(true), interval=0.1, repeat=True, )

#     # ani.save('model.gif', writer='Pillow', fps=40)
plt.show()
# # print(true)
