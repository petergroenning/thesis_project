import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


p = pd.read_csv('p3.csv')
# pw = pd.read_csv('p3w.csv')
sim = pd.read_csv('sim.csv')
true = pd.read_csv('data.csv')

d = np.linspace(0, 15, 16)

def get_pred(idx):
    return p.iloc[idx,1:].values
def get_predw(idx):
    return pw.iloc[idx,1:].values
def get_true(idx):
    return true.iloc[idx, 2:18].values
def get_time(idx):
    return true.X[idx]
def get_sim(idx):
    return sim.iloc[idx, 1:].values



fig, ax = plt.subplots(1,1,figsize = (10, 6))
def animate(i):
    p = get_pred(i)
    # pw = get_predw(i)
    s = get_sim(i)
    # s, sw = get_sd(i)
    ax.clear()
    ax.set_xlim(-0.5, 16.5)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.set_xlabel('Depth (m)')
    ax.set_ylabel('Temperature (C)')
    ax.set_title(f'Water temperature at different depths {get_time(i)}')
    ax.plot(d, get_true(i), label='True', marker = 'o')
    ax.plot(d, p, label='Predicted (24h)', marker = 'o')
    # ax.plot(d, pw, label='Predicted (1w)', marker = 'o')
    ax.plot(d, s, label='Simulated', marker = 'o')
    ax.legend()



# print(simw)
print(sim)
ani = FuncAnimation(fig, animate, frames=len(sim), interval=1, repeat=True, )

    # ani.save('model.gif', writer='Pillow', fps=40)
plt.show()
# print(true)
