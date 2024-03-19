import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
path = 'data.csv'
data = pd.read_csv(path)
# data.rename(columns={'Unnamed: 0': 't'}, inplace=True)
# data.set_index('X', inplace=True)
data.t = data.t - data.t.min()
data.set_index('t', inplace=True)
def rk4(f, x0, t0, tf, dt):
    t = np.arange(t0, tf, dt)
    n = len(t)
    x = np.zeros((n, len(x0)))
    x[0] = x0
    
    for i in tqdm(range(n-1)):
        inputs = data.loc[i, ['Fbot', 'Fmid', 'Ftop', 'Tbot', 'Tmid', 'Ttop', 'FbotIn', 'FbotOut', 'FmidIn', 'FmidOut', 'FtopIn', 'FtopOut']]
        k1 = dt * f(t[i], x[i], **inputs)
        k2 = dt * f(t[i] + 0.5*dt, x[i] + 0.5*k1, **inputs)
        k3 = dt * f(t[i] + 0.5*dt, x[i] + 0.5*k2, **inputs)
        k4 = dt * f(t[i] + dt, x[i] + k3, **inputs)
        x[i+1] = x[i] + (k1 + 2*k2 + 2*k3 + k4) / 6
    return t, x


def f(t, x):
    return 

x15 = data.X15.values

plt.plot(x15)
plt.show()

