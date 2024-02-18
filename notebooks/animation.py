import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from functools import partial
import pandas as pd

df = pd.read_csv('data/processed/data_depths.csv', index_col=0)
df = df[df.EVENT_TIME != '2023-09-04 14:00:00.000']

df2 = df[df.DEPTHS > 0.5]
times = df2.EVENT_TIME.sort_values().unique()



D = []
T = []

for t in times:
    temp = df2[df2.EVENT_TIME == t]
    temp = temp.sort_values('DEPTHS')

    D.append(temp.DEPTHS.values)
    T.append(temp.VALUE_CUR)



fig, ax = plt.subplots()


def animate(i):
    ax.clear()
    ax.set_xlim(0, 14)
    ax.set_ylim(50, 90)
    ax.grid()
    ax.set_ylabel('Temperature (°C)', fontsize=12)
    ax.set_xlabel('Depth (m)', fontsize=12)
    ax.scatter(D[i], T[i], c='r', label='Temperature (°C) at sensor', alpha = .7)
    ax.set_title(f'Time: {times[i]}', fontsize=15)
    ax.legend()


ani = FuncAnimation(fig, animate, frames=len(D), interval=2, repeat=True)
# Save as gif
# ani.save('animation.gif', writer='Pillow', fps=40)

plt.show()