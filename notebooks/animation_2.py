import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import scipy.optimize as opt
from tqdm import tqdm

water_sensors = pd.read_csv('data/processed/dronninglund/water_sensors.csv', parse_dates=True, index_col=0)
inputs = pd.read_csv('data/processed/dronninglund/inputs.csv', parse_dates=True, index_col=0)
times = water_sensors.index


d = np.linspace(0, 1, 16)

def get_water_temp(idx):
    return water_sensors.iloc[idx]

def get_input(idx):
    return inputs.iloc[idx]

def logistic(d: float, t: float = 85, k: float = 10, a: float = 0.5):
    ''' Computes the temperature as a function of depth d using a logistic function
        Note: k is the steepness of the curve, a is the inflection point (k >= 0)
    '''
    return t / (1 + np.exp( (a-d) / k)) - t / (1 + np.exp(a/k))

def logistic_temp(d: float, t0: float = 50, t1: float = 30, k1: float = 0.5, a1: float = 0.5, t2: float = 10, k2: float = 0.5, a2: float = 0.5) -> float:
    ''' Computes the temperature as a function of depth d using a logistic function
        Note: k is the steepness of the curve, a is the inflection point (k >= 0)
    '''

    l1 = logistic(d, t1, k1, a1)
    l2 = logistic(d, t2, k2, a2)

    return t0 + l1 + l2

def MSE(temp: np.ndarray, pred: np.ndarray, *args) -> float:
    ''' Computes the mean squared error between the true and predicted values'''
    return np.mean((temp - pred)**2)



def l1_loss(temp: np.ndarray, pred: np.ndarray, *args) -> float:
    ''' Computes the loss function'''
    mse = MSE(temp, pred)
    params = np.array(args)
    return mse + 0.001 * np.sum(np.abs(params))


def loss(params: np.ndarray, temp: np.ndarray, d: np.ndarray, criterion: callable = MSE) -> float:
    ''' Computes the loss function'''
    return criterion(temp, logistic_temp(d, *params))


def minimize_loss(temp: np.ndarray, d: np.ndarray, criterion: callable = MSE, loss: callable = loss, initial_guess: np.ndarray = np.array([0, 0.2, 0.5, 0.5, 0.1, 0.5, 0.5])) -> np.ndarray:
    ''' Minimizes the loss function using the scipy.optimize.minimize function'''
    return opt.minimize(loss, initial_guess, args=(temp, d, criterion),
                        bounds = [(0, 1), (0, 1), (1e-5, 1), (0,1), (0, 1), (1e-5,1), (0,1)]
                        ).x


def get_params(temp: np.ndarray, d: np.ndarray, criterion: callable = MSE, loss: callable = loss, initial_guess: np.ndarray = np.array([0, 0.2, 0.5, 0.5, 0.1, 0.5, 0.5])) -> np.ndarray:
    ''' Minimizes the loss function using the scipy.optimize.minimize function'''
    P = []

    for i in tqdm(range(start,start + N_frames)):

        temp = get_water_temp(i)
        params = minimize_loss(temp, d, criterion, loss, initial_guess)

        P.append(params)
        initial_guess = params
    return P


def temp_profile(ax, idx):
    idx += start
    t = get_water_temp(idx)
    i = get_input(idx)

    params = P[idx - start]

    log_temp = logistic_temp(d, *params)
    l1 = logistic(d, params[1], params[2], params[3])
    l2 = logistic(d, params[4], params[5], params[6])
    
    res = t - log_temp

    ax[0,0].scatter(d, t, c = 'C1')
    ax[0,0].plot(d, log_temp)
    ax[0,0].plot(d, l1, '--', label = 'l1')
    ax[0,0].plot(d, l2, '--', label = 'l2')

    ax[0,0].set_ylim(0, 1)
    ax[0,0].grid('--')
    ax[0,0].set_title(t.name)
    ax[0,0].legend()


    ax[1,0].set_xticks(np.arange(3))
    ax[1,0].set_xticklabels(['F_bot', 'F_mid', 'F_top'])
    ax[1,0].set_ylim(-1, 1)
    ax[1,0].grid('--')
    ax[1,0].bar(np.arange(3), i.values[:3])
    ax[1,0].set_title('Input')

    ax[0,1].scatter(d, res)
    ax[0,1].grid('--')
    ax[0,1].set_ylim(-0.02, 0.02)
    ax[0,1].set_title('Residuals')

    params_list = ['t0', 't1', 'k1', 'a1', 't2', 'k2', 'a2']
    for i in range(1,7):
        ax[1, 1].plot(np.arange(idx - start), np.array(P)[:(idx-start),i], label = params_list[i], c = f'C{i}', linestyle = '-.')
    ax[1, 1].set_xlim(0, N_frames)
    ax[1, 1].set_ylim(0, 1)
    ax[1, 1].grid('--')
    ax[1, 1].legend(loc = 'upper right', ncol = 2)
    # ax[1, 1].scatter([idx - start], params[1], label = 't1', c = 'C0', s = 1)
    # ax[1, 1].scatter([idx - start], params[2], label = 'k1', c = 'C1', s = 1)
    # ax[1, 1].scatter([idx - start], params[3], label = 'a1', c = 'C2', s = 1)
    # ax[1, 1].scatter([idx - start], params[4], label = 't2', c = 'C3', s = 1)
    # ax[1, 1].scatter([idx - start], params[5], label = 'k2', c = 'C4', s = 1)
    # ax[1, 1].scatter([idx - start], params[6], label = 'a2', c = 'C5', s = 1)

    
    


start = 300
N_frames = 500

P = get_params(get_water_temp(0), d)

fig, ax = plt.subplots(2,2, figsize = (10,8))

# ax[1, 1].legend(loc = 'upper right')

def animate(i):
    ax[0,0].clear()
    ax[0,1].clear()
    ax[1,0].clear()
    ax[1,1].clear()

    temp_profile(ax, i)
ani = FuncAnimation(fig, animate, frames=N_frames, interval=1, repeat=True, )
# Save as gif
# ani.save('test.gif', writer='Pillow', fps=40)
plt.show()
    
