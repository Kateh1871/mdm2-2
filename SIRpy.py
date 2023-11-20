import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

beta = 0.005
gamma = 0.0005
delta = 0.001

def model(y, t):
    s = y[0]
    i = y[1]
    r = y[2]
    dsdt = -beta * s * i / (s + i + r) + delta * r / (s + i + r)
    didt = beta * s * i / (s + i + r) - gamma * i / (s + i + r)
    drdt = gamma * i / (s + i + r) - delta * r / (s + i + r)
    return [dsdt, didt, drdt]

y0 = [0.99, 0.01, 0]
tSIR = np.linspace(0, 10000, 100)
t = np.linspace(0, 10000, 100)

ySIR = odeint(model, y0, t)

y = odeint(model, y0, t)
#print(y[1000])

plt.plot(t, y[:, 0], 'r', label='s(t)')
plt.plot(t, y[:, 1], 'b', label='i(t)')
plt.plot(t, y[:, 2], 'g', label='r(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.savefig("IdealGraph")
plt.show()
