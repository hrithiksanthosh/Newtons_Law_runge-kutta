import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable

time_value = np.array([])
temp_value = np.array([])
tempTable = PrettyTable(["Time", "K1", "K2", "K3", "K4", "Temperature"])
time_value_AD = np.array([])
temp_value_AD = np.array([])
tempTable_AD = PrettyTable(["Time", "Temperature"])


def f(time, Temp):
    K = 0.04
    amp_temp = 23
    return -1 * K * (Temp - amp_temp)


def rk4(x0, y0, xn, h):
    n = (int)((xn - x0) / h)
    global time_value, temp_value
    time_value = np.append(time_value, x0)
    temp_value = np.append(temp_value, y0)
    tempTable.add_row([x0, "", "", "", "", y0])
    for i in range(n):
        k1 = h * (f(x0, y0))
        k2 = h * (f((x0 + h / 2), (y0 + k1 / 2)))
        k3 = h * (f((x0 + h / 2), (y0 + k2 / 2)))
        k4 = h * (f((x0 + h), (y0 + k3)))
        k = (k1 + 2 * k2 + 2 * k3 + k4) / 6
        yn = y0 + k
        y0 = yn
        x0 = x0 + h
        time_value = np.append(time_value, x0)
        temp_value = np.append(temp_value, y0)
        tempTable.add_row([x0, k1, k2, k3, k4, y0])
    print(tempTable)


def AD(x0, xn, y0, h):
    N = (int)((xn - x0) / h)
    global time_value_AD, temp_value_AD
    time_value_AD = np.append(time_value_AD, x0)
    temp_value_AD = np.append(temp_value_AD, y0)
    tempTable_AD.add_row([x0, y0])
    for i in range(0, N):
        if i in range(0, 3):
            k1 = h * (f(x0, y0))
            k2 = h * (f((x0 + h / 2), (y0 + k1 / 2)))
            k3 = h * (f((x0 + h / 2), (y0 + k2 / 2)))
            k4 = h * (f((x0 + h), (y0 + k3)))
            k = (k1 + 2 * k2 + 2 * k3 + k4) / 6
            yn = y0 + k
            y0 = yn
            x0 = x0 + h
            time_value_AD = np.append(time_value_AD, x0)
            temp_value_AD = np.append(temp_value_AD, y0)
            tempTable_AD.add_row([x0, "Runge Kutta Method = " + str(y0)])
        else:
            yn = y0 + h * (
                    55.0 * f(time_value_AD[i], temp_value_AD[i]) - 59.0 *
                    f(time_value_AD[i - 1],

                      temp_value_AD[i - 1]) + 37.0 * f(
                time_value_AD[i - 2], temp_value_AD[i - 2]) - 9.0 * f(
                time_value_AD[i - 3], temp_value_AD[i - 3])) / 24.0
            y0 = yn
            x0 = x0 + h
            time_value_AD = np.append(time_value_AD, x0)
            temp_value_AD = np.append(temp_value_AD, y0)
            tempTable_AD.add_row([x0, y0])
    print(tempTable_AD)
    return (time_value_AD, temp_value_AD)


print('Enter initial conditions:')
x0 = float(input('Initial time = '))
y0 = float(input('Initial temperature = '))
print('Enter calculation point: ')
xn = float(input('Final time = '))
print('Enter number of steps:')
step = float(input('Step size need to consider = '))
rk4(x0, y0, xn, step)
AD(x0, xn, y0, step)
fig, ax = plt.subplots(2, figsize=(10, 10))
ax[0].plot(time_value, temp_value, color="green")
ax[0].set_xlabel("Time", color="Blue")
ax[0].set_ylabel("Temperature - Runge Kutta Method", color="Brown")
ax[0].axvline(0, c='black', ls='--')
ax[0].axhline(0, c='black', ls='--')
ax[1].plot(time_value_AD, temp_value_AD, color="red")
ax[1].set_xlabel("Time", color="Blue")
ax[1].set_ylabel("Temperature - Adam Bashforth", color="Orange")
ax[1].axvline(0, c='black', ls='--')
ax[1].axhline(0, c='black', ls='--')
plt.show()
