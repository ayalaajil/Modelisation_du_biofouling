from scipy.integrate import solve_ivp
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def solve_edp_Tof(n, Text, Rtuyau, lba, c, qm, e, d):
    dz = 3.1/n
    T = np.empty(n)
    T[0] = 290.15

    for i in range(n-1):
        Rth = Rtuyau + e[i]/lba
        dTdz = (Text-T[i])*(math.pi)*d/(Rth*c*qm)
        T[i+1] = T[i]+dTdz*dz

    return T


def solve_edp_Tvst(einf, lba, n, E, R, k25):
    t = np.linspace(0, 100, n)
    dt = 100/n
    z = np.linspace(0, 3.1, n)
    e = np.empty((n, n))
    T = np.empty((n, n))
    T[:, 0] = solve_edp_Tof(n, 308.15, 5e-4, 0.6, 4180,
                            216, np.zeros(n), 12e-3)
    e[:, 0] = 0.000001
    for i in range(n-1):
        k = k25*np.exp((-E/R)*(1/T[:, i] - 1/298.15))
        dedt = (k/lba)*(einf-e[:, i])*e[:, i]
        T1 = solve_edp_Tof(n, 308.15, 5e-4, 0.6, 4180, 0.208, e[:, i], 12e-3)
        T[:, i+1] = T1
        e[:, i+1] = e[:, i]+dedt*dt
        if i == 50:
            e_test = e[:, i]
    Ts = solve_edp_Tof(1000, 308.15, 1.6e-4, 0.6, 4180,
                       0.208, np.ones(1000)*0.58e-3, 12e-3)
    fig, axs = plt.subplots(2, 2)
    T[:, 0] = np.nan
    axs[0, 0].plot(t, e[-1, :], '-r')
    axs[0, 0].set_title("espessure vs temp à z=3.1")
    axs[0, 1].plot(t, T[-1, :], '-b')
    axs[0, 1].set_title("Temperature vs temp à z=3.1")
    axs[1, 0].plot(z, Ts, '-g')
    axs[1, 0].set_title("Temperature vs z à t infinite")
    axs[1, 1].plot(z, e_test, '-k')
    axs[1, 1].set_title("espessure vs z à t=5 days")
    plt.show()


solve_edp_Tvst((np.ones(1000)*0.2e-3), 0.6, 1000, 25000, 8.314, 715)
