#*****************************************************************************#
from math import *
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

if __name__ == "__main__":

    print('')
    print('')
    print('....')


#*****************************************************************************#
#**********************   DÉCLARATION DE VARIABLES  **************************#
#*****************************************************************************#

    n = 300  # Le nombre d itération en temps
    nz = 300  # Le nombre d'itération en longeur
    j_max = 50  # Durée de l'expérience (en s)
    t_max = j_max*24*3600  # Durée de l'expérience (en s)
    l = 3.1  # La longeur du tuyau (en m)
    einf = 0.0002  # epaisseur maximale (en m)
    Rpropre = 5e-4  # résistance thermique propre du tuyau (m2.K/W)
    # résistance thermique d'une section élémentaire du tuyau (m.K/W)
    Rtuyau = Rpropre
    lambdafilm = 0.6  # Conductivité thermique du biofilm (en W/m/K)
    c = 4186  # capacité thermique massique de l'eau (en J/kg/K)
    qm = 0.2086  # débit massique en kg/s
    Texté = 35+273.15  # Température du fluide extérieur
    d = 0.012  # diamètre tuyau
    k25c = 1000/(3600*24)
    k25b = 716/(3600*24)  # W/m2/K/s
    k25a = 500/(3600*24)
    E = 40000              # J/mol


#*****************************************************************************#
#*************************   DÉCLARATION DE TABLES  **************************#
#*****************************************************************************#

    t = linspace(0, j_max, n)  # Déclaration du tableau temps à n éléments
    # Déclaration du tableau température à n*nz éléments
    tempa = np.zeros((n, nz), float)
    tempb = np.zeros((n, nz), float)
    tempc = np.zeros((n, nz), float)
    # Déclaration du tableau épaisseur à n*nz éléments
    ea = np.zeros((n, nz), float)
    eb = np.zeros((n, nz), float)
    ec = np.zeros((n, nz), float)
    k = np.zeros((n, nz), float)  # Déclaration du tableau k à n*nz éléments

    # Définition du pas de temps de calcul à l'aide du temps de calcul et du nombre d'itération
    dt = t_max/(n-1)
    # Définition du pas de longeur à l'aide de la longeur totale et du nombre d'itération en longeur
    dz = l/(nz-1)

#*****************************************************************************#
#***************************    INITIALISATION    ****************************#
#*****************************************************************************#

    # Initialisation épaisseur(t=0)
    ea[0, :] = (6*(10**-6))*np.ones((1, nz), float)
    eb[0, :] = (6*(10**-6))*np.ones((1, nz), float)
    ec[0, :] = (6*(10**-6))*np.ones((1, nz), float)  # print(e[0,0])
    # initialisation de la température :
    tempa[0, 0] = 17+273.15
    tempb[0, 0] = 17+273.15
    tempc[0, 0] = 17+273.15
    Rtha = Rtuyau + ea[0, 0]/lambdafilm
    Rthb = Rtuyau + eb[0, 0]/lambdafilm
    Rthc = Rtuyau + ec[0, 0]/lambdafilm
    for iz in range(nz-1):
        tempa[0, iz+1] = tempa[0, iz]+dz*np.pi * \
            d*(Texté-tempa[0, iz])/(Rtha*qm*c)
        tempb[0, iz+1] = tempb[0, iz]+dz*np.pi * \
            d*(Texté-tempb[0, iz])/(Rthb*qm*c)
        tempc[0, iz+1] = tempc[0, iz]+dz*np.pi * \
            d*(Texté-tempc[0, iz])/(Rthc*qm*c)


#*****************************************************************************#
#***************************        CALCUL        ****************************#
#*****************************************************************************#

    # Loop over integration steps:

    # Boucle d'itération allant de 0 à n-1 car on calcule temp(t+dt), k(t+dt) et e(t+dt)
    for it in range(n-1):
        # initialisation de la température d'entré de l'eau
        tempa[it+1, 0] = tempa[it, 0]
        tempb[it+1, 0] = tempb[it, 0]
        tempc[it+1, 0] = tempc[it, 0]
        for iz in range(nz-1):
            Rtha = Rtuyau + ea[it, iz+1]/lambdafilm
            tempa[it+1, iz+1] = tempa[it+1, iz]+dz * \
                (np.pi*d*(Texté-tempa[it, iz+1])/(Rtha*qm*c))

            ka = k25a*np.exp(-(E/8.314)*((1/tempa[it+1, iz])-(1/298.15)))
            ea[it+1, iz+1] = ea[it, iz+1] + dt * ka * \
                (einf - ea[it, iz+1]) * ea[it, iz+1] / lambdafilm

            Rthb = Rtuyau + eb[it, iz+1]/lambdafilm
            tempb[it+1, iz+1] = tempb[it+1, iz]+dz * \
                (np.pi*d*(Texté-tempb[it, iz+1])/(Rthb*qm*c))

            kb = k25b*np.exp(-(E/8.314)*((1/tempb[it+1, iz])-(1/298.15)))
            eb[it+1, iz+1] = eb[it, iz+1] + dt * kb * \
                (einf - eb[it, iz+1]) * eb[it, iz+1] / lambdafilm

            Rthc = Rtuyau + ec[it, iz+1]/lambdafilm
            tempc[it+1, iz+1] = tempc[it+1, iz]+dz * \
                (np.pi*d*(Texté-tempc[it, iz+1])/(Rthc*qm*c))

            kc = k25c*np.exp(-(E/8.314)*((1/tempc[it+1, iz])-(1/298.15)))
            ec[it+1, iz+1] = ec[it, iz+1] + dt * kc * \
                (einf - ec[it, iz+1]) * ec[it, iz+1] / lambdafilm

#*****************************************************************************#
#************************        AFFICHAGE        ****************************#
#*****************************************************************************#

    figure(1)
    plot(t, eb[:, nz-1]*1e6, 'r', label="longeur totale")
    plot(t, eb[:, 3*nz//4]*1e6, 'y', label="75%")
    plot(t, eb[:, nz//2]*1e6, 'b', label="50%")
    plot(t, eb[:, nz//4]*1e6, 'g', label="25%")
    legend()
    xlabel("temps (en jours)")
    ylabel("épaisseur du biofilm (en µm)")

    figure(2)
    plot(linspace(0, 3.1, nz-1), eb[n//10, 1:]*1e6, 'b', label="t=5 jours")
    plot(linspace(0, 3.1, nz-1), eb[n//5, 1:]*1e6, 'r', label="t=10jours")
    plot(linspace(0, 3.1, nz-1), eb[n*3//10, 1:]*1e6, 'c', label="t=15 jours")
    plot(linspace(0, 3.1, nz-1), eb[n*2//5, 1:]*1e6, 'm', label="t=20 jours")
    plot(linspace(0, 3.1, nz-1), eb[n*5//10, 1:]*1e6, 'y', label="t=25 jours")
    plot(linspace(0, 3.1, nz-1), eb[n*3//5, 1:]
         * 1e6, 'dodgerblue', label="t=30 jours")
    plot(linspace(0, 3.1, nz-1), eb[n*7//10, 1:]
         * 1e6, 'orange', label="t=35 jours")
    plot(linspace(0, 3.1, nz-1), eb[n*4//5, 1:]*1e6, 'k', label="t=40 jours")
    plot(linspace(0, 3.1, nz-1), eb[n-1, 1:] *
         1e6, 'slategray', label="t=50 jours")
    legend()
    xlabel("longeur (en m)")
    ylabel("épaisseur de biofilm (en µm)")

    figure(3)
    plot(t, tempb[:, nz-1], 'r', label="longeur totale")
    plot(t, tempb[:, 3*nz//4], 'y', label="75%")
    plot(t, tempb[:, nz//2], 'b', label="50%")
    plot(t, tempb[:, nz//4], 'g', label="25%")
    legend()
    xlabel("temps (en jours)")
    ylabel("Température (en K)")

    figure(4)
    ax = plt.axes(projection='3d')
    img_3d = ax.contour3D(linspace(0, 3.1, nz), t, tempb, 90,
                          cmap='viridis', edgecolor='none')
    ax.set_xlabel('longeur (en m)')
    ax.set_ylabel('temps (en jours)')
    ax.set_zlabel('température (en K)')
    ax.view_init(40, -140)
    colorbar(img_3d)

    Rfa = np.zeros(n-1)
    Rfb = np.zeros(n-1)
    Rfc = np.zeros(n-1)
    for it in range(n-1):
        m = 0
        for iz in range(nz-1):
            m += ea[it+1, iz+1]
        Rfa[it] = m
    for it in range(n-1):
        m = 0
        for iz in range(nz-1):
            m += eb[it+1, iz+1]
        Rfb[it] = m
    for it in range(n-1):
        m = 0
        for iz in range(nz-1):
            m += ec[it+1, iz+1]
        Rfc[it] = m
    Rfa = Rfa/(lambdafilm*(nz-1))
    Rfb = Rfb/(lambdafilm*(nz-1))
    Rfc = Rfc/(lambdafilm*(nz-1))

    figure(5)
    plot(t[1:], Rfa, '-r', label='k25=500')
    plot(t[1:], Rfb, '-g', label='k25=716')
    plot(t[1:], Rfc, '-b', label='k25=1000')
    legend()
    plot(np.array([0.19, 0.76, 1.33, 1.71, 2.29, 3.05, 3.82, 4.77, 5.92, 7.06, 7.83, 8.97, 9.93, 11.08, 12.03, 13.37, 14.15, 14.91, 15.88, 17.40, 17.98, 19.12, 20.08, 20.08, 21.04, 22.00, 23.13, 23.90, 25.04, 26.18, 26.95, 28.09, 29.05, 30.20, 31.15, 31.92, 34.20, 34.98, 35.93, 37.27, 38.04, 38.98, 39.95, 40.91, 42.40, 43.35, 44.12, 44.69, 46.41, 47.16, 47.76, 49.08, 50.24, 51.19, 52.14, 52.89, 53.85, 54.99, 56.12, 56.89, 58.00, 58.93]), np.array([1.27e-05, 5.08e-05, 8.89e-05, 1.14e-04, 4.15e-03, 1.62e-02, 3.63e-02, 1.63e-02, 3.64e-02, 4.05e-02, 4.45e-02,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                    4.86e-02, 7.27e-02, 7.67e-02, 8.48e-02, 9.29e-02, 1.37e-01, 1.33e-01, 1.77e-01, 1.65e-01, 1.97e-01, 1.89e-01, 2.01e-01, 2.09e-01, 2.37e-01, 2.41e-01, 2.26e-01, 2.42e-01, 2.42e-01, 2.34e-01, 2.38e-01, 2.42e-01, 2.62e-01, 2.78e-01, 2.70e-01, 2.82e-01, 2.86e-01, 3.30e-01, 3.06e-01, 3.30e-01, 3.59e-01, 3.31e-01, 3.71e-01, 3.91e-01, 2.91e-01, 2.99e-01, 3.19e-01, 3.03e-01, 3.19e-01, 2.91e-01, 3.67e-01, 3.39e-01, 3.71e-01, 3.67e-01, 3.55e-01, 3.24e-01, 3.48e-01, 3.52e-01, 3.04e-01, 3.24e-01, 2.40e-01, 2.00e-01])*1e-3, 's', color='orange')
    xlabel("temps (en jours)")
    ylabel("Rf")
    plt.show()
