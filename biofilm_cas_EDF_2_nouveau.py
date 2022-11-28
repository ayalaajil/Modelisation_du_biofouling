#*****************************************************************************#
from math import *
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

if __name__ == "__main__":

    print('')
    print('')
    print('....')


#*****************************************************************************#
#**********************   DÉCLARATION DE VARIABLES  **************************#
#*****************************************************************************#

    # Le nombre d itération en temps (offre un modèle cohérent dès 10 points et pour un temps correct même au delà de 300 points)
    n = 300
    # Le nombre d'itération en longeur (offre un modèle cohérent dès 10 points et pour un temps correct même au delà de 300 points)
    nz = 100
    j_max = 50              # Durée de l'expérience (en jours)
    t_max = j_max*24*3600   # Durée de l'expérience (en s)
    j_max2 = 150
    t_max2 = j_max2*24*3600

    l = 10                  # La longeur du tuyau (en m)
    einf = 2e-4             # epaisseur maximale (en m)
    Rpropre = 5e-4          # résistance thermique propre du tuyau (m2.K/W)
    lambdafilm = 0.6        # Conductivité thermique du biofilm (en W/m/K)
    # capacité thermique massique de l'Eeuropeu (en J/kg/K)
    c = 4184
    qm = 0.5                # débit massique en kg/s
    Texté = 60+273.15       # Température du fluide extérieur
    TinEurope = 12+273.15   # Température d'entrée de l'eau en Europe
    TinAsie = 30+273.15     # Température d'entrée de l'eau en Asie
    d = 0.015               # diamètre interne du tuyau (en m)
    dext = 0.02             # diamètre externe du tuyau (en m)
    k25 = 650/(3600*24)     # W/(m2*K*s)
    E = 15000               # J/mol


#*****************************************************************************#
#*************************   DÉCLARATION DE TABLES  **************************#
#*****************************************************************************#

    t = linspace(0, j_max, n)  # Déclaration de la table temps à n éléments
    t2 = linspace(0, j_max2, n)
    # Déclaration des tables température à n*nz éléments
    tempeurope = np.zeros((n, nz), float)
    tempasie = np.zeros((n, nz), float)
    # Déclaration des tables épaisseur à n*nz éléments
    Eeurope = np.zeros((n, nz), float)
    Easie = np.zeros((n, nz), float)

    # Définition du pas de temps de calcul à l'aide du temps de calcul et du nombre d'itération
    dt = t_max/(n-1)
    dt2 = t_max2/(n-1)
    # Définition du pas de longeur à l'aide de la longeur totale et du nombre d'itération en longeur
    dz = l/(nz-1)

#*****************************************************************************#
#***************************    INITIALISATION    ****************************#
#*****************************************************************************#

    # Initialisation épaisseur(t=0)
    Eeurope[0, :] = (6*(10**-6))*np.ones((1, nz), float)
    Easie[0, :] = (6*(10**-6))*np.ones((1, nz), float)
    # initialisation de la température :
    tempeurope[0, 0] = TinEurope
    tempasie[0, 0] = TinAsie
    Rtheurope = Rpropre + Eeurope[0, 0]/lambdafilm
    Rtheasie = Rpropre + Easie[0, 0]/lambdafilm
    for iz in range(nz-1):
        tempeurope[0, iz+1] = tempeurope[0, iz]+dz*np.pi * \
            d*(Texté-tempeurope[0, iz])/(Rtheurope*qm*c)
        tempasie[0, iz+1] = tempasie[0, iz]+dz*np.pi * \
            d*(Texté-tempasie[0, iz])/(Rtheasie*qm*c)


#*****************************************************************************#
#***************************        CALCUL        ****************************#
#*****************************************************************************#

    # Loop over integration steps:

    # Boucle d'itération allant de 0 à n-1 car on calcule temp(t+dt), k(t+dt) et e(t+dt)
    for it in range(n-1):
        # initialisation de la température d'entré de l'Eeuropeu
        tempeurope[it+1, 0] = tempeurope[it, 0]
        tempasie[it+1, 0] = tempasie[it, 0]
        for iz in range(nz-1):

            # Cas de l'Europe
            Rtheurope = Rpropre + Eeurope[it, iz+1]/lambdafilm
            tempeurope[it+1, iz+1] = tempeurope[it+1, iz]+dz * \
                (np.pi*d*(Texté-tempeurope[it, iz+1])/(Rtheurope*qm*c))

            keurope = k25 * \
                np.exp(-(E/8.314)*((1/tempeurope[it+1, iz])-(1/298.15)))
            Eeurope[it+1, iz+1] = Eeurope[it, iz+1] + dt * keurope * \
                (einf - Eeurope[it, iz+1]) * Eeurope[it, iz+1] / lambdafilm

            # Cas de l'Asie
            Rtheasie = Rpropre + Easie[it, iz+1]/lambdafilm
            tempasie[it+1, iz+1] = tempasie[it+1, iz]+dz * \
                (np.pi*d*(Texté-tempasie[it, iz+1])/(Rtheasie*qm*c))

            kasie = k25*np.exp(-(E/8.314) *
                               ((1/tempasie[it+1, iz])-(1/298.15)))
            Easie[it+1, iz+1] = Easie[it, iz+1] + dt * kasie * \
                (einf - Easie[it, iz+1]) * Easie[it, iz+1] / lambdafilm

#*****************************************************************************#
#************************        AFFICHAGE        ****************************#
#*****************************************************************************#
    figure(3)
    plot(t, tempeurope[:, nz-1]-(273.15)*np.ones(len(tempeurope[:, nz-1])),
         'y', label="Temperature dans la sortie")
    plot(t, (27)*np.ones(len(t)), 'k--', label='Limit Temperature Europe')
    legend()
    title("Évolution de la température de l'eau dans la sortie CAS 1")
    xlabel("temps (en jours)")
    ylabel("Température (en °C)")

#*****************************************************************************#
#********************        CAS 1 : RÉFÉRENCE        ************************#
#*****************************************************************************#

    print("CAS 1")

    print("Europe")
    Tfinaleu = tempeurope[-1, -1]
    Pfinaleu = c*qm*(Tfinaleu-TinEurope)
    Neurope = 1+(20e6//Pfinaleu)
    CAPEX = 1e-3 * 20234 * (np.pi * l * dext * Neurope)**0.61  # en k€
    TAC = CAPEX * 0.094  # en k€/an
    print("N :", Neurope)
    print("CAPEX :", CAPEX)
    print("TAC :", TAC)

    print("Asie")
    Tfinalas = tempasie[-1, -1]
    Pfinalas = c*qm*(Tfinalas-TinAsie)
    Nasie = 1+(20e6//Pfinalas)
    CAPEX = 1e-3 * 20234 * (np.pi * l * dext * Nasie)**0.61  # en k€
    TAC = CAPEX * 0.094  # en k€/an
    print("N :", Nasie)
    print("CAPEX :", CAPEX)
    print("TAC :", TAC)

    figure(4)
    # plot de puissance pendant le temp Europe et Asia CAS 1
    puissance_europe = qm*c*(tempeurope[:, -1]-tempeurope[:, 0])  # en W
    puissance_asia = qm*c*(tempasie[:, -1]-tempasie[:, 0])  # en W
    plot(t, puissance_europe/1000, '-c', label='Europe')
    plot(t, puissance_asia/1000, '-y', label='Asie et Sud-Est')
    legend()
    title("Puissance vs Temp CAS 1")
    xlabel("temp (en jours)")
    ylabel("puissance (en kW)")
    #*****************************************************************************#
    #*************        CAS 2 : USAGE DE BIOCIDE CHLORÉ        *****************#
    #*****************************************************************************#

    # On a considéré que l'usage de biocide à la concentration de 5ppm permet une réduction de 70% de einf
    # Dans ce cas on va actualiser notre modèle avec cette nouvelle valeur
    ebiocide = einf*0.3

    # Déclaration des tables température à n*nz éléments
    tempeurope2 = np.zeros((n, nz), float)
    tempasie2 = np.zeros((n, nz), float)
    # Déclaration des tables épaisseur à n*nz éléments
    Eeurope2 = np.zeros((n, nz), float)
    Easie2 = np.zeros((n, nz), float)
    # Initialisation épaisseur(t=0)
    Eeurope2[0, :] = (6*(10**-6))*np.ones((1, nz), float)
    Easie2[0, :] = (6*(10**-6))*np.ones((1, nz), float)
    # initialisation de la température :
    tempeurope2[0, 0] = TinEurope
    tempasie2[0, 0] = TinAsie
    Rtheurope2 = Rpropre + Eeurope2[0, 0]/lambdafilm
    Rtheasie2 = Rpropre + Easie2[0, 0]/lambdafilm

    for iz in range(nz-1):
        tempeurope2[0, iz+1] = tempeurope2[0, iz]+dz*np.pi * \
            d*(Texté-tempeurope2[0, iz])/(Rtheurope2*qm*c)
        tempasie2[0, iz+1] = tempasie2[0, iz]+dz*np.pi * \
            d*(Texté-tempasie2[0, iz])/(Rtheasie2*qm*c)

    # Boucle d'itération allant de 0 à n-1 car on calcule temp(t+dt), k(t+dt) et e(t+dt)
    for it in range(n-1):
        # initialisation de la température d'entré de l'Eeuropeu
        tempeurope2[it+1, 0] = tempeurope2[it, 0]
        tempasie2[it+1, 0] = tempasie2[it, 0]
        for iz in range(nz-1):
            # Cas de l'Europe
            Rtheurope2 = Rpropre + Eeurope2[it, iz+1]/lambdafilm
            tempeurope2[it+1, iz+1] = tempeurope2[it+1, iz]+dz * \
                (np.pi*d*(Texté-tempeurope2[it, iz+1])/(Rtheurope2*qm*c))

            keurope2 = k25 * \
                np.exp(-(E/8.314)*((1/tempeurope2[it+1, iz])-(1/298.15)))
            Eeurope2[it+1, iz+1] = Eeurope2[it, iz+1] + dt2 * keurope2 * \
                (ebiocide - Eeurope2[it, iz+1]) * \
                Eeurope2[it, iz+1] / lambdafilm

            # Cas de l'Asie
            Rtheasie2 = Rpropre + Easie2[it, iz+1]/lambdafilm
            tempasie2[it+1, iz+1] = tempasie2[it+1, iz]+dz * \
                (np.pi*d*(Texté-tempasie2[it, iz+1])/(Rtheasie2*qm*c))

            kasie2 = k25 * \
                np.exp(-(E/8.314)*((1/tempasie2[it+1, iz])-(1/298.15)))
            Easie2[it+1, iz+1] = Easie2[it, iz+1] + dt2 * kasie2 * \
                (ebiocide - Easie2[it, iz+1]) * Easie2[it, iz+1] / lambdafilm

    print("")
    print("CAS 2")

    print("Europe")
    Tfinaleu = tempeurope2[-1, -1]
    Pfinaleu = c*qm*(Tfinaleu-TinEurope)
    Neurope = 1+(20e6//Pfinaleu)
    CAPEX = 1e-3 * 20234 * (np.pi * l * dext * Neurope)**0.61  # en k€
    OPEX = Neurope * 25.704e-3  # en k€/an
    TAC = CAPEX * 0.094 + OPEX  # en k€/an
    print("N :", Neurope)
    print("CAPEX :", CAPEX)
    print("OPEX :", OPEX)
    print("TAC :", TAC)

    print("Asie")
    Tfinalas = tempasie2[-1, -1]
    Pfinalas = c*qm*(Tfinalas-TinAsie)
    Nasie = 1+(20e6//Pfinalas)
    CAPEX = 1e-3 * 20234 * (np.pi * l * dext * Nasie)**0.61  # en k€
    OPEX = Nasie * 25.704e-3  # en k€/an
    TAC = CAPEX * 0.094 + OPEX  # en k€/an
    print("N :", Nasie)
    print("CAPEX :", CAPEX)
    print("OPEX :", OPEX)
    print("TAC :", TAC)

    figure(5)
    # plot de la puissance pendant le temp Europe et Asia CAS 2
    puissance_europe2 = qm*c*(tempeurope2[:, -1]-tempeurope2[:, 0])  # en W
    puissance_asia2 = qm*c*(tempasie2[:, -1]-tempasie2[:, 0])  # en W
    plot(t2, puissance_europe2/1000, '-g', label="Europe")
    plot(t2, puissance_asia2/1000, '-k', label="Asia et Sud-Est")
    xlabel("temp(en jours)")
    ylabel("puissance(en kW)")
    title("Puissance vs Temp CAS 2")
    legend()

    #*****************************************************************************#
    #***************        CAS 3 : SYSTÈME AUXILIAIRE        ********************#
    #*****************************************************************************#

    # On peut fixer une durée de cycle au bout de laquelle on lançera le netttoyage

    PuissEur = c*qm*(tempeurope[:, nz-1]-TinEurope)
    PuissAsie = c*qm*(tempasie[:, nz-1]-TinAsie)

    # array de durée d'un cycle d'utilisation du refroidisseur
    dcycle = np.arange(1, 50, 1)
    print("")
    print("CAS 3")

    print("Europe")
    P3Europe = np.zeros(len(dcycle))
    for i in range(len(dcycle)):
        # Puissance à partir de laquelle on considèrera le refroidisseur comme encrassé
        P3Europe[i] = PuissEur[n*dcycle[i]//j_max]
    Neurope = np.ones(49)+(20e6//P3Europe)
    CAPEXEur = 2*1e-3 * 20234 * (np.pi * l * dext * Neurope)**0.61  # en k€
    OPEXEur = (340e-3/40)*((0.35*3600e-9*(np.ones(len(dcycle)) *
                            PuissEur[-1]-P3Europe)) + 75*(np.pi * l * dext * Neurope)**0.6)  # en k€/an
    TACEur = CAPEXEur * 0.094 + OPEXEur  # en k€/an
    print("N :", Neurope[14])  # cas especific de 15 jours
    print("CAPEX :", CAPEXEur[14])  # cas especific de 15 jours
    print("OPEX :", OPEXEur[14])  # cas especific de 15 jours
    print("TAC :", TACEur[14])  # cas especific de 15 jours

    print("Asie")
    # Puissance à partir de laquelle on considèrera le refroidisseur comme encrassé
    P3Asie = PuissAsie[n*dcycle[14]//j_max]
    Nasie = 1+(20e6//P3Asie)
    CAPEX = 2*1e-3 * 20234 * (np.pi * l * dext * Nasie)**0.61  # en k€
    OPEX = (340e-3/40)*(0.35*3600e-9 *
                        (PuissAsie[-1]-P3Asie) + 75*(np.pi * l * dext * Nasie)**0.6)  # en k€/an
    TAC = CAPEX * 0.094 + OPEX  # en k€/an
    print("N :", Nasie)
    print("CAPEX :", CAPEX)
    print("OPEX :", OPEX)
    print("TAC :", TAC)

    figure(6)
    # plot de la puissance pendant le temp pour le CAS 3
    puissance_europe3 = qm*c * \
        (tempeurope[:n*dcycle[14]//j_max, -1] -
         tempeurope[:n*dcycle[14]//j_max, 0])  # en W
    puissance_asia3 = qm*c * \
        (tempasie[:n*dcycle[14]//j_max, -1] -
         tempasie[:n*dcycle[14]//j_max, 0])  # en W
    list_europe = puissance_europe3.tolist()
    list_asia = puissance_asia3.tolist()
    list_europe = list_europe + list_europe[:n//j_max]+list_europe + list_europe[:n //
                                                                                 j_max] + list_europe + list_europe[:n//j_max] + list_europe[:n*2//j_max]
    list_asia = list_asia + list_asia[:n//j_max]+list_asia + list_asia[:n //
                                                                       j_max] + list_asia+list_asia[:n//j_max] + list_asia[:n*2//j_max]

    plot(t, np.array(list_europe)/1000, '-r', label='Europe')
    plot(t, np.array(list_asia)/1000, '-b', label="Asie et Sud-Est")
    legend()
    title("Puissance vs Temp CAS 3 pour cycle de 15 jours")
    xlabel("temp(en jours)")
    ylabel("puissance(en kW)")

    #*****************************************************************************#
    #***********        CAS 4 : SYSTÈME AUXILIAIRE ET BIOCIDE        *************#
    #*****************************************************************************#

    # On peut fixer une durée de cycle au bout de laquelle on lançera le netttoyage

    PuissEur2 = c*qm*(tempeurope2[:, nz-1]-TinEurope)
    PuissAsie2 = c*qm*(tempasie2[:, nz-1]-TinAsie)

    # array de durée d'un cycle d'utilisation du refroidisseur
    dcycle = np.arange(1, 50, 1)
    print("")
    print("CAS 4")

    print("Europe")
    P4Europe = np.zeros(len(dcycle))
    for i in range(len(dcycle)):
        # Puissance à partir de laquelle on considèrera le refroidisseur comme encrassé
        P4Europe[i] = PuissEur2[n*dcycle[i]//j_max2]
    Neurope4 = np.ones(49)+(20e6//P4Europe)
    CAPEXEur4 = 2*1e-3 * 20234 * (np.pi * l * dext * Neurope4)**0.61  # en k€
    OPEXEur4 = (340e-3/40)*((0.35*3600e-9*(np.ones(len(dcycle)) *
                            PuissEur[-1]-P4Europe)) + 75*(np.pi * l * dext * Neurope)**0.6) + Neurope4*25.704e-3  # en k€/an
    TACEur4 = CAPEXEur4 * 0.094 + OPEXEur4  # en k€/an
    print("N :", Neurope4[14])  # cas spécifique de 15 jours
    print("CAPEX :", CAPEXEur4[14])  # cas spécifique de 15 jours
    print("OPEX :", OPEXEur4[14])  # cas spécifique de 15 jours
    print("TAC :", TACEur4[14])  # cas spécifique de 15 jours

    print("Asie")
    # Puissance à partir de laquelle on considèrera le refroidisseur comme encrassé
    P4Asie = PuissAsie[n*dcycle[14]//j_max2]
    Nasie4 = 1+(20e6//P4Asie)
    CAPEX4 = 2*1e-3 * 20234 * (np.pi * l * dext * Nasie4)**0.61  # en k€
    OPEX4 = (340e-3/40)*(0.35*3600e-9 *
                         (PuissAsie[-1]-P4Asie) + 75*(np.pi * l * dext * Nasie)**0.6)+Nasie*25.704e-3  # en k€/an
    TAC4 = CAPEX4 * 0.094 + OPEX4  # en k€/an
    print("N :", Nasie4)
    print("CAPEX :", CAPEX4)
    print("OPEX :", OPEX4)
    print("TAC :", TAC4)

    figure(7)
    # plot de la puissance pendant le temp pour le CAS 4
    puissance_europe4 = qm*c * \
        (tempeurope2[:n*dcycle[14]//j_max2, -1] -
         tempeurope2[:n*dcycle[14]//j_max2, 0])  # en W
    puissance_asia4 = qm*c * \
        (tempasie2[:n*dcycle[14]//j_max2, -1] -
         tempasie2[:n*dcycle[14]//j_max2, 0])  # en W
    list_europe4 = puissance_europe4.tolist()
    list_asia4 = puissance_asia4.tolist()
    list_europe4 = (
        list_europe4 + list_europe4[:n//j_max2])*9+list_europe4[:n*6//j_max2]
    list_asia4 = (
        list_asia4 + list_asia4[:n//j_max2])*9 + list_asia4[:n*6//j_max2]

    plot(t2, np.array(list_europe4)/1000, '-r', label='Europe')
    plot(t2, np.array(list_asia4)/1000, '-b', label="Asie et Sud-Est")
    legend()
    title("Puissance vs Temp CAS 4 pour cycle de 15 jours")
    xlabel("temp(en jours)")
    ylabel("puissance(en kW)")

    figure(8)
    # plot de CAPEX, OPEX et TAC en fonction du nombre de cycles
    plot(dcycle, CAPEXEur, 'b.', label="CAS 3")
    plot(dcycle, CAPEXEur4, 'r.', label="CAS 4")
    legend()
    title("CAPEX en fonction de nombre de cycles")
    xlabel("durée de un cycle (en jours)")
    ylabel("Value (en k€)")
    figure(9)
    plot(dcycle, OPEXEur, 'bx', label="CAS 3")
    plot(dcycle, OPEXEur4, 'rx', label="CAS 4")
    legend()
    title("OPEX en fonction de nombre de cycles")
    xlabel("durée de un cycle (en jours)")
    ylabel("Value (en k€)")
    figure(10)
    plot(dcycle, TACEur, 'k^', label="CAS 3")
    plot(dcycle, TACEur4, 'c^', label="CAS 4")
    legend()
    title("TAC en fonction de nombre de cycles")
    xlabel("durée de un cycle (en jours)")
    ylabel("Value (en k€)")

    # plot animation de la épaissur Asia et Europe CAS 1
    fig12, ax12 = plt.subplots()
    line, = ax12.plot(linspace(0, 10, nz-1),
                      Eeurope[0, 1:]*1e6, '-b', label='Europe')
    lineA, = ax12.plot(linspace(0, 10, nz-1),
                       Easie[0, 1:]*1e6, '-g', label='Asie et Sud-Est')
    ax12.legend(loc='upper right')
    ax12.set_xlabel("longeur (en m)")
    ax12.set_ylabel("épaisseur de biofilm (en µm)")
    ax12.set_title("Croissance épaisseur de biofilm comparaison CAS 1")
    ax12.set_ylim([min(Eeurope[0, 1:])*1e6 - 1, max(Eeurope[-1, 1:])*1e6+1])
    time_template = 'time = %.1f jours'
    time_text = ax12.text(0.05, 0.9, '', transform=ax12.transAxes)

    def animate(i):
        line.set_ydata(Eeurope[i, 1:]*1e6)
        lineA.set_ydata(Easie[i, 1:]*1e6)
        time_text.set_text(time_template % ((i/86400)*dt))
        return line, lineA, time_text

    ani = animation.FuncAnimation(
        fig12, animate, range(1, 300), interval=50, blit=False, save_count=50)

    # plot animation de la épaisseur Asia et Europe CAS 2
    fig11, ax11 = plt.subplots()
    line1, = ax11.plot(linspace(0, 10, nz-1),
                       Eeurope2[0, 1:]*1e6, '-b', label='Europe')
    line1A, = ax11.plot(linspace(0, 10, nz-1),
                        Easie2[0, 1:]*1e6, '-g', label='Asie et Sud-Est')
    ax11.legend(loc='upper right')
    ax11.set_xlabel("longeur (en m)")
    ax11.set_ylabel("épaisseur de biofilm (en µm)")
    ax11.set_title("Croissance épaisseur de biofilm comparaison CAS 2")
    ax11.set_ylim([min(Eeurope2[0, 1:])*1e6 - 1, max(Easie2[-1, 1:])*1e6+1])
    time_template1 = 'time = %.1f jours'
    time_text1 = ax11.text(0.05, 0.9, '', transform=ax11.transAxes)

    def animate2(i):
        line1.set_ydata(Eeurope2[i, 1:]*1e6)
        line1A.set_ydata(Easie2[i, 1:]*1e6)
        time_text1.set_text(time_template1 % ((i/86400)*dt2))
        return line1, line1A, time_text1

    ani2 = animation.FuncAnimation(
        fig11, animate2, range(1, 300), interval=50, blit=False, save_count=50)
    figure(13)
    # plot de la temperature de sortie pendant le temp pour CAS 2
    plot(t2, tempeurope2[:, nz-1]-(273.15)*np.ones(len(tempeurope2[:, nz-1])),
         'b', label="Temperature de l'eau dans la sortie")
    plot(t2, (27)*np.ones(len(t2)), 'k--', label='Limit Temperature Europe')
    legend()
    title("Évolution de la température de l'eau dans la sortie CAS 2")
    xlabel("temps (en jours)")
    ylabel("Température (en °C)")

    figure(14)
    # plot de la temperature de sortie pendant le temp pour CAS 3
    plot(t, np.array(list_europe)/(qm*c) + TinEurope *
         np.ones(len(list_europe)) - (273.15)*np.ones(len(list_europe)), 'r', label="Temperature dans la sortie")
    plot(t, (27)*np.ones(len(t)), 'k--', label='Limit Temperature Europe')
    legend()
    title("Évolution de la température de l'eau dans la sortie CAS 3 pour cycle de 15 jours")
    xlabel("temps (en jours)")
    ylabel("Température (en °C)")

    figure(15)
    # plot de la temperature de sortie pendant le temp pour CAS 4
    plot(t2, np.array(list_europe4)/(qm*c) + TinEurope *
         np.ones(len(list_europe4)) - (273.15)*np.ones(len(list_europe4)), 'g', label="Temperature dans la sortie")
    plot(t2, (27)*np.ones(len(t)), 'k--', label='Limit Temperature Europe')
    legend()
    title("Évolution de la température de l'eau dans la sortie CAS 4 pour cycle de 15 jours")
    xlabel("temps (en jours)")
    ylabel("Température (en °C)")

    figure(16)
    plot(t, tempasie[:, nz-1]-(273.15)*np.ones(len(tempasie[:, nz-1])),
         'y', label="Temperature dans la sortie")
    plot(t, (45)*np.ones(len(t)), 'k--', label='Limit Temperature Asie')
    legend()
    title("Évolution de la température de l'eau dans la sortie CAS 1")
    xlabel("temps (en jours)")
    ylabel("Température (en °C)")

    figure(17)
    plot(t2, tempasie2[:, nz-1]-(273.15)*np.ones(len(tempasie2[:, nz-1])),
         'b', label="Temperature de l'eau dans la sortie")
    plot(t2, (45)*np.ones(len(t2)), 'k--', label='Limit Temperature Asie')
    legend()
    title("Évolution de la température de l'eau dans la sortie CAS 2")
    xlabel("temps (en jours)")
    ylabel("Température (en °C)")

    figure(18)
    # plot de la temperature de sortie pendant le temp pour CAS 3
    plot(t, np.array(list_asia)/(qm*c) + TinAsie *
         np.ones(len(list_europe)) - (273.15)*np.ones(len(list_asia)), 'r', label="Temperature dans la sortie")
    plot(t, (45)*np.ones(len(t)), 'k--', label='Limit Temperature Asie')
    legend()
    title("Évolution de la température de l'eau dans la sortie CAS 3 pour cycle de 15 jours")
    xlabel("temps (en jours)")
    ylabel("Température (en °C)")

    figure(19)
    # plot de la temperature de sortie pendant le temp pour CAS 4
    plot(t2, np.array(list_asia4)/(qm*c) + TinAsie *
         np.ones(len(list_asia4)) - (273.15)*np.ones(len(list_asia4)), 'g', label="Temperature dans la sortie")
    plot(t2, (45)*np.ones(len(t)), 'k--', label='Limit Temperature Asie')
    legend()
    title("Évolution de la température de l'eau dans la sortie CAS 4 pour cycle de 15 jours")
    xlabel("temps (en jours)")
    ylabel("Température (en °C)")
    plt.show()
