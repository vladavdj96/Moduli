import numpy as np
import pandas as pd
import math


def nosivsotAB(b, h, d1, As1, As2, fck, fyk):
    d = h - d1
    # Beton
    gamma_c = 1.5
    fcd = 0.85 * fck / gamma_c
    eps_c2 = -0.002
    eps_cu2 = -0.0035

    # Celik
    gamma_s = 1.15
    fyd = fyk / gamma_s
    Es = 20000  # kN/cm2
    eps_su = 0.01
    eps_y = fyd / Es

    # Radni dijagram betona
    def sig_eps_beton(eps):
        n = len(eps)
        sig = np.zeros(n)
        for i in range(n):
            if eps[i] >= 0:
                sig[i] = 0
            elif 0 > eps[i] >= eps_c2:
                sig[i] = -fcd * (1 - (1 - eps[i] / eps_c2) ** 2)
            elif eps[i] < eps_c2 and eps[i] >= eps_cu2:
                sig[i] = -fcd
            else:
                print('Loše uneto')
        return (sig)

    # Radni dijagram celika
    def sig_eps_celik(eps):
        n = len(eps)
        sig = np.zeros(n)
        for i in range(n):
            if eps[i] >= -eps_su and eps[i] < -eps_y:
                sig[i] = -fyd
            elif eps[i] >= -eps_y and eps[i] <= eps_y:
                sig[i] = eps[i] * Es
            elif eps[i] > eps_y and eps[i] <= eps_su:
                sig[i] = fyd
            else:
                print('Lose uneto')
        return (sig)

    n_step = 50

    # Domen 1
    yn_sup = -9999
    yn_inf = 0
    yn = np.linspace(yn_sup, yn_inf, n_step)
    psi = eps_su / (d - yn)
    eps_s1 = np.full(n_step, eps_su)
    eps_s2 = -psi*(yn - d1)
    sig_s1 = sig_eps_celik(eps_s1)
    sig_s2 = sig_eps_celik(eps_s2)

    Nrd_1 = As1 * sig_s1 + As2 * sig_s2
    Mrd_1 = As1 * sig_s1 * (h / 2 - d1) - As2 * sig_s2 * (h / 2 - d1)

    # DOmen 2
    yn_sup = 0.1
    yn_inf = d * eps_cu2 / (eps_su - eps_cu2)
    yn = np.linspace(yn_sup, yn_inf, n_step)
    psi = eps_su / (d - yn)
    eps_s1 = np.full(n_step, eps_su)
    eps_s2 = -psi*(yn - d1)
    sig_s1 = sig_eps_celik(eps_s1)
    sig_s2 = sig_eps_celik(eps_s2)

    Nc = np.zeros(n_step)
    Mc = np.zeros(n_step)

    for i in range(n_step):
        y = np.linspace(0, yn[i], 10)
        eps_c = -(yn[i] - y) * psi[i]
        eps_c = np.round(eps_c, 7)
        sig_c = sig_eps_beton(eps_c)
        Nc[i] = b * np.trapz(sig_c, y)
        ygc = np.nan_to_num(np.trapz(sig_c*y, y) / np.trapz(sig_c, y))
        Mc[i] = Nc[i] * (h / 2 - ygc)

    Nrd_2 = Nc + As1 * sig_s1 + As2 * sig_s2
    Mrd_2 = -Mc + As1 * sig_s1 * (h / 2 - d1) - As2 * sig_s2 * (h / 2 - d1)

    # Domen 3
    yn_sup = -d * eps_cu2 / (eps_su - eps_cu2)
    yn_inf = -d * eps_cu2 / (eps_y - eps_cu2)
    yn = np.linspace(yn_sup, yn_inf, n_step)
    psi = -eps_cu2 / yn
    eps_s1 = -psi * (yn - d)
    eps_s2 = -psi * (yn - d1)
    eps_s1 = np.round(eps_s1, 5)
    eps_s2 = np.round(eps_s2, 5)
    sig_s1 = sig_eps_celik(eps_s1)
    sig_s2 = sig_eps_celik(eps_s2)
    Nc = np.zeros(n_step)
    Mc = np.zeros(n_step)
    for i in range(n_step):
        y = np.linspace(0, yn[i], 50)
        eps_c = -(yn[i] - y) * psi[i]
        eps_c = np.round(eps_c, 7)
        sig_c = sig_eps_beton(eps_c)
        Nc[i] = b * np.trapz(sig_c, y)
        ygc = np.trapz(sig_c * y, y) / np.trapz(sig_c, y)
        Mc[i] = Nc[i] * (h / 2 - ygc)

    Nrd_3 = Nc + As1 * sig_s1 + As2 * sig_s2
    Mrd_3 = -Mc + As1 * sig_s1 * (h / 2 - d1) - As2 * sig_s2 * (h / 2 - d1)

    # Domen 4
    yn_sup = -d * eps_cu2 / (eps_y - eps_cu2)
    yn_inf = h
    yn = np.linspace(yn_sup, yn_inf, n_step)
    psi = -eps_cu2 / yn
    eps_s1 = -psi * (yn - d)
    eps_s2 = -psi * (yn - d1)
    sig_s1 = sig_eps_celik(eps_s1)
    sig_s2 = sig_eps_celik(eps_s2)
    Nc = np.zeros(n_step)
    Mc = np.zeros(n_step)
    for i in range(n_step):
        y = np.linspace(0, yn[i], 50)
        eps_c = -(yn[i] - y) * psi[i]
        eps_c = np.round(eps_c, 7)
        sig_c = sig_eps_beton(eps_c)
        Nc[i] = b * np.trapz(sig_c, y)
        ygc = np.trapz(sig_c * y, y) / np.trapz(sig_c, y)
        Mc[i] = Nc[i] * (h / 2 - ygc)

    Nrd_4 = Nc + As1 * sig_s1 + As2 * sig_s2
    Mrd_4 = -Mc + As1 * sig_s1 * (h / 2 - d1) - As2 * sig_s2 * (h / 2 - d1)

    # Domen 5
    yn_sup = h
    yn_inf = h + 9999
    yn = np.linspace(yn_sup, yn_inf, n_step)
    t = 3 / 7 * h
    psi = -eps_c2 / (yn - t)
    eps_s1 = -psi * (yn - d)
    eps_s2 = -psi * (yn - d1)
    sig_s1 = sig_eps_celik(eps_s1)
    sig_s2 = sig_eps_celik(eps_s2)
    Nc = np.zeros(n_step)
    Mc = np.zeros(n_step)
    for i in range(n_step):
        y = np.linspace(0, h, 50)
        eps_c = -(yn[i] - y) * psi[i]
        eps_c = np.round(eps_c, 7)
        sig_c = sig_eps_beton(eps_c)
        Nc[i] = b * np.trapz(sig_c, y)
        ygc = np.trapz(sig_c * y, y) / np.trapz(sig_c, y)
        Mc[i] = Nc[i] * (h / 2 - ygc)

    Nrd_5 = Nc + As1 * sig_s1 + As2 * sig_s2
    Mrd_5 = -Mc + As1 * sig_s1 * (h / 2 - d1) - As2 * sig_s2 * (h / 2 - d1)

    Mrd = np.hstack((Mrd_1 * 1E-2, Mrd_2 * 1E-2, Mrd_3 * 1E-2, Mrd_4 * 1E-2, Mrd_5 * 1E-2))
    Nrd = np.hstack((Nrd_1, Nrd_2, Nrd_3, Nrd_4, Nrd_5))

    return ((Nrd, Mrd))


def klasabetona(klasa):
    klasa = klasa.capitalize()
    if klasa[0] == 'C':
        klasa = klasa
    else:
        klasa = 'C' + klasa
    klasa = klasa.replace('/', '')
    klasa_lista = ['C1215', 'C1620', 'C2025', 'C2530', 'C3037', 'C3545', 'C4050', 'C4555', 'C5060']
    fck_lista = [1.2, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]  # kN/cm2
    fctm_lista = [0.16, 0.19, 0.22, 0.26, 0.29, 0.32, 0.35, 0.38, 0.41]  # kN/cm2
    Ecm_lista = [2700, 2900, 3000, 3100, 3300, 3400, 3500, 3600, 3700]  # kN/cm2
    tabela = pd.DataFrame({'fck': fck_lista, 'fctm': fctm_lista, 'Ecm': Ecm_lista}, index=klasa_lista)
    lista = tabela.loc[klasa]
    fck = lista.loc['fck']
    fctm = lista.loc['fctm']
    Ecm = lista.loc['Ecm']
    return(fck, fctm, Ecm)


def poduzna_armatura(h, b, Med, Ned, fck, fyk, Cnom):
    sipke_precnik = [8, 10, 12, 14, 16, 20, 22, 25]  # mm
    tabela_n = np.zeros(len(sipke_precnik))
    tabela_fi = np.zeros(len(sipke_precnik))
    tabela_As = np.zeros(len(sipke_precnik))
    tabela_ba = [''] * len(sipke_precnik)
    i = 0
    for fi in sipke_precnik:
        a1 = fi**2*math.pi/4
        s_max = b-2*Cnom-fi/10
        s_min = max(3.2 + 5 + fi/10, fi/10 + 2, fi/10 + 2*fi/10)
        n_red = math.floor(s_max/s_min)+1
        n_mog = np.arange(2, 2*n_red, 1)
        for n in n_mog:
            As1 = n*a1/100
            As2 = 1*As1/100
            if n < n_red:
                d1 = Cnom+(8+fi/2)/10
            else:
                d1 = Cnom+(8+fi/2+fi)/10
            Nrd, Mrd = nosivsotAB(b, h, d1, As1, As2, fck, fyk)
            Mrd1 = np.interp(Ned, Nrd[::-1], Mrd[::-1])
            if Mrd1 > Med:
                tabela_n[i] = int(n)
                tabela_fi[i] = int(fi)
                tabela_As[i] = np.round(n*(fi/10)**2*math.pi/4, 2)
                if n <= n_red:
                    tabela_ba[i] = 'jedan red armature'
                else:
                    tabela_ba[i] = 'dva reda armature'
                break
        i += 1

    lista_precnik = ['8', '10', '12', '14', '16', '20', '22', '25', ]
    tabela = pd.DataFrame({'broj sipki': tabela_n, 'precnik[mm]': tabela_fi,
                           'povrsina armature[cm2]': tabela_As,
                           'broj redova armature': tabela_ba},
                          index=lista_precnik)
    return tabela
