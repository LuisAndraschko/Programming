import pandas as pd
import numpy as np
from math import radians as rd
from math import degrees as deg
from cmath import polar as polar
from cmath import rect as ret
from cmath import phase
from cmath import acos, cos, sin as acos, cos, sin


class AdmittanceMatrix:
    def __init__(self, df):
        self.df = df

    def create(self):
        termos = []
        cnxs = {}
        for terminais in self.df['CONEXOES']:
            for terminal in terminais:
                if terminal not in cnxs.keys() and terminal != terminais[1] and terminal != '0':

                    cnxs[terminal] = []
                    cnxs[terminal].append(int(terminais[0]))
                    cnxs[terminal].append(int(terminais[-1]))
                elif terminal in cnxs.keys() and terminais[-1] != '0':
                    cnxs[terminal].append(int(terminais[0]))
                    cnxs[terminal].append(int(terminais[-1]))
                try:
                    termos.append(int(terminal))
                except:
                    pass
        for bar in cnxs.keys():
            cnxs[bar] = list(set(cnxs[bar]))
        d = max(termos)
        Y = np.zeros((d, d), dtype=complex)#(max(termos)) = Número de barras
        for zre, zimg, cnx in zip(self.df['ZRE'], self.df['ZIMG'], self.df['CONEXOES']):
            for i in range(len(Y)):
                if cnx[0] == str(i+1) or cnx[-1] == str(i+1):
                    Y[i][i] += 1/(complex(zre, zimg))
                for j in range(len(Y)):
                    if cnx[0] == str(i+1) and cnx[-1] == str(j+1):
                        Y[i][j] -= 1/complex(zre, zimg)
                        Y[j][i] -= 1/complex(zre, zimg)
        for i in range(len(Y)):
            for j in range(len(Y)):
                Y[i][j] = complex(round(Y[i][j].real, 4), round(Y[i][j].imag, 4))

        return [Y, cnxs]


class SplitComplexMatrix:
    def __init__(self, m):
        self.m = m

    def split(self):
        n = len(self.m)
        g = np.zeros((n, n))
        b = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                g[i][j] = self.m[i][j].real
                b[i][j] = self.m[i][j].imag
        return g, b


class InjectedPower:
    def __init__(self, sg, sc):
        self.sg = sg
        self.sc = sc

    def evaluate(self):
        n = len(self.sg)
        pcalc = np.zeros(n)
        qcalc = np.zeros(n)
        for i in range(n):
            pcalc[i] = apparent_power_generated[i].real - apparent_power_consumed[i].real
            qcalc[i] = apparent_power_generated[i].imag - apparent_power_consumed[i].imag
        return pcalc, qcalc



class Bars:
    def __init__(self, df):
        self.df = df

    def info(self):
        size = np.shape(np.array(self.df["BARRAS"]))[0]
        bars = np.zeros(size, dtype=int)
        types = []
        scs = np.zeros(size, dtype=complex)
        sgs = np.zeros(size, dtype=complex)
        vs = np.zeros(size, dtype=complex)
        for bar, type, pc, qc, pg, qg, v, o in zip(self.df["BARRAS"], self.df["TIPO DE BARRAS"],
                                        self.df["Pc (pu)"], self.df["Qc (pu)"],
                                        self.df["Pg (pu)"], self.df["Qg (pu)"],
                                        self.df["|V| (pu)"], self.df["Ov (°)"]):
            i = int(bar) - 1
            bars[i] = int(bar)
            types.append(type)
            scs[i] = complex(pc, qc)
            sgs[i] = complex(pg, qg)
            vs[i] = ret(v, rd(o))

        nb = max(bars)
        npv, npq = 0, 0
        for bt in types:
            if bt == "PV":
                npv += 1
            elif bt == "PQ":
                npq += 1
        return bars, types, sgs, scs, vs, npq, npv


class PowerEquation:
    def __init__(self, bar, bar_type, connections, g, b, powerc, powerg, voltages, npq, npv):
        self.bar = bar
        self.bt = bar_type
        self.connections = connections
        self.g = g
        self.b = b
        self.powerc = powerc
        self.powerg = powerg
        self.voltages = voltages
        self.npq = npq
        self.npv = npv

    def evaluate(self):
        dim = 2 * self.npq + self.npv
        f = np.zeros(dim)
        """bar começa do 1 (barra slack) e vai até o 9, porém F tem tamanho 2NPQ+NPV e no python seu " \
           "primeiro elemento tem índice 0. Como desejo escrever as equações das barras PV e a única ac" \
           "ima da 1ª barra PV na minha base de dados é a slack, contabilizo -1 para ajustar os índices" \
           " das barras com o índice de F no python e -1 para contabilizar a barra slack, portanto o ín-" \
           "dice de F é bar - 2."""
        for bt, bar in zip(self.bt, self.bar):
            for cnx in self.connections[str(bar)]:
                i = bar - 1
                j = cnx - 1
                g = self.g[i][j].real
                b = self.b[i][j].real
                v = self.voltages[j].__abs__().real
                o = phase(self.voltages[i]).real - phase(self.voltages[j]).real
                seno, cosseno = sin(o).real, cos(o).real
                if bt == "PV":
                    f[bar - 2] -= v * (g*cosseno + b*seno)
                elif bt == "PQ":
                    f[bar - 2] -= v * (g*cosseno + b*seno)
                    f[(bar + self.npq) - 2] -= v * (g * seno - b * cosseno)
            f[bar - 2] *= self.voltages[i].__abs__()
            f[bar - 2] += self.powerg[i].real - self.powerc[i].real
            f[(bar + self.npq) - 2] *= self.voltages[i].__abs__()
            f[(bar + self.npq) - 2] += self.powerg[i].imag - self.powerc[i].imag
            f[bar - 2] = round(f[bar - 2], 6)
            f[(bar + self.npq) - 2] = round(f[(bar + self.npq) - 2], 6)

        return f



class Jacobian:
    def __init__(self, bar, bar_type, connections, g, b, voltages, npq, npv):
        self.bar = bar
        self.bt = bar_type
        self.connections = connections
        self.g = g
        self.b = b
        self.voltages = voltages
        self.npq = npq
        self.npv = npv

    def terms(self, i, j):
        gkm = self.g[i][j].real
        bkm = self.b[i][j].real
        vk = self.voltages[i].__abs__().real
        vm = self.voltages[j].__abs__().real
        o = phase(self.voltages[i]).real - phase(self.voltages[j]).real
        seno, cosseno = sin(o), cos(o)
        return gkm, bkm, vk, vm, seno.real, cosseno.real


    def h_matrix(self):
        n_columns = self.npq + self.npv
        m_rows = n_columns
        h_matrix = np.zeros((m_rows, n_columns))
        for row in range(m_rows):
            for col in range(n_columns):
                if row == col:
                    for cnx in self.connections[str(row + 2)]:
                        i = row + 1
                        j = cnx - 1
                        gkm, bkm, vk, vm, seno, cosseno = self.terms(i, j)

                        h_matrix[row][col] -= vm * (-gkm*seno + bkm*cosseno)
                    h_matrix[row][col] *= vk
                if row != col:
                    if (col + 2) in self.connections[str(row + 2)]:
                        i = row + 1
                        j = (col + 2) - 1
                        gkm, bkm, vk, vm, seno, cosseno = self.terms(i, j)
                        h_matrix[row][col] -= vm * (gkm * seno - bkm * cosseno)
                        h_matrix[row][col] *= vk
                h_matrix[row][col] = round(h_matrix[row][col], 6)
        return h_matrix

    def n_matrix(self):
        n_columns = self.npq
        m_rows = self.npq + self.npv
        n_matrix = np.zeros((m_rows, n_columns))
        for row in range(m_rows):
            for col in range(n_columns):
                i = row + 1
                j = col + (self.npv + 1)
                terminal = j + 1
                barra = i + 1
                if i == j:
                    for cnx in self.connections[str(j + 1)]:
                        gkm, bkm, vk, vm, seno, cosseno = self.terms(i, (cnx - 1))
                        if cnx - 1 == i:
                            n_matrix[row][col] -= 2 * vk * gkm
                        else:
                            n_matrix[row][col] -= vm * (gkm * cosseno - bkm * seno)

                elif terminal in self.connections[str(barra)]:
                    gkm, bkm, vk, vm, seno, cosseno = self.terms(i, (terminal - 1))
                    n_matrix[row][col] -= vm * (gkm * cosseno - bkm * seno)
                n_matrix[row][col] = round(n_matrix[row][col], 6)

        return n_matrix

    def m_matrix(self):
        n_columns = self.npq + self.npv
        m_rows = self.npq
        m_matrix = np.zeros((m_rows, n_columns))
        for row in range(m_rows):
            for col in range(n_columns):
                i = row + (self.npv + 1)
                jgb = col + 1
                terminal = jgb + 1
                barra = i + 1
                if i == jgb:
                    for cnx in self.connections[str(terminal)]:
                        gkm, bkm, vk, vm, seno, cosseno = self.terms(i, (cnx - 1))
                        m_matrix[row][col] -= vm * (gkm * cosseno + bkm * seno)
                    m_matrix[row][col] *= vk

                elif terminal in self.connections[str(barra)]:
                    gkm, bkm, vk, vm, seno, cosseno = self.terms(i, (terminal - 1))
                    m_matrix[row][col] += vk * vm * (gkm * cosseno + bkm * seno)
                m_matrix[row][col] = round(m_matrix[row][col], 6)
        return m_matrix

    def l_matrix(self):
        n_columns = self.npq
        m_rows = self.npq
        l_matrix = np.zeros((m_rows, n_columns))
        for row in range(m_rows):
            for col in range(n_columns):
                i = row + (self.npv + 1)
                jgb = col + (self.npv + 1)
                terminal = jgb + 1
                barra = i + 1
                if i == jgb:
                    for cnx in self.connections[str(terminal)]:
                        gkm, bkm, vk, vm, seno, cosseno = self.terms(i, (cnx - 1))
                        if (cnx - 1) == jgb:
                            l_matrix[row][col] += 2 * vk * bkm
                        else:
                            l_matrix[row][col] -= vm * (gkm * seno - bkm * cosseno)
                elif terminal in self.connections[str(barra)]:
                    gkm, bkm, vk, vm, seno, cosseno = self.terms(i, (terminal - 1))
                    l_matrix[row][col] -= vk * (gkm * seno - bkm * cosseno)
                l_matrix[row][col] = round(l_matrix[row][col], 6)
        return l_matrix

    def matrix(self):
        hm = self.h_matrix()
        nm = self.n_matrix()
        mm = self.m_matrix()
        lm = self.l_matrix()
        jacobian = np.bmat([[hm, nm], [mm, lm]])
        return jacobian


class ReformatToCompare:
    def __init__(self, bar, bar_type, npv, npq, voltages):
        self.bar = bar
        self.bt = bar_type
        self.npq = npq
        self.npv = npv
        self.v = voltages

    def compare(self):
        x = np.zeros(self.npv + (2 * self.npq))
        for bt, bar in zip(self.bt, self.bar):
            i = int(bar) - 1
            o = phase(self.v[i])
            v = self.v[i].__abs__()
            if bt != "Slack":
                x[i - 1] = o
            if bt == "PQ":
                x[(i - 1) + self.npq] = v
        return x


class ReformatToIter:
    def __init__(self, bar, bar_type, xk1, npv, npq, voltages):
        self.bar = bar
        self.bt = bar_type
        self.xk1 = xk1
        self.npq = npq
        self.npv = npv
        self.v = voltages

    def iter(self):
        size = np.shape(np.array(self.bar))[0]
        nv = np.zeros(size, dtype=complex)
        for bt, bar in zip(self.bt, self.bar):
            i = int(bar) - 1
            ipqv = int(bar) - 1
            ipq = ((int(bar) - 2) + self.npq)
            if bt == "Slack":
                nv[i] = self.v[ipqv]
            elif bt == "PV":
                nv[i] = ret(self.v[ipqv].__abs__(), self.xk1[ipqv - 1])
            elif bt == "PQ":
                nv[i] = ret(self.xk1[ipq], self.xk1[ipqv - 1])
        return nv


class NewtonRaphson:
    def __init__(self, bar, bar_type, voltages, npv, npq, F, J):
        self.bar = bar
        self.bt = bar_type
        self.v = voltages
        self.npq = npq
        self.npv = npv
        self.f = F
        self.j = J

    def iterate(self):
        self.f = np.array(self.f)
        self.j = np.array(self.j)
        #print(np.linalg.det(self.j))
        x = ReformatToCompare(self.bar, self.bt, self.npv, self.npq, self.v).compare()
        xk1 = x - np.dot(np.linalg.inv(self.j), self.f)
        return xk1


df = pd.read_excel(r"C:\Users\luis_\PycharmProjects\pythonProject\ISEP\ISEP.xlsx")
df = df.fillna(0)
tol = 0.001
exp = {}
fnzrows = (df["BARRAS"] != 0) # Filtering not zero rows
#print(df.loc[fnzrows, ["BARRAS", "TIPO DE BARRAS"]])

bars, bar_types, apparent_power_generated, apparent_power_consumed, phasor_voltages, npq, npv = Bars(
    df.loc[fnzrows, ["BARRAS", "TIPO DE BARRAS",
                     "Pc (pu)", "Qc (pu)", "Pg (pu)", "Qg (pu)",
                     "|V| (pu)", "Ov (°)"]]).info()
am, cnxs = AdmittanceMatrix(df).create()[0], AdmittanceMatrix(df).create()[1]
gb_matrices = SplitComplexMatrix(am).split()
f_x = PowerEquation(bars, bar_types, cnxs, gb_matrices[0], gb_matrices[1],
                    apparent_power_consumed, apparent_power_generated, phasor_voltages,
                    npq, npv).evaluate()
j = Jacobian(bars, bar_types, cnxs, gb_matrices[0], gb_matrices[1], phasor_voltages, npq, npv).matrix()
x = ReformatToCompare(bars, bar_types, npv, npq, phasor_voltages).compare()
xk1 = NewtonRaphson(bars, bar_types, phasor_voltages, npv, npq, f_x, j).iterate()
bool_array = (xk1 - x).__abs__() < 1
xk1 = pd.DataFrame(bool_array)

teste = ReformatToCompare(bars, bar_types, npv, npq, phasor_voltages).compare()
teste2 = ReformatToIter(bars, bar_types, teste, npv, npq, phasor_voltages).iter()
teste = pd.DataFrame(teste)
teste2 = pd.DataFrame(teste2)
print(teste.to_string())
print(teste2.to_string())
"""
Yd = pd.DataFrame(Y)
print(Yd.to_string())
"""