from pprint import pprint

import pandas as pd
import numpy as np
from math import radians as rd
from math import degrees as deg
from cmath import polar as polar
from cmath import rect as ret
from cmath import phase
from cmath import acos, cos, sin as acos, cos, sin


class PuConversions:
    def __init__(self, basis, elements, identifier):
        self.basis = basis
        self.elements = elements
        self.identifier = identifier
        self.pu = {identifier: []}

    def convert_to_pu(self):
        for item in self.elements:
            self.pu[self.identifier].append((item / self.basis))
        return self.pu


class BasisMono:
    def __init__(self, base_voltage, base_power):
        self.b_V = base_voltage
        self.b_S = base_power
        self.pu_I = []
        self.pu_Z = []
        self.pu_Y = []

    def base_currents(self):
        for item in self.b_V:
            self.pu_I.append(self.b_S / item)
        return self.pu_I

    def base_impedances(self):
        for item in self.b_V:
            self.pu_Z.append((item ** 2) / self.b_S)
        return self.pu_Z

    def base_admittances(self):
        for item in self.b_V:
            self.pu_Y.append(self.b_S / (item ** 2))
        return self.pu_Y


class BasisTri:
    def __init__(self, voltage, power):
        self.voltage = voltage
        self.power = power
        self.currents_pu = []
        self.impedances_pu = []
        self.admittances_pu = []

    def base_currents(self):
        for item in self.voltage:
            self.currents_pu.append(self.power / (item * (3 ** 0.5)))
        return self.currents_pu

    def base_impedances(self):
        return BasisMono(self.voltage, self.power).base_impedances()

    def base_admittances(self):
        return BasisMono(self.voltage, self.power).base_admittances()


class PhaseLineConversionEE:
    def __init__(self, element):
        self.element = element
        self.converted = []

    def convert_voltage(self):
        if self.element[0] == 'fase':
            for element in self.element[1]:
                self.converted.append((element * (3 ** 0.5)) * ret(1, rd(30)))
        else:
            for element in self.element[1]:
                self.converted.append((element / (3 ** 0.5)) * ret(1, rd(-30)))
        return self.converted

    def convert_current(self):
        return self.element[1]


class PhaseLineConversionET:
    def __init__(self, element):
        self.element = element
        self.converted = []

    def convert_voltage(self):
        return PhaseLineConversionEE(self.element).convert_voltage()

    def convert_current(self):
        if self.element[0] == 'fase':
            for element in self.element[1]:
                self.converted.append((element * (3 ** 0.5)) * ret(1, rd(-30)))
        else:
            for element in self.element[1]:
                self.converted.append((element / (3 ** 0.5)) * ret(1, rd(30)))
        return self.converted


class PhaseLineConversionTT:
    def __init__(self, element):
        self.element = element
        self.converted = []

    def convert_voltage(self):
        return self.element

    def convert_current(self):
        return PhaseLineConversionET(self.element).convert_current()


class PhaseLineConversionTE:
    def __init__(self, element):
        self.element = element
        self.converted = []

    def convert_voltage(self):
        return self.element

    def convert_current(self):
        return self.element


class SystemBasis:  # All transformer relations must be given like: {'trasformer 1': [V1, V2]}
    def __init__(self, base_voltage, transformers_relations):
        self.base_voltage = base_voltage
        self.transformers_relations = transformers_relations
        self.sys_voltages = {'Vbs': [self.base_voltage]}

    def source_to_load(self):
        for relations in self.transformers_relations.values():
            self.sys_voltages['Vbs'].append(self.sys_voltages['Vbs'][-1] * (relations[1] / relations[0]))
        return self.sys_voltages

    def load_to_source(self):
        for relations in reversed(self.transformers_relations.values()):
            self.sys_voltages['Vbs'].append(self.sys_voltages['Vbs'][-1] * (relations[0] / relations[1]))
        return self.sys_voltages


class ImpedanceBCG:
    def __init__(self, pu_impedance, nominal_voltage, nominal_power, base_voltage, base_power):
        self.pu_Z = pu_impedance
        self.n_V = nominal_voltage
        self.n_S = nominal_power
        self.b_V = base_voltage
        self.b_S = base_power
        self.r_Z, self.pu_Znew = 0, 0

    def generator(self):
        self.r_Z = (self.pu_Z * (self.n_V ** 2)) / self.n_S
        self.pu_Znew = (self.pu_Z * (self.n_V ** 2) * self.b_S) / (self.n_S * (self.b_V ** 2))
        return {'puZg': self.pu_Znew, 'ohmZg': self.r_Z}


class ImpedanceBCT:  # To get the right answer i need the correct voltage of the transformer.
    def __init__(self, pu_impedance, nominal_voltages, nominal_power, base_voltages, base_power):
        self.pu_Z = pu_impedance
        self.n_Vs = nominal_voltages
        self.n_S = nominal_power
        self.b_Vs = base_voltages
        self.b_S = base_power
        self.r_Z = []
        self.pu_Znew = []

    def transformer(self):
        for n_Vs, b_Vs in zip(self.n_Vs, self.b_Vs):
            self.r_Z.append((self.pu_Z * (n_Vs ** 2)) / self.n_S)
            self.pu_Znew.append((self.pu_Z * (n_Vs ** 2) * self.b_S) / (self.n_S * (b_Vs ** 2)))
        return {'puZt': self.pu_Znew, 'ohmZt': self.r_Z}


class DistributedParametersConversion:
    def __init__(self, d_parameter, tl_length):
        self.d_parameter = d_parameter
        self.tl_length = tl_length
        self.c_parameter = {}

    def inductance(self):
        for dp, lt in zip(self.d_parameter, self.tl_length):
            self.c_parameter['Zl'] = dp * lt
        return self.c_parameter

    def capacitance(self):
        for dp, lt in zip(self.d_parameter, self.tl_length):
            self.c_parameter['Zc'] = dp / lt
        return self.c_parameter


class BankTransformers:
    def __init__(self, mono_voltages, mono_power, pu_impedance, connections):
        self.mono_voltages = mono_voltages
        self.mono_power = mono_power
        self.pu_impedance = pu_impedance
        self.connections = connections
        self.three_phase_specs = {}

    def three_phase(self):
        if self.connections == 'EE':
            self.three_phase_specs['p_V'] = polar(PhaseLineConversionEE(
                ['fase', [ret(self.mono_voltages[0], rd(0))]]).convert_voltage()[0])[0]
            self.three_phase_specs['s_V'] = polar(PhaseLineConversionEE(
                ['fase', [ret(self.mono_voltages[1], rd(0))]]).convert_voltage()[0])[0]
            self.three_phase_specs['3S'] = self.mono_power * 3
            self.three_phase_specs['pu_3S'] = self.pu_impedance
        elif self.connections == 'ET':
            self.three_phase_specs['p_V'] = polar(PhaseLineConversionEE(
                ['fase', [ret(self.mono_voltages[0], rd(0))]]).convert_voltage()[0])[0]
            self.three_phase_specs['s_V'] = self.mono_voltages[1]
            self.three_phase_specs['3S'] = self.mono_power * 3
            self.three_phase_specs['pu_3Z'] = self.pu_impedance
        elif self.connections == 'TT':
            self.three_phase_specs['p_V'] = self.mono_voltages[0]
            self.three_phase_specs['s_V'] = self.mono_voltages[1]
            self.three_phase_specs['3S'] = self.mono_power * 3
            self.three_phase_specs['pu_3Z'] = self.pu_impedance
        elif self.connections == 'TE':
            self.three_phase_specs['p_V'] = self.mono_voltages[0]
            self.three_phase_specs['s_V'] = polar(PhaseLineConversionEE(
                ['fase', [ret(self.mono_voltages[1], rd(0))]]).convert_voltage()[0])[0]
            self.three_phase_specs['3S'] = self.mono_power * 3
            self.three_phase_specs['pu_3Z'] = self.pu_impedance
        else: print('Please, check the connection type given')
        return self.three_phase_specs


class PuConstantPowerLoad:
    def __init__(self, power_factor, apparent_power, base_power, identifier):
        self.fp = power_factor
        self.aS = apparent_power
        self.bS = base_power
        self.identifier = identifier
        self.cpower_pu = {identifier: []}

    def convert(self):
        if self.fp[1] == 'atraso':
            self.cpower_pu[self.identifier].append(complex((self.aS * self.fp[0]), ((self.aS ** 2) -
                                                            ((self.aS * self.fp[0]) ** 2)) ** 0.5) /
                                                            self.bS)
        elif self.fp[1] == 'avanço':
            self.cpower_pu[self.identifier].append(complex((self.aS * self.fp[0]), -((self.aS ** 2) -
                                                            ((self.aS * self.fp[0]) ** 2)) ** 0.5) /
                                                            self.bS)
        return self.cpower_pu


class SeriesImpedance:
    def __init__(self, impedances, identifier):
        self.Zs = impedances
        self.identifier = identifier
        self.Z_equivalent = {identifier: 0}

    def calculate(self):
        for Z in self.Zs:
            self.Z_equivalent[self.identifier] += Z
        return self.Z_equivalent


class ParallelImpedance:
    def __init__(self, impedances, identifier):
        self.Zs = impedances
        self.identifier = identifier
        self.Z_equivalent = {identifier: 0}

    def calculate(self):
        for Z in self.Zs:
            self.Z_equivalent[self.identifier] += 1 / Z
        self.Z_equivalent[self.identifier] = self.Z_equivalent[self.identifier] ** (- 1)
        return self.Z_equivalent


class SourceTransform:
    def __init__(self, source_voltages, series_impedances):
        """
        :param source_voltages: [V1, V2,  ..., Vn]
        :param series_impedance: [Z1, Z2, ..., Zn]
        """
        self.sVs = source_voltages
        self.sZs = series_impedances
        self.newSources = {'I': []}

    def calculate(self):
        """
        :return: {'I': [I1, I2, ..., In]}
        """
        for V, Z in zip(self.sVs, self.sZs):
            self.newSources['I'].append(V / Z)
        return self.newSources


class CreateVector:
    def __init__(self, N):
        self.n = N

    def create(self):
        v = []
        for size in range(self.n):
            v.append(0)
        return v


class CreateMatrix:
    def __init__(self, m, n):
        self.m = m
        self.n = n

    def create(self):
        m_rows = []
        matrix = []
        for r in range(self.m):
            for c in range(self.n):
                m_rows.append(0)
            matrix.append(m_rows)
            m_rows = []
        return matrix


class LU:
    def __init__(self, M):
        self.m = M

    def decompose(self):
        n = len(self.m)
        L = CreateMatrix(n, n).create()
        U = self.m
        for i in range(n):
            for j in range(n):
                if j == i:
                    L[i][i] = 1
                elif j > i:
                    L[j][i] = U[j][i]/U[i][i]
                    for k in range(n):
                        U[j][k] = U[j][k] - L[j][i] * U[i][k]

        return L, U


class LULinearSystem:
    def __init__(self, L, U, vr):
        self.l = L
        self.u = U
        self.vr = vr

    def solve(self):
        n = len(self.l)
        y = CreateVector(n).create()
        x = y[:]
        ant = y[:]
        for i in range(n):
            for j in range(n):
                y[i] += - self.l[i][j] * ant[j] + (self.vr[i]/n)
            ant = y[:]
        ant = CreateVector(n).create()
        for i in range(n):
            for j in range(n):
                x[(n - 1) - i] += (- self.u[(n - 1) - i][j] * ant[j] + (y[(n - 1) - i] / n))\
                                  / self.u[(n - 1) - i][(n - 1) - i]
            ant = x[:]
        return y, x


'''L, U = LU([[3, 2, 4], [1, 1, 2], [4, 3, -2]]).decompose()
print(L)
print(U)
y, x = LULinearSystem(L, U, [1, 2, 3]).solve()
print(y, x)
'''
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
        pkcalc = np.zeros(self.npq + self.npv)
        qkcalc = np.zeros(self.npq + self.npv)

        for item in self.connections.items():
            if int(item[0]) in item[1]:
                del (item[1][item[1].index(int(item[0]))])

        for bt, bar in zip(self.bt, self.bar):
            if bt != "Slack":
                for cnx in self.connections[str(bar)]:
                    i = bar - 1
                    j = cnx - 1
                    gkm, gkk, bkm, bkk, vk, vm, seno, cosseno = self.terms(i, j)

                    if bt == "PV":
                        f[bar - 2] -= vm * (gkm * cosseno + bkm * seno)
                        pkcalc[bar - 2] += vm * (gkm*cosseno + bkm*seno)
                        qkcalc[bar - 2] += vm * (gkm*seno - bkm*cosseno)
                    elif bt == "PQ":
                        f[bar - 2] -= vm * (gkm*cosseno + bkm*seno)
                        f[(bar + self.npq) - 2] -= vm * (gkm * seno - bkm * cosseno)
                        pkcalc[bar - 2] += vm * (gkm * cosseno + bkm * seno)
                        qkcalc[bar - 2] += vm * (gkm * seno - bkm * cosseno)

                deltaPk = self.powerg[i].real - self.powerc[i].real
                deltaQk = self.powerg[i].imag - self.powerc[i].imag

                f[bar - 2] *= vk
                pkcalc[bar - 2] *= vk
                qkcalc[bar - 2] *= vk

                f[bar - 2] -= pow(vk, 2) * gkk
                pkcalc[bar - 2] += pow(vk, 2) * gkk
                qkcalc[bar - 2] -= pow(vk, 2) * bkk

                f[bar - 2] += deltaPk

                if bt == "PQ":
                    f[(bar + self.npq) - 2] *= vk
                    f[(bar + self.npq) - 2] += pow(vk, 2) * bkk
                    f[(bar + self.npq) - 2] += deltaQk

                f[bar - 2] = round(f[bar - 2], 6)
                f[(bar + self.npq) - 2] = round(f[(bar + self.npq) - 2], 6)
        return f, pkcalc, qkcalc

    def terms(self, i, j):
        gkm = self.g[i][j].real
        gkk = self.g[i][i].real
        bkm = self.b[i][j].real
        bkk = self.b[i][i].real
        vk = self.voltages[i].__abs__().real
        vm = self.voltages[j].__abs__().real
        o = phase(self.voltages[i]).real - phase(self.voltages[j]).real
        seno, cosseno = sin(o), cos(o)
        return gkm, gkk, bkm, bkk, vk, vm, seno.real, cosseno.real




class Jacobian:
    def __init__(self, bar, bar_type, connections, g, b, voltages, npq, npv, pkcalc, qkcalc):
        self.bar = bar
        self.bt = bar_type
        self.connections = connections
        self.g = g
        self.b = b
        self.voltages = voltages
        self.npq = npq
        self.npv = npv
        self.pkcalc = pkcalc
        self.qkcalc = qkcalc

    def terms(self, i, j):
        gkm = self.g[i][j].real
        gkk = self.g[i][i].real
        bkm = self.b[i][j].real
        bkk = self.b[i][i].real
        vk = self.voltages[i].__abs__().real
        vm = self.voltages[j].__abs__().real
        o = phase(self.voltages[i]).real - phase(self.voltages[j]).real
        seno, cosseno = sin(o), cos(o)
        return gkm, gkk, bkm, bkk, vk, vm, seno.real, cosseno.real

    def h_matrix(self):
        n_columns = self.npq + self.npv
        m_rows = n_columns
        h_matrix = np.zeros((m_rows, n_columns))
        for row in range(m_rows):
            for col in range(n_columns):
                if row == col:
                    i = row + 1
                    gkm, gkk, bkm, bkk, vk, vm, seno, cosseno = self.terms(i, i)
                    h_matrix[row][col] = - self.qkcalc[row] - pow(vk, 2) * bkk
                if row != col:
                    if (col + 2) in self.connections[str(row + 2)]:
                        i = row + 1
                        j = (col + 2) - 1
                        gkm, gkk, bkm, bkk, vk, vm, seno, cosseno = self.terms(i, j)
                        h_matrix[row][col] = vk * vm * (gkm * seno - bkm * cosseno)
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
                    gkm, gkk, bkm, bkk, vk, vm, seno, cosseno = self.terms(i, i)
                    n_matrix[row][col] = (self.pkcalc[row] + (pow(vk, 2) * gkk)) / vk
                elif terminal in self.connections[str(barra)]:
                    gkm, gkk, bkm, bkk, vk, vm, seno, cosseno = self.terms(i, (terminal - 1))
                    n_matrix[row][col] = vm * (gkm * cosseno + bkm * seno)
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
                    gkm, gkk, bkm, bkk, vk, vm, seno, cosseno = self.terms(i, i)
                    m_matrix[row][col] = self.pkcalc[row] - (pow(vk, 2) * gkk)
                elif terminal in self.connections[str(barra)]:
                    gkm, gkk, bkm, bkk, vk, vm, seno, cosseno = self.terms(i, (terminal - 1))
                    m_matrix[row][col] = - (vk * vm * (gkm * cosseno + bkm * seno))
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
                    gkm, gkk, bkm, bkk, vk, vm, seno, cosseno = self.terms(i, i)
                    l_matrix[row][col] = (self.qkcalc[row] - (pow(vk, 2) * bkk)) / vk
                elif terminal in self.connections[str(barra)]:
                    gkm, gkk, bkm, bkk, vk, vm, seno, cosseno = self.terms(i, (terminal - 1))
                    l_matrix[row][col] = vk * (gkm * seno - bkm * cosseno)
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
        # self.xk1 é sempre o mesmo valor aqui, antes de passar pelo reformat to compare
        if np.shape(self.xk1)[0] != 14:
            self.xk1 = ReformatToCompare(self.bar, self.bt, self.npv, self.npq, self.xk1).compare()
        p = pd.DataFrame(self.xk1)
        # print(p.to_string())
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

        # Arredonda os valores para não ter diferença entre os valores em uma casa decimal sem importância
        return np.around(nv, 6)

class NewtonRaphson:
    def __init__(self, bar, bar_type, voltages, npv, npq, F, J):
        self.bar = bar
        self.bt = bar_type
        self.v = voltages  # muda a cada iteração
        self.npq = npq
        self.npv = npv
        self.f = F  # muda a cada iteração
        self.j = J  # muda a cada iteração

    def iterate(self):
        self.f = np.array(self.f)
        self.j = np.array(self.j)
        x = ReformatToCompare(self.bar, self.bt, self.npv, self.npq, self.v).compare()
        try:
            xk1 = x + np.dot(np.linalg.inv(self.j), self.f)
        except:
            print(f'Determinante = {np.linalg.det(self.j)}')
        return xk1


df = pd.read_excel(r"ISEP.xlsx")
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

##########################################   X0  ################################################################
f_x0, pkcalc0, qkcalc0 = PowerEquation(bars, bar_types, cnxs, gb_matrices[0], gb_matrices[1],
                                       apparent_power_consumed, apparent_power_generated, phasor_voltages,
                                       npq, npv).evaluate()
j0 = Jacobian(bars, bar_types, cnxs, gb_matrices[0], gb_matrices[1], phasor_voltages, npq, npv, pkcalc0, qkcalc0).matrix()

# Um jeito melhor era reformatar e já retornar no formato certo, pra não precisar ficar reformatando toda hora.
xk1 = [NewtonRaphson(bars, bar_types, phasor_voltages, npv, npq, f_x0, j0).iterate()]
xk1[0] = ReformatToIter(bars, bar_types, xk1[0], npv, npq, phasor_voltages).iter()
tolerance_bool = f_x0 < tol


for i in range(10):
    f_x1, pkcalc1, qkcalc1 = PowerEquation(bars, bar_types, cnxs, gb_matrices[0], gb_matrices[1],
                                           apparent_power_consumed, apparent_power_generated, xk1[i],
                                           npq, npv).evaluate()
    if np.all(f_x1 < tol):
        print(i)
        print(f_x1)
        break
    j1 = Jacobian(bars, bar_types, cnxs, gb_matrices[0], gb_matrices[1], xk1[i], npq, npv, pkcalc1,
                  qkcalc1).matrix()
    # Calcula o primeiro valor para xk1[i+1]
    xk1.append(NewtonRaphson(bars, bar_types, xk1[i], npv, npq, f_x1, j1).iterate())
    # Reformata o valor para xk1[i+1]
    xk1[i + 1] = ReformatToIter(bars, bar_types, xk1[i+1], npv, npq, phasor_voltages).iter()
    assert not np.array_equal(xk1[i], xk1[i + 1])

# while not np.all(tolerance_bool):
#     f_x1, pkcalc1, qkcalc1 = PowerEquation(bars, bar_types, cnxs, gb_matrices[0], gb_matrices[1],
#                                            apparent_power_consumed, apparent_power_generated, xk1[i],
#                                            npq, npv).evaluate()
#     j1 = Jacobian(bars, bar_types, cnxs, gb_matrices[0], gb_matrices[1], xk1[i], npq, npv, pkcalc1,
#                   qkcalc1).matrix()
#     # Calcula o primeiro valor para xk1[i+1]
#     xk1.append(NewtonRaphson(bars, bar_types, xk1[i], npv, npq, f_x1, j1).iterate())
#     # Reformata o valor para xk1[i+1]
#     xk1[i+1] = ReformatToIter(bars, bar_types, xk1[i], npv, npq, phasor_voltages).iter()
#     # assert not np.array_equal(xk1[i], xk1[i + 1])
#     f_x1, pkcalc1, qkcalc1 = PowerEquation(bars, bar_types, cnxs, gb_matrices[0], gb_matrices[1],
#                                            apparent_power_consumed, apparent_power_generated, xk1[i+1],
#                                            npq, npv).evaluate()
#     i += 1
#     tolerance_bool = f_x1 < tol

#xk2p = pd.DataFrame(xk1[0])
#print(xk2p.to_string())
#f_xp = pd.DataFrame(f_x0)
#print(f_xp.to_string())



"""
##########################################   X1  ################################################################
f_x1, pkcalc1, qkcalc1 = PowerEquation(bars, bar_types, cnxs, gb_matrices[0], gb_matrices[1],
                                    apparent_power_consumed, apparent_power_generated, xk1,
                                    npq, npv).evaluate()
j1 = Jacobian(bars, bar_types, cnxs, gb_matrices[0], gb_matrices[1], xk1, npq, npv, pkcalc1,
                 qkcalc1).matrix()
xk1 = ReformatToCompare(bars, bar_types, npv, npq, xk1).compare()
xk2 = NewtonRaphson(bars, bar_types, xk1, npv, npq, f_x1, j1).iterate()
#xk2 = ReformatToIter(bars, bar_types, xk2, npv, npq, phasor_voltages).iter()
tolerance_bool = f_x1 < tol
"""







"""
while not np.all(tolerance_bool):
    j = Jacobian(bars, bar_types, cnxs, gb_matrices[0], gb_matrices[1], xk1[0], npq, npv, pkcalc,
                 qkcalc).matrix()
    xk1.append(NewtonRaphson(bars, bar_types, xk1[0], npv, npq, f_x, j).iterate())
    xk1[i] = ReformatToCompare(bars, bar_types, npv, npq, phasor_voltages).compare()
    f_x, pkcalc, qkcalc = PowerEquation(bars, bar_types, cnxs, gb_matrices[0], gb_matrices[1],
                                        apparent_power_consumed, apparent_power_generated, xk1[i],
                                        npq, npv).evaluate()
    i += 1
    tolerance_bool = f_x < tol
    print(np.all(tolerance_bool))
print(f'Finalizado na {i} ª iteração')"""
#bool_array = (xk1 - x).__abs__() < 1
#bool_array_data = pd.DataFrame(bool_array)
#print(xk1.to_string())




"""expy = pd.DataFrame(am)
expj = pd.DataFrame(j)
expf = pd.DataFrame(f_x)
expx = pd.DataFrame(x)
expxk1 = pd.DataFrame(xk1)

with pd.ExcelWriter('Output.xlsx') as writer:
    expy.to_excel(writer, sheet_name='Matriz de Admitância')
    expj.to_excel(writer, sheet_name='Jacobiana.0')
    expf.to_excel(writer, sheet_name='F_x.0')
    expx.to_excel(writer, sheet_name='x.0')
    expxk1.to_excel(writer, sheet_name='xk1.0')"""
#print(xk1p.to_string())
#print(tolerance_boolp.to_string())
#print(f_xp.to_string())
#print(pkcalcp.to_string())
#print(qkcalcp.to_string())
#print(jp.to_string())
"""
L com sinal trocado
N com sinal trocado
M 

f_x = pd.DataFrame(f_x)
print(f_x.to_string())
pkcalc = pd.DataFrame(pkcalc)
print(pkcalc.to_string())
qkcalc = pd.DataFrame(qkcalc)
print(qkcalc.to_string())
"""








"""################################################################################################
expy = pd.DataFrame(am)
expj = pd.DataFrame(j)
expf = pd.DataFrame(f_x)
expx = pd.DataFrame(x)
expxk1 = pd.DataFrame(xk1)

with pd.ExcelWriter('Output.xlsx') as writer:
    expy.to_excel(writer, sheet_name='Matriz de Admitância')
    expj.to_excel(writer, sheet_name='Jacobiana.0')
    expf.to_excel(writer, sheet_name='F_x.0')
    expx.to_excel(writer, sheet_name='x.0')
    expxk1.to_excel(writer, sheet_name='xk1.0')

voltages = ReformatToIter(bars, bar_types, xk1, npv, npq, phasor_voltages).iter()
f_x = PowerEquation(bars, bar_types, cnxs, gb_matrices[0], gb_matrices[1],
                    apparent_power_consumed, apparent_power_generated, voltages,
                    npq, npv).evaluate()
j = Jacobian(bars, bar_types, cnxs, gb_matrices[0], gb_matrices[1], voltages, npq, npv).matrix()
x = ReformatToCompare(bars, bar_types, npv, npq, voltages).compare()
xk1 = NewtonRaphson(bars, bar_types, voltages, npv, npq, f_x, j).iterate()
bool_array = (xk1 - x).__abs__() < tol


expy = pd.DataFrame(am)
expj = pd.DataFrame(j)
expf = pd.DataFrame(f_x)
expx = pd.DataFrame(x)
expxk1 = pd.DataFrame(xk1)

with pd.ExcelWriter('Output.xlsx') as writer:
    expy.to_excel(writer, sheet_name='Matriz de Admitância')
    expj.to_excel(writer, sheet_name='Jacobiana.0')
    expf.to_excel(writer, sheet_name='F_x.0')
    expx.to_excel(writer, sheet_name='x.0')
    expxk1.to_excel(writer, sheet_name='xk1.0')



while not np.all((xk1 - x) < tol):
    voltages = ReformatToIter(bars, bar_types, xk1, npv, npq, phasor_voltages).iter()
    f_x = PowerEquation(bars, bar_types, cnxs, gb_matrices[0], gb_matrices[1],
                        apparent_power_consumed, apparent_power_generated, voltages,
                        npq, npv).evaluate()
    j = Jacobian(bars, bar_types, cnxs, gb_matrices[0], gb_matrices[1], voltages, npq, npv).matrix()
    x = ReformatToCompare(bars, bar_types, npv, npq, voltages).compare()
    xk1 = NewtonRaphson(bars, bar_types, voltages, npv, npq, f_x, j).iterate()
    bool_array = (xk1 - x).__abs__() < tol


#print(xk1)

#print(expj.to_string())
#print(j)
#print('------------------------------------')
f_x = np.array(f_x)
j = np.array(j)


for row, row_i in zip(am, range(len(am))):
    exp.update({row_i: row})
#print(exp[4][3].imag)

"""
"""
teste = ReformatToCompare(bars, bar_types, npv, npq, phasor_voltages).compare()
teste2 = ReformatToIter(bars, bar_types, teste, npv, npq, phasor_voltages).iter()
teste = pd.DataFrame(teste)
teste2 = pd.DataFrame(teste2)
print(teste.to_string())
print(teste2.to_string())


print(df['CONEXOES'])               #Getting Columns
print(df[['CONEXOES', 'BARRAS']])   #Getting Columns
print(df.iloc[1])                   #Getting Rows
print(df.iloc[[1,2]])               #Getting Rows"""

"""
FILTERING: The rows that we want are the first argument and the columns that we want are the second argument
First method: iloc function which stands for integer location
print(df.iloc[[1, 2], [1, 0]])
print(df.iloc[[1, 2], 1])
Second method: loc function which stands for integer location
print(df.loc[[1, 2], "BARRAS"])
print(df.loc[[1, 2], ["BARRAS", "CONEXOES"]])
"""






