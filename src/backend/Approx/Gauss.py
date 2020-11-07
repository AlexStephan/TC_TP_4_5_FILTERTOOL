from src.backend.Filter.Filter import *
from scipy import signal
from scipy import special
import numpy as np

class Gauss(object):
    def __init__(self, filter: Filter):
        self.filter = filter
        self.type = self.filter.get_type(self.filter)
        self.order = None
        self.gdo = None
        self.z = None
        self.p = None
        self.k = None
        self.num = None
        self.den = None
        self.sos = None

    def calc_Order(self):
        if self.type is "Group Delay":
            self.gdo = self.filter.reqData[FilterData.GD.value]
            ft = self.filter.reqData[FilterData.ft.value]
            tol = self.filter.get_reqData[FilterData.tolerance.value]
            Nmin = self.filter.reqData[FilterData.Nmin.value]
            Nmax = self.filter.reqData[FilterData.Nmax.value]
            wtn = 2 * np.pi * ft * self.gdo
            order = 0
            good_enough = False
            while good_enough is False:
                order += 1
                w, gd = self.get_Group_Delay(self, order)
                for i in range(0, len(w)):
                    if w[i] >= wtn and gd[i] >= 1 - tol:
                        good_enough = True
                        break
            if Nmin > order:
                self.order = Nmin
            elif Nmax < order:
                self.order = Nmin
            else:
                self.order = order

    def get_Normalized_zpk(self):
        z, p, k = self.get_Gauss_Exp_zpk(self, self.order)
        w, gd = self.get_Group_Delay(self, self.order)
        k = 1
        for i in range(0, len(p)):
            p[i] *= gd[0]
            k *= p[i]
        return z, p, k

    def denormalize(self):
        z, p, k = self.get_Normalized_zpk(self)
        self.z, self.p, self.k = signal.lp2lp_zpk(z, p, k, wo=1 / self.gdo)

    def get_zpk(self):
        return self.z, self.p, self.k

    def check_Q(self) -> bool:
        Qmax = self.filter.reqData[FilterData.Qmax.value]
        if self.order > 1 and Qmax is not None:
            z, p, k = self.get_zpk(self)
            q_arr = []
            for pole in p:
                q = abs(abs(pole) / (2 * pole.real))
                q_arr.append(q)
            q_sys = np.max(q_arr)
            if q_sys > Qmax and self.order > 1:
                self.order = self.order - 1
                return False
            else:
                return True
        else:
            return True

    def calculate(self):
        self.calc_Order(self)
        self.denormalize(self)
        while self.check_Q(self) is False:
            self.denormalize(self)



    #########################
    #       Gauss Calc      #
    #########################

    def get_Gauss_Exp_Poly(self, n):
        gamma = 1
        coefs = []
        for i in range(0, n + 1):
            coefs.append(0)
            coefs.append(gamma ** (n - i) / special.factorial(n - i))
        return np.poly1d(coefs)

    def get_Gauss_Exp_zpk(self, n):
        z = []
        r = 1j * np.roots(self.get_Gauss_Exp_Poly(self, n))
        p = []
        for i in range(0, len(r)):
            if r[i].real < 0:
                p.append(r[i])
        k = np.polyval(self.get_Gauss_Exp_Poly(self, n), 0)
        for pole in p:
            k *= pole
        return z, p, k

    def get_Gauss_Exp_System(self, n):
        z, p, k = self.get_Gauss_Exp_zpk(self, n)
        return signal.lti(z, p, k)

    def get_Normalized_Group_Delay(self, n):
        w, mag, phase = signal.bode(self.get_Gauss_Exp_System(self, n), w=np.logspace(-3, 3, num=10000) , n=10000)
        gd = - np.diff(phase) / np.diff(w)
        gdn = gd / gd[0]
        w = w.tolist()
        gdn = gdn.tolist()
        gdn.append(gdn[len(gd) - 1])
        return w, gdn

    def get_Group_Delay(self, n):
        w, mag, pha = signal.bode(self.get_Gauss_Exp_System(self, n), w=np.logspace(-3, 3, num=10000), n=10000)
        gd = - np.diff(pha) / np.diff(w)
        w = w.tolist()
        gd = gd.tolist()
        gd.append(gd[len(gd) - 1])
        return w, gd
