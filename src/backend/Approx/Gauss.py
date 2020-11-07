from src.backend.Filter.Filter import *
from scipy import signal
from scipy import special
import numpy as np

class Gauss(object):
    def __init__(self, filter: Filter):
        self.filter = filter
        self.type = self.filter.get_type()
        self.order = None
        self.gdo = None
        self.z = None
        self.p = None
        self.k = None
        self.num = None
        self.den = None
        self.sos = None
        self.w_bode = None
        self.mag = None
        self.pha = None
        self.A = None
        self.w_tf = None
        self.h = None
        self.wgd = None
        self.GroupDelay = None
        self.calculate()

    def calc_Order(self):
        if self.type == "Group Delay":
            self.gdo = self.filter.reqData[FilterData.GD]
            ft = self.filter.reqData[FilterData.ft]
            tol = self.filter.reqData[FilterData.tolerance] / 100
            Nmin = self.filter.reqData[FilterData.Nmin]
            Nmax = self.filter.reqData[FilterData.Nmax]
            wtn = 2 * np.pi * ft * self.gdo
            order = 0
            good_enough = False
            while good_enough is False:
                order += 1
                w, gd = self.get_UnNormalized_Group_Delay(order)
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
        z, p, k = self.get_Gauss_Exp_zpk(self.order)
        w, gd = self.get_UnNormalized_Group_Delay(self.order)
        k = 1
        for i in range(0, len(p)):
            p[i] *= gd[0] / 100
            k *= p[i]
        return z, p, k

    def denormalize(self):
        z, p, k = self.get_Normalized_zpk()
        self.z, self.p, self.k = signal.lp2lp_zpk(z, p, k, wo=1 / self.gdo)

    def check_Q(self) -> bool:
        Qmax = self.filter.reqData[FilterData.Qmax]
        if self.order > 1 and Qmax is not None:
            z, p, k = self.get_zpk()
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

    def get_Gain(self):
        return self.filter.reqData[FilterData.gain]

    def calc_TransFunc(self):
        val, msg = self.filter.validate()
        if val is False:
            return msg
        z, p, k = self.get_zpk()
        if self.type is "Group Delay":
            sys = signal.ZerosPolesGain(z, p, k)
            self.w_tf, self.h = signal.TransferFunction(sys)
        else:
            message = "Error: Enter Filter Type."
            return message

    def calc_MagAndPhase(self):  # return angular frequency, Mag and Phase
        val, msg = self.filter.validate()
        if val is False:
            return msg
        z, p, k = self.get_zpk()
        if self.type is "Group Delay":
            sys = signal.ZerosPolesGain(z, p, k)
            self.w_bode, self.mag, self.pha = signal.bode(sys)
        else:
            message = "Error: Enter Filter Type."
            return message

    def calc_Attenuation(self):
        A = self.mag
        for i in range(0, len(A)):
            A[i] = 1 / A[i]
        self.A = A

    def calc_Group_Delay(self):
        w, mag, pha = self.get_MagAndPhaseWithoutGain()
        gd = - np.diff(pha) / np.diff(w)
        gd = gd.tolist()
        gd.append(gd[len(gd) - 1])
        self.GroupDelay = gd
        w = w.tolist()
        self.wgd = w

    def calculate(self):
        self.calc_Order()
        self.denormalize()
        while self.check_Q() is False:
            self.denormalize()
        self.calc_TransFunc()
        self.calc_MagAndPhase()
        self.calc_Group_Delay()

    #####################
    #       ALEX        #
    #####################

    def get_Attenuation(self):
        return self.w_bode, self.A

    def get_Group_Delay(self):
        return self.wgd, self.GroupDelay

    def get_Order(self):
        return self.order

    def get_TransFuncWithGain(self):
        hg = self.h
        gain = self.get_Gain()
        for i in range(0, len(hg)):
            hg[i] = hg[i] * 10 ** (gain / 20)
        return self.w_tf, hg

    def get_TransFuncWithoutGain(self):
        return self.w_tf, self.h

    def get_MagAndPhaseWithGain(self):
        magg = self.mag
        gain = self.get_Gain()
        for i in range(0, len(magg)):
            magg[i] = magg[i] + gain
        return self.w_bode, magg, self.pha

    def get_MagAndPhaseWithoutGain(self):
        return self.w_bode, self.mag, self.pha

    def get_zpk(self):
        return self.z, self.p, self.k

    def get_zpGk(self):
        Gk = self.k * 10 ** (self.get_Gain() / 20)
        return self.z, self.p, Gk

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
        r = 1j * np.roots(self.get_Gauss_Exp_Poly(n))
        p = []
        for i in range(0, len(r)):
            if r[i].real < 0:
                p.append(r[i])
        k = np.polyval(self.get_Gauss_Exp_Poly(n), 0)
        for pole in p:
            k *= pole
        return z, p, k

    def get_Gauss_Exp_System(self, n):
        z, p, k = self.get_Gauss_Exp_zpk(n)
        return signal.lti(z, p, k)

    def get_Normalized_Group_Delay(self, n):
        w, mag, phase = signal.bode(self.get_Gauss_Exp_System(n), w=np.logspace(-3, 3, num=10000) , n=10000)
        gd = - np.diff(phase) / np.diff(w)
        gdn = gd / gd[0]
        w = w.tolist()
        gdn = gdn.tolist()
        gdn.append(gdn[len(gd) - 1])
        return w, gdn

    def get_UnNormalized_Group_Delay(self, n):
        w, mag, pha = signal.bode(self.get_Gauss_Exp_System(n), w=np.logspace(-3, 3, num=10000), n=10000)
        gd = - np.diff(pha) / np.diff(w)
        w = w.tolist()
        gd = gd.tolist()
        gd.append(gd[len(gd) - 1])
        return w, gd
