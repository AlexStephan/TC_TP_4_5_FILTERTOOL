from src.backend.Filter.Filter import *
from scipy import signal
from scipy import special
import numpy as np

class Legendre(object):
    def __init__(self, filter: Filter):
        self.filter = filter
        self.type = self.filter.get_type(self.filter)
        self.order = None
        self.fo = None
        self.z = None
        self.p = None
        self.k = None
        self.num = None
        self.den = None
        self.sos = None
        self.wan = None
        self.epsilon2 = None
        self.w_bode = None
        self.mag = None
        self.pha = None
        self.A = None
        self.w_tf = None
        self.h = None
        self.wgd = None
        self.GroupDelay = None
        self.calculate(self)

    def calc_Order(self):
        val, msg = self.filter.validate(self.filter)
        if val is False:
            return msg

        fpMin = self.filter.reqData[FilterData.fpMin.value]
        fpMax = self.filter.reqData[FilterData.fpMax.value]
        Ap = self.filter.reqData[FilterData.Ap.value]
        if Ap is 0:
            Ap = 1e-9

        faMin = self.filter.reqData[FilterData.faMin.value]
        faMax = self.filter.reqData[FilterData.fpMax.value]
        Aa = self.filter.reqData[FilterData.Aa.value]

        Nmin = self.filter.reqData[FilterData.Nmin.value]
        Nmax = self.filter.reqData[FilterData.Nmax.value]

        self.epsilon2 = 10 ** (Ap / 10) - 1
        order = 1

        if self.type is "Low Pass":
            self.wan = faMin / fpMin
            while self.get_L_Poly_Value(self, order, self.wan) < (10 ** (Aa / 10) - 1) / self.epsilon2:
                order = order + 1
            if Nmin > order:
                self.order = Nmin
            elif Nmax < order:
                self.order = Nmax
            else:
                self.order = order
        elif self.type is "High Pass":
            self.wan = fpMin / faMin
            while self.get_L_Poly_Value(self, order, self.wan) < (10 ** (Aa / 10) - 1) / self.epsilon2:
                order = order + 1
            if Nmin > order:
                self.order = Nmin
            elif Nmax < order:
                self.order = Nmax
            else:
                self.order = order
        elif self.type is "Band Pass":
            self.wan = (faMax - faMin) / (fpMax - fpMin)
            while self.get_L_Poly_Value(self, order, self.wan) < (10 ** (Aa / 10) - 1) / self.epsilon2:
                order = order + 1
            if Nmin > order:
                self.order = Nmin
            elif Nmax < order:
                self.order = Nmax
            else:
                self.order = order
        elif self.type is "Band Reject":
            self.wan = (fpMax - fpMin) / (faMax - faMin)
            while self.get_L_Poly_Value(self, order, self.wan) < (10 ** (Aa / 10) - 1) / self.epsilon2:
                order = order + 1
            if Nmin > order:
                self.order = Nmin
            elif Nmax < order:
                self.order = Nmax
            else:
                self.order = order
        else:
            message = "Error: Enter Filter Type."
            return message

    def calc_fo(self):
        val, msg = self.filter.validate(self.filter)
        if val is False:
            return msg

        fpMin = self.filter.reqData[FilterData.fpMin.value]
        fpMax = self.filter.reqData[FilterData.fpMax.value]

        Denorm = self.filter.reqData[FilterData.Denorm.value]

        if self.type is "Low Pass":
            wo1 = self.get_L_wo(self)
            wo2 = wo1 * self.wan / self.get_L_wa(self)
            fod = 10 ** (np.log10(wo1) * (1 - Denorm / 100) + np.log10(wo2) * Denorm / 100) / (2 * np.pi)
            self.fo = fod * fpMin
        elif self.type is "High Pass":
            wo1 = self.get_L_wo(self)
            wo2 = wo1 * self.wan / self.get_L_wa(self)
            fod = 10 ** (np.log10(wo1) * (1 - Denorm / 100) + np.log10(wo2) * Denorm / 100) / (2 * np.pi)
            self.fo = fpMin / fod
        elif self.type is "Band Pass":
            Bw = 2 * np.pi * fpMax - fpMin
            wp = 2 * np.pi * np.sqrt(fpMin * fpMax)
            wo1 = self.get_L_wo(self)
            wo2 = wo1 * self.wan / self.get_L_wa(self)
            wod = 10 ** (np.log10(wo1) * (1 - Denorm / 100) + np.log10(wo2) * Denorm / 100)
            wo1 = (np.sqrt((wod * Bw) ** 2 + 4 * wp ** 2) + wod * Bw) / 2
            wo2 = (np.sqrt((wod * Bw) ** 2 + 4 * wp ** 2) - wod * Bw) / 2
            self.fo = np.sqrt(wo1 * wo2) / (2 * np.pi)
        elif self.type is "Band Reject":
            Bw = 2 * np.pi * fpMax - fpMin
            wp = 2 * np.pi * np.sqrt(fpMin * fpMax)
            wo1 = self.get_L_wo(self)
            wo2 = wo1 * self.wan / self.get_L_wa(self)
            wod = 10 ** (np.log10(wo1) * (1 - Denorm / 100) + np.log10(wo2) * Denorm / 100)
            wo1 = (np.sqrt((Bw / wod) ** 2 + 4 * wp ** 2) + Bw / wod) / 2
            wo2 = (np.sqrt((Bw / wod) ** 2 + 4 * wp ** 2) - Bw / wod ) / 2
            self.fo = np.sqrt(wo1 * wo2) / (2 * np.pi)
        else:
            message = "Error: Enter Filter Type."
            return message

    def calc_Denormalization_zpk(self):
        if self.type is "Low Pass":
            z, p, k = self.get_L_zpk(self)
            self.z, self.p, self.k = signal.lp2lp_zpk(z, p, k, wo=self.fo / (2 * np.pi))
        elif self.type is "High Pass":
            z, p, k = self.get_L_zpk(self)
            self.z, self.p, self.k = signal.lp2hp_zpk(z, p, k, wo=self.fo / (2 * np.pi))
        elif self.type is "Band Pass":
            z, p, k = self.get_L_zpk(self)
            self.z, self.p, self.k = signal.lp2bp_zpk(z, p, k, wo=self.fo / (2 * np.pi))
        elif self.type is "Band Reject":
            z, p, k = self.get_L_zpk(self)
            self.z, self.p, self.k = signal.lp2bs_zpk(z, p, k, wo=self.fo / (2 * np.pi))

    def get_zpk(self):
        return self.z, self.p, self.k

    def get_zpGk(self):
        Gk = self.k * 10 ** (self.get_Gain(self) / 20)
        return self.z, self.p, Gk

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

    def get_Gain(self):
        return self.filter.reqData[FilterData.gain.value]

    def calc_TransFunc(self):
        val, msg = self.filter.validate(self.filter)
        if val is False:
            return msg
        z, p, k = self.get_zpk(self)
        if self.type is "Low Pass":
            sys = signal.ZerosPolesGain(z, p, k)
            self.w_tf, self.h = signal.TransferFunction(sys)
        elif self.type is "High Pass":
            sys = signal.ZerosPolesGain(z, p, k)
            self.w_tf, self.h = signal.TransferFunction(sys)
        elif self.type is "Band Pass":
            sys = signal.ZerosPolesGain(z, p, k)
            self.w_tf, self.h = signal.TransferFunction(sys)
        elif self.type is "Band Reject":
            sys = signal.ZerosPolesGain(z, p, k)
            self.w_tf, self.h = signal.TransferFunction(sys)
        else:
            message = "Error: Enter Filter Type."
            return message

    def calc_MagAndPhase(self):                                     # return angular frequency, Mag and Phase
        val, msg = self.filter.validate(self.filter)
        if val is False:
            return msg
        z, p, k = self.get_zpk(self)
        if self.type is "Low Pass":
            sys = signal.ZerosPolesGain(z, p, k)
            self.w_bode, self.mag, self.pha = signal.bode(sys)
        elif self.type is "High Pass":
            sys = signal.ZerosPolesGain(z, p, k)
            self.w_bode, self.mag, self.pha = signal.bode(sys)
        elif self.type is "Band Pass":
            sys = signal.ZerosPolesGain(z, p, k)
            self.w_bode, self.mag, self.pha = signal.bode(sys)
        elif self.type is "Band Reject":
            sys = signal.ZerosPolesGain(z, p, k)
            self.w_bode, self.mag, self.pha = signal.bode(sys)
        else:
            message = "Error: Enter Filter Type."
            return message

    def get_TransFuncWithGain(self):
        hg = self.h
        gain = self.get_Gain(self)
        for i in range(0, len(hg)):
            hg[i] = hg[i] * 10 ** (gain / 20)
        return self.w_tf, hg

    def get_TransFuncWithoutGain(self):
        return self.w_tf, self.h

    def get_MagAndPhaseWithGain(self):
        magg = self.mag
        gain = self.get_Gain(self)
        for i in range(0, len(magg)):
            magg[i] = magg[i] + gain
        return self.w_bode, magg, self.pha

    def get_MagAndPhaseWithoutGain(self):
        return self.w_bode, self.mag, self.pha

    def calc_Attenuation(self):
        A = self.mag
        for i in range(0, len(A)):
            A[i] = 1 / A[i]
        self.A = A

    def get_Attenuation(self):
        return self.w_bode, self.A

    def calc_Group_Delay(self):
        w, mag, pha = self.get_MagAndPhaseWithoutGain(self)
        gd = - np.diff(pha) / np.diff(w)
        gd = gd.tolist()
        gd.append(gd[len(gd) - 1])
        self.GroupDelay = gd
        w = w.tolist()
        self.wgd = w

    def get_Group_Delay(self):
        return self.wgd, self.GroupDelay

    def calculate(self):
        self.calc_Order(self)
        self.calc_fo(self)
        self.calc_Denormalization_zpk(self)
        while self.check_Q(self) is False:
            self.calc_fo(self)
            self.calc_Denormalization_zpk(self)
        self.calc_TransFunc(self)
        self.calc_MagAndPhase(self)
        self.calc_Group_Delay(self)

    #############################
    #       Legendre Calc       #
    #############################

    def get_Legendre_Poly(self, order):
        return special.legendre(order)

    def get_even_L_Optimum(self, n):
        k = (n - 2) / 2
        if k % 2 == 0:
            a0 = 1 / (np.sqrt((k + 1) * (k + 2)))
            coefs = [a0]
            i = 1
            while i <= k:
                if i % 2 == 0:
                    coefs.append((2 * i + 1) * a0)

                else:
                    coefs.append(0)
                i += 1

            poly_sum = coefs[0] * self.get_Legendre_Poly(self, 0)
            i = 1
            while i <= k:
                poly_sum = np.polyadd(poly_sum, coefs[i] * self.get_Legendre_Poly(self, i))
                i += 1

            dpoly = np.polymul(poly_sum, poly_sum)
            dpoly = np.polymul(dpoly, np.poly1d([1, 1]))
            poly = np.polyint(dpoly)
            uplim = np.poly1d([2, 0, -1])
            lowlim = -1
            poly = np.polysub(np.polyval(poly, uplim), np.polyval(poly, lowlim))

        else:
            a1 = 3 / (np.sqrt((k + 1) * (k + 2)))
            coefs = [0, a1]
            i = 2
            while i <= k:
                if i % 2 == 0:
                    coefs.append(0)
                else:
                    coefs.append((2 * i + 1) * a1 / 3)
                i += 1

            poly_sum = coefs[0] * self.get_Legendre_Poly(self, 0)
            i = 1
            while i <= k:
                poly_sum = np.polyadd(poly_sum, coefs[i] * self.get_Legendre_Poly(self, i))
                i += 1

            dpoly = np.polymul(poly_sum, poly_sum)
            dpoly = np.polymul(dpoly, np.poly1d([1, 1]))
            poly = np.polyint(dpoly)
            uplim = np.poly1d([2, 0, -1])
            lowlim = -1
            poly = np.polysub(np.polyval(poly, uplim), np.polyval(poly, lowlim))
        return poly

    def get_odd_L_Optimum(self, n):
        k = (n - 1) / 2
        ao = 1 / (np.sqrt(2) * (k + 1))
        coefs = [ao]
        i = 1
        while i <= k:
            coefs.append((2 * i + 1) * ao)
            i += 1

        poly_sum = coefs[0] * self.get_Legendre_Poly(0)
        i = 1
        while i <= k:
            poly_sum = np.polyadd(poly_sum, coefs[i] * self.get_Legendre_Poly(i))
            i += 1

        dpoly = np.polymul(poly_sum, poly_sum)
        poly = np.polyint(dpoly)
        uplim = np.poly1d([2, 0, -1])
        lowlim = -1
        poly = np.polysub(np.polyval(poly, uplim), np.polyval(poly, lowlim))
        return poly

    def get_L_Poly(self, n):
        if n % 2 == 0:
            Ln = self.get_even_L_Optimum(self, n)
        else:
            Ln = self.get_odd_L_Optimum(self, n)
        return Ln

    def get_L_Poly_Value(self, n, wan):
        value = np.polyval(self.get_L_Poly(self, n), wan)
        return value

    def get_L_Filter_Poly(self):
        Ln = self.get_L_Poly(self, self.order)
        e2Ln = np.polymul(np.poly1d([self.epsilon2]), Ln)
        A = np.polyadd(np.poly1d([1]), e2Ln)
        return A

    def get_L_zpk(self):
        z = []
        r = 1j * np.roots(self.get_L_Filter_Poly(self))
        p = []
        for i in range(0, len(r)):
            if r[i].real < 0:
                p.append(r[i])
        k = np.polyval(self.get_L_Filter_Poly(self), 0)
        for pole in p:
            k *= pole
        return z, p, k

    def get_L_System(self):
        z, p, k = self.get_L_zpk(self)
        return signal.lti(z, p, k)

    def get_L_wo(self):
        sys = self.get_L_System(self)
        w, mag, pha = signal.bode(sys, w=np.logspace(-1, 1, num=100000))
        for i in range(0, len(mag)):
            if mag[i] <= -3:
                wo = w[i]
                return wo

    def get_L_wa(self):
        sys = self.get_L_System(self)
        w, mag, pha = signal.bode(sys, w=np.logspace(-1, 1, num=100000))
        for i in range(0, len(mag)):
            if mag[i] <= -self.Aa:
                wa = w[i]
                return wa
