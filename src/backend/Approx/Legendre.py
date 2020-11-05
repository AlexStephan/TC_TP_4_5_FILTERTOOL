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
        self.epsilon2 = None

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
            wan = faMin / fpMin
            while self.get_L_Poly_Value(self, order, wan) < (10 ** (Aa / 10) - 1) / epsilon2:
                order = order + 1
            if Nmin > order:
                order = Nmin
                self.order = Nmin
            elif Nmax < order:
                order = Nmax
                self.order = Nmax
            else:
                self.order = order
        elif self.type is "High Pass":
            wan = fpMin / faMin
            while self.get_L_Poly_Value(self, order, wan) < (10 ** (Aa / 10) - 1) / epsilon2:
                order = order + 1
            if Nmin > order:
                order = Nmin
                self.order = Nmin
            elif Nmax < order:
                order = Nmax
                self.order = Nmax
            else:
                self.order = order
        elif self.type is "Band Pass":
            wan = (faMax - faMin) / (fpMax - fpMin)
            while self.get_L_Poly_Value(self, order, wan) < (10 ** (Aa / 10) - 1) / epsilon2:
                order = order + 1
            if Nmin > order:
                order = Nmin
                self.order = Nmin
            elif Nmax < order:
                order = Nmax
                self.order = Nmax
            else:
                self.order = order
        elif self.type is "Band Reject":
            wan = (fpMax - fpMin) / (faMax - faMin)
            while self.get_L_Poly_Value(self, order, wan) < (10 ** (Aa / 10) - 1) / epsilon2:
                order = order + 1
            if Nmin > order:
                order = Nmin
                self.order = Nmin
            elif Nmax < order:
                order = Nmax
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
        Ap = self.filter.reqData[FilterData.Ap.value]

        faMin = self.filter.reqData[FilterData.faMin.value]
        faMax = self.filter.reqData[FilterData.fpMax.value]
        Aa = self.filter.reqData[FilterData.Aa.value]

        Denorm = self.filter.reqData[FilterData.Denorm.value]

        if self.type is "Low Pass":
            fo1 = fpMin / ((10 ** (Ap / 10) - 1) ** (1 / (2 * self.order)))
            fo2 = faMin / ((10 ** (Aa / 10) - 1) ** (1 / (2 * self.order)))
            self.fo = 10 ** (np.log10(fo1) * (1 - Denorm / 100) + np.log10(fo2) * Denorm / 100)
        elif self.type is "High Pass":
            fo1 = fpMin * ((10 ** (Ap / 10) - 1) ** (1 / (2 * self.order)))
            fo2 = faMin * ((10 ** (Aa / 10) - 1) ** (1 / (2 * self.order)))
            self.fo = 10 ** (np.log10(fo1) * (1 - Denorm / 100) + np.log10(fo2) * Denorm / 100)
        elif self.type is "Band Pass":
            fop1 = fpMin * (10 ** (Ap / 10) - 1) ** (1 / self.order)
            foa1 = faMin * (10 ** (Aa / 10) - 1) ** (1 / self.order)
            fop2 = fpMax / (10 ** (Ap / 10) - 1) ** (1 / self.order)
            foa2 = faMax / (10 ** (Aa / 10) - 1) ** (1 / self.order)
            fo1 = 10 ** (np.log10(fop1) * (1 - Denorm / 100) + np.log10(foa1) * Denorm / 100)
            fo2 = 10 ** (np.log10(fop2) * (1 - Denorm / 100) + np.log10(foa2) * Denorm / 100)
            self.fo = np.sqrt(fo1 * fo2)
        elif self.type is "Band Reject":
            fop1 = fpMin / (10 ** (Ap / 10) - 1) ** (1 / self.order)
            foa1 = faMin / (10 ** (Aa / 10) - 1) ** (1 / self.order)
            fop2 = fpMax * (10 ** (Ap / 10) - 1) ** (1 / self.order)
            foa2 = faMax * (10 ** (Aa / 10) - 1) ** (1 / self.order)
            fo1 = 10 ** (np.log10(fop1) * (1 - Denorm / 100) + np.log10(foa1) * Denorm / 100)
            fo2 = 10 ** (np.log10(fop2) * (1 - Denorm / 100) + np.log10(foa2) * Denorm / 100)
            self.fo = np.sqrt(fo1 * fo2)
        else:
            message = "Error: Enter Filter Type."
            return message


    def get_Legendre_Poly(self, order):
        return special.legendre(order)

    def get_even_L_Optimum(self, n):
        k = (n - 2) / 2
        if k % 2 == 0:
            a0 = 1 / (np.sqrt((k + 1) * (k + 2)))
            coefs = []
            coefs.append(a0)
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
            coefs = []
            coefs.append(0)
            coefs.append(a1)
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
        coefs = []
        coefs.append(ao)
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

    def get_L_Attenuation(self, order):
        Ln = self.get_L_Poly(self, order)
        Aw = np.polyadd(np.poly1d([1]), self.epsilon2 * Ln)
        return Aw

    def get_L_zpk(self, order):
        poles = np.roots(self.get_L_Attenuation(self, order))
