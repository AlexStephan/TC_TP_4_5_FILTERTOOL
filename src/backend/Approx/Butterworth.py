from src.backend.Filter.Filter import *
from scipy import signal
import numpy as np

class Butterworth(object):
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

    def get_Order(self):
        val, msg = self.filter.validate(self.filter)
        if val is False:
            return msg
        if self.type is "Low Pass":
            data = {FilterData.Aa: None, FilterData.faMin: None,
                    FilterData.Ap: None, FilterData.fpMin: None,
                    FilterData.gain: None, FilterData.Nmax: None,
                    FilterData.Nmin: None, FilterData.Qmax: None,
                    FilterData.Denorm: None}
            data = self.filter.get_reqData(self.filter)
            self.order, wo = signal.buttord(2 * np.pi * data.fpMin, 2 * np.pi * data.faMin, data.Ap, data.Aa, analog=True)
            if data.Nmin > self.order:
                self.order = data.Nmin
                return data.Nmin
            elif data.Nmax < self.order:
                self.order = data.Nmax
                return data.Nmax
            else:
                return self.order
        elif self.type is "High Pass":
            data = {FilterData.Aa: None, FilterData.faMin: None,
                    FilterData.Ap: None, FilterData.fpMin: None,
                    FilterData.gain: None, FilterData.Nmax: None,
                    FilterData.Nmin: None, FilterData.Qmax: None,
                    FilterData.Denorm: None}
            data = self.filter.get_reqData(self.filter)
            self.order, wo = signal.buttord(2 * np.pi * data.fpMin, 2 * np.pi * data.faMin, data.Ap, data.Aa, analog=True)
            if data.Nmin > self.order:
                self.order = data.Nmin
                return data.Nmin
            elif data.Nmax < self.order:
                self.order = data.Nmax
                return data.Nmax
            else:
                return self.order

    def get_fo(self):
        val, msg = self.filter.validate(self.filter)
        if val is False:
            return msg
        if self.type is "Low Pass":
            data = {FilterData.Aa: None, FilterData.faMin: None,
                    FilterData.Ap: None, FilterData.fpMin: None,
                    FilterData.gain: None, FilterData.Nmax: None,
                    FilterData.Nmin: None, FilterData.Qmax: None,
                    FilterData.Denorm: None}
            data = self.filter.get_reqData(self.filter)
            fo1 = data.fpMin / ((10 ** (data.Ap / 10) - 1) ** (1 / (2 * self.get_Order(self))))
            fo2 = data.faMin / ((10 ** (data.Aa / 10) - 1) ** (1 / (2 * self.get_Order(self))))
            self.fo = 10 ** (np.log10(fo1) * (1 - data.Denorm) + np.log10(fo2) * data.Denorm)
        elif self.type is "High Pass":
            data = {FilterData.Aa: None, FilterData.faMin: None,
                    FilterData.Ap: None, FilterData.fpMin: None,
                    FilterData.gain: None, FilterData.Nmax: None,
                    FilterData.Nmin: None, FilterData.Qmax: None,
                    FilterData.Denorm: None}
            data = self.filter.get_reqData(self.filter)
            fo1 = data.fpMin / ((10 ** (data.Ap / 10) - 1) ** (1 / (2 * self.get_Order(self))))
            fo2 = data.faMin / ((10 ** (data.Aa / 10) - 1) ** (1 / (2 * self.get_Order(self))))
            self.fo = 10 ** (np.log10(fo1) * (1 - data.Denorm) + np.log10(fo2) * data.Denorm)

    def get_NumDen(self):
        val, msg = self.filter.validate(self.filter)
        if val is False:
            return msg
        if self.type is "Low Pass":
            self.num, self.den = signal.butter(self.get_Order(self), 2 * np.pi * self.get_fo(self), btype='low', analog=True, output='ba')
            return self.num, self.den
        elif self.type is "High Pass":
            self.num, self.den = signal.butter(self.get_Order(self), 2 * np.pi * self.get_fo(self), btype='high', analog=True, output='ba')
            return self.num, self.den

    def get_zpk(self):
        val, msg = self.filter.validate(self.filter)
        if val is False:
            return msg
        if self.type is "Low Pass":
            self.z, self.p, self.k = signal.butter(self.get_Order(self), 2*np.pi*self.get_fo(self), btype='low', analog=True, output='zpk')
            return self.z, self.p, self.k
        elif self.type is "High Pass":
            self.z, self.p, self.k = signal.butter(self.get_Order(self), 2*np.pi*self.get_fo(self), btype='high', analog=True, output='zpk')
            return self.z, self.p, self.k

    def get_SOS(self):
        val, msg = self.filter.validate(self.filter)
        if val is False:
            return msg
        if self.type is "Low Pass":
            self.sos = signal.butter(self.n, 2 * np.pi * self.fo, btype='low', analog=True, output='sos')
            return self.sos

    def get_TransFunc(self):                                        #return angular frequency and H(s)
        val, msg = self.filter.validate(self.filter)
        if val is False:
            return msg
        if self.type is "Low Pass":
            data = {FilterData.Aa: None, FilterData.faMin: None,
                    FilterData.Ap: None, FilterData.fpMin: None,
                    FilterData.gain: None, FilterData.Nmax: None,
                    FilterData.Nmin: None, FilterData.Qmax: None,
                    FilterData.Denorm: None}
            data = self.filter.get_reqData(self.filter)
            z, p, k = self.get_zpk(self)
            w, h = signal.freqs_zpk(z, p, k)
            h = h * 10 ** (data.gain / 20)
            return w, h