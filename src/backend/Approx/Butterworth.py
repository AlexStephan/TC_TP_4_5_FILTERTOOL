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

    def get_Data(self):
        data = {FilterData.Aa: None, FilterData.faMin: None,
                FilterData.Ap: None, FilterData.fpMin: None,
                FilterData.gain: None, FilterData.Nmax: None,
                FilterData.Nmin: None, FilterData.Qmax: None,
                FilterData.Denorm: None}
        data = self.filter.get_reqData(self.filter)
        return data

    def calc_Order(self):
        val, msg = self.filter.validate(self.filter)
        if val is False:
            return msg
        if self.type is "Low Pass":
            data = self.get_Data(self)
            order, wo = signal.buttord(2 * np.pi * data.fpMin, 2 * np.pi * data.faMin, data.Ap, data.Aa, analog=True)
            if data.Nmin > order:
                order = data.Nmin
                self.order = data.Nmin
            elif data.Nmax < order:
                order = data.Nmax
                self.order = data.Nmax
            else:
                self.order = order
        elif self.type is "High Pass":
            data = self.get_Data(self)
            order, wo = signal.buttord(2 * np.pi * data.fpMin, 2 * np.pi * data.faMin, data.Ap, data.Aa, analog=True)
            if data.Nmin > order:
                order = data.Nmin
                self.order = data.Nmin
            elif data.Nmax < order:
                order = data.Nmax
                self.order = data.Nmax
            else:
                self.order = order
        elif self.type is "Band Pass":
            data = self.get_Data(self)
            order, wo = signal.buttord([2 * np.pi * data.fpMin, 2 * np.pi * data.fpMax],
                                            [2 * np.pi * data.faMin, 2 * np.pi * data.faMax],
                                            data.Ap, data.Aa, analog=True)
            if data.Nmin > order:
                order = data.Nmin
                self.order = data.Nmin
            elif data.Nmax < order:
                order = data.Nmax
                self.order = data.Nmax
            else:
                self.order = order
        elif self.type is "Band Reject":
            data = self.get_Data(self)
            order, wo = signal.buttord([2 * np.pi * data.fpMin, 2 * np.pi * data.fpMax],
                                            [2 * np.pi * data.faMin, 2 * np.pi * data.faMax],
                                            data.Ap, data.Aa, analog=True)
            if data.Nmin > order:
                order = data.Nmin
                self.order = data.Nmin
            elif data.Nmax < order:
                order = data.Nmax
                self.order = data.Nmax
            else:
                self.order = order
        else:
            message = "Error: Enter Filter Type."
            return message

    def calc_fo(self):
        val, msg = self.filter.validate(self.filter)
        if val is False:
            return msg
        if self.type is "Low Pass":
            data = self.get_Data(self)
            fo1 = data.fpMin / ((10 ** (data.Ap / 10) - 1) ** (1 / (2 * self.order)))
            fo2 = data.faMin / ((10 ** (data.Aa / 10) - 1) ** (1 / (2 * self.order)))
            self.fo = 10 ** (np.log10(fo1) * (1 - data.Denorm / 100) + np.log10(fo2) * data.Denorm / 100)
        elif self.type is "High Pass":
            data = self.get_Data(self)
            fo1 = data.fpMin * ((10 ** (data.Ap / 10) - 1) ** (1 / (2 * self.order)))
            fo2 = data.faMin * ((10 ** (data.Aa / 10) - 1) ** (1 / (2 * self.order)))
            self.fo = 10 ** (np.log10(fo1) * (1 - data.Denorm / 100) + np.log10(fo2) * data.Denorm / 100)
        elif self.type is "Band Pass":
            data = self.get_Data(self)
            fop1 = data.fpMin * (10 ** (data.Ap / 10) - 1) ** (1 / self.order)
            foa1 = data.faMin * (10 ** (data.Aa / 10) - 1) ** (1 / self.order)
            fop2 = data.fpMax / (10 ** (data.Ap / 10) - 1) ** (1 / self.order)
            foa2 = data.faMax / (10 ** (data.Aa / 10) - 1) ** (1 / self.order)
            fo1 = 10 ** (np.log10(fop1) * (1 - data.Denorm / 100) + np.log10(foa1) * data.Denorm / 100)
            fo2 = 10 ** (np.log10(fop2) * (1 - data.Denorm / 100) + np.log10(foa2) * data.Denorm / 100)
            self.fo = np.sqrt(fo1 * fo2)
        elif self.type is "Band Reject":
            data = self.get_Data(self)
            fop1 = data.fpMin / (10 ** (data.Ap / 10) - 1) ** (1 / self.order)
            foa1 = data.faMin / (10 ** (data.Aa / 10) - 1) ** (1 / self.order)
            fop2 = data.fpMax * (10 ** (data.Ap / 10) - 1) ** (1 / self.order)
            foa2 = data.faMax * (10 ** (data.Aa / 10) - 1) ** (1 / self.order)
            fo1 = 10 ** (np.log10(fop1) * (1 - data.Denorm / 100) + np.log10(foa1) * data.Denorm / 100)
            fo2 = 10 ** (np.log10(fop2) * (1 - data.Denorm / 100) + np.log10(foa2) * data.Denorm / 100)
            self.fo = np.sqrt(fo1 * fo2)
        else:
            message = "Error: Enter Filter Type."
            return message


    def calc_NumDen(self):
        val, msg = self.filter.validate(self.filter)
        if val is False:
            return msg
        if self.type is "Low Pass":
            self.num, self.den = signal.butter(self.order, 2 * np.pi * self.fo,
                                               btype='lowpass', analog=True, output='ba')
        elif self.type is "High Pass":
            self.num, self.den = signal.butter(self.order, 2 * np.pi * self.fo,
                                               btype='highpass', analog=True, output='ba')
        elif self.type is "Band Pass":
            self.num, self.den = signal.butter(self.order, 2 * np.pi * self.fo,
                                               btype='bandpass', analog=True, output='ba')
        elif self.type is "Band Reject":
            self.num, self.den = signal.butter(self.order, 2 * np.pi * self.fo,
                                               btype='bandstop', analog=True, output='ba')
        else:
            message = "Error: Enter Filter Type."
            return message

    def calc_zpk(self):
        val, msg = self.filter.validate(self.filter)
        if val is False:
            return msg
        if self.type is "Low Pass":
            self.z, self.p, self.k = signal.butter(self.order, 2 * np.pi * self.fo,
                                                   btype='lowpass', analog=True, output='zpk')
        elif self.type is "High Pass":
            self.z, self.p, self.k = signal.butter(self.order, 2 * np.pi * self.fo,
                                                   btype='highpass', analog=True, output='zpk')
        elif self.type is "Band Pass":
            self.z, self.p, self.k = signal.butter(self.order, 2 * np.pi * self.fo,
                                                   btype='bandpass', analog=True, output='zpk')
        elif self.type is "Band Reject":
            self.z, self.p, self.k = signal.butter(self.order, 2 * np.pi * self.fo,
                                                   btype='bandstop', analog=True, output='zpk')
        else:
            message = "Error: Enter Filter Type."
            return message

    def get_zpk(self):
        return self.z, self.p, self.k

    def get_SOS(self):
        val, msg = self.filter.validate(self.filter)
        if val is False:
            return msg
        if self.type is "Low Pass":
            self.sos = signal.butter(self.order, 2 * np.pi * self.fo,
                                     btype='lowpass', analog=True, output='sos')
        elif self.type is "High Pass":
            self.sos = signal.butter(self.order, 2 * np.pi * self.fo,
                                     btype='highpass', analog=True, output='sos')
        elif self.type is "Band Pass":
            self.sos = signal.butter(self.order, 2 * np.pi * self.fo,
                                     btype='bandpass', analog=True, output='sos')
        elif self.type is "Band Stop":
            self.sos = signal.butter(self.order, 2 * np.pi * self.fo,
                                     btype='bandstop', analog=True, output='sos')
        else:
            message = "Error: Enter Filter Type."
            return message

    def get_TransFunc(self):                                        #return angular frequency and H(s)
        val, msg = self.filter.validate(self.filter)
        if val is False:
            return msg
        if self.type is "Low Pass":
            data = self.get_Data(self)
            z, p, k = self.get_zpk(self)
            sys = signal.ZerosPolesGain(z, p, k)
            w, mag, pha = signal.bode(sys)
            mag = mag * 10 ** (data.gain / 20)
            return w, mag, pha
        elif self.type is "High Pass":
            data = self.get_Data(self)
            z, p, k = self.get_zpk(self)
            sys = signal.ZerosPolesGain(z, p, k)
            w, mag, pha = signal.bode(sys)
            mag = mag * 10 ** (data.gain / 20)
            return w, mag, pha
        elif self.type is "Band Pass":
            data = self.get_Data(self)
            z, p, k = self.get_zpk(self)
            sys = signal.ZerosPolesGain(z, p, k)
            w, mag, pha = signal.bode(sys)
            mag = mag * 10 ** (data.gain / 20)
            return w, mag, pha
        elif self.type is "Band Reject":
            data = self.get_Data(self)
            z, p, k = self.get_zpk(self)
            sys = signal.ZerosPolesGain(z, p, k)
            w, mag, pha = signal.bode(sys)
            mag = mag * 10 ** (data.gain / 20)
            return w, mag, pha
        else:
            message = "Error: Enter Filter Type."
            return message

    def check_Q(self) -> bool:
        if self.order > 1:
            data = self.get_Data(self)
            z, p, k = self.get_zpk(self)
            q_arr = []
            for pole in p:
                q = abs(abs(pole) / (2 * pole.real))
                q_arr.append(q)
            q_sys = np.max(q_arr)
            if q_sys > data.Qmax and self.order > 1:
                self.order = self.order - 1
                return 1
            else:
                return 0
        else:
            return 0

    def calculate(self):
        self.calc_Order(self)
        self.calc_fo(self)
        self.calc_zpk(self)
        while self.check_Q(self) is 1:
            self.calc_fo(self)
            self.calc_zpk(self)
        return self.get_TransFunc(self)