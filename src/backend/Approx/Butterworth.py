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

    def calc_Order(self):
        val, msg = self.filter.validate(self.filter)
        if val is False:
            return msg

        fpMin = self.filter.reqData[FilterData.fpMin.value]
        fpMax = self.filter.reqData[FilterData.fpMax.value]
        Ap = self.filter.reqData[FilterData.Ap.value]

        faMin = self.filter.reqData[FilterData.faMin.value]
        faMax = self.filter.reqData[FilterData.fpMax.value]
        Aa = self.filter.reqData[FilterData.Aa.value]

        Nmin = self.filter.reqData[FilterData.Nmin.value]
        Nmax = self.filter.reqData[FilterData.Nmax.value]

        if self.type is "Low Pass":
            order, wo = signal.buttord(2 * np.pi * fpMin, 2 * np.pi * faMin, Ap, Aa, analog=True)
            if Nmin > order:
                order = Nmin
                self.order = Nmin
            elif Nmax < order:
                order = Nmax
                self.order = Nmax
            else:
                self.order = order
        elif self.type is "High Pass":
            order, wo = signal.buttord(2 * np.pi * fpMin, 2 * np.pi * faMin, Ap, Aa, analog=True)
            if Nmin > order:
                order = Nmin
                self.order = Nmin
            elif Nmax < order:
                order = Nmax
                self.order = Nmax
            else:
                self.order = order
        elif self.type is "Band Pass":
            order, wo = signal.buttord([2 * np.pi * fpMin, 2 * np.pi * fpMax],
                                            [2 * np.pi * faMin, 2 * np.pi * faMax],
                                            Ap, Aa, analog=True)
            if Nmin > order:
                order = Nmin
                self.order = Nmin
            elif Nmax < order:
                order = Nmax
                self.order = Nmax
            else:
                self.order = order
        elif self.type is "Band Reject":
            order, wo = signal.buttord([2 * np.pi * fpMin, 2 * np.pi * fpMax],
                                            [2 * np.pi * faMin, 2 * np.pi * faMax],
                                            Ap, Aa, analog=True)
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

    def get_NumDen(self):
        return self.num, self.den

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

    def get_zpGk(self):
        gain = self.filter.reqData[FilterData.gain.value]
        Gk = self.k * 10 ** (gain / 20)
        return self.z, self.p, Gk

    def get_Gain(self):
        return self.filter.reqData[FilterData.gain.value]

    def calc_sos(self):
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

    def get_TransFuncWithGain(self):                                        #return angular frequency and H(s)
        val, msg = self.filter.validate(self.filter)
        if val is False:
            return msg

        gain = self.filter.reqData[FilterData.gain.value]
        z, p, k = self.get_zpk(self)

        if self.type is "Low Pass":
            sys = signal.ZerosPolesGain(z, p, k)
            w, mag, pha = signal.bode(sys)
            for i in range(len(mag)):
                mag[i] = mag[i] + gain
            return w, mag, pha
        elif self.type is "High Pass":
            sys = signal.ZerosPolesGain(z, p, k)
            w, mag, pha = signal.bode(sys)
            for i in range(len(mag)):
                mag[i] = mag[i] + gain
            return w, mag, pha
        elif self.type is "Band Pass":
            sys = signal.ZerosPolesGain(z, p, k)
            w, mag, pha = signal.bode(sys)
            for i in range(len(mag)):
                mag[i] = mag[i] + gain
            return w, mag, pha
        elif self.type is "Band Reject":
            sys = signal.ZerosPolesGain(z, p, k)
            w, mag, pha = signal.bode(sys)
            for i in range(len(mag)):
                mag[i] = mag[i] + gain
            return w, mag, pha
        else:
            message = "Error: Enter Filter Type."
            return message

    def get_TransFuncWithoutGain(self):  # return angular frequency and H(s)
        val, msg = self.filter.validate(self.filter)
        if val is False:
            return msg

        z, p, k = self.get_zpk(self)

        if self.type is "Low Pass":
            sys = signal.ZerosPolesGain(z, p, k)
            w, mag, pha = signal.bode(sys)
            return w, mag, pha
        elif self.type is "High Pass":
            sys = signal.ZerosPolesGain(z, p, k)
            w, mag, pha = signal.bode(sys)
            return w, mag, pha
        elif self.type is "Band Pass":
            sys = signal.ZerosPolesGain(z, p, k)
            w, mag, pha = signal.bode(sys)
            return w, mag, pha
        elif self.type is "Band Reject":
            sys = signal.ZerosPolesGain(z, p, k)
            w, mag, pha = signal.bode(sys)
            return w, mag, pha
        else:
            message = "Error: Enter Filter Type."
            return message

    #TODO: MATI: FIJARTE SI ESTO ESTA BIEN ESCRITO
    def get_ssTransferFunction(self):
        return signal.TransferFunction(self.num,self.den)

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
        self.calc_fo(self)
        self.calc_zpk(self)
        while self.check_Q(self) is False:
            self.calc_fo(self)
            self.calc_zpk(self)
        return self.get_TransFuncWithGain(self)