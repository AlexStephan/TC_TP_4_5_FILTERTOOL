from src.backend.Filter.Filter import *
from scipy import signal
import numpy as np

class Butterworth(object):
    def __init__(self, filter: Filter):
        self.filter = filter
        self.type = self.filter.get_type()
        self.order = None
        self.fo = None
        self.z = None
        self.p = None
        self.k = None
        self.num = None
        self.den = None
        self.sos = None
        self.w_tf = None
        self.h = None
        self.w_bode = None
        self.mag = None
        self.pha = None
        self.A = None
        self.wgd = None
        self.GroupDelay = None
        self.calculate()


    def calc_Order(self):
        val, msg = self.filter.validate()
        if val is False:
            return msg

        fpMin = self.filter.reqData[FilterData.fpMin]
        fpMax = self.filter.reqData[FilterData.fpMax]
        Ap = self.filter.reqData[FilterData.Ap]

        faMin = self.filter.reqData[FilterData.faMin]
        faMax = self.filter.reqData[FilterData.faMax]
        Aa = self.filter.reqData[FilterData.Aa]

        Nmin = self.filter.reqData[FilterData.Nmin]
        Nmax = self.filter.reqData[FilterData.Nmax]

        if self.type == "Low Pass":
            order, wo = signal.buttord(2 * np.pi * fpMin, 2 * np.pi * faMin, Ap, Aa, analog=True)
            if Nmin is not None and Nmin > order:
                self.order = Nmin
            elif Nmax is not None and Nmax < order:
                self.order = Nmax
            else:
                self.order = order
        elif self.type == "High Pass":
            order, wo = signal.buttord(2 * np.pi * fpMin, 2 * np.pi * faMin, Ap, Aa, analog=True)
            if Nmin is not None and Nmin > order:
                self.order = Nmin
            elif Nmax is not None and Nmax < order:
                self.order = Nmax
            else:
                self.order = order
        elif self.type == "Band Pass":
            order, wo = signal.buttord([2 * np.pi * fpMin, 2 * np.pi * fpMax],
                                            [2 * np.pi * faMin,2 * np.pi * faMax],
                                            Ap, Aa, analog=True)
            if Nmin is not None and Nmin > order and Nmin >= 2:
                self.order = Nmin
            elif Nmax is not None and order > Nmax >= 2:
                self.order = Nmax
            elif order >= 2:
                self.order = order
        elif self.type == "Band Reject":
            order, wo = signal.buttord([2 * np.pi * fpMin, 2 * np.pi * fpMax],
                                            [2 * np.pi * faMin, 2 * np.pi * faMax],
                                            Ap, Aa, analog=True)
            if Nmin is not None and Nmin > order and Nmin >= 2:
                self.order = Nmin
            elif Nmax is not None and order > Nmax >= 2:
                self.order = Nmax
            elif order >= 2:
                self.order = order
        else:
            message = "Error: Enter Filter Type."
            return message

    def calc_fo(self):
        val, msg = self.filter.validate()
        if val is False:
            return msg

        fpMin = self.filter.reqData[FilterData.fpMin]
        fpMax = self.filter.reqData[FilterData.fpMax]
        Ap = self.filter.reqData[FilterData.Ap]

        faMin = self.filter.reqData[FilterData.faMin]
        faMax = self.filter.reqData[FilterData.fpMax]
        Aa = self.filter.reqData[FilterData.Aa]

        Denorm = self.filter.reqData[FilterData.Denorm]

        if self.type == "Low Pass":
            fo1 = fpMin / ((10 ** (Ap / 10) - 1) ** (1 / (2 * self.order)))
            fo2 = faMin / ((10 ** (Aa / 10) - 1) ** (1 / (2 * self.order)))
            self.fo = 10 ** (np.log10(fo1) * (1 - Denorm / 100) + np.log10(fo2) * Denorm / 100)
        elif self.type == "High Pass":
            fo1 = fpMin * ((10 ** (Ap / 10) - 1) ** (1 / (2 * self.order)))
            fo2 = faMin * ((10 ** (Aa / 10) - 1) ** (1 / (2 * self.order)))
            self.fo = 10 ** (np.log10(fo1) * (1 - Denorm / 100) + np.log10(fo2) * Denorm / 100)
        elif self.type == "Band Pass":
            fop1 = fpMin * (10 ** (Ap / 10) - 1) ** (1 / self.order)
            foa1 = faMin * (10 ** (Aa / 10) - 1) ** (1 / self.order)
            fop2 = fpMax / (10 ** (Ap / 10) - 1) ** (1 / self.order)
            foa2 = faMax / (10 ** (Aa / 10) - 1) ** (1 / self.order)
            fo1 = 10 ** (np.log10(fop1) * (1 - Denorm / 100) + np.log10(foa1) * Denorm / 100)
            fo2 = 10 ** (np.log10(fop2) * (1 - Denorm / 100) + np.log10(foa2) * Denorm / 100)
            self.fo = [fo1, fo2]
        elif self.type == "Band Reject":
            fop1 = fpMin / (10 ** (Ap / 10) - 1) ** (1 / self.order)
            foa1 = faMin / (10 ** (Aa / 10) - 1) ** (1 / self.order)
            fop2 = fpMax * (10 ** (Ap / 10) - 1) ** (1 / self.order)
            foa2 = faMax * (10 ** (Aa / 10) - 1) ** (1 / self.order)
            fo1 = 10 ** (np.log10(fop1) * (1 - Denorm / 100) + np.log10(foa1) * Denorm / 100)
            fo2 = 10 ** (np.log10(fop2) * (1 - Denorm / 100) + np.log10(foa2) * Denorm / 100)
            self.fo = [fo1, fo2]
        else:
            message = "Error: Enter Filter Type."
            return message


    def calc_NumDen(self):
        val, msg = self.filter.validate()
        if val is False:
            return msg
        if self.type == "Low Pass":
            self.num, self.den = signal.butter(self.order, 2 * np.pi * self.fo,
                                               btype='lowpass', analog=True, output='ba')
        elif self.type == "High Pass":
            self.num, self.den = signal.butter(self.order, 2 * np.pi * self.fo,
                                               btype='highpass', analog=True, output='ba')
        elif self.type == "Band Pass":
            self.num, self.den = signal.butter(self.order, np.multiply(self.fo, 2*np.pi),
                                               btype='bandpass', analog=True, output='ba')
        elif self.type == "Band Reject":
            self.num, self.den = signal.butter(self.order, np.multiply(self.fo, 2*np.pi),
                                               btype='bandstop', analog=True, output='ba')
        else:
            message = "Error: Enter Filter Type."
            return message

    def calc_zpk(self):
        val, msg = self.filter.validate()
        if val is False:
            return msg
        if self.type == "Low Pass":
            self.z, self.p, self.k = signal.butter(self.order, 2 * np.pi * self.fo,
                                                   btype='lowpass', analog=True, output='zpk')
        elif self.type == "High Pass":
            self.z, self.p, self.k = signal.butter(self.order, 2 * np.pi * self.fo,
                                                   btype='highpass', analog=True, output='zpk')
        elif self.type == "Band Pass":
            self.z, self.p, self.k = signal.butter(self.order, np.multiply(self.fo, 2*np.pi),
                                                   btype='bandpass', analog=True, output='zpk')
        elif self.type == "Band Reject":
            self.z, self.p, self.k = signal.butter(self.order, np.multiply(self.fo, 2*np.pi),
                                                   btype='bandstop', analog=True, output='zpk')
        else:
            message = "Error: Enter Filter Type."
            return message

    def calc_sos(self):
        val, msg = self.filter.validate()
        if val is False:
            return msg
        if self.type == "Low Pass":
            self.sos = signal.butter(self.order, 2 * np.pi * self.fo,
                                     btype='lowpass', analog=True, output='sos')
        elif self.type == "High Pass":
            self.sos = signal.butter(self.order, 2 * np.pi * self.fo,
                                     btype='highpass', analog=True, output='sos')
        elif self.type == "Band Pass":
            self.sos = signal.butter(self.order, np.multiply(self.fo, 2*np.pi),
                                     btype='bandpass', analog=True, output='sos')
        elif self.type == "Band Stop":
            self.sos = signal.butter(self.order, np.multiply(self.fo, 2*np.pi),
                                     btype='bandstop', analog=True, output='sos')
        else:
            message = "Error: Enter Filter Type."
            return message

    def calc_TransFunc(self):
        val, msg = self.filter.validate()
        if val is False:
            return msg
        z, p, k = self.get_zpk()
        if self.type == "Low Pass":
            sys = signal.lti(z, p, k)
            self.w_tf, self.h = sys.freqresp()
        elif self.type == "High Pass":
            sys = signal.lti(z, p, k)
            self.w_tf, self.h = sys.freqresp()
        elif self.type == "Band Pass":
            sys = signal.lti(z, p, k)
            self.w_tf, self.h = sys.freqresp()
        elif self.type == "Band Reject":
            sys = signal.lti(z, p, k)
            self.w_tf, self.h = sys.freqresp()
        else:
            message = "Error: Enter Filter Type."
            return message


    def calc_MagAndPhase(self):              # return angular frequency, Mag and Phase
        val, msg = self.filter.validate()
        if val is False:
            return msg
        z, p, k = self.get_zpk()
        if self.type == "Low Pass":
            sys = signal.ZerosPolesGain(z, p, k)
            self.w_bode, self.mag, self.pha = signal.bode(sys)
        elif self.type == "High Pass":
            sys = signal.ZerosPolesGain(z, p, k)
            self.w_bode, self.mag, self.pha = signal.bode(sys)
        elif self.type == "Band Pass":
            sys = signal.ZerosPolesGain(z, p, k)
            self.w_bode, self.mag, self.pha = signal.bode(sys)
        elif self.type == "Band Reject":
            sys = signal.ZerosPolesGain(z, p, k)
            self.w_bode, self.mag, self.pha = signal.bode(sys)
        else:
            message = "Error: Enter Filter Type."
            return message

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

    def calculate(self):
        self.calc_Order()
        self.calc_fo()
        self.calc_NumDen()
        self.calc_zpk()
        while self.check_Q() is False:
            self.calc_fo()
            self.calc_NumDen()
            self.calc_zpk()
        self.calc_TransFunc()
        self.calc_MagAndPhase()
        self.calc_Group_Delay()

    def calc_Attenuation(self):
        A = self.mag
        for i in range(0, len(A)):
            A[i] = 1 / A[i]
        self.A = A

    def calc_Group_Delay(self):
        w, mag, pha = self.get_MagAndPhaseWithoutGain()
        gd = np.divide(- np.diff(pha), np.diff(w))
        gd = gd.tolist()
        gd.append(gd[len(gd) - 1])
        self.GroupDelay = gd
        w = w.tolist()
        self.wgd = w

    def get_Gain(self):
        return self.filter.reqData[FilterData.gain]

    #####################
    #       ALEX        #
    #####################

    def get_NumDen(self):
        return self.num, self.den

    def get_zpk(self):
        return self.z, self.p, self.k

    def get_zpGk(self):
        gain = self.get_Gain()
        Gk = self.k * 10 ** (gain / 20)
        return self.z, self.p, Gk

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

    def get_Attenuation(self):
        return self.w_bode, self.A

    def get_Norm_Attenuation(self):
        num, den = signal.normalize(self.num, self.den)
        sys = signal.lti(num, den)
        w, mag, pha = signal.bode(sys)
        return w, self.A

    def get_Group_Delay(self):
        return self.wgd, self.GroupDelay

    def get_Order(self):
        return self.order
