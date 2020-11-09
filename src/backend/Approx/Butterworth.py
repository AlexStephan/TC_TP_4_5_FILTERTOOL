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
        self.w_tfn = None
        self.h_n = None
        self.w_tf = None
        self.h = None
        self.w_bode = None
        self.mag = None
        self.pha = None
        self.w_att = None
        self.w_anorm = None
        self.A = None
        self.A_norm = None
        self.wgd = None
        self.GroupDelay = None
        self.timp = None
        self.impresp = None
        self.tstep = None
        self.stepresp = None

        #TODO CALCULAR WAN
        self.wan = 1
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
        faMax = self.filter.reqData[FilterData.faMax]
        Aa = self.filter.reqData[FilterData.Aa]

        Denorm = self.filter.reqData[FilterData.Denorm]
        self.d = Denorm

        if self.type == "Low Pass":
            self.f1 = fpMin / ((10 ** (Ap / 10) - 1) ** (1 / (2 * self.order)))
            self.f2 = faMin / ((10 ** (Aa / 10) - 1) ** (1 / (2 * self.order)))
            self.fo = 10 ** (np.log10(self.f1) * (1 - Denorm / 100) + np.log10(self.f2) * Denorm / 100)
        elif self.type == "High Pass":
            self.f1 = fpMin * ((10 ** (Ap / 10) - 1) ** (1 / (2 * self.order)))
            self.f2 = faMin * ((10 ** (Aa / 10) - 1) ** (1 / (2 * self.order)))
            self.fo = 10 ** (np.log10(self.f1) * (1 - Denorm / 100) + np.log10(self.f2) * Denorm / 100)
        elif self.type == "Band Pass":
            fop1 = fpMin * ((10 ** (Ap / 10) - 1) ** (1 / (2 * self.order)))
            foa1 = faMin * ((10 ** (Aa / 10) - 1) ** (1 / (2 * self.order)))
            fop2 = fpMax / ((10 ** (Ap / 10) - 1) ** (1 / (2 * self.order)))
            foa2 = faMax / ((10 ** (Aa / 10) - 1) ** (1 / (2 * self.order)))
            self.f1 = 10 ** (np.log10(fop1) * (1 - (Denorm + 50) / 100) + np.log10(foa1) * (Denorm + 50) / 100)
            self.f2 = 10 ** (np.log10(fop2) * (1 - (Denorm + 50) / 100) + np.log10(foa2) * (Denorm + 50) / 100)
            self.fo = [self.f1, self.f2]
        elif self.type == "Band Reject":
            fop1 = fpMin / ((10 ** (Ap / 10) - 1) ** (1 / (2 * self.order)))
            foa1 = faMin / ((10 ** (Aa / 10) - 1) ** (1 / (2 * self.order)))
            fop2 = fpMax * ((10 ** (Ap / 10) - 1) ** (1 / (2 * self.order)))
            foa2 = faMax * ((10 ** (Aa / 10) - 1) ** (1 / (2 * self.order)))
            self.f1 = 10 ** (np.log10(fop1) * (1 - (Denorm + 50) / 100) + np.log10(foa1) * (Denorm + 50) / 100)
            self.f2 = 10 ** (np.log10(fop2) * (1 - (Denorm + 50) / 100) + np.log10(foa2) * (Denorm + 50) / 100)
            self.fo = [self.f1, self.f2]
        else:
            message = "Error: Enter Filter Type."
            return message

    def calc_wan(self):
        val, msg = self.filter.validate()
        if val is False:
            return msg

        fpMin = self.filter.reqData[FilterData.fpMin]
        fpMax = self.filter.reqData[FilterData.fpMax]

        faMin = self.filter.reqData[FilterData.faMin]
        faMax = self.filter.reqData[FilterData.fpMax]

        if self.type == "Low Pass":
            self.wan = faMin / fpMin
        elif self.type == "High Pass":
            self.wan = fpMin / faMin
        elif self.type == "Band Pass":
            self.wan = np.minimum(fpMin / faMin, faMax / fpMax)
        elif self.type == "Band Reject":
            self.wan = np.minimum(faMin / fpMin, fpMax / faMax)
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
            self.w_tf, self.h = sys.freqresp(w=np.logspace(-3, 9, num=100000))
        elif self.type == "High Pass":
            sys = signal.lti(z, p, k)
            self.w_tf, self.h = sys.freqresp(w=np.logspace(-3, 9, num=100000))
        elif self.type == "Band Pass":
            sys = signal.lti(z, p, k)
            self.w_tf, self.h = sys.freqresp(w=np.logspace(-3, 9, num=100000))
        elif self.type == "Band Reject":
            sys = signal.lti(z, p, k)
            self.w_tf, self.h = sys.freqresp(w=np.logspace(-3, 9, num=100000))
        else:
            message = "Error: Enter Filter Type."
            return message

    def calc_Norm_TransFunc(self):
        val, msg = self.filter.validate()
        if val is False:
            return msg
        if self.type == "Low Pass":
            z, p, k = signal.butter(self.order, 2 * np.pi * self.fo,
                                    btype='lowpass', analog=True, output='zpk')
            sys = signal.lti(z, p, k)
            self.w_tfn, self.h_n = sys.freqresp(w=np.logspace(-3, 9, num=100000))
        elif self.type == "High Pass":
            z, p, k = signal.butter(self.order, 2 * np.pi * self.fo,
                                    btype='lowpass', analog=True, output='zpk')
            sys = signal.lti(z, p, k)
            self.w_tfn, self.h_n = sys.freqresp(w=np.logspace(-3, 9, num=100000))
        elif self.type == "Band Pass":
            z, p, k = signal.butter(self.order, 2 * np.pi * self.fo,
                                    btype='lowpass', analog=True, output='zpk')
            sys = signal.lti(z, p, k)
            self.w_tfn, self.h_n = sys.freqresp(w=np.logspace(-3, 9, num=100000))
        elif self.type == "Band Reject":
            z, p, k = signal.butter(self.order, 2 * np.pi * self.fo,
                                    btype='lowpass', analog=True, output='zpk')
            sys = signal.lti(z, p, k)
            self.w_tfn, self.h_n = sys.freqresp(w=np.logspace(-3, 9, num=100000))
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
            self.w_bode, self.mag, self.pha = signal.bode(sys, w=np.logspace(-3, 9, num=100000))
        elif self.type == "High Pass":
            sys = signal.ZerosPolesGain(z, p, k)
            self.w_bode, self.mag, self.pha = signal.bode(sys, w=np.logspace(-3, 9, num=100000))
        elif self.type == "Band Pass":
            sys = signal.ZerosPolesGain(z, p, k)
            self.w_bode, self.mag, self.pha = signal.bode(sys, w=np.logspace(-3, 9, num=100000))
        elif self.type == "Band Reject":
            sys = signal.ZerosPolesGain(z, p, k)
            self.w_bode, self.mag, self.pha = signal.bode(sys, w=np.logspace(-3, 9, num=100000))
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
        self.calc_wan()
        self.calc_TransFunc()
        self.calc_Norm_TransFunc()
        self.calc_MagAndPhase()
        self.calc_Group_Delay()
        self.calc_Impulse_Response()
        self.calc_Step_Response()
        self.calc_Attenuation()
        self.calc_Norm_Attenuation()

    def calc_Attenuation(self):
        w, h = self.get_TransFuncWithoutGain()
        A = []
        for i in range(len(h)):
            A.append(20 * log10(abs(1 / h[i])))
        self.w_att = w
        self.A = A

    def calc_Norm_Attenuation(self):
        w, h = self.get_Norm_TransFunc()
        wn = np.divide(w, 2 * np.pi * self.fo)
        An = []
        for i in range(len(h)):
            An.append(20 * log10(abs(1 / h[i])))
        self.w_anorm = wn
        self.A_norm = An

    def get_Norm_TransFunc(self):
        return self.w_tfn, self.h_n

    def calc_Group_Delay(self):
        w, mag, pha = self.get_MagAndPhaseWithoutGain()
        gd = np.divide(- np.diff(pha), np.diff(w))
        gd = gd.tolist()
        gd.append(gd[len(gd) - 1])
        self.GroupDelay = gd
        w = w.tolist()
        self.wgd = w

    def calc_Impulse_Response(self):
        t, out = signal.impulse(self.get_lti(), T=np.linspace(0, 1, num=100000))
        self.timp = t
        self.impresp = out

    def calc_Step_Response(self):
        t, out = signal.step(self.get_lti(), T=np.linspace(0, 1, num=100000))
        self.tstep = t
        self.stepresp = out

    def get_Gain(self):
        return self.filter.reqData[FilterData.gain]

    def get_lti(self):
        z, p, k = self.get_zpGk()
        return signal.lti(z, p, k)

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
        return self.w_att, self.A

    def get_Norm_Attenuation(self):
        return self.w_anorm, self.A_norm

    def get_Group_Delay(self):
        return self.wgd, self.GroupDelay

    def get_Order(self):
        return self.order

    def get_Impulse_Response(self):
        return self.timp, self.impresp

    def get_Step_Response(self):
        return self.tstep, self.stepresp

    def get_Qs(self):
        z, p, k = self.get_zpk()
        q_arr = []
        for pole in p:
            q = abs(abs(pole) / (2 * pole.real))
            q_arr.append(q)
        return q_arr

    #####################################
    # ALEX LO HIZO
    #####################################

    def get_very_useful_data(self): #colocar TAL CUAL en las otras aprox
        fpMin = self.filter.reqData[FilterData.fpMin]
        fpMax = self.filter.reqData[FilterData.fpMax]
        Ap = self.filter.reqData[FilterData.Ap]

        faMin = self.filter.reqData[FilterData.faMin]
        faMax = self.filter.reqData[FilterData.faMax]
        Aa = self.filter.reqData[FilterData.Aa]

        #Nmin = self.filter.reqData[FilterData.Nmin]
        #Nmax = self.filter.reqData[FilterData.Nmax]
        return fpMin,fpMax,Ap,faMin,faMax,Aa

    def get_wan(self):
        return self.wan
