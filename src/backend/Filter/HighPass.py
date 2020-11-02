from src.backend.Filter.Filter import *
from src.backend.Filter.TemplateLimit import *


class HighPass(Filter):
    def __init__(self):
        self.type = FilterType.HP
        self.reqData = {FilterData.Aa: None, FilterData.faMin: None,
                        FilterData.Ap: None, FilterData.fpMin: None,
                        FilterData.gain: None,
                        FilterData.Nmax: None, FilterData.Nmin: None, FilterData.Qmax: None,
                        FilterData.Denorm: None}
        self.default = {FilterData.Aa: 30, FilterData.faMin: 9e3,
                        FilterData.Ap: 5, FilterData.fpMin: 10e3,
                        FilterData.gain: 0,
                        FilterData.Nmax: None, FilterData.Nmin: None, FilterData.Qmax: None,
                        FilterData.Denorm: 0}

    def validate(self) -> (bool, str):  # Returns true and "Ok" if everything is fine, false and "ErrMsg" if not
        valid = False
        message = "Ok"
        if self.reqData[FilterData.Aa.value] < 0 or self.reqData[FilterData.Ap.value] < 0:
            message = "Error: Enter positive values for Aa and Ap."
        elif self.reqData[FilterData.Aa] < self.reqData[FilterData.Ap]:
            message = "Error: Aa must be greater than Ap."
        elif self.reqData[FilterData.faMin] > self.reqData[FilterData.fpMin]:
            message = "Error: fa must be lesser than fp."
        else:
            valid = True
        return valid, message

    def get_template_limits(self):  # Create one set of squares for denormalized graph, and one set for
        # normalized graph
        Ap = self.reqData[FilterData.Ap.value]
        Aa = self.reqData[FilterData.Aa.value]
        fpMin = self.reqData[FilterData.fpMin.value]
        faMin = self.reqData[FilterData.faMin.value]

        denormLimit1 = Limit(Dot(0, Aa), Dot(faMin, Aa), Dot(0, 0), Dot(faMin, 0))
        denormLimit2 = Limit(Dot(fpMin, 1e9), Dot(1e12, 1e9), Dot(fpMin, Ap), Dot(1e12, Ap))
        denormLimit = [denormLimit1, denormLimit2]

        selectivity = faMin / fpMin  # K = fa/fp
        normalizedF1 = 1 / (2 * pi)
        normalizedF2 = 1 / (2 * pi * selectivity)

        normLimit1 = Limit(Dot(0, 1e9), Dot(normalizedF1, 1e9), Dot(0, Ap), Dot(normalizedF1, Ap))
        normLimit2 = Limit(Dot(normalizedF2, Aa), Dot(1e12, Aa), Dot(normalizedF2, 0), Dot(1e12, 0))
        normLimit = [normLimit1, normLimit2]

        return [denormLimit, normLimit]
