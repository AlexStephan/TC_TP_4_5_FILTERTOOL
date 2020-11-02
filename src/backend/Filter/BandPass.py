from src.backend.Filter.Filter import *


class LowPass(Filter):
    def __init__(self):
        self.type = FilterType.LP
        self.reqData = {FilterData.Aa: None, FilterData.faMin: None,
                        FilterData.Ap: None, FilterData.fpMin: None,
                        FilterData.gain: None,
                        FilterData.Nmax: None, FilterData.Nmin: None, FilterData.Qmax: None}
        self.default = {FilterData.Aa: 30, FilterData.faMin: 10e3,
                        FilterData.Ap: 5, FilterData.fpMin: 9e3,
                        FilterData.gain: 0,
                        FilterData.Nmax: None, FilterData.Nmin: None, FilterData.Qmax: None}

    def validate(self) -> (bool, str):  # Returns true and "Ok" if everything is fine, false and "ErrMsg" if not
        valid = False
        message = "Ok"
        if self.reqData[FilterData.Aa.value] < 0 or self.reqData[FilterData.Ap.value] < 0:
            message = "Error: Enter positive values for Aa and Ap."
        elif self.reqData[FilterData.Aa] < self.reqData[FilterData.Ap]:
            message = "Error: Aa must be greater than Ap."
        elif self.reqData[FilterData.faMin] < self.reqData[FilterData.fpMin]:
            message = "Error: fa must be greater than fp."
        else:
            valid = True
        return valid, message

    def get_template_limits(self):  # Create one set of squares for denormalized graph, and one set for
                                    # normalized graph

        return
