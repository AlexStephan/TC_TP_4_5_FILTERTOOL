from src.backend.Filter.Filter import *
from src.backend.Filter.TemplateLimit import *


class GroupDelay(Filter):
    def __init__(self):
        self.type = FilterType.GD
        self.reqData = {FilterData.ft: None,
                        FilterData.GD: None,
                        FilterData.tolerance: None,
                        FilterData.gain: None}
        self.default = {FilterData.ft: 10e3,            #In Hz
                        FilterData.GD: 100,             #In useg
                        FilterData.tolerance: 5,        #In %
                        FilterData.gain: 0}             #In dB

    def validate(self) -> (bool, str):  # Returns true and "Ok" if everything is fine, false and "ErrMsg" if not
        valid = False
        message = "Ok"
        if self.reqData[FilterData.tolerance.value] < 0 or self.reqData[FilterData.tolerance.value] > 100:
            message = "Error: Tolerance must be between 0 and 100."
        elif self.reqData[FilterData.GD.value] < 0:
            message = "Error: Group Delay must be positive."
        else:
            valid = True
        return valid, message

    def get_template_limits(self):  # Create one set of squares for denormalized graph, and one set for
                                    # normalized graph
        ft = self.reqData[FilterData.ft.value]
        maxGD = self.reqData[FilterData.GD.value]*(1-tol)
        GD = self.reqData[FilterData.GD.value]
        tol = self.reqData[FilterData.tolerance.value]

        denormLimit1 = Limit(Dot(0, maxGD), Dot(ft, maxGD), Dot(0, -1e9), Dot(ft, -1e9))
        denormLimit = [denormLimit1]

        normLimit1 = Limit(Dot(0, 1-tol), Dot(ft*GD*1e-6, 1-tol), Dot(0, -1e9), Dot(ft*GD*1e-6, -1e9))
        normLimit = [normLimit1]

        return [denormLimit, normLimit]
