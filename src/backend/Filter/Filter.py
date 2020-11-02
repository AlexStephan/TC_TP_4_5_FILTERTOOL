from enum import Enum


class FilterData(Enum):
    Ap = "Ap"
    fpMin = "fp-"
    fpMax = "fp+"

    Aa = "Aa"
    faMin = "fa-"
    faMax = "fa+"

    GD = "Group delay"
    ft = "ft"
    tolerance = "Tolerance"
    gain = "Gain"


class FilterType(Enum):
    LP = "Low Pass"
    HP = "High Pass"
    BP = "Band Pass"
    BR = "Band Reject"
    GD = "Group Delay"


class Filter(object):

    def __init__(self, type: FilterType):
        self.type = type
        self.reqData = {}  # Data for implementation
        self.norm = {"Poles": [], "Zeros": [], "GainDB": None,
                     "Order": None}  # Order and Gain are floats, Poles and Zeros are Array#
        self.denorm = {"Poles": [], "Zeros": [], "GainDB": None, "Order": None,
                       "Q": None}  # Q is inherent to denormalized transfer
        self.default = []  # Default values are given at assignment of filter type

    def get_type(self):
        return self.type

    def get_reqData(self):
        return self.reqData

    def get_norm(self):
        return self.norm

    def get_denorm(self):
        return self.denorm

    def get_default(self):
        return self.denorm
