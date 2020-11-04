from Lib.random import random

# Qt Modules
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget
from PyQt5 import QtWidgets, QtGui
from PyQt5.QtWidgets import QFileDialog
from src.ui.filterToolGUI_v2 import Ui_Form

# Matplotlib Modules
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

# Python Modules
import numpy as np
import scipy.signal as ss
from enum import Enum
import array as arr

# SymPy modules

import sympy as sp
from sympy.parsing.sympy_parser import parse_expr
from sympy import I

class SimpleHs(object):
    def __init__(self, zeros, poles, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if len(zeros) == 0:
            self.numerator=[1]
        elif len(zeros) == 1:
            self.numerator=[1, np.real(-zeros[0])]
        else:
            self.numerator=[1, np.real(-zeros[0]-zeros[1]), np.real(zeros[0]*zeros[1])]

        if len(poles) == 0:
            self.denominator=[1]
        elif len(poles) == 1:
            self.denominator=[1, np.real(-poles[0])]
        else:
            self.denominator=[1, np.real(-poles[0]-poles[1]), np.real(poles[0]*poles[1])]

        #print(self.numerator)
        #print(self.denominator)
        #print(np.roots(self.numerator))
        #print(np.roots(self.denominator))

        self.K = 1
        self.visible = True
        self.order = self.__order(self.denominator)
        if self.order == 2:
            if self.denominator[2] != 0:
                self.w0 = np.sqrt(np.absolute(self.denominator[2]))
                if self.denominator[1] != 0:
                    self.Q = np.absolute(self.w0/self.denominator[1])
                else:
                    self.Q = None
            else:
                self.w0 = None
                self.Q = None
        elif self.order == 1:
            self.Q = None
            if self.denominator[1] != 0:
                self.w0 = np.absolute(self.denominator[1])
            else:
                self.w0 = None
        else:
            self.Q = None
            self.w0 = None

    def validate(self):
        if self.__order(self.numerator) > self.__order(self.denominator):
            return False
        elif self.__order(self.denominator) < 1:
            return False
        elif self.__order(self.numerator) > 2 or self.__order(self.denominator) > 2:
            return False
        else:
            return True

    def getData(self):
        if self.validate() == True:
            return self.K,self.visible,self.order,self.w0/(2*np.pi),self.Q
        else:
            return None,None,None,None,None

    def updateGain(self,gain):
        self.K = np.power(10,gain/20)

    def __order(self, var):
        return len(var)-1
