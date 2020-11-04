# Imports
#

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

# My Own Modules
from src.backend.Stages.SimpleHs import SimpleHs
from src.backend.Filter.Filter import FilterData

class filterType(Enum):
    LP=0
    HP=1
    BP=2
    BS=3
    GD=4

class approxTypeALL(Enum):
    Butterworth = 0
    Gauss = 1
    Legendre = 2

class FilterTool(QWidget,Ui_Form):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setupUi(self)
        self.setWindowTitle("TP GRUPAL 1 - TEORÍA DE CIRCUITOS")
        self.setWindowIcon(QtGui.QIcon('py.png'))

        self.__setConstants()
        self.__init_graphs()
        self.__setCallbacks()

        self.__showHideState_CreateNewStage()

    def __createFilter(self):
        lol = SimpleHs([2+1j,2-1j],[1,4])

    def __setCallbacks(self):
        self.pushButton_CanceNewStage.clicked.connect(self.__cancelNewStage)
        self.pushButton_CreateNewStage.clicked.connect(self.__crateNewStage)
        self.pushButton_1stOrderPole.clicked.connect(self.__clicked_1stOrderPole)
        self.pushButton_SecondOrderPoles.clicked.connect(self.__clicked_secondOrderPoles)
        self.pushButton_ComplexPoles.clicked.connect(self.__clicked_complexPoles)
        self.pushButton_RealPoles.clicked.connect(self.__clicked_realPoles)
        self.pushButton_SelectRealPole.clicked.connect(self.__selectRealPole)
        self.pushButton_SelectRealPoles.clicked.connect(self.__selectRealPoles)
        self.pushButton_SelectComplexPoles.clicked.connect(self.__selectComplexPoles)
        self.pushButton_FINISH1stOrder.clicked.connect(self.__finish1stOrderPole)
        self.pushButton_AddZero_1stOrder.clicked.connect(self.__clicked_addZero)
        self.pushButton_FINISH2ndOrder.clicked.connect(self.__finish2ndOrderPole)
        self.pushButton_Add1Zero.clicked.connect(self.__clicked_add1Zero)
        self.pushButton_Add2Zeros.clicked.connect(self.__clicked_add2Zeros)
        self.pushButton_ComplexConjZeros.clicked.connect(self.__clicked_complexZeros)
        self.pushButton_RealZeros.clicked.connect(self.__clicked_realZeros)
        self.pushButton_SelectRealZero.clicked.connect(self.__finishRealZero)
        self.pushButton_SelectRealZeros.clicked.connect(self.__finishRealZeros)
        self.pushButton_SelectComplexZeros.clicked.connect(self.__finishComplexZeros)

        self.comboBox_filterType.currentIndexChanged.connect(self.__showAndHideParameters)
        self.__showAndHideParameters()

        self.pushButton_createFilter.clicked.connect(self.__createFilter)

    def __cancelNewStage(self):
        self.__showHideState_CreateNewStage()

    def __crateNewStage(self):
        self.__showHideState_SelectOrder()

    def __clicked_1stOrderPole(self):
        self.__showHideState_SelectRealPole()

    def __clicked_secondOrderPoles(self):
        self.__showHideState_SelectTypeOfPoles()

    def __clicked_complexPoles(self):
        self.__showHideState_SelectComplexPoles()

    def __clicked_realPoles(self):
        self.__showHideState_SelectRealPoles()

    def __init_graphs(self):
        self.figure_Magnitude = Figure()
        self.canvas_Magnitude = FigureCanvas(self.figure_Magnitude)
        self.index_Magnitude = self.stackedWidget_Magnitude.addWidget(self.canvas_Magnitude)
        self.stackedWidget_Magnitude.setCurrentIndex(self.index_Magnitude)
        self.toolbar_Magnitude = NavigationToolbar(self.canvas_Magnitude,self)
        self.horizontalLayout_Magnitude.addWidget(self.toolbar_Magnitude)
        self.axis_Magnitude = self.figure_Magnitude.add_subplot()

        self.figure_Attenuation = Figure()
        self.canvas_Attenuation = FigureCanvas(self.figure_Attenuation)
        self.index_Attenuation = self.stackedWidget_Attenuation.addWidget(self.canvas_Attenuation)
        self.stackedWidget_Attenuation.setCurrentIndex(self.index_Attenuation)
        self.toolbar_Attenuation = NavigationToolbar(self.canvas_Attenuation,self)
        self.horizontalLayout_Attenuation.addWidget(self.toolbar_Attenuation)
        self.axis_Attenuation = self.figure_Attenuation.add_subplot()

        self.figure_NormalizedAttenuation = Figure()
        self.canvas_NormalizedAttenuation = FigureCanvas(self.figure_NormalizedAttenuation)
        self.index_NormalizedAttenuation = self.stackedWidget_NormalizedAttenuation.addWidget(self.canvas_NormalizedAttenuation)
        self.stackedWidget_NormalizedAttenuation.setCurrentIndex(self.index_NormalizedAttenuation)
        self.toolbar_NormalizedAttenuation = NavigationToolbar(self.canvas_NormalizedAttenuation,self)
        self.horizontalLayout_NormalizedAttenuation.addWidget(self.toolbar_NormalizedAttenuation)
        self.axis_NormalizedAttenuation = self.figure_NormalizedAttenuation.add_subplot()

        self.figure_Phase = Figure()
        self.canvas_Phase = FigureCanvas(self.figure_Phase)
        self.index_Phase = self.stackedWidget_Phase.addWidget(self.canvas_Phase)
        self.stackedWidget_Phase.setCurrentIndex(self.index_Phase)
        self.toolbar_Phase = NavigationToolbar(self.canvas_Phase,self)
        self.horizontalLayout_Phase.addWidget(self.toolbar_Phase)
        self.axis_Phase = self.figure_Phase.add_subplot()

        self.figure_GroupDelay = Figure()
        self.canvas_GroupDelay = FigureCanvas(self.figure_GroupDelay)
        self.index_GroupDelay = self.stackedWidget_GroupDelay.addWidget(self.canvas_GroupDelay)
        self.stackedWidget_GroupDelay.setCurrentIndex(self.index_GroupDelay)
        self.toolbar_GroupDelay = NavigationToolbar(self.canvas_GroupDelay,self)
        self.horizontalLayout_GroupDelay.addWidget(self.toolbar_GroupDelay)
        self.axis_GroupDelay = self.figure_GroupDelay.add_subplot()

        self.figure_ZerosAndPoles = Figure()
        self.canvas_ZerosAndPoles = FigureCanvas(self.figure_ZerosAndPoles)
        self.index_ZerosAndPoles = self.stackedWidget_ZerosAndPoles.addWidget(self.canvas_ZerosAndPoles)
        self.stackedWidget_ZerosAndPoles.setCurrentIndex(self.index_ZerosAndPoles)
        self.toolbar_ZerosAndPoles = NavigationToolbar(self.canvas_ZerosAndPoles,self)
        self.horizontalLayout_ZerosAndPoles.addWidget(self.toolbar_ZerosAndPoles)
        self.axis_ZerosAndPoles = self.figure_ZerosAndPoles.add_subplot()

        self.figure_ImpulseResponse = Figure()
        self.canvas_ImpulseResponse = FigureCanvas(self.figure_ImpulseResponse)
        self.index_ImpulseResponse = self.stackedWidget_ImpulseResponse.addWidget(self.canvas_ImpulseResponse)
        self.stackedWidget_ImpulseResponse.setCurrentIndex(self.index_ImpulseResponse)
        self.toolbar_ImpulseResponse = NavigationToolbar(self.canvas_ImpulseResponse,self)
        self.horizontalLayout_ImpulseResponse.addWidget(self.toolbar_ImpulseResponse)
        self.axis_ImpulseResponse = self.figure_ImpulseResponse.add_subplot()

        self.figure_StepResponse = Figure()
        self.canvas_StepResponse = FigureCanvas(self.figure_StepResponse)
        self.index_StepResponse = self.stackedWidget_StepResponse.addWidget(self.canvas_StepResponse)
        self.stackedWidget_StepResponse.setCurrentIndex(self.index_StepResponse)
        self.toolbar_StepResponse = NavigationToolbar(self.canvas_StepResponse,self)
        self.horizontalLayout_StepResponse.addWidget(self.toolbar_StepResponse)
        self.axis_StepResponse = self.figure_StepResponse.add_subplot()

        self.figure_Q = Figure()
        self.canvas_Q = FigureCanvas(self.figure_Q)
        self.index_Q = self.stackedWidget_Q.addWidget(self.canvas_Q)
        self.stackedWidget_Q.setCurrentIndex(self.index_Q)
        self.toolbar_Q = NavigationToolbar(self.canvas_Q,self)
        self.horizontalLayout_Q.addWidget(self.toolbar_Q)
        self.axis_Q = self.figure_Q.add_subplot()

        self.figure_StagesGain = Figure()
        self.canvas_StagesGain = FigureCanvas(self.figure_StagesGain)
        self.index_StagesGain = self.stackedWidget_StagesGain.addWidget(self.canvas_StagesGain)
        self.stackedWidget_StagesGain.setCurrentIndex(self.index_StagesGain)
        self.toolbar_StagesGain = NavigationToolbar(self.canvas_StagesGain,self)
        self.horizontalLayout_StagesGain.addWidget(self.toolbar_StagesGain)
        self.axis_StagesGain = self.figure_StagesGain.add_subplot()

        self.figure_StagesPhase = Figure()
        self.canvas_StagesPhase = FigureCanvas(self.figure_StagesPhase)
        self.index_StagesPhase = self.stackedWidget_StagesPhase.addWidget(self.canvas_StagesPhase)
        self.stackedWidget_StagesPhase.setCurrentIndex(self.index_StagesPhase)
        self.toolbar_StagesPhase = NavigationToolbar(self.canvas_StagesPhase,self)
        self.horizontalLayout_StagesPhase.addWidget(self.toolbar_StagesPhase)
        self.axis_StagesPhase = self.figure_StagesPhase.add_subplot()

    def __selectRealPole(self):
        self.__showHideState_1stOrderPoleReady()

    def __selectRealPoles(self):
        self.__showHideState_2ndOrderPolesReady()

    def __selectComplexPoles(self):
        self.__showHideState_2ndOrderPolesReady()

    def __finish1stOrderPole(self):
        self.__showHideState_CreateNewStage()

    def __clicked_addZero(self):
        self.__showHideState_SelectRealZero()

    def __finish2ndOrderPole(self):
        self.__showHideState_CreateNewStage()

    def __clicked_add1Zero(self):
        self.__showHideState_SelectRealZero()

    def __clicked_add2Zeros(self):
        self.__showHideState_SelectTypeOfZeros()

    def __clicked_complexZeros(self):
        self.__showHideState_SelectComplexZeros()

    def __clicked_realZeros(self):
        self.__showHideState_SelectRealZeros()

    def __finishRealZero(self):
        self.__showHideState_CreateNewStage()

    def __finishRealZeros(self):
        self.__showHideState_CreateNewStage()

    def __finishComplexZeros(self):
        self.__showHideState_CreateNewStage()

    def __showHideState_CreateNewStage(self):
        self.__showAndHideButtons([1, 0,0,0, 0,0,0, 0,0, 0,0,0, 0,0, 0,0,0, 0,0,0,0, 0,0,0, 0,0, 0,0,0, 0,0])

    def __showHideState_SelectOrder(self):
        self.__showAndHideButtons([0, 1,1,1, 0,0,0, 0,0, 0,0,0, 0,0, 0,0,0, 0,0,0,0, 0,0,0, 0,0, 0,0,0, 0,0])

    def __showHideState_SelectTypeOfPoles(self):
        self.__showAndHideButtons([0, 0,0,0, 1,1,1, 0,0, 0,0,0, 0,0, 0,0,0, 0,0,0,0, 0,0,0, 0,0, 0,0,0, 0,0])

    def __showHideState_SelectRealPole(self):
        self.__showAndHideButtons([0, 0,0,0, 0,0,0, 1,1, 0,0,0, 0,0, 0,0,0, 0,0,0,0, 0,0,0, 0,0, 0,0,0, 0,0])

    def __showHideState_SelectRealPoles(self):
        self.__showAndHideButtons([0, 0,0,0, 0,0,0, 0,0, 1,1,1, 0,0, 0,0,0, 0,0,0,0, 0,0,0, 0,0, 0,0,0, 0,0])

    def __showHideState_SelectComplexPoles(self):
        self.__showAndHideButtons([0, 0,0,0, 0,0,0, 0,0, 0,0,0, 1,1, 0,0,0, 0,0,0,0, 0,0,0, 0,0, 0,0,0, 0,0])

    def __showHideState_1stOrderPoleReady(self):
        self.__showAndHideButtons([0, 0,0,0, 0,0,0, 0,0, 0,0,0, 0,0, 1,1,1, 0,0,0,0, 0,0,0, 0,0, 0,0,0, 0,0])

    def __showHideState_2ndOrderPolesReady(self):
        self.__showAndHideButtons([0, 0,0,0, 0,0,0, 0,0, 0,0,0, 0,0, 0,0,0, 1,1,1,1, 0,0,0, 0,0, 0,0,0, 0,0])

    def __showHideState_SelectTypeOfZeros(self):
        self.__showAndHideButtons([0, 0,0,0, 0,0,0, 0,0, 0,0,0, 0,0, 0,0,0, 0,0,0,0, 1,1,1, 0,0, 0,0,0, 0,0])

    def __showHideState_SelectRealZero(self):
        self.__showAndHideButtons([0, 0,0,0, 0,0,0, 0,0, 0,0,0, 0,0, 0,0,0, 0,0,0,0, 0,0,0, 1,1, 0,0,0, 0,0])

    def __showHideState_SelectRealZeros(self):
        self.__showAndHideButtons([0, 0,0,0, 0,0,0, 0,0, 0,0,0, 0,0, 0,0,0, 0,0,0,0, 0,0,0, 0,0, 1,1,1, 0,0])

    def __showHideState_SelectComplexZeros(self):
        self.__showAndHideButtons([0, 0,0,0, 0,0,0, 0,0, 0,0,0, 0,0, 0,0,0, 0,0,0,0, 0,0,0, 0,0, 0,0,0, 1,1])

    def __showAndHideButtons(self, arr):
        length = len(self.itemList)
        for i in range(length):
            self.__showOrHideSingle(self.itemList[i],arr[i])

    def __showAndHideParameters(self):
        filterT = self.comboBox_filterType.currentIndex()
        if filterT is filterType.LP.value or filterT is filterType.HP.value:
            self.label_5.show()
            self.label_6.show()
            self.label_7.show()
            self.label_8.show()
            self.doubleSpinBox_Aa.show()
            self.doubleSpinBox_Ap.show()
            self.doubleSpinBox_fa_minus.show()
            self.doubleSpinBox_fp_minus.show()

            self.label_9.hide()
            self.label_10.hide()
            self.doubleSpinBox_fa_plus.hide()
            self.doubleSpinBox_fp_plus.hide()

            self.label_47.hide()
            self.label_49.hide()
            self.label_51.hide()
            self.doubleSpinBox_ft.hide()
            self.doubleSpinBox_Tolerance.hide()
            self.doubleSpinBox_GroupDelay.hide()
        elif filterT is filterType.BP.value or filterT is filterType.BS.value:
            self.label_5.show()
            self.label_6.show()
            self.label_7.show()
            self.label_8.show()
            self.doubleSpinBox_Aa.show()
            self.doubleSpinBox_Ap.show()
            self.doubleSpinBox_fa_minus.show()
            self.doubleSpinBox_fp_minus.show()

            self.label_9.show()
            self.label_10.show()
            self.doubleSpinBox_fa_plus.show()
            self.doubleSpinBox_fp_plus.show()

            self.label_47.hide()
            self.label_49.hide()
            self.label_51.hide()
            self.label_47.hide()
            self.label_49.hide()
            self.label_51.hide()
            self.doubleSpinBox_ft.hide()
            self.doubleSpinBox_Tolerance.hide()
            self.doubleSpinBox_GroupDelay.hide()
        else:
            self.label_5.hide()
            self.label_6.hide()
            self.label_7.hide()
            self.label_8.hide()
            self.doubleSpinBox_Aa.hide()
            self.doubleSpinBox_Ap.hide()
            self.doubleSpinBox_fa_minus.hide()
            self.doubleSpinBox_fp_minus.hide()

            self.label_9.hide()
            self.label_10.hide()
            self.doubleSpinBox_fa_plus.hide()
            self.doubleSpinBox_fp_plus.hide()

            self.label_47.show()
            self.label_49.show()
            self.label_51.show()
            self.doubleSpinBox_ft.show()
            self.doubleSpinBox_Tolerance.show()
            self.doubleSpinBox_GroupDelay.show()

    def __showOrHideSingle(self, item, show):
        if show:
            item.show()
        else:
            item.hide()

    def __setConstants(self):
        self.itemList=[
            self.pushButton_CreateNewStage,

            self.label_33,
            self.pushButton_1stOrderPole,
            self.pushButton_SecondOrderPoles,

            self.label_34,
            self.pushButton_ComplexPoles,
            self.pushButton_RealPoles,

            self.comboBox_SelectRealPole,
            self.pushButton_SelectRealPole,

            self.comboBox_Select1stRealPole,
            self.comboBox_Select2ndRealPole,
            self.pushButton_SelectRealPoles,

            self.comboBox_SelectComplexPoles,
            self.pushButton_SelectComplexPoles,

            self.label_37,
            self.pushButton_FINISH1stOrder,
            self.pushButton_AddZero_1stOrder,

            self.label_38,
            self.pushButton_FINISH2ndOrder,
            self.pushButton_Add1Zero,
            self.pushButton_Add2Zeros,

            self.label_39,
            self.pushButton_ComplexConjZeros,
            self.pushButton_RealZeros,

            self.comboBox_SelectRealZero,
            self.pushButton_SelectRealZero,

            self.comboBox_Select1stRealZero,
            self.comboBox_Select2ndRealZero,
            self.pushButton_SelectRealZeros,

            self.comboBox_SelectComplexZeros,
            self.pushButton_SelectComplexZeros
        ]
        self.errorBox = QtWidgets.QMessageBox()

    def __error_message(self, description):
        self.errorBox.setWindowTitle("Error")
        self.errorBox.setIcon(self.errorBox.Information)
        self.errorBox.setText(description)
        self.errorBox.exec()