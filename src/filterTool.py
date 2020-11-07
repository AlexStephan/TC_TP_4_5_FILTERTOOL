# Imports
#

from Lib.random import random, seed

# Qt Modules
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget
from PyQt5 import QtWidgets, QtGui
from PyQt5.QtWidgets import QFileDialog,QWidget, QGridLayout,QPushButton, QApplication, QLabel, QCheckBox
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

from src.backend.Filter.LowPass import LowPass
from src.backend.Filter.HighPass import HighPass
from src.backend.Filter.BandPass import BandPass
from src.backend.Filter.BandReject import BandReject
from src.backend.Filter.GroupDelay import GroupDelay

#from src.backend.Approx.Gauss import Gauss
from src.backend.Approx.Legendre import Legendre
from src.backend.Approx.Butterworth import Butterworth

DEBUG = True

t_domain = np.linspace(0,0.001,10000)

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

        seed()

        self.__setConstants()
        self.__init_graphs()
        self.__setCallbacks()
        self.__manageDebug()

        self.__showHideState_CreateNewStage()
        self.__refreshFilterMakerGraphs()
        self.__refreshStagesGraphs()

    def __createFilter(self):
        if self.checkBox_Test.isChecked() == False:
            gain = self.doubleSpinBox_Gain.value()
            Aa = self.doubleSpinBox_Aa.value()
            Ap = self.doubleSpinBox_Ap.value()
            fa_minus = self.doubleSpinBox_fa_minus.value()
            fp_minus = self.doubleSpinBox_fp_minus.value()
            fa_plus = self.doubleSpinBox_fa_plus.value()
            fp_plus = self.doubleSpinBox_fp_plus.value()
            ft = self.doubleSpinBox_ft.value()
            tol = self.doubleSpinBox_Tolerance.value()
            gDelay = self.doubleSpinBox_GroupDelay.value()

            denorm = self.doubleSpinBox_denorm.value()
            if self.checkBox_Nmin.isChecked():
                Nmin = self.spinBox_Nmin.value()
            else:
                Nmin = None
            if self.checkBox_Nmax.isChecked():
                Nmax = self.spinBox_Nmax.value()
            else:
                Nmax = None
            if self.checkBox_Qmax.isChecked():
                Qmax = self.doubleSpinBox_Qmax.value()
            else:
                Qmax = None

            Name = self.lineEdit_name.text()

            filtertype = self.comboBox_filterType.currentIndex()
            if filtertype == filterType.LP.value:
                name_filterType = "LowPass"
                newFilter = LowPass(Aa,fa_minus,Ap,fp_minus,gain,Nmax,Nmin,Qmax,denorm)
            elif filtertype == filterType.HP.value:
                name_filterType = "HighPass"
                newFilter = HighPass(Aa,fa_minus,Ap,fp_minus,gain,Nmax,Nmin,Qmax,denorm)
            elif filtertype == filterType.BP.value:
                name_filterType = "BandPass"
                newFilter = BandPass(Aa,fa_minus,fa_plus,Ap,fp_minus,fp_plus,gain,Nmax,Nmin,Qmax,denorm)
            elif filtertype == filterType.BS.value:
                name_filterType = "BandStop"
                newFilter = BandReject(Aa,fa_minus,fa_plus,Ap,fp_minus,fp_plus,gain,Nmax,Nmin,Qmax,denorm)
            else:
                name_filterType = "GroupDelay"
                newFilter = GroupDelay(ft,gDelay,tol,gain)

            valid,message = newFilter.validate()
            if not valid:
                self.__error_message(message)
                return

            approxtype = self.comboBox_approximation.currentIndex()
            if approxtype == approxTypeALL.Gauss.value:
                name_approxType = "Gauss"
                print("lol")
                return
            elif approxtype == approxTypeALL.Butterworth.value:
                name_approxType = "Butterworth"
                newApprox = Butterworth(newFilter)
            elif approxtype == approxTypeALL.Legendre.value:
                name_approxType = "Legendre"
                newApprox = Legendre(newFilter)
            else:
                self.__error_message("Invalid Approximation Type")
                return

            #TODO VERIFICAR Q LA APROX SEA VALIDA
            w,mag,pha = newApprox.calculate()
            newTransFunc = [w,mag,pha]
            fullname = Name + " - " + name_filterType + " - " + name_approxType
            self.myFilters.append([fullname,newApprox,newTransFunc,True])
            self.comboBox_YourFilters.addItem(fullname)
            self.comboBox_SelectYourFilter.addItem(fullname)
        else:
            newApprox = myFilterTest()
            fullname = "TEST" + str(self.testvar1)
            self.testvar1 += 1
            w,mag,pha = newApprox.calculate()
            newTransFunc = [w,mag,pha]
            self.myFilters.append([fullname,newApprox,newTransFunc,True])
            self.comboBox_YourFilters.addItem(fullname)
            self.comboBox_SelectYourFilter.addItem(fullname)

        self.__refreshFilterMakerGraphs()

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

        self.comboBox_YourFilters.currentIndexChanged.connect(self.__indexChanged_yourFilters)
        self.checkBox_VisibleFilter.clicked.connect(self.__clicked_visibleFilter)
        self.__indexChanged_yourFilters()
        self.pushButton_EraseFilter.clicked.connect(self.__clicked_EraseFilter)

        #####################################################################

        self.comboBox_SelectYourFilter.currentIndexChanged.connect(self.__indexChanged_SelectYourFilter)
        self.__indexChanged_SelectYourFilter()

        #####################################################################
        self.pushButton_TEST.clicked.connect(self.__test)
        self.pushButton_TEST_2.clicked.connect(self.__test2)

    def __indexChanged_yourFilters(self):
        if self.comboBox_YourFilters.currentIndex() == 0:
            self.checkBox_VisibleFilter.setDisabled(True)
        else:
            self.checkBox_VisibleFilter.setDisabled(False)
            i = self.comboBox_YourFilters.currentIndex()-1
            self.checkBox_VisibleFilter.setChecked(self.myFilters[i][3])

    def __indexChanged_SelectYourFilter(self):
        self.__cleanStagesWidgets()

        if self.comboBox_SelectYourFilter.currentIndex() == 0:
            self.checkBox_SelectedFilterVisible.setDisabled(True)
            self.checkBox_SelectedFilterVisible.setChecked(False)
            self.label_SelectedFilterGain.setText("")
            self.label_SelectedFilterK.setText("")
        else:
            self.checkBox_SelectedFilterVisible.setDisabled(False)
            self.checkBox_SelectedFilterVisible.setChecked(True)
            i = self.comboBox_SelectYourFilter.currentIndex()-1
            self.__getAndOrderPolesAndZeros(i)
            self.label_SelectedFilterGain.setText(str(self.myFilters[i][1].get_Gain()))
            self.__refreshStagesGraphs()
            self.__updateStagesAvailable()
            #DO THINGS

    def __cleanStagesWidgets(self):
        self.ComplexPoles = []
        self.RealPoles = []
        self.ComplexZeros = []
        self.RealZeros = []
        self.ComplexZerosAvailable = []
        self.ComplexZerosAvailable = []
        self.ComplexZerosAvailable = []
        self.ComplexZerosAvailable = []
        self.__cleanThisGridLayout(self.gridLayout_ComplexPoles,self.ComplexPolesWidgets)
        self.__cleanThisGridLayout(self.gridLayout_RealPoles, self.RealPolesWidgets)
        self.__cleanThisGridLayout(self.gridLayout_ComplesZeros, self.ComplexZerosWidgets)
        self.__cleanThisGridLayout(self.gridLayout_RealZeros, self.RealZerosWidgets)

    def __cleanThisGridLayout(self,layout: QGridLayout,widgets):
        if DEBUG:
            print(layout.columnCount())
            print(layout.rowCount())
        #x = layout.columnCount()
        #y = layout.rowCount()
        for i in widgets:
            layout.removeWidget(i)
            i.deleteLater()
            del i
        widgets.clear()

    def __updateStagesAvailable(self):
        for i in range(len(self.ComplexPoles)):
            reZ = QLabel('{:.2f}'.format(np.real(self.ComplexPoles[i])))
            #reZ.setMaximumWidth(50)
            self.gridLayout_ComplexPoles.addWidget(reZ,i+1,0)
            self.ComplexPolesWidgets.append(reZ)

            imZ = QLabel('{:.2f}'.format(np.imag(self.ComplexPoles[i])))
            #imZ.setMaximumWidth(50)
            self.gridLayout_ComplexPoles.addWidget(imZ,i+1,1)
            self.ComplexPolesWidgets.append(imZ)

            foZ = QLabel('{:.2f}'.format(np.absolute(self.ComplexPoles[i])/(2*np.pi)))
            #foZ.setMaximumWidth(50)
            self.gridLayout_ComplexPoles.addWidget(foZ,i+1,2)
            self.ComplexPolesWidgets.append(foZ)

            if np.real(self.ComplexPoles[i]) == 0:
                Q = "inf"
            else:
                Q = '{:.2f}'.format(np.abs(self.ComplexPoles[i])/(2*np.abs(np.real(self.ComplexPoles[i]))))
            QZ = QLabel(Q)
            self.gridLayout_ComplexPoles.addWidget(QZ,i+1,3)
            self.ComplexPolesWidgets.append(QZ)

            visible = QCheckBox()
            visible.setMaximumWidth(16)
            visible.setChecked(False)
            visible.setDisabled(True)
            self.gridLayout_ComplexPoles.addWidget(visible,i+1,4)
            self.ComplexPolesWidgets.append(visible)

        for i in range(len(self.RealPoles)):
            reZ = QLabel('{:.2f}'.format(np.real(self.RealPoles[i])))
            #reZ.setMaximumWidth(50)
            self.gridLayout_RealPoles.addWidget(reZ,i+1,0)
            self.RealPolesWidgets.append(reZ)

            foZ = QLabel('{:.2f}'.format(np.absolute(self.RealPoles[i])/(2*np.pi)))
            #foZ.setMaximumWidth(50)
            self.gridLayout_RealPoles.addWidget(foZ,i+1,1)
            self.RealPolesWidgets.append(foZ)

            visible = QCheckBox()
            visible.setMaximumWidth(16)
            visible.setChecked(False)
            visible.setDisabled(True)
            self.gridLayout_RealPoles.addWidget(visible,i+1,2)
            self.RealPolesWidgets.append(visible)

        for i in range(len(self.ComplexZeros)):
            reZ = QLabel('{:.2f}'.format(np.real(self.ComplexZeros[i])))
            #reZ.setMaximumWidth(50)
            self.gridLayout_ComplesZeros.addWidget(reZ,i+1,0)
            self.ComplexZerosWidgets.append(reZ)

            imZ = QLabel('{:.2f}'.format(np.imag(self.ComplexZeros[i])))
            #imZ.setMaximumWidth(50)
            self.gridLayout_ComplesZeros.addWidget(imZ,i+1,1)
            self.ComplexZerosWidgets.append(imZ)

            visible = QCheckBox()
            visible.setMaximumWidth(16)
            visible.setChecked(False)
            visible.setDisabled(True)
            self.gridLayout_ComplesZeros.addWidget(visible,i+1,2)
            self.ComplexZerosWidgets.append(visible)

        for i in range(len(self.RealZeros)):
            reZ = QLabel('{:.2f}'.format(np.real(self.RealZeros[i])))
            #reZ.setMaximumWidth(50)
            self.gridLayout_RealZeros.addWidget(reZ,i+1,0)
            self.RealZerosWidgets.append(reZ)

            visible = QCheckBox()
            visible.setMaximumWidth(16)
            visible.setChecked(False)
            visible.setDisabled(True)
            self.gridLayout_RealZeros.addWidget(visible,i+1,1)
            self.RealZerosWidgets.append(visible)


    def __getAndOrderPolesAndZeros(self,i):
        z, p, Gk = self.myFilters[i][1].get_zpGk()
        self.label_SelectedFilterK.setText(str(Gk))
        for i in z:
            if np.imag(i) == 0:
                self.RealZeros.append(i)
                self.ComplexZerosAvailable.append(False)
            else:
                if (i in self.ComplexZeros) or (np.conjugate(i) in self.ComplexZeros):
                    print("Found z="+str(i))
                else:
                    self.ComplexZeros.append(i)
                    self.ComplexZerosAvailable.append(False)

        for i in p:
            if np.imag(i) == 0:
                self.RealPoles.append(i)
                self.ComplexZerosAvailable.append(False)
            else:
                if i in self.ComplexPoles or np.conjugate(i) in self.ComplexPoles:
                    print("Found p="+str(i))
                else:
                    self.ComplexPoles.append(i)
                    self.ComplexZerosAvailable.append(False)

    def __refreshStagesGraphs(self):
        self.__cleanStagesGraphs()
        if self.checkBox_SelectedFilterVisible.isChecked() and self.comboBox_SelectYourFilter.currentIndex() > 0:
            i = self.comboBox_SelectYourFilter.currentIndex() - 1
            self.__addGraphicsForStages(self.myFilters[i][0]+" - Desired Filter",self.myFilters[i][2])

    def __cleanStagesGraphs(self):
        self.axis_StagesGain.clear()
        self.axis_StagesGain.grid()
        self.canvas_StagesGain.draw()

        self.axis_StagesPhase.clear()
        self.axis_StagesPhase.grid()
        self.canvas_StagesPhase.draw()

    def __addGraphicsForStages(self,name, transfunc):
        self.axis_StagesGain.semilogx(transfunc[0]/(2*np.pi),transfunc[1],label=name)
        self.axis_StagesGain.legend()
        self.canvas_StagesGain.draw()

        self.axis_StagesPhase.semilogx(transfunc[0]/(2*np.pi),transfunc[2],label=name)
        self.axis_StagesPhase.legend()
        self.canvas_StagesPhase.draw()

    def __clicked_visibleFilter(self):
        if self.comboBox_YourFilters.currentIndex() > 0:
            i = self.comboBox_YourFilters.currentIndex()-1
            self.myFilters[i][3] = self.checkBox_VisibleFilter.isChecked()
            self.__refreshFilterMakerGraphs()

    def __clicked_EraseFilter(self):
        if self.comboBox_YourFilters.currentIndex() > 0:
            i = self.comboBox_YourFilters.currentIndex()-1
            self.myFilters.pop(i)
            self.comboBox_YourFilters.removeItem(i+1)
            self.comboBox_SelectYourFilter.removeItem(i+1)
            self.__refreshFilterMakerGraphs()
            self.comboBox_YourFilters.setCurrentIndex(0)
            self.__indexChanged_yourFilters()
            #self.comboBox_SelectYourFilter.setCurrentIndex(0)

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

    def __refreshFilterMakerGraphs(self):
        self.__cleanFilterMakerGraphs()
        for i in self.myFilters:
            if i[3]:
                self.__addGraphicsForFilterMaker(i[0],i[1],i[2])

    def __addGraphicsForFilterMaker(self,name,filter,transfunc):
        self.axis_Magnitude.semilogx(transfunc[0]/(2*np.pi),transfunc[1],label=name)
        self.axis_Magnitude.legend()
        self.canvas_Magnitude.draw()

        self.axis_Phase.semilogx(transfunc[0]/(2*np.pi),transfunc[2],label=name)
        self.axis_Phase.legend()
        self.canvas_Phase.draw()

        z,p,gk=filter.get_zpGk()
        zReal = []
        zImag = []
        for i in z:
            zReal.append(np.real(i))
            zImag.append(np.imag(i))
        pReal = []
        pImag = []
        for i in p:
            pReal.append(np.real(i))
            pImag.append(np.imag(i))

        temp = self.axis_ZerosAndPoles.scatter(pReal,pImag,marker="x",label=name+' - Poles')
        color = temp.get_facecolor()[0]
        if len(zReal) != 0:
            self.axis_ZerosAndPoles.scatter(zReal,zImag,marker="o",label=name+' - Zeros',color=color)
        self.axis_ZerosAndPoles.legend()
        self.canvas_ZerosAndPoles.draw()

        Hs = filter.get_ssTransferFunction()
        impulseResponse = ss.impulse(Hs)
        stepResponse = ss.step(Hs)

        self.axis_ImpulseResponse.plot(impulseResponse[0],impulseResponse[1],label=name)
        self.axis_ImpulseResponse.legend()
        self.canvas_ImpulseResponse.draw()

        self.axis_StepResponse.plot(stepResponse[0],stepResponse[1],label=name)
        self.axis_StepResponse.legend()
        self.canvas_StepResponse.draw()

    def __cleanFilterMakerGraphs(self):
        self.axis_Magnitude.clear()
        self.axis_Magnitude.grid()
        self.canvas_Magnitude.draw()

        self.axis_Attenuation.clear()
        self.axis_Attenuation.grid()
        self.canvas_Attenuation.draw()

        self.axis_NormalizedAttenuation.clear()
        self.axis_NormalizedAttenuation.grid()
        self.canvas_NormalizedAttenuation.draw()

        self.axis_Phase.clear()
        self.axis_Phase.grid()
        self.canvas_Phase.draw()

        self.axis_GroupDelay.clear()
        self.axis_GroupDelay.grid()
        self.canvas_GroupDelay.draw()

        self.axis_ZerosAndPoles.clear()
        self.axis_ZerosAndPoles.grid()
        self.canvas_ZerosAndPoles.draw()

        self.axis_ImpulseResponse.clear()
        self.axis_ImpulseResponse.grid()
        self.canvas_ImpulseResponse.draw()

        self.axis_StepResponse.clear()
        self.axis_StepResponse.grid()
        self.canvas_StepResponse.draw()

        self.axis_Q.clear()
        self.axis_Q.grid()
        self.canvas_Q.draw()

    ##################################################################################
    # CREACION DE STAGES
    ##################################################################################

    def __cancelNewStage(self):
        self.__showHideState_CreateNewStage()
        self.tempPolos = []
        self.tempZeros = []

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

    ###

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

    #####################################################################33

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
        self.myFilters = []
        self.myStages = []
        self.currentFilter = None
        self.currentStage = None

        self.ComplexPoles = []
        self.RealPoles = []
        self.ComplexZeros = []
        self.RealZeros = []
        self.ComplexZerosAvailable = []
        self.ComplexZerosAvailable = []
        self.ComplexZerosAvailable = []
        self.ComplexZerosAvailable = []
        self.ComplexPolesWidgets = []
        self.RealPolesWidgets = []
        self.ComplexZerosWidgets = []
        self.RealZerosWidgets = []

        self.tempPolos = []
        self.tempZeros = []
        ########################################################
        self.testvar1 = 0
        self.testvar2 = 0

    def __error_message(self, description):
        self.errorBox.setWindowTitle("Error")
        self.errorBox.setIcon(self.errorBox.Information)
        self.errorBox.setText(description)
        self.errorBox.exec()

    def __test(self):
        print("Test")
        #self.comboBox_YourFilters.addItem("LOL"+str(self.testvar1))
        #self.comboBox_SelectYourFilter.addItem("LOL"+str(self.testvar1))
        #self.testvar1+=1
        #self.myFilters.append([str(self.testvar1),None,None,True])
        #mfl = myFilterTest()

        for y in reversed(range(self.gridLayout_TEST.rowCount()-1)):
            for x in reversed(range(self.gridLayout_TEST.columnCount()-1)):
                temp = self.gridLayout_TEST.itemAt(x + 3*y)
                self.gridLayout_TEST.removeItem(temp)
                del temp

    def __test2(self):
        #print("Test2")
        #if self.comboBox_YourFilters.count() >= 2 + 1:
        #    self.comboBox_YourFilters.removeItem(2)
        #    self.comboBox_SelectYourFilter.removeItem(2)
        #    self.myFilters.pop(2)
        #else:
        #    print("YOLO")

        #x = [0,1,0,2,3,4,5]
        #print(str(0 in x))

        for x in range(3):
            for y in range(3):
                button = QPushButton(str(str(3*x+y)))
                button.setMaximumWidth(30)
                self.gridLayout_TEST.addWidget(button, x, y)

        print("TEST2")
        print(str(sp.oo))
        print(self.gridLayout_ComplexPoles.rowCount())

    def __manageDebug(self):
        if DEBUG:
            print("DEBUG")
        else:
            self.pushButton_TEST.hide()
            self.pushButton_TEST_2.hide()
            self.checkBox_Test.hide()
            self.line_16.hide()

class myFilterTest(object):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #self.z = [-15000j,15000j,-5000j,5000j,0,-3000,6000]
        #self.p = [-200000,-150000+150000j,-150000-150000j,
        #    -50000+200000j,-50000-200000j,
        #     -20000, -15000 + 15000j, -15000 - 15000j,
        #     -5000 + 20000j, -5000 - 20000j]
        r1 = random()
        r2 = random()
        r3 = random()
        r4 = random()
        r5 = random()
        r6 = random()
        r7 = random()
        self.z = [0, -4000*r2, 5000j*r1, -5000j*r1]
        self.p = [-4000+7000j*r3,-4000-7000j*r3,-6000*r4, -2000*r5]
        self.num = np.polynomial.polynomial.polyfromroots(self.z).tolist()
        self.den = np.polynomial.polynomial.polyfromroots(self.p).tolist()
        self.num.reverse()
        self.den.reverse()
        print(self.num)
        print(self.den)
        self.k = 1*r6
        self.gain = 150*r7

        self.Hs=ss.TransferFunction(np.real(self.num)*self.k*10**(self.gain/20),np.real(self.den))
        print(self.Hs)
        print(self.Hs.zeros)
        print(self.Hs.poles)

    def calculate(self):
        self.sys = ss.ZerosPolesGain(self.z,self.p,self.k)
        self.w,self.mag,self.pha = ss.bode(self.sys,np.logspace(0,9,10000)*(2*np.pi))
        self.mag += self.gain
        return self.w,self.mag,self.pha

    def get_zpGk(self):
        return self.z,self.p, self.k * np.power(10,self.gain/20)

    def get_ssTransferFunction(self):
        return self.Hs

    def get_Gain(self):
        return self.gain
