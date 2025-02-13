# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import xlrd #reading from excel
import xlwt
from xlutils.copy import copy
import os
from MenuFunctions import ViewFileName

global read_from
#read_from = 'System_data3.xls'
#file = GetFileName(filext="xls")
file = ViewFileName(filext="xls")
# file = OpenFile()
read_from = file

#-----------------------------------------------------------------------------
# Parameters ----------------------------------------------------------------
#-----------------------------------------------------------------------------

class Parameters:
    'Common base class for all Generators'
    def __init__(self): 
        self.Sbase = 0      # System reference value
        self.Vbase = 0      # Voltage reference level
        self.NumBus = 0     # Number of buses in the system
        self.NumLoad = 0    # Number of load buses in the system
        self.NumGen = 0     # Number of generating buses in the system
        self.NumLin = 0     # Number of transmission lines
  

    def displayParameters(self):
        print('\n','System parameters :')
        print ('Sbase    :','{:6.1f}'.format(self.Sbase),' Vbase   :', '{:6.1f}'.format(self.Vbase))
        print ('NumBus   :', '{:4.0f}'.format(self.NumBus),'   NumLoad : ', '{:4.0f}'.format(self.NumLoad),
             '    NumGen  :', '{:4.0f}'.format(self.NumGen), ' NumLin : ', '{:4.0f}'.format(self.NumLin))
        print(' ')
      
    def read_from_sheet(self):
        wb = xlrd.open_workbook(read_from)
        wr = wb.sheet_by_name('Parameters')
        self.Sbase = float(wr.cell(1, 0).value)
        self.Vbase = float(wr.cell(1, 1).value)
        self.NumBus = int(wr.cell(1, 2).value)
        self.NumLoad = int(wr.cell(1, 3).value)
        self.NumGen = int(wr.cell(1, 4).value)
        self.NumLin = int(wr.cell(1, 5).value)
   

# Buses ----------------------------------------------------------------
class Bus:
    'Common base class for all Buses'
    BusCount = 0
    def __init__(self):
        Bus.BusCount += 1
        self.IntNum = 0
        self.BusNo = 0
        self.Name = ' '
        self.BaskV = 0
        self.BusCode =0.0
        self.GL = 0.0
        self.BL = 0.0
        self.Area = 0.0
        self.Zone = 0.0
        self.Vmag = 0.0
        self.Vang = 0.0
        self.Ipload = -1
        self.Ipgen = -1
        self.plambda = 0
        self.qlambda = 0
   
    def displayCount(self):
        print(' ')
        print ("Total number of Buses %d" % Bus.BusCount)

    def displayBus(self):
        print ('Bus Name :', self.Name,' No :', self.BusNo, " Code : ", self.BusCode)
      
    def read_from_sheet(self,inum):
        wb = xlrd.open_workbook(read_from)
        wr = wb.sheet_by_name('Buses')
        self.BusNo = int(wr.cell(inum, 0).value)
        self.Name = 'Bus'+chr(self.BusNo)
        self.BusCode = int(wr.cell(inum, 1).value)
        self.GL = float(wr.cell(inum, 4).value)
        self.BL = float(wr.cell(inum, 5).value)
        self.Vmag = float(wr.cell(inum, 7).value)
        self.Vang = float(wr.cell(inum, 8).value)
        self.BaskV = float(wr.cell(inum, 9).value)

#-----------------------------------------------------------------------------
# Loads ----------------------------------------------------------------
#-----------------------------------------------------------------------------        

class Load:
    'Common base class for all Loads'
    LoadCount = 0
    def __init__(self):
        Load.LoadCount += 1
        self.IntNum = 0
        self.BusNo = 0 # Note possibly character
        self.LoadId = ' '
        self.LoadStat = 0
        self.Area =0
        self.Zone = 0
        self.Pload = 0.0            # Active power load in pu
        self.Qload = 0.0            # Reactive power load in pu
        self.IP = 0.0
        self.QP = 0.0
        self.YP = 0.0
        self.QP = 0.0
        self.AlfaLoad = 0.0         # Relative increase in reactive load
        self.BetaLoad = 0.0         # Relative increase in active load
        self.AccPload = 0.0         # Accumulated Pload in CPF sequence
        self.AccQload = 0.0         # Accumulated Qload in CPF sequence
   
    def displayCount(self):
     print ("Total number of Loads %d" % Load.LoadCount)

    def displayLoad(self):
      print ('Load No :', '{:4.0f}'.format(self.BusNo),'   Stat(0/1):', '{:4.0f}'.format(self.LoadStat),
             '   Pload : ', '{:6.2f}'.format(self.Pload), ' Qload : ', '{:6.2f}'.format(self.Qload))
      
    def read_from_sheet(self,inum):
        wb = xlrd.open_workbook(read_from)
        wr = wb.sheet_by_name('Loads')
        self.IntNum = int(wr.cell(inum, 0).value)
        self.BusNo = int(wr.cell(inum, 1).value)
        self.LoadId = int(wr.cell(inum, 2).value)
        self.LoadStat = int(wr.cell(inum, 3).value)
        self.Area = int(wr.cell(inum, 4).value)
        self.Zone = int(wr.cell(inum, 5).value)
        self.Pload = float(wr.cell(inum, 6).value)
        self.Qload = float(wr.cell(inum, 7).value)

#-----------------------------------------------------------------------------
# Generators ----------------------------------------------------------------
#-----------------------------------------------------------------------------

class Generator:
    'Common base class for all Generators'
    GenCount = 0
    def __init__(self):
        Generator.GenCount += 1
        self.IntNum = 0
        self.BusNo = 0
        self.GenId = ' '
        self.Pgen = 0.0
        self.Qgen = 0.0
        self.Qmax = 0.0
        self.Qmin = 0.0
        self.Vset = 0.0
        self.mBase = 0.0
        self.istat = 0
        self.Pmax = 0.0
        self.Pmin = 0.0
        self.Cost1 = 0.0
        self.Cost2 = 0.0

        self.BetaGen = 0.0      # Relative increase in generation
        self.AccPgen = 0.0           # Accumulated increase in generation
   
    def displayCount(self):
     print ("Total number of Generators %d" % Generator.GenCount)

    def displayGenerator(self):
      print ('Gen no : ', '{:4.0f}'.format(self.BusNo), ' Pgen :', '{:5.2f}'.format(self.Pgen),
             ' Qgen : ', '{:5.2f}'.format(self.Qgen), ' Qmax :', '{:5.2f}'.format(self.Qmax),
             ' Qmin : ', '{:5.2f}'.format(self.Qmin), ' Vset : ', '{:5.2f}'.format(self.Vset))

      
    def read_from_sheet(self,inum):
        wb = xlrd.open_workbook(read_from)
        wr = wb.sheet_by_name('Generators')
        self.IntNum = int(wr.cell(inum, 0).value)
        self.BusNo = int(wr.cell(inum, 1).value)
        self.GenId = 1
        self.Pgen = float(wr.cell(inum, 2).value)
        self.Qgen = float(wr.cell(inum, 3).value)
        self.Qmax = float(wr.cell(inum, 4).value)
        self.Qmin = float(wr.cell(inum, 5).value)
        self.Vset = float(wr.cell(inum, 6).value)
        self.mBase = float(wr.cell(inum, 7).value)
        self.istat = int(wr.cell(inum, 8).value)
        self.Pmax = float(wr.cell(inum, 9).value)
        self.Pmin = float(wr.cell(inum, 10).value)
        self.Cost1 = float(wr.cell(inum,11).value)
        self.Cost2 =  float(wr.cell(inum,12).value)


   
#-----------------------------------------------------------------------------
# Transmission ----------------------------------------------------------------
#-----------------------------------------------------------------------------       

class Transmission:
    'Common base class for all Transmission lines'
    TransCount = 0
    def __init__(self):
        Transmission.TransCount += 1
        self.IntNum = 0
        self.FromBus = 0
        self.ToBus = 0
        self.R =0.0
        self.X = 0.0
        self.Bc2 = 0.0
        self.RateA = 0.0
        self.RateB = 0.0
        self.RateC = 0.0
        self.Ratio = 0.0
        self.G = 0.0
        self.B = 0.0
        self.Ibstat = 1
   
    def displayCount(self):
     print ("Total number of Transmission Lines %d" % Transmission.TransCount)

    def bij(self, R, X):
        return (1.0/complex(R, X)).imag

    def gij(self, R, X):
        return (1.0/complex(R, X)).real

    def displayTransmission(self):
      print ('Line :', '{:4.0f}'.format(self.IntNum), ' From :', '{:4.0f}'.format(self.FromBus),
             ' To :', '{:4.0f}'.format(self.ToBus), " Res (R) :", '{:7.4f}'.format(self.R),
             " React (X) :", '{:7.4f}'.format(self.X))
      
    def read_from_sheet(self,inum):
        wb = xlrd.open_workbook(read_from)
        wr = wb.sheet_by_name('TR_Lines')
        self.IntNum = int(wr.cell(inum, 0).value)
        self.FromBus = int(wr.cell(inum, 1).value)
        self.ToBus = int(wr.cell(inum, 2).value)
        self.R = float(wr.cell(inum, 3).value)
        self.X = float(wr.cell(inum, 4).value)
        self.Bc2 = float(wr.cell(inum, 5).value)
        self.RateA = float(wr.cell(inum, 6).value)
        self.RateB = float(wr.cell(inum, 7).value)
        self.RateC = float(wr.cell(inum, 8).value)
        self.Ratio = float(wr.cell(inum, 9).value)
        self.Ibstat = int(wr.cell(inum, 11).value)
        self.G = self.gij(self.R, self.X)
        self.B = self.bij(self.R, self.X)

    def jac11DC(self, X, FromBus, ToBus):
        loclist = []
        loclist.append([FromBus, FromBus, 1.0/X])
        loclist.append([ToBus, ToBus, 1.0/X])
        loclist.append([FromBus, ToBus, -1.0/X])
        loclist.append([ToBus, FromBus, -1.0/X])
        return loclist
    
    def jac11DC2(self):
        loclist = []
        loclist.append([self.FromBus, self.FromBus, 1.0/self.X])
        loclist.append([self.ToBus, self.ToBus, 1.0/self.X])
        loclist.append([self.FromBus, self.ToBus, -1.0/self.X])
        loclist.append([self.ToBus, self.FromBus, -1.0/self.X])
        return loclist
