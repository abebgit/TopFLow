# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import xlrd #reading from excel
import xlwt
from xlutils.copy import copy
import os

global read_from
#read_from = 'System_data3.xls'

# Lists on the global level keeping the objects
BusList = []
GenList = []
LoadList = []
TransList = []

inum = 0

# All object definations of the Power System

from PowerSystemObjects_psse import *

from AC_CPF_Power_Flows import ACLoadFlow
from BuildSystem import *



#
# Local functions
#

def compvar(var1, var2,num1,num2):
    iloop = num1
    while iloop < num2:
        print('{:6.0f}'.format(iloop), '   Var1 :','{:10.5f}'.format(var1[iloop]), '   Var2 :','{:10.5f}'.format(var2[iloop]),
              '   Diff :','{:10.7f}'.format(var1[iloop]-var2[iloop]))
        iloop += 1
        

# ----------------------------------------------------------
#
# Main program ---------------------------------------------
#
# -----------------------------------------------------------

print('\n','******** System setup ********')

param, BusList, GenList, LoadList, TransList = BuildSystem(DispData=False, BusList = BusList, GenList=GenList, LoadList= LoadList, TransList=TransList)

print('\n','******** System setup - completed ********')  

acflow = ACLoadFlow(BusList,GenList,LoadList,TransList)
jacobi, busjac = acflow.SolveAC(sparse=False,itmax=4)
acflow.UpdateGen()      # Update generation on slackbus

acflow.dispVolt()   # Display voltages
acflow.dispVolt(fromBus=10, toBus=22, tpres=True)       # Display in table format (only 13 buses at a time)

# There are a lot of functions for displaying results. Most of them are reached with acflow.dispxxxx (look in AC_CPF_Power_Flow.py






