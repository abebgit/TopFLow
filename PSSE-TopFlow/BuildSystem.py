# -*- coding: utf-8 -*-

from PowerSystemObjects_psse import *
# ---------------------------------------------------------------
# Read data from the specified file in csv-format
#
# Build the system description
# -------------------------------------------------------------------

def BuildSystem(DispData=False,BusList=[], GenList=[], LoadList=[], TransList=[]):
    param = Parameters()
    param.read_from_sheet()
    param.displayParameters()

    # Read Bus data  --------------------------------------------
    iloop = 0
    print(' ')
    while iloop < param.NumBus:
        BusList.append(Bus())
        BusList[iloop].read_from_sheet(iloop + 1)
        if DispData:
            BusList[iloop].displayBus()
        iloop = iloop + 1

    # Read Generator data  ---------------------------------------
    iloop = 0
    #    print(' ')
    while iloop < param.NumGen:
        GenList.append(Generator())
        GenList[iloop].read_from_sheet(iloop + 1)
        if DispData:
            GenList[iloop].displayGenerator()
        iloop = iloop + 1

    # Read Load data  ---------------------------------------------
    iloop = 0
    #    print(' ')
    while iloop < param.NumLoad:
        LoadList.append(Load())
        LoadList[iloop].read_from_sheet(iloop + 1)
        if DispData:
            LoadList[iloop].displayLoad()
        iloop = iloop + 1

    # Read Transmission line data -----------------------------------
    iloop = 0
    #    print(' ')
    while iloop < param.NumLin:
        TransList.append(Transmission())
        TransList[iloop].read_from_sheet(iloop + 1)
        if DispData:
            TransList[iloop].displayTransmission()
        iloop = iloop + 1

    #   print(BusList, GenList, LoadList, TransList)

    # Set up internal pointers

    iloop = 0
    while iloop < len(LoadList):
        ibus = LoadList[iloop].BusNo
        BusList[ibus - 1].Ipload = iloop
        #        print('Ipload', BusList[ibus - 1].Ipload)
        iloop += 1

    iloop = 0
    while iloop < len(GenList):
        ibus = GenList[iloop].BusNo
        BusList[ibus - 1].Ipgen = iloop
        #        print('Ipgeb', BusList[ibus - 1].Ipgen)
        iloop += 1

    return param, BusList, GenList, LoadList, TransList

