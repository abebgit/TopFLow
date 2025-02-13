#!/usr/bin/python
# Copyright (c) 2021, Olav B. Fosso, NTNU
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
#     * Redistributions of source code must retain the above copyright notice,
#       this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright notice,
#       this list of conditions and the following disclaimer in the documentation
#       and/or other materials provided with the distribution.

import numpy as np
#from DistribObjects_v2 import *
import pandas as pd
import pandapower as pp
from pandapower.plotting import create_generic_coordinates

#create empty net
net = pp.create_empty_network()

from MenuFunctions_v2 import ViewFileName

def BuildSystem3():

# Read the data from an xls-fil and reformats it into Pandapower-format

    BusList = []
    LineList = []
    LoadList = []
    BusType = []
 #   basv = 22.0         # Must be adapted to the real system
 #   sbase = 10.0         # Must be adapted to the real system
    file = ViewFileName(filext="xls")
    xls = pd.ExcelFile(file)
# Read parameters
    df2 = pd.read_excel(xls, 'Parameters', usecols="A:F")
    values = df2.values
    basv = float(values[0, 1])
    sbase = float(values[0, 0])
    nbus = int(values[0,2])
    nload = int(values[0,3])
    ngen = int(values[0,4])
    nlin = int(values[0,5])
    zref = basv * basv / sbase


# Read bus data
    df2 = pd.read_excel(xls, 'Buses', usecols="A:M")
    values = df2.values

# Read Bus data  --------------------------------------------
    iloop = 0
    # print(' ')
    if nbus != len(values):
        print("Matmitch between specified and entered number of buses")
    while iloop < nbus:
#        Enter all bus-data
        BusList.append(pp.create_bus(net, vn_kv=float(values[iloop, 9]), name='Bus'+str(int(values[iloop, 0])), index=None, geodata=None, type='b', zone=None,
                              in_service=True, max_vm_pu=float(values[iloop, 11]), min_vm_pu=float(values[iloop, 12]), coords=None))
        BusType.append(int(values[iloop, 1]))
        iloop += 1

    iloop = 0
# Read bus data
    df2 = pd.read_excel(xls, 'Loads', usecols="B:K")
    values = df2.values


    if nload != len(values):
       print("Matmitch between specified and entered number of loads")
    while iloop < nload:
# Enter Load-data
        print(values[iloop, 0]-1, len(values))
        LoadList.append(pp.create_load(net, bus=int(values[iloop, 0])-1, p_mw= float(values[iloop, 5])*sbase, q_mvar=float(values[iloop, 6])*sbase,
                       const_z_percent=0, const_i_percent=0, sn_mva=sbase, name="Load" + str(int(values[iloop, 0])), scaling=1.0,
                       index=None, in_service=True, type='wye'))
        iloop += 1

# Define an external grid
#    pp.create_ext_grid(net, bus=BusList[0], vm_pu=1.0, va_degree=0.0)

    iloop = 0
# Read Generator data
    df2 = pd.read_excel(xls, 'Generators', usecols="B:R")
    values = df2.values


    if ngen != len(values):
      print("Matmitch between specified and entered number of generators")
    while iloop < ngen:
        busref = int(values[iloop,0]) - 1
        print(BusType[int(values[iloop, 0]) - 1])
        if BusType[int(values[iloop,0]) - 1] == 3:
            pp.create_gen(net, bus=int(values[iloop,0]) - 1, p_mw = float(values[iloop, 1])* sbase, vm_pu= float(values[iloop, 5]), slack=True, name= "Gen" + str(int(values[iloop, 0])))  # Create a generator (Type 3)
        else:
            pp.create_gen(net, bus=int(values[iloop, 0]) - 1, p_mw = float(values[iloop, 1])* sbase, vm_pu= float(values[iloop, 5]), slack=False, name= "Gen" + str(int(values[iloop, 0])))   # Create a generator (Type 2)
        iloop += 1

# Transmission lines
    df2 = pd.read_excel(xls, 'TR_Lines', usecols="B:L")
    values = df2.values
    # Read Line data  --------------------------------------------
    iloop = 0
    # print(' ')


    if nlin != len(values):
        print("Matmitch between specified and entered number of buses")
    while iloop < nlin:

        pp.create_line_from_parameters(net, from_bus = int(values[iloop, 0])-1, to_bus = int(values[iloop, 1]) -1, length_km=1.0,
                                       r_ohm_per_km = float(values[iloop, 2])*zref , x_ohm_per_km = float(values[iloop, 3])*zref, c_nf_per_km = float(values[iloop, 4])*10**9/(314*zref), max_i_ka=10,
                                        name="Line" + str(iloop), index=None, type=None, geodata=None, in_service=True, df=1.0, parallel=1,
                                        g_us_per_km=0.0, max_loading_percent=80)
        iloop += 1



    return BusList, LineList, LoadList

BusList, LineList, LoadList = BuildSystem3()

#pp.replace_impedance_by_line(net)      # Alternative representation
#net = create_generic_coordinates(net)  # For creating coordinates.

# Small test
# Do a loadflow and save the results
# This could have been done a a separate file calling the others.
#

# Run power flow
pp.runpp(net)

# Export results
pp.to_excel(net, filename="testpandapower.xlsx", include_empty_tables=False, include_results=True)
