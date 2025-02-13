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
    basv = 22.0         # Must be adapted to the real system
    sbase = 10.0         # Must be adapted to the real system
    zref = basv*basv/sbase
    file = ViewFileName(filext="xls")
    xls = pd.ExcelFile(file)
    df2 = pd.read_excel(xls, 'bus', usecols="A:M")
    values = df2.values
    # Read Bus data  --------------------------------------------
    iloop = 0
    # print(' ')
    while iloop < len(values):
#        Enter all bus-data
        BusList.append(pp.create_bus(net, vn_kv=float(values[iloop, 9]), name='Bus'+str(iloop), index=None, geodata=None, type='b', zone=None,
                              in_service=True, max_vm_pu=float(values[iloop, 11]), min_vm_pu=float(values[iloop, 12]), coords=None))
# Enter Load-data
        LoadList.append(pp.create_load(net, bus=BusList[iloop], p_mw= float(values[iloop, 2])*sbase, q_mvar=float(values[iloop, 3])*sbase,
                       const_z_percent=0, const_i_percent=0, sn_mva=10, name="Load" + str(iloop), scaling=1.0,
                       index=None, in_service=True, type='wye'))
        iloop += 1

# Define an external grid
    pp.create_ext_grid(net, bus=BusList[0], vm_pu=1.0, va_degree=0.0)
#    pp.create_gen(net, bus=BusList[0], p_mw = 0.0, vm_pu=1.0, slack=True)  # Could create a generator
    df2 = pd.read_excel(xls, 'branch', usecols="A:L")
    values = df2.values
    # Read Line data  --------------------------------------------
    iloop = 0
    # print(' ')
    while iloop < len(values):

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
pp.to_excel(net, filename="test5.xlsx", include_empty_tables=False, include_results=True)
