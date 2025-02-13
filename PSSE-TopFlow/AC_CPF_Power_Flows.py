# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix, csr_matrix


# ---------------------------------------------------------------------
# AC Load Flow - Regular Newton Rhapson-------------------------------------------------
# ----------------------------------------------------------------------

class ACLoadFlow:
    """
    Common base class on AC  Load Flow
    Input:
        BusList      - List off all Bus ojects
        GenList      - List of all Generator objects
        LoadList     - List of all Load objects
        TransList    - List of all transmission lines objets
    Returns: None
        
    """
    
    def __init__(self,BusList,GenList,LoadList,TransList):
        self.BusList = BusList
        self.GenList = GenList
        self.LoadList = LoadList
        self.TransList = TransList
        self.voang = np.zeros(len(self.BusList))
        self.vomag = np.ones(len(self.BusList))
        self.plambda = np.ones(len(self.BusList))
        self.qlambda = np.zeros(len(self.BusList))
        



    def BuildJacobi(self):   # Build the Jacobian Matrix (complete) ---
        """
        Desc:    Build the complete Jacobian Matrix.
                 The elements are stored in a list containing lists as (x, y, elem)
        Input:   None
        Returns: elem - List of sub lists
        """         
        
        def uij(gij,bij,tetai,tetaj):
            return (gij*np.sin(tetai-tetaj)-bij*np.cos(tetai-tetaj))

        def tij(gij,bij,tetai,tetaj):
            return (gij*np.cos(tetai-tetaj)+bij*np.sin(tetai-tetaj))

        bsize = len(self.BusList)
        elem = []
        iloop = 0
        while iloop < len(self.BusList):            # Add bus shunt elements (positiv in)
            if self.BusList[iloop].BL > 0.0:
                elem.append([iloop + bsize, iloop + bsize, -2*self.BusList[iloop].BL*self.vomag[iloop]])  # dQ1/dV1 - Bus shunt
            iloop += 1       

        iloop = 0
        while iloop < len(self.TransList):
            if self.TransList[iloop].Ibstat:
                ifr = self.TransList[iloop].FromBus
                itr = self.TransList[iloop].ToBus
                teta1 = self.voang[ifr-1]
                teta2 = self.voang[itr-1]
                v1 = self.vomag[ifr-1]
                v2 = self.vomag[itr-1]
                b = self.TransList[iloop].B
                g = self.TransList[iloop].G
                bsh = self.TransList[iloop].Bc2/2.0
        # dP/dTeta
                elem.append([ifr, ifr, v1*v2*uij(g,b,teta1,teta2)])  # dP1/dT1
                elem.append([ifr, itr, -v1*v2*uij(g,b,teta1,teta2)]) # dP1/dT1
                elem.append([itr, itr, v1*v2*uij(g,b,teta2,teta1)])  # dP2/dT2
                elem.append([itr, ifr, -v1*v2*uij(g,b,teta2,teta1)]) # dP2/dT1
        # dP/dV               
                elem.append([ifr, ifr + bsize, 2*g*v1 - v2*tij(g,b,teta1,teta2)])   # dP1/dV1
                elem.append([ifr, itr + bsize, -v1*tij(g,b,teta1,teta2)])           # dP1/dV2
                elem.append([itr, ifr + bsize, -v2*tij(g,b,teta2,teta1)])           # dP2/dV1
                elem.append([itr, itr + bsize, 2*v2*g - v1*tij(g,b,teta2,teta1)])   # dP2/dV2
        # dQ/dTeta
                elem.append([ifr + bsize, ifr, -v1*v2*tij(g,b,teta1,teta2)])    # dQ1/dT1
                elem.append([ifr + bsize, itr, v1*v2*tij(g,b,teta1,teta2)])     # dQ1/dT2
                elem.append([itr + bsize, itr, -v1*v2*tij(g,b,teta2,teta1)])    # dQ2/dT2
                elem.append([itr + bsize, ifr, v1*v2*tij(g,b,teta2,teta1)])     # dQ2/dT1
        # dQ/dV               
                elem.append([ifr + bsize, ifr + bsize, -2*(b+bsh)*v1 - v2*uij(g,b,teta1,teta2)])  # dQ1/dV1
                elem.append([ifr + bsize, itr + bsize, -v1*uij(g,b,teta1,teta2)])           # dQ1/dV2
                elem.append([itr + bsize, ifr + bsize, -v2*uij(g,b,teta2,teta1)])           # dQ2/dV1
                elem.append([itr + bsize, itr + bsize, -2*v2*(b+bsh) - v1*uij(g,b,teta2,teta1)])  # dQ2/dV2
                
            iloop = iloop + 1

        return elem





    def dispGen(self,tpres=False):
        """
        Desc:    Display generation at all gen buses
        Input:   tpres= False (Display in tableformat if True
        Returns: None
        """                    
        mainlist = []
        rowno = []
        
        iloop = 0
        print(' ')
        while iloop < len(self.GenList):
            oref = self.GenList[iloop]
            if tpres == False:
                print(' Bus no :', '{:4.0f}'.format(oref.BusNo),
                      ' Pgen :',  '{:7.4f}'.format(oref.Pgen),
                      ' AccPgen :',  '{:7.4f}'.format(oref.AccPgen),
                      ' Qgen :',  '{:7.4f}'.format(oref.Qgen))
# Prepare for graphics presentation                
            sublist = ['{:4.0f}'.format(oref.BusNo), '{:7.4f}'.format(oref.Pgen),
                       '{:7.4f}'.format(oref.AccPgen),
                   '{:7.4f}'.format(oref.Qgen)]
            mainlist.append(sublist)
            rowno.append('G '+ str(iloop+1))
            iloop += 1
# Present table
        if tpres:
            title = 'Individual Generators'
            colind = ['Bus no', 'Pgen', 'AccPgen','Qgen']
            colw = [0.12,0.22,0.22,0.22,0.22]
            self.tableplot(mainlist, title, colind, rowno, columncol=[], rowcol=[], colw=colw)

    def loadScale(self, pscale=1.0, qscale=1.0):
        """
        Scale the loads of the system. Idependent scaling of P and Q
        :param pscale: Default 1.0
        :param qscale: Default 1.0
        :return: none
        """
        iloop = 0
        while iloop < len(self.LoadList):
            oref = self.LoadList[iloop]
            oref.Pload *= pscale
            oref.Qload *= qscale
            iloop += 1

    def dispLoad(self, tpres=False):
        """
        Desc:    Display load at all load buses
        Input:   tpres= False (Display in tableformat if True
        Returns: None
        """

        mainlist = []
        rowno = []
        
        iloop = 0
        print(' ')
        while iloop < len(self.LoadList):
            oref = self.LoadList[iloop]
            if tpres == False:
                print(' Bus no :', '{:4.0f}'.format(oref.BusNo),
                      ' Pload :',  '{:7.4f}'.format(oref.Pload),
                      ' Qload :',  '{:7.4f}'.format(oref.Qload),
                      ' AccPload :','{:7.4f}'.format(oref.AccPload),
                      ' AccQload :',  '{:7.4f}'.format(oref.AccQload))
                
# Prepare for graphics presentation                
            sublist = ['{:4.0f}'.format(oref.BusNo), '{:7.4f}'.format(oref.Pload),
                       '{:7.4f}'.format(oref.Qload),
                   '{:7.4f}'.format(oref.AccPload), '{:7.4f}'.format(oref.AccQload)]
            mainlist.append(sublist)
            rowno.append('L '+ str(iloop+1))
            
            iloop += 1
# Present table
        if tpres:
            title = 'Individual Bus Loads'
            colind = ['Bus no', 'Pload', 'Qload','AccPload','AccQload']
            colw = [0.12,0.22,0.22,0.22,0.22]
            self.tableplot(mainlist,title,colind,rowno,columncol=[],rowcol=[],colw=colw)

#
#
    def dispNonZeroLoad(self, fromBus=0, toBus=0, tpres=False):
        """
        Desc:    Display buses with nin-zero load  for a bus range
        Input:   tpres= False (Display in tableformat if True
        Returns: None
        """
        mainlist = []
        rowno = []

        if toBus == 0:
            toBus = len(self.BusList)
        if tpres:
            toBus = np.minimum(fromBus + 13, toBus)
        iloop = fromBus
        while iloop < toBus:
            oref = self.LoadList[iloop]
            if oref.Pload != 0.0 or oref.Qload != 0.0:
                if tpres == False:
                    print(' Bus no :', '{:4.0f}'.format(oref.BusNo),
                          ' Pload :', '{:7.4f}'.format(oref.Pload),
                          ' Qload :', '{:7.4f}'.format(oref.Qload))

            # Prepare for graphics presentation
                sublist = ['{:4.0f}'.format(oref.BusNo), '{:7.4f}'.format(oref.Pload),
                       '{:7.4f}'.format(oref.Qload)]
                mainlist.append(sublist)
                rowno.append('L ' + str(iloop+1))
            iloop += 1
# Present table
        if tpres:
            title = 'Individual Bus Loads'
            colind = ['Bus no', 'Pload', 'Qload']
            colw = [0.12,0.22,0.22]
            self.tableplot(mainlist,title,colind,rowno,columncol=[],rowcol=[], colw=colw)

#
#
    def dispLoadAcc(self, fromBus=0, toBus=0, tpres=False):
        """
        Desc:    Display bus loads and accumulated loads  for a bus range
        Input:   tpres= False (Display in tableformat if True
        Returns: None
        """
        mainlist = []
        rowno = []

        if toBus == 0:
            toBus = len(self.BusList)
        if tpres:
            toBus = np.minimum(fromBus + 13, toBus)
        iloop = fromBus
        while iloop < toBus:
            oref = self.LoadList[iloop]
            if tpres == False:
                print(' Bus no :', '{:4.0f}'.format(oref.BusNo),
                      ' Pload :', '{:7.4f}'.format(oref.Pload),
                      ' Qload :', '{:7.4f}'.format(oref.Qload),
                      ' AccPload :', '{:7.4f}'.format(oref.AccPload),
                      ' AccQload :', '{:7.4f}'.format(oref.AccQload))

                # Prepare for graphics presentation
            sublist = ['{:4.0f}'.format(oref.BusNo), '{:7.4f}'.format(oref.Pload),
                       '{:7.4f}'.format(oref.Qload),
                       '{:7.4f}'.format(oref.AccPload), '{:7.4f}'.format(oref.AccQload)]
            mainlist.append(sublist)
            rowno.append('L' + str(iloop+1))
            iloop += 1
# Present table
        if tpres:
            title = 'Individual Bus Loads'
            colind = ['Bus no', 'Pload', 'Qload','AccPload','AccQload']
            colw = [0.12,0.22,0.22,0.22,0.22]
            self.tableplot(mainlist,title,colind,rowno,columncol=[],rowcol=[], colw=colw)


    def dispLoss(self):
        """
        Desc:    Display total load, generation and loss
        Input:   None
        Returns: None
        """

        pltot = 0.0
        qltot = 0.0
        pgtot = 0.0
        qgtot = 0.0

# Find total load        
        iloop = 0
        print(' ')
        while iloop < len(self.LoadList):
            oref = self.LoadList[iloop]
            pltot += oref.Pload
            qltot += oref.Qload
            iloop += 1

            
# Find total generation
        iloop = 0
        while iloop < len(self.GenList):
            oref = self.GenList[iloop]
            pgtot += oref.Pgen
            qgtot += oref.Qgen
            iloop += 1
        
# Display results
        print(' Pload :',  '{:7.4f}'.format(pltot),
              ' Qload :',  '{:7.4f}'.format(qltot),
              ' Pgen  :','{:7.4f}'.format(pgtot),
              ' Qgen  :',  '{:7.4f}'.format(qgtot),
              ' Ploss :','{:7.5f}'.format(pgtot-pltot),
              ' Qloss :','{:7.5f}'.format(qgtot-qltot))
                
            

    def dispBusData(self, fromBus=0, toBus = 0, tpres=True):
        """
        Desc: Displays all the data associated with the buses.
              Extension for larger cases to define them into blocks
        """

        rowno = []
        mainlist = []

        if toBus == 0:
            toBus = len(self.BusList)
        if tpres:
            toBus = np.minimum(fromBus + 13, toBus)

        iloop = fromBus
        while iloop < toBus:
            oref = self.BusList[iloop]
            sublist = ['{:4.0f}'.format(oref.BusNo),
                       '{:7.4f}'.format(oref.BaskV),'{:4.0f}'.format(oref.BusCode),
                       '{:4.0f}'.format(oref.Ipload),'{:4.0f}'.format(oref.Ipgen)]
            mainlist.append(sublist)
            rowno.append(' Bus ' + str(iloop + 1))
            iloop += 1

        if tpres:
            title = 'Bus Data'
            colind = ['Bus no', 'Base kV',
                      'BusCode', 'Load No', 'Gen No']
            colw = [0.12,0.22,0.22,0.22,0.22]

            self.tableplot(mainlist,title,colind,rowno,columncol=[],rowcol=[], colw=colw)
        
    #
    # Display the voltages.
    #
    def dispVolt(self, fromBus=0, toBus = 0, tpres=False):
        """
        Desc:    Display voltages at all buses
        Input:   tpres= False (Display in tableformat if True)
                 fromBus and toBus defines the block, If tpres=True, it will display 13 lines from fromBus
        Returns: None
        """
        mainlist = []
        rowno = []
        if toBus == 0:
            toBus = len(self.BusList)
        if tpres:
            toBus = np.minimum(fromBus + 13, toBus)

        iloop = fromBus
        print(' ')
        while iloop < toBus:
            oref = self.BusList[iloop]
            if tpres == False:
                print(' Bus no :', '{:4.0f}'.format(oref.BusNo),
                      ' Vmag :', '{:7.5f}'.format(self.vomag[iloop]),
                      ' Theta :', '{:7.5f}'.format(self.voang[iloop] * 180 / np.pi))
            # Prepare for graphics presentation
            sublist = ['{:4.0f}'.format(oref.BusNo),
                       '{:7.5f}'.format(self.vomag[iloop]),
                       '{:7.5f}'.format(self.voang[iloop] * 180 / np.pi)]

            mainlist.append(sublist)
            rowno.append('Bus ' + str(iloop+1))
            iloop += 1
        # Present table
        if tpres:
            title = 'Bus Voltages'
            colind = ['Bus no', 'Vmag', 'Theta']
            colw = [0.12,0.22,0.22]
            self.tableplot(mainlist, title, colind, rowno, columncol=[], rowcol=[],colw=colw)

   #
    # Display the voltages.
    #
    def dispLowVolt(self, fromBus=0, toBus = 0, tpres=False, vmax = 1.1):
        """
        Desc:    Display voltages at all buses below or equal to the limit vmax
        Input:   tpres= False (Display in tableformat if True)
                 fromBus and toBus defines the block, If tpres=True, it will display 13 lines from fromBus
                 vmax = Upper voltage limit (default 1.1 pu)
        Returns: None
        """
        mainlist = []
        rowno = []
        if toBus == 0:
            toBus = len(self.BusList)
#        if tpres:
#            toBus = np.minimum(fromBus + 13, toBus)

        iloop = fromBus
        print(' ')
        while iloop < toBus:
            oref = self.BusList[iloop]
            if self.vomag[iloop] <= vmax:
                if tpres == False:
                    print(' Bus no :', '{:4.0f}'.format(oref.BusNo),
                          ' Vmag :', '{:7.5f}'.format(self.vomag[iloop]),
                          ' Theta :', '{:7.5f}'.format(self.voang[iloop] * 180 / np.pi))
                # Prepare for graphics presentation
                sublist = ['{:4.0f}'.format(oref.BusNo),
                           '{:7.5f}'.format(self.vomag[iloop]),
                           '{:7.5f}'.format(self.voang[iloop] * 180 / np.pi)]

                mainlist.append(sublist)
                rowno.append('Bus ' + str(iloop + 1))
            iloop += 1
        # Present table
        if tpres:
            title = 'Bus Voltages'
            colind = ['Bus no', 'Vmag', 'Theta']
            colw = [0.12,0.22,0.22,0.22,0.22]
            self.tableplot(mainlist, title, colind, rowno, columncol=[], rowcol=[],colw=colw)

#
    # Display the voltages.
    #
    def dispVoltRange(self, fromBus=0, toBus = 0, tpres=False, vmin=0.9, vmax = 1.1):
        """
        Desc:    Display voltages at all buses below or equal to the limit vmax
        Input:   tpres= False (Display in tableformat if True)
                 fromBus and toBus defines the block, If tpres=True, it will display 13 lines from fromBus
                 vmax = Upper voltage limit (default 1.1 pu)
                 vmin = Lower voltage limit (defualt 0.9 pu)
        Returns: None
        """
        mainlist = []
        rowno = []
        if toBus == 0:
            toBus = len(self.BusList)
#        if tpres:
#            toBus = np.minimum(fromBus + 13, toBus)

        iloop = fromBus
        print(' ')
        while iloop < toBus:
            oref = self.BusList[iloop]
            if self.vomag[iloop] <= vmax and self.vomag[iloop] >= vmin:
                if tpres == False:
                    print(' Bus no :', '{:4.0f}'.format(oref.BusNo),
                          ' Vmag :', '{:7.5f}'.format(self.vomag[iloop]),
                          ' Theta :', '{:7.5f}'.format(self.voang[iloop] * 180 / np.pi))
                # Prepare for graphics presentation
                sublist = ['{:4.0f}'.format(oref.BusNo),
                           '{:7.5f}'.format(self.vomag[iloop]),
                           '{:7.5f}'.format(self.voang[iloop] * 180 / np.pi)]

                mainlist.append(sublist)
                rowno.append('Bus ' + str(iloop + 1))
            iloop += 1
        # Present table
        if tpres:
            title = 'Bus Voltages'
            colind = ['Bus no', 'Vmag', 'Theta']
            colw = [0.12,0.22,0.22,0.22,0.22]
            self.tableplot(mainlist, title, colind, rowno, columncol=[], rowcol=[], colw=colw)



    def dispVolt2(self,tpres=False):
        """
        Desc:    Display voltages at all buses
        Input:   tpres= False (Display in tableformat if True)
        Returns: None
        """
        mainlist = []
        rowno = []
        
        iloop = 0
        print(' ')
        while iloop < len(self.BusList):
            oref = self.BusList[iloop]
            if tpres == False:
                print(' Bus no :', '{:4.0f}'.format(oref.BusNo),
                      ' Vmag :',  '{:7.5f}'.format(self.vomag[iloop]),
                      ' Theta :',  '{:7.5f}'.format(self.voang[iloop]*180/np.pi))
# Prepare for graphics presentation                
            sublist = ['{:4.0f}'.format(oref.BusNo),
                       '{:7.5f}'.format(self.vomag[iloop]),
                       '{:7.5f}'.format(self.voang[iloop]*180/np.pi)]
            
            mainlist.append(sublist)
            rowno.append('Bus '+ str(iloop+1))
            iloop += 1
# Present table
        if tpres:
            title = 'Bus Voltages'
            colind = ['Bus no', 'Vmag', 'Theta']
            colw = [0.12,0.22,0.22,0.22,0.22]
            self.tableplot(mainlist,title,colind,rowno,columncol=[],rowcol=[], colw=colw)

    def dispDistFact(self, jacobi, busjac, tpres=False):  # Define function for calculating flow on a transmission line
        """
        Desc:    Display transmission line flows
        Input:   tpres= False (Display in tableformat if True)
        Returns: None
        """

        def uij(gij, bij, tetai, tetaj):
            return (gij * np.sin(tetai - tetaj) - bij * np.cos(tetai - tetaj))

        def tij(gij, bij, tetai, tetaj):
            return (gij * np.cos(tetai - tetaj) + bij * np.sin(tetai - tetaj))

        mainlist = []
        rowno = []
        iloop = 0

        print(' ')
        #while iloop < len(self.TransList):
        while iloop < 5:
            ifr = self.TransList[iloop].FromBus
            itr = self.TransList[iloop].ToBus
            if self.TransList[iloop].Ibstat:
                teta1 = self.voang[ifr - 1]
                teta2 = self.voang[itr - 1]
                v1 = self.vomag[ifr - 1]
                v2 = self.vomag[itr - 1]
                b = self.TransList[iloop].B
                g = self.TransList[iloop].G
                bsh = self.TransList[iloop].Bc2 / 2.0
                taprat = self.TransList[iloop].Ratio
                if taprat > 0.0:
                    taprat = 1.0 / taprat
                    bsh1 = -taprat * (1.0 - taprat) * b
                    bsh2 = (1.0 - taprat) * b
                    b = taprat * b
                    dPfrdTi = v1 * v2 * uij(g, b, teta1, teta2)
                    dPfrdTj = - v1 * v2 * uij(g, b, teta1, teta2)
                    dPfrdVi = 2*g * v1 - v2 * tij(g, b, teta1, teta2)
                    dPfrdVj = -v1 * tij(g, b, teta1, teta2)
                else:
                    # dPfrdTi = v1 * v2 * uij(g, b, teta1, teta2)
                    # dPfrdTj = - v1 * v2 * uij(g, b, teta1, teta2)
                    # dPfrdVi = 2*g * v1 - v2 * tij(g, b, teta1, teta2)
                    # dPfrdVj = -v1 * tij(g, b, teta1, teta2)

                    dPfrdTj = -v1 * v2 * uij(g, b, teta2, teta1)
                    dPfrdTi =  v1 * v2 * uij(g, b, teta2, teta1)
                    dPfrdVi = -2*g * v1 + v2 * tij(g, b, teta2, teta1)  # Switched i and j to get impact on injection on to bus
                    dPfrdVj = v1 * tij(g, b, teta2, teta1)  # Switching v1 and v2 will give the flow in the end of the line (by bus j)

                  #  dPfrdVi = 0
                 #   dPfrdVj = 0


                    iloop2 = 0  # Identify the right rows in active and reactive power blocks
                    rhs = np.zeros(len(jacobi))
                    while iloop2 < len(self.BusList) - 1:
                        if ifr == self.BusList[busjac[iloop2]].BusNo:
                            rhs[iloop2] = dPfrdTi
                        if itr == self.BusList[busjac[iloop2]].BusNo:
                            rhs[iloop2] = dPfrdTj
                        iloop2 += 1

                    iloop2 = 0
                    inext = len(self.BusList) - 1
                    while iloop2 < len(self.BusList) - 1:
                        if ifr == self.BusList[busjac[iloop2 + inext]].BusNo:
                            rhs[iloop2 + inext] = dPfrdVi
                        if itr == self.BusList[busjac[iloop2 + inext]].BusNo:
                            rhs[iloop2 + inext] = dPfrdVj
                        iloop2 += 1

                   # print(rhs)
                    x = np.linalg.solve(jacobi, rhs)

                    print('Line:', ifr, itr)
                    print(x)



            # if tpres == False:
            #     print(' FromBus :', '{:4.0f}'.format(ifr), ' ToBus :', '{:4.0f}'.format(itr),
            #           ' Pfrom :', '{:7.4f}'.format(Pfrom), ' Qfrom : ', '{:7.4f}'.format(Qfrom),
            #           ' Pto :', '{:7.4f}'.format(Pto), ' Qto :', '{:7.4f}'.format(Qto))
            #
            # sublist = [ifr, itr, '{:7.4f}'.format(Pfrom), '{:7.4f}'.format(Qfrom),
            #            '{:7.4f}'.format(Pto), '{:7.4f}'.format(Qfrom)]
            # mainlist.append(sublist)
            # rowno.append('L ' + str(iloop + 1))

            #        print(pinj, qinj)
            #        print(iloop)
            iloop = iloop + 1

        if tpres:
            title = 'Transmission line flow'
            colind = ['FromBus', 'ToBus', 'Pfrom', ' Qfrom', 'Pto', 'Qto']
            self.tableplot(mainlist, title, colind, rowno, columncol=[], rowcol=[])

    def dispDistFactList(self, jacobi, busjac, distlinelist = [], tpres=False):  # Define function for calculating flow on a transmission line
        """
        Desc:    Display transmission line flows
        Input:   tpres= False (Display in tableformat if True)
        Returns: None
        """

        def uij(gij, bij, tetai, tetaj):
            return (gij * np.sin(tetai - tetaj) - bij * np.cos(tetai - tetaj))

        def tij(gij, bij, tetai, tetaj):
            return (gij * np.cos(tetai - tetaj) + bij * np.sin(tetai - tetaj))

        objectlist = []
        for i in distlinelist:
            objectlist.append(self.TransList[i-1])

        mainlist = []
        rowno = []
        iloop = 0

        print(' ')
        while iloop < len(objectlist):
            ifr = objectlist[iloop].FromBus
            itr = objectlist[iloop].ToBus
            if objectlist[iloop].Ibstat:
                teta1 = self.voang[ifr - 1]
                teta2 = self.voang[itr - 1]
                v1 = self.vomag[ifr - 1]
                v2 = self.vomag[itr - 1]
                b = objectlist[iloop].B
                g = objectlist[iloop].G
                bsh = objectlist[iloop].Bc2 / 2.0
                taprat = objectlist[iloop].Ratio
                if taprat > 0.0:
                    taprat = 1.0 / taprat
                    bsh1 = -taprat * (1.0 - taprat) * b
                    bsh2 = (1.0 - taprat) * b
                    b = taprat * b
                    dPfrdTi = v1 * v2 * uij(g, b, teta1, teta2)
                    dPfrdTj = - v1 * v2 * uij(g, b, teta1, teta2)
                    dPfrdVi = 2*g * v1 - v2 * tij(g, b, teta1, teta2)
                    dPfrdVj = -v1 * tij(g, b, teta1, teta2)
                else:
                    # dPfrdTi = v1 * v2 * uij(g, b, teta1, teta2)
                    # dPfrdTj = - v1 * v2 * uij(g, b, teta1, teta2)
                    # dPfrdVi = 2*g * v1 - v2 * tij(g, b, teta1, teta2)
                    # dPfrdVj = -v1 * tij(g, b, teta1, teta2)

                    dPfrdTj = -v1 * v2 * uij(g, b, teta2, teta1)
                    dPfrdTi =  v1 * v2 * uij(g, b, teta2, teta1)
                    dPfrdVi = -2*g * v1 + v2 * tij(g, b, teta2, teta1)  # Switched i and j to get impact on injection on to bus
                    dPfrdVj = v1 * tij(g, b, teta2, teta1)  # Switching v1 and v2 will give the flow in the end of the line (by bus j)

                  #  dPfrdVi = 0
                 #   dPfrdVj = 0


                    iloop2 = 0  # Identify the right rows in active and reactive power blocks
                    rhs = np.zeros(len(jacobi))
                    while iloop2 < len(self.BusList) - 1:
                        if ifr == self.BusList[busjac[iloop2]].BusNo:
                            rhs[iloop2] = dPfrdTi
                        if itr == self.BusList[busjac[iloop2]].BusNo:
                            rhs[iloop2] = dPfrdTj
                        iloop2 += 1

                    iloop2 = 0
                    inext = len(self.BusList) - 1
                    while iloop2 < len(self.BusList) - 1:
                        if ifr == self.BusList[busjac[iloop2 + inext]].BusNo:
                            rhs[iloop2 + inext] = dPfrdVi
                        if itr == self.BusList[busjac[iloop2 + inext]].BusNo:
                            rhs[iloop2 + inext] = dPfrdVj
                        iloop2 += 1

                   # print(rhs)
                    x = np.linalg.solve(jacobi, rhs)

                    print('Line:', ifr, itr)
                    print(x)



            # if tpres == False:
            #     print(' FromBus :', '{:4.0f}'.format(ifr), ' ToBus :', '{:4.0f}'.format(itr),
            #           ' Pfrom :', '{:7.4f}'.format(Pfrom), ' Qfrom : ', '{:7.4f}'.format(Qfrom),
            #           ' Pto :', '{:7.4f}'.format(Pto), ' Qto :', '{:7.4f}'.format(Qto))
            #
            # sublist = [ifr, itr, '{:7.4f}'.format(Pfrom), '{:7.4f}'.format(Qfrom),
            #            '{:7.4f}'.format(Pto), '{:7.4f}'.format(Qfrom)]
            # mainlist.append(sublist)
            # rowno.append('L ' + str(iloop + 1))

            #        print(pinj, qinj)
            #        print(iloop)
            iloop += 1

        if tpres:
            title = 'Transmission line flow'
            colind = ['FromBus', 'ToBus', 'Pfrom', ' Qfrom', 'Pto', 'Qto']
            self.tableplot(mainlist, title, colind, rowno, columncol=[], rowcol=[])


    def dispFlow(self,tpres=False): # Define function for calculating flow on a transmission line
        """
        Desc:    Display transmission line flows
        Input:   tpres= False (Display in tableformat if True)
        Returns: None
        """            
        def uij(gij,bij,tetai,tetaj):
            return (gij*np.sin(tetai-tetaj)-bij*np.cos(tetai-tetaj))

        def tij(gij,bij,tetai,tetaj):
            return (gij*np.cos(tetai-tetaj)+bij*np.sin(tetai-tetaj))
        mainlist = []
        rowno = []
        iloop = 0
        print(' ')
        while iloop < len(self.TransList):
            ifr = self.TransList[iloop].FromBus
            itr = self.TransList[iloop].ToBus            
            if self.TransList[iloop].Ibstat:
                teta1 = self.voang[ifr-1]
                teta2 = self.voang[itr-1]
                v1 = self.vomag[ifr-1]
                v2 = self.vomag[itr-1]
                b = self.TransList[iloop].B
                g = self.TransList[iloop].G
                bsh = self.TransList[iloop].Bc2/2.0
                taprat = self.TransList[iloop].Ratio
                if  taprat > 0.0:
                    taprat = 1.0/taprat
                    bsh1 = -taprat*(1.0-taprat)*b
                    bsh2 = (1.0-taprat)*b
                    b = taprat*b
                    Pfrom = g*v1*v1 - v1*v2*tij(g,b,teta1,teta2)
                    Pto = g*v2*v2 - v1*v2*tij(g,b,teta2,teta1)                    
                    Qfrom = -(b+bsh1)*v1*v1 - v1*v2*uij(g,b,teta1,teta2)
                    Qto = -(b+bsh2)*v2*v2 - v1*v2*uij(g,b,teta2,teta1)
                else:
                    Pfrom = g*v1*v1 - v1*v2*tij(g,b,teta1,teta2)
                    Pto = g*v2*v2 - v1*v2*tij(g,b,teta2,teta1)
                    Qfrom = -(b+bsh)*v1*v1 - v1*v2*uij(g,b,teta1,teta2)
                    Qto = -(b+bsh)*v2*v2 - v1*v2*uij(g,b,teta2,teta1)                

            else:
                Pfrom = 0.0
                Pto = 0.0
                Qfrom = 0.0
                Qto = 0.0
            if tpres == False:
                print(' FromBus :', '{:4.0f}'.format(ifr), ' ToBus :','{:4.0f}'.format(itr),
                      ' Pfrom :', '{:7.4f}'.format(Pfrom), ' Qfrom : ', '{:7.4f}'.format(Qfrom),
                      ' Pto :', '{:7.4f}'.format(Pto), ' Qto :', '{:7.4f}'.format(Qto))
            
            sublist = [ifr, itr, '{:7.4f}'.format(Pfrom), '{:7.4f}'.format(Qfrom),
                       '{:7.4f}'.format(Pto), '{:7.4f}'.format(Qfrom)]
            mainlist.append(sublist)
            rowno.append('L '+ str(iloop+1))

            
    #        print(pinj, qinj)
    #        print(iloop)
            iloop = iloop + 1

        if tpres:
            title = 'Transmission line flow'
            colind = ['FromBus', 'ToBus', 'Pfrom',' Qfrom','Pto','Qto']
            self.tableplot(mainlist,title,colind,rowno,columncol=[],rowcol=[])



    def dispFlowViol(self,rating='RateA',rmult= 1.0,tpres=False): # Define function for calculating flow on a transmission line
        """
        Desc:    Display transmission lines with flow violation
        Input:   rating = 'RateA' (Alternatives RateB, RateC)
                 tpres= False (Display in tableformat if True)
        Returns: None
        """
        self.rating = rating
        self.rmult = rmult
        
        def uij(gij,bij,tetai,tetaj):
            return (gij*np.sin(tetai-tetaj)-bij*np.cos(tetai-tetaj))

        def tij(gij,bij,tetai,tetaj):
            return (gij*np.cos(tetai-tetaj)+bij*np.sin(tetai-tetaj))

        mainlist = []
        rowno = []
        
        iloop = 0
        print(' ')
        while iloop < len(self.TransList):
            ifr = self.TransList[iloop].FromBus
            itr = self.TransList[iloop].ToBus
            if self.TransList[iloop].Ibstat:
                teta1 = self.voang[ifr-1]
                teta2 = self.voang[itr-1]
                v1 = self.vomag[ifr-1]
                v2 = self.vomag[itr-1]
                b = self.TransList[iloop].B
                g = self.TransList[iloop].G
                bsh = self.TransList[iloop].Bc2/2.0
                taprat = self.TransList[iloop].Ratio
                if  taprat > 0.0:
                    taprat = 1.0/taprat
                    bsh1 = -taprat*(1.0-taprat)*b
                    bsh2 = (1.0-taprat)*b
                    b = taprat*b
                    Pfrom = g*v1*v1 - v1*v2*tij(g,b,teta1,teta2)
                    Pto = g*v2*v2 - v1*v2*tij(g,b,teta2,teta1)                    
                    Qfrom = -(b+bsh1)*v1*v1 - v1*v2*uij(g,b,teta1,teta2)
                    Qto = -(b+bsh2)*v2*v2 - v1*v2*uij(g,b,teta2,teta1)
                else:
                    Pfrom = g*v1*v1 - v1*v2*tij(g,b,teta1,teta2)
                    Pto = g*v2*v2 - v1*v2*tij(g,b,teta2,teta1)
                    Qfrom = -(b+bsh)*v1*v1 - v1*v2*uij(g,b,teta1,teta2)
                    Qto = -(b+bsh)*v2*v2 - v1*v2*uij(g,b,teta2,teta1)                      
            else:
                Pfrom = 0.0
                Pto = 0.0
                Qfrom = 0.0
                Qto = 0.0
            s1 = np.sqrt(Pfrom**2 + Qfrom**2)
            s2 = np.sqrt(Pto**2 + Qto**2)
            if self.rating == 'RateA':
                rate = self.TransList[iloop].RateA*self.rmult
            elif self.rating == 'RateB':
                rate = self.TransList[iloop].RateB*self.rmult
            elif self.rating == 'RateC':
                rate = self.TransList[iloop].RateC*self.rmult
# Print results to screen                
            if s1 > rate or s2 > rate:
                if tpres == False:
                    print(' FromBus :', '{:4.0f}'.format(ifr), ' ToBus :','{:4.0f}'.format(itr),
                          ' Pfrom :', '{:7.4f}'.format(Pfrom), ' Qfrom : ', '{:7.4f}'.format(Qfrom),
                          ' Pto :', '{:7.4f}'.format(Pto), ' Qto :', '{:7.4f}'.format(Qto),
                          self.rating, '{:7.4f}'.format(rate))
# Prepare for graphics presentation                
                sublist = [ifr, itr, '{:7.4f}'.format(Pfrom), '{:7.4f}'.format(Qfrom),
                       '{:7.4f}'.format(Pto), '{:7.4f}'.format(Qfrom), '{:7.4f}'.format(rate)]
                mainlist.append(sublist)
                rowno.append('Line '+ str(iloop+1))

            iloop = iloop + 1
            
# Present the results in a table format
        if tpres:
            title = 'Overloaded transmission lines'
            colind = ['FromBus', 'ToBus', 'Pfrom','Qfrom','Pto','Qto', self.rating]
            self.tableplot(mainlist,title,colind,rowno,columncol=[],rowcol=[])

    def dispMarginalCost(self, fromBus=0, toBus = 0, tpres=False):
        """
        Desc:    Display marginal costs at all buses at all buses
        Input:   tpres= False (Display in tableformat if True)
        Returns: None
        """
        mainlist = []
        rowno = []
        if toBus == 0:
            toBus = len(self.BusList)
        if tpres:
            toBus = np.minimum(fromBus + 13, toBus)

        iloop = fromBus
        print(' ')
        while iloop < toBus:
            oref = self.BusList[iloop]
            if tpres == False:
                print(' Bus no :', '{:4.0f}'.format(oref.BusNo),
                      ' Pmc :', '{:7.5f}'.format(self.plambda[iloop]),
                      ' Qmc :', '{:7.5f}'.format(self.qlambda[iloop]))
            # Prepare for graphics presentation
            sublist = ['{:4.0f}'.format(oref.BusNo),
                       '{:7.5f}'.format(self.plambda[iloop]),
                       '{:7.5f}'.format(self.qlambda[iloop])]

            mainlist.append(sublist)
            rowno.append('Bus ' + str(iloop + 1))
            iloop += 1
        # Present table
        if tpres:
            title = 'Bus Marginal Costs'
            colind = [' Bus no ', ' Pmc ', ' Qmc ']
            colw = [0.12,0.22,0.22,0.22,0.22]
            self.tableplot(mainlist, title, colind, rowno, columncol=[], rowcol=[], colw=colw)



    def ExportValues(self):                 # Make the voltages available on top level
        """
        Desc:    Make the locally stored voltages available on top level
        Input:   None
        Returns: voang, vomag - Voltage angles and magnitudes
                   
        """         
        return self.voang, self.vomag



    def ExtrNZ(self,amat):
        """
        Desc:    extracts non-zero elements from a matrix
        Input:   amat - square matrix assumed so far
        Returns: listnz - List of lists of non-zero elements (i, j elem)
        """
        id1 = len(amat)
        ir1 = 0
        listnz = []
        while ir1 < id1:
            ic1 = 0
            while ic1 < id1:
                if amat[ir1,ic1] != 0.0:
                    listnz.append([ir1, ic1,'{:10.5f}'.format(amat[ir1,ic1])])
                ic1 += 1
            ir1 += 1
        return listnz


    

    def InitVolt(self):    # Initialize the voltage angles and magnitude - Flat start so far
        """
        Desc:    Initializes the bus voltages according to gen Vset
                 Flat start assumed so far
        Input:   None
        Returns: None
                   
        """
        
        iloop = 0
        while iloop < len(self.GenList):
            self.vomag[self.GenList[iloop].BusNo -1] = self.GenList[iloop].Vset
            iloop += 1



    def MakeMatrix(self, elem): # Make the complete matrix (full representation) --
        """
        Desc:    Makes a full matrix (temporary solution during testing)
        Input:   elem - List of lists (Jacobi elements)
        Returns: amat - full matrix 
                  
        """        
        iloop = 0
        amat = np.zeros([2*len(self.BusList),2*len(self.BusList)])
        while iloop < len(elem):
            sub = elem[iloop]
            amat[sub[0]-1][sub[1]-1] += sub[2]
            iloop += 1
        return amat

 


    def Netinj(self):   # Define function for calculating net injection at nodes ---
        """
        Desc:    Calculates net injection of active and reactive power
        Input:   None
        Returns: pinj, qinj - net injection at all nodes
                   
        """            
        def uij(gij,bij,tetai,tetaj):
            return (gij*np.sin(tetai-tetaj)-bij*np.cos(tetai-tetaj))

        def tij(gij,bij,tetai,tetaj):
            return (gij*np.cos(tetai-tetaj)+bij*np.sin(tetai-tetaj))

        pinj = np.zeros(len(self.BusList))
        qinj = np.zeros(len(self.BusList))
        iloop = 0
        while iloop < len(self.TransList):
            if self.TransList[iloop].Ibstat:
                ifr = self.TransList[iloop].FromBus
                itr = self.TransList[iloop].ToBus
                teta1 = self.voang[ifr-1]
                teta2 = self.voang[itr-1]
                v1 = self.vomag[ifr-1]
                v2 = self.vomag[itr-1]
                b = self.TransList[iloop].B
                g = self.TransList[iloop].G
                bsh = self.TransList[iloop].Bc2/2.0
                taprat = self.TransList[iloop].Ratio
                if  taprat > 0.0:
                    taprat = 1.0/taprat
                    bsh1 = -taprat*(1.0-taprat)*b
                    bsh2 = (1.0-taprat)*b
                    b = taprat*b
                    pinj[ifr-1] += g*v1*v1 - v1*v2*tij(g,b,teta1,teta2)
                    pinj[itr-1] += g*v2*v2 - v1*v2*tij(g,b,teta2,teta1)                    
                    qinj[ifr-1] += -(b+bsh1)*v1*v1 - v1*v2*uij(g,b,teta1,teta2)
                    qinj[itr-1] += -(b+bsh2)*v2*v2 - v1*v2*uij(g,b,teta2,teta1)
                else:
                    pinj[ifr-1] += g*v1*v1 - v1*v2*tij(g,b,teta1,teta2)
                    pinj[itr-1] += g*v2*v2 - v1*v2*tij(g,b,teta2,teta1)
                    qinj[ifr-1] += -(b+bsh)*v1*v1 - v1*v2*uij(g,b,teta1,teta2)
                    qinj[itr-1] += -(b+bsh)*v2*v2 - v1*v2*uij(g,b,teta2,teta1)
    #        print(pinj, qinj)
    #        print(iloop)
            iloop = iloop + 1
        return pinj, qinj



    def PQmism(self, pinj, qinj): # Calculate the mismatch of the AC Load Flow
        """
        Desc:    Calculate the mismatch of the AC Load Flow
        Input:   pinj, qinj - net injection
        Returns: rhs - mismatch vector
        """               
        
        rhs = np.zeros(2*len(self.BusList))
        iloop = 0
        while iloop < len(self.BusList):
            rhs[iloop] = -pinj[iloop] 
            rhs[iloop + len(self.BusList)] = -qinj[iloop] + self.BusList[iloop].BL*self.vomag[iloop]**2
            iloop += 1
    # Modify for specified values for injections
        iloop = 0
        while iloop < len(self.LoadList):
            oref = self.LoadList[iloop]
            rhs[oref.BusNo -1] -= oref.Pload + oref.AccPload    # Update for loads
            rhs[oref.BusNo -1 + len(self.BusList)] -= oref.Qload + oref.AccQload
            iloop += 1
        iloop = 0
        while iloop < len(self.GenList):
            oref = self.GenList[iloop]
            rhs[oref.BusNo -1] += oref.Pgen  + oref.AccPgen  # Update for generators
            if self.BusList[oref.BusNo-1].BusCode == -2:
               rhs[oref.BusNo -1 + len(self.BusList)] += oref.Qgen 
            iloop += 1
        rhs = self.TrimRHS(rhs)
        return rhs



    def print_f105(self,alist, icount=0,varname=[]):
        """
        Desc:    Print a list in one row.
        Input:   alist, icount=0 (optional), varname=[] (optional)
        Returns: None
        """          
        self.icount = icount
        self.name = varname
        self.printlist = alist
        if self.icount > 0:
            print('Count:',self.icount, end=" ")
        if self.name != []:
            print(' ',self.name, end=" ")
        a1 = [print('{:10.5f}'.format(i), end=" ") for i in alist]
        print(end="\n")



    def print_mat(self, amat):
        """
        Desc: Prints a square matrix to screen in a given format

        """
        self.amat = amat
        id1 = len(amat)
        ic1 = 0
        ir1 = 0
        while ir1 < id1:
            alist = amat[ir1,:]
            a1 = [print('{:10.5f}'.format(i), end=" ") for i in alist[:]]
            print('\n')
            ir1 += 1

    def testSym(self,amat):
        id1 = len(amat)
        ic1 = 0
        ir1 = 0
        isym = True
        while ir1 < id1:
            diff1 = amat[ic1,ir1] - amat[ir1,ic1]
            if np.abs(diff1) > 0.000001:
                isym = False
                print('Element : ', ir1, ic1, '  different from', ic1, ir1)
            ir1 += 1
        print( '\n',' Matrix is symmetric ','\n')
        

    def QlimAct(self):
        """
        Desc:    Activates Q limiations on violated generator buses,
                 Sets limits and updates BodCode
        Input:   None
        Returns: None
        """               
        iloop = 0
        while iloop < len(self.GenList):
            igen = self.GenList[iloop]
            if igen.Qgen > igen.Qmax:
                igen.Qgen = igen.Qmax
                self.BusList[igen.BusNo-1].BusCode = -2
            elif igen.Qgen < igen.Qmin:
                igen.Qgen = igen.Qmin
                self.BusList[igen.BusNo-1].BusCode = -2
            iloop += 1




    def QlimDeact(self):            # Check this logic
        """
        Desc:    Deactivates Q limiations on generator buses,
                 Sets Vset and updates BodCode
        Input:   None
        Returns: None
        """                    
        iloop = 0
        while iloop < len(self.GenList):
            igen = self.GenList[iloop]
            if self.vomag[igen.BusNo-1] > igen.Vset:
                self.BusList[igen.BusNo-1].BusCode = 2
            elif self.vomag[igen.BusNo-1] < igen.Vset and  np.abs(igen.Qgen - igen.Qmin) < 0.005:       
                self.BusList[igen.BusNo-1].BusCode = 2
            iloop += 1



    def SolveAC(self, flatStart = True, enfqlim=False, itmax=4,sparse=False):
        """
        Desc:    Solve the AC load flow case
        Input:   None
        Returns: None
        """
        self.enfqlim = enfqlim
        self.itmax = itmax
        self.sparse = sparse
        self.flatStart = flatStart

        print('\n','******** Running AC Power Flow ********','\n')
        if self.flatStart:
            self.InitVolt()
#        print('Initial vomag :', self.vomag)
        iloop = 0
        while iloop < self.itmax:
            elem = self.BuildJacobi()
            amat = self.MakeMatrix(elem)
#            print(amat)
            jacobi, busjac = self.TrimMatrix(amat)
#            print('jacobi :', len(jacobi))
            
            pinj, qinj = self.Netinj()
            rhs1 = self.PQmism(pinj, qinj)

            if self.sparse:
                jacobi1 = csc_matrix(jacobi, dtype=float)
                x = spsolve(jacobi1,rhs1)
            else:
                x = np.linalg.solve(jacobi,rhs1)
            
            iloop2 = 0
            while iloop2 < len(x):
                if iloop2 < len(self.BusList) -1:    # Reduce due to slack bus
                    if self.BusList[busjac[iloop2]].BusCode < 3:
                        self.voang[busjac[iloop2]] += x[iloop2]
                elif self.BusList[busjac[iloop2]].BusCode < 2:
                     self.vomag[busjac[iloop2]] += x[iloop2]
                iloop2 += 1
            iloop += 1
        return jacobi, busjac

    def PQmismOPF(self, pinj, qinj):  # Calculate the mismatch of the AC Load Flow
        """
        Desc:    Calculate the mismatch of the Newton OPF Load Flow (all buses)
        Input:   pinj, qinj - net injection
        Returns: rhs - mismatch vector
        """

        rhs = np.zeros(2 * len(self.BusList))
        iloop = 0
        while iloop < len(self.BusList):
            rhs[iloop] = -pinj[iloop]
            rhs[iloop + len(self.BusList)] = -qinj[iloop] + self.BusList[iloop].BL * self.vomag[iloop] ** 2
            iloop += 1
        # Modify for specified values for injections
        iloop = 0
        while iloop < len(self.LoadList):
            oref = self.LoadList[iloop]
            rhs[oref.BusNo - 1] -= oref.Pload + oref.AccPload  # Update for loads
            rhs[oref.BusNo - 1 + len(self.BusList)] -= oref.Qload + oref.AccQload
            iloop += 1
        iloop = 0
        while iloop < len(self.GenList):
            oref = self.GenList[iloop]
            rhs[oref.BusNo - 1] += oref.Pgen + oref.AccPgen  # Update for generators
            rhs[oref.BusNo - 1 + len(self.BusList)] += oref.Qgen    # All buses kept in the Newton OPF
            if self.BusList[oref.BusNo - 1].BusCode == -2:
                rhs[oref.BusNo - 1 + len(self.BusList)] += oref.Qgen
            iloop += 1
#        rhs = self.TrimRHS(rhs)
        return rhs







    def tableplot(self,table_data,title,columns,rows,columncol=[],rowcol=[],colw=None):
        """
        Desc:   Make a table of the provided data. There must be a row and a column
                data correpsonding to the table
        Input:  table_data  - np.array
                title - string
                columns - string vector
                rows    - string vector
                columncol - colors of each column label (default [])
                rowcol - colors of each row lable
        """
        
        fig = plt.figure(dpi=150)
        
        ax = fig.add_subplot(1,1,1)
        
        tdim = np.shape(table_data)
        iloop = 0
        if rowcol==[]:
            while iloop < tdim[0]:
                rowcol.append('cyan')
                iloop += 1
        iloop = 0        
        if columncol == []:
            while iloop < tdim[1]:
                columncol.append('cyan')
                iloop += 1            

        table = ax.table(cellText=table_data, rowLabels=rows, colColours=columncol,
                         rowColours=rowcol,colLabels=columns, colWidths=colw, loc='center',cellLoc='center')
        
        table.set_animated = True
#        table.scale(1,1.5)
        table.scale(1,1)
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        ax.set_title(title, fontsize=14)
        ax.axis('off')
        plt.show()
      
 
    def TrimRHS(self, rhs):   # Remove rows corresponsing to slack and voltage controlled buses
        """
        Desc:    Remove row corresponding to slack bus and gen under voltage control
        Input:   rhs - mismatch vector
        Returns: rhs - mismatch vector without slack bus
        """             
        iloop = 0
        rembus = []
        while iloop < len(self.BusList):
            bcode = self.BusList[iloop].BusCode
            if bcode == 3:
               rembus.append(iloop)
               rembus.append(iloop + len(self.BusList))
            elif bcode == 2:
               rembus.append(iloop + len(self.BusList))
            iloop +=1

        rembus.sort(reverse=True)   # Sort buses in decending order

        iloop = 0
        while iloop < len(rembus):  # Remove Rows and Columns (the last first)
            rhs = np.delete(rhs,rembus[iloop],0)
            iloop += 1
        return rhs




    def TrimMatrix(self, amat):  # Remove columns and rows corr. to slack and voltage contr. buses
        """
        Desc:    Remove columns and rows corr. to slack and voltage contr. buses
        Input:   amat - Full Jacobian matrix
        Returns: amat - reduced matrix
                 busjac - Vector of external bus references
                  
        """  
        iloop = 0
        rembus = []
        busjac = np.zeros(2*len(self.BusList), dtype=int)
        while iloop < len(self.BusList):
            bcode = self.BusList[iloop].BusCode
            busjac[iloop] = iloop               # Points at the internal numbering in self.BusList
            busjac[iloop + len(self.BusList)] = iloop
            if bcode == 3:
               rembus.append(iloop)
               rembus.append(iloop + len(self.BusList))
            elif bcode == 2:
               rembus.append(iloop + len(self.BusList))
            iloop +=1

        rembus.sort(reverse=True)   # Sort buses in decending order
   
        iloop = 0
        while iloop < len(rembus):       # Remove Rows and Columns (the last first)
            amat = np.delete(amat,rembus[iloop],0)
            amat = np.delete(amat,rembus[iloop],1)
            busjac = np.delete(busjac, rembus[iloop], 0)
            iloop += 1
        return amat, busjac

    

    # def UpdateGen(self):
    #     """
    #     Desc:    Updates active generation at slack bus and reactive generation at gen buses
    #     Input:   None
    #     Returns: None
    #     """
    #     pinj, qinj = self.Netinj()
    #     iloop = 0
    #  #   print('Pinj :',pinj)
    # #   print('Qinj :',qinj)
    #
    #     while iloop < len(self.GenList):
    #         ibus = self.GenList[iloop].BusNo -1
    #         iload = self.LoadList[self.BusList[ibus].Ipload] # Update reactive power gen
    #         self.GenList[iloop].Qgen = (qinj[ibus] +iload.Qload +  iload.AccQload
    #                                     - self.BusList[ibus].BL**self.vomag[ibus]**2)
    #         self.GenList[iloop].Pgen = pinj[ibus] + iload.Pload + iload.AccPload    # Pgen will be updated in case of mismatch
    #         if self.BusList[ibus].BusCode == 3:             # Update slack bus generation
    #             iload = self.LoadList[self.BusList[ibus].Ipload]
    #             self.GenList[iloop].Pgen = pinj[ibus] +iload.Pload  +  iload.AccPload       # Was commented 30102021
    #         iloop += 1


    def UpdateGen(self):
        """
        Desc:    Updates active generation at slack bus and reactive generation at gen buses
        Input:   None
        Returns: None

        May need some update in the DPF-case if load is accumulated at non-load buses 17.07.2021
        """
        pinj, qinj = self.Netinj()
        iloop = 0
        #   print('Pinj :',pinj)
        #   print('Qinj :',qinj)

        while iloop < len(self.GenList):
            ibus = self.GenList[iloop].BusNo - 1
            if self.BusList[ibus].Ipload >= 0:
                iload = self.LoadList[self.BusList[ibus].Ipload]  # Update reactive power gen
                self.GenList[iloop].Qgen = (qinj[ibus] + iload.Qload + iload.AccQload
                                            - self.BusList[ibus].BL ** self.vomag[ibus] ** 2)
                self.GenList[iloop].Pgen = pinj[
                                               ibus] + iload.Pload + iload.AccPload  # Pgen will be updated in case of mismatch
            else:
                self.GenList[iloop].Qgen = (qinj[ibus]
                                            - self.BusList[ibus].BL ** self.vomag[ibus] ** 2)
                self.GenList[iloop].Pgen = pinj[ibus]  # Pgen will be updated in case of mismatch

            if self.BusList[ibus].BusCode == 3:  # Update slack bus generation
                iload = self.LoadList[self.BusList[ibus].Ipload]
                if self.BusList[ibus].Ipload >= 0:
                    self.GenList[iloop].Pgen = pinj[ibus] + iload.Pload + iload.AccPload  # Was commented 30102021
                else:
                    self.GenList[iloop].Pgen = pinj[ibus]
            iloop += 1

# End of file
