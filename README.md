# MCR-SVC-Based Voltage Regulation Study
This repository contains the code and data used to evaluate and optimize reactive power compensation in a wind farm equipped with a Magnetically Controlled Reactor-Static Var Compensator (MCR-SVC) system. The analysis is based on steady-state power flow simulations using operational SCADA data collected between 2023 and 2024.
The repository supports the findings presented in the associated research study and includes tools for compensation assessment, voltage profile analysis, and sensitivity analysis.
## Requirements
The code was developed and tested using Python 3.12 and requires the following packages:
- pandapower
- pandas
- numpy
- matplotlib
- scipy
- xlrd
- xlwt
- xlutils
### Installation
Install the required packages using:
```bash
pip install pandapower pandas numpy matplotlib scipy xlrd xlwt xlutils
