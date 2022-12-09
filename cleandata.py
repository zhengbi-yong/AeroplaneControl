import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt


LiftFilePath = "./Lift.csv"
DragFilePath = "./Drag.csv"
PTFilePath = "./PT.csv"
if __name__ == "__main__":
    LiftData = genfromtxt(LiftFilePath, delimiter=',')
    DragData = genfromtxt(DragFilePath, delimiter=',')
    PTData = genfromtxt(PTFilePath, delimiter=',')
    ulist = LiftData[1:, 0]
    vlist = LiftData[1:, 1]
    CL0list = LiftData[1:, 2]
    CLalist = LiftData[1:, 3]
    CD0list = DragData[1:, 2]
    CDalist = DragData[1:, 3]
    CDa2list = DragData[1:, 4]
    Cm0list = PTData[1:, 2]
    Cmalist = PTData[1:, 3]
    # print(f"ulist:{ulist}")
    # print(f"vlist:{vlist}")
    # print(f"CL0list:{CL0list}")
    # print(f"CLalist:{CLalist}")
    # print(f"CD0list:{CD0list}")
    # print(f"CDalist:{CDalist}")
    # print(f"CDa2list:{CDa2list}")
    # print(f"Cm0list:{Cm0list}")
    # print(f"Cmalist:{Cmalist}")
