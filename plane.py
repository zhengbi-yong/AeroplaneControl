class Plane(object):
    def __init__(self, m, l, Sw, Ix, Iy, Iz, Ixy, Iyz, Izx, cA):
        # 飞机物理参数
        self._m = m  # 质量
        self._l = l  # 翼展
        self._Sw = Sw  # 翼面积
        self._Ix = Ix  # 滚转转动惯量
        self._Iy = Iy  # 俯仰转动惯量
        self._Iz = Iz  # 偏航转动惯量
        self._Ixy = Ixy  # 惯性积
        self._Iyz = Iyz
        self._Izx = Izx
        self._cA = cA  # 平均气动弦长
        self._u = 0.0  # 翼展变形率
        self._v = 0.0  # 翼角变形率
        self._zt = 0.0
        self._TV = 0.0
        self._gamma0 = 0.0
        # 飞机气动参数
        self._CL0 = 0.2
        self._CLa = 3.0
        self._CLM = 0.2
        self._CLdeltae = 0.0
        self._CD0 = 0.03
        self._CDa = 0.2
        self._CDM = 0.05
        self._Cm0 = 0.0
        self._Cma = -0.5
        self._Cmq = -7.0
        self._Cmadot = -3.0
        self._CmM = -0.06
        self._Cmdeltae = 0.0
        # 大导数
        self._XV = 0.0
        self._Xa = 0.0
        self._Xtheta = 0.0
        self._ZV = 0.0
        self._Za = 0.0
        self._MV = 0.0
        self._Ma = 0.0
        self._Madot = 0.0
        self._Mq = 0.0
        # 中间量
        self._DV = 0.0
        self._Da = 0.0
        self._LV = 0.0
        self._La = 0.0
        self._Ldeltae = 0.0
        self._MaV = 0.0
        self._Maa = 0.0
        self._Maadot = 0.0
        self._Maq = 0.0
        self._Madeltae = 0.0
        # 操纵力矩
        self._Xdeltat = 0.0
        self._Zdeltae = 0.0
        self._Mdeltae = 0.0
        self._Mdeltat = 0.0
        self._Tdeltat = 0.0

    def printm(self):
        print(f"质量：  {self._m:>10} kg")

    def printl(self):
        print(f"翼展：  {self._l:>10} m")

    def printSw(self):
        print(f"翼面积：{self._Sw:>10} m^2")

    def printIx(self):
        print(
            f"Ix ：{self._Ix:>10} kg*m^2 Iy：{self._Iy:>10} kg*m^2 Iz：{self._Iz:>10} kg*m^2")

    def printIxx(self):
        print(
            f"Ixy：{self._Ixy:>13} m^4 Iyz：{self._Iyz:>12} m^4 Izx：{self._Izx:>12} m^4")

    def printcA(self):
        print(f"平均气动弦长：{self._cA} m")

    def printuv(self):
        print(f"翼展变形率：{self._u}\n翼角变形率：{self._v}")

    def printzt(self):
        print(f"zt：{self._zt:>10}")

    def printTV(self):
        print(f"TV：{self._TV:>10}")

    def setCL0(self, CL0):
        self._CL0 = CL0

    def printCL0(self):
        print(f"CL0：{self._CL0:>10}")

    def setCLa(self, CLa):
        self._CLa = CLa

    def printCLa(self):
        print(f"CLa：{self._CLa:>10}")

    def setCLM(self, CLM):
        self._CLM = CLM

    def printCLM(self):
        print(f"CLM：{self._CLM:>10}")

    def setCD0(self, CD0):
        self._CD0 = CD0

    def printCD0(self):
        print(f"CD0：{self._CD0:>10}")

    def setCDa(self, CDa):
        self._CDa = CDa

    def printCDa(self):
        print(f"CDa：{self._CDa:>10}")

    def setCDM(self, CDM):
        self._CDM = CDM

    def printCDM(self):
        print(f"CDM：{self._CDM:>10}")

    def setCma(self, Cma):
        self._Cma = Cma

    def printCma(self):
        print(f"Cma：{self._Cma:>10}")

    def setCmq(self, Cmq):
        self._Cmq = Cmq

    def printCmq(self):
        print(f"Cmq：{self._Cmq:>10}")

    def setCmadot(self, Cmadot):
        self._Cmadot = Cmadot

    def printCmadot(self):
        print(f"Cmadot：{self._Cmadot:>10}")

    def setCmM(self, CmM):
        self._CmM = CmM

    def printCmM(self):
        print(f"CmM：{self._CmM:>10}")

    def printInfo(self):
        print('*'*80)
        print("飞机参数")
        print('-'*80)
        self.printm()
        self.printl()
        self.printSw()
        self.printIx()
        self.printIxx()
        self.printcA()
        self.printuv()
        self.printzt()
        self.printTV()
        self.printCL0()
        self.printCLa()
        self.printCLM()
        self.printCD0()
        self.printCDa()
        self.printCDM()
        self.printCma()
        self.printCmq()
        self.printCmadot()
        self.printCmM()
        print('-'*80)
        print('*'*80)

    def computeDV(self, V0, rho, M0):
        self._DV = (1/V0)*((1/2)*rho*(V0**2)) * \
            self._Sw*(2*self._CD0+M0*self._CDM)

    def computeDa(self, V0, rho):
        self._Da = self._CDa*(1/2)*rho*V0**2*self._Sw
        # print('*'*20)
        # print("computeDa")
        # print('*'*20)
        # print(f"CDa:{self._CDa}")
        # print(f"rho:{rho}")
        # print(f"V0:{V0}")
        # print(f"Sw:{self._Sw}")
        # print(f"Da:{self._Da}")

    def computeLV(self, V0, rho, M0):
        self._LV = (1/V0)*((1/2)*rho*(V0**2)) * \
            self._Sw*(2*self._CL0+M0*self._CLM)
        # print('*'*20)
        # print("computeLV")
        # print('*'*20)
        # print(f"V0:{V0}")
        # print(f"rho:{rho}")
        # print(f"Sw:{self._Sw}")
        # print(f"CL0:{self._CL0}")
        # print(f"M0:{M0}")
        # print(f"CLM:{self._CLM}")
        # print(f"LV:{self._LV}")

    def computeLa(self, V0, rho):
        self._La = self._CLa*(1/2)*rho*V0**2*self._Sw

    def computeLdeltae(self, V0, rho):
        self._Ldeltae = self._CLdeltae*(1/2)*rho*V0**2*self._Sw

    def computeMaV(self, V0, rho, M0):
        self._MaV = ((2*self._Cm0+M0*self._CmM)/V0) * \
            ((1/2)*rho*V0**2)*self._cA*self._Sw

    def computeMaa(self, V0, rho):
        self._Maa = ((1/2)*rho*V0**2)*self._cA*self._Sw*self._Cma

    def computeMaadot(self, V0, rho):
        self._Maadot = ((1/2)*rho*V0**2)*(self._cA **
                                          2/(2*V0))*self._Sw*self._Cmadot

    def computeMaq(self, V0, rho):
        self._Maq = ((1/2)*rho*V0**2)*(self._cA **
                                       2/(2*V0))*self._Sw*self._Cmq
        # print('*'*20)
        # print("computeMaq")
        # print('*'*20)
        # print(f"rho:{rho}")
        # print(f"V0:{V0}")
        # print(f"cA:{self._cA}")
        # print(f"Sw:{self._Sw}")
        # print(f"Cmq:{self._Cmq}")

    def computeMadeltae(self, V0, rho):
        self._Madeltae = ((1/2)*rho*V0**2)*(self._cA**2 /
                                            (2*V0))*self._Sw*self._Cmdeltae

    def computeXV(self):
        self._XV = (self._DV-self._TV)/self._m

    def computeXa(self, V0, g):
        self._Xa = self._Da/(self._m*V0)-g/V0
        # print('*'*20)
        # print("computeXa")
        # print('*'*20)
        # print(f"Da:{self._Da}")
        # print(f"m:{self._m}")
        # print(f"V0:{V0}")
        # print(f"g:{g}")
        # print(f"Xa:{self._Xa}")

    def computeXtheta(self, V0, g):
        self._Xtheta = g/V0

    def computeXdeltat(self, V0):
        self._Xdeltat = -(self._Tdeltat/(self._m*V0))

    def computeZV(self, V0):
        self._ZV = self._LV/(self._m*V0)
        # print('*'*20)
        # print("computeZV")
        # print('*'*20)
        # print(f"LV:{self._LV}")
        # print(f"m:{self._m}")
        # print(f"V0:{V0}")

    def computeZa(self, V0):
        self._Za = self._La/(self._m*V0)

    def computeZdeltae(self, V0):
        self._Zdeltae = self._Ldeltae/(self._m*V0)

    def computeMV(self, V0):
        self._MV = -(V0*(self._MaV+self._TV*self._zt)/self._Iy)

    def computeMadot(self):
        self._Madot = -self._Maadot/self._Iy

    def computeMa(self):
        self._Ma = -self._Maa/self._Iy

    def computeMq(self):
        self._Mq = -self._Maq/self._Iy
        # print('*'*20)
        # print("computeMq")
        # print('*'*20)
        # print(f"Maq:{self._Maq}")
        # print(f"Iy:{self._Iy}")

    def computeMdeltae(self):
        self._Mdeltae = -self._Madeltae/self._Iy

    def computeMdeltat(self):
        self._Mdeltat = -(self._Tdeltat*self._zt)/self._Iy

    def computeBigDerivative(self, V0, rho, M0, g):
        self.computeDV(V0, rho, M0)
        self.computeDa(V0, rho)
        self.computeLV(V0, rho, M0)
        self.computeLa(V0, rho)
        self.computeLdeltae(V0, rho)
        self.computeMaV(V0, rho, M0)
        self.computeMaa(V0, rho)
        self.computeMaadot(V0, rho)
        self.computeMaq(V0, rho)
        self.computeMadeltae(V0, rho)
        self.computeXV()
        self.computeXa(V0, g)
        self.computeXtheta(V0, g)
        self.computeXdeltat(V0)
        self.computeZV(V0)
        self.computeZa(V0)
        self.computeZdeltae(V0)
        self.computeMV(V0)
        self.computeMadot()
        self.computeMa()
        self.computeMq()
        self.computeMdeltae()
        self.computeMdeltat()

    def printBigDerivative(self):
        print('*'*80)
        print("气动导数")
        print('-'*80)
        print(f"XV：{self._XV:<25} Xa：{self._Xa:<25} Xtheta：{self._Xtheta:<25}")
        print(f"ZV：{self._ZV:<25} Za：{self._Za:<25} MV：{self._MV:<25}")
        print(f"Ma：{self._Ma:<25} Madot：{self._Madot:<22} Mq：{self._Mq:<25}")
        print('-'*80)
        print('*'*80)


g = 9.8
M = 343.0
M0 = 0.9
V0 = M0*M
h = 11000.0

if __name__ == "__main__":
    print('*'*80)
    print("Plane类测试")
    print('*'*80)
    # 设置参数
    m = 900.0
    l = 10.0
    Sw = 27.95
    cA = 3.097
    Ix = 0.0
    Iy = 7447.0
    Iz = 0.0
    V0 = 266.0
    rho = 0.0371
    a = 285.0
    a0 = 3.62
    CL0 = 0.246
    CLa = 3.9
    CLM = 0.23
    CD0 = 0.0306
    CDa = 0.284
    CDM = 0.055
    Cma = -0.562
    Cmq = -7.06
    Cmadot = -2.792
    CmM = -0.0654
    TV = 0.0
    gamma0 = 0.0
    zt = 0.0
    Ixy = 0.0
    Iyz = 0.0
    Izx = 0.0
    # 初始化飞机
    testPlane = Plane(m, l, Sw, Ix, Iy, Iz, Ixy, Iyz, Izx, cA)
    testPlane.setCD0(CD0=CD0)
    testPlane.setCDa(CDa=CDa)
    testPlane.setCDM(CDM=CDM)
    testPlane.setCL0(CL0=CL0)
    testPlane.setCLa(CLa=CLa)
    testPlane.setCLM(CLM=CLM)
    testPlane.setCma(Cma=Cma)
    testPlane.setCmadot(Cmadot=Cmadot)
    testPlane.setCmq(Cmq=Cmq)
    testPlane.setCmM(CmM=CmM)
    # 打印飞机信息
    testPlane.printInfo()
    # 求大导数
    testPlane.computeBigDerivative(V0, rho, M0, g)
    # 打印大导数
    testPlane.printBigDerivative()
