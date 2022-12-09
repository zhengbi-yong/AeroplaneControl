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
        self._zt = 0.0  # 发动机推力向量到机体 x 轴的距离
        self._TV = 0.0  # 发动机推力 T 对飞行速度 V 的偏导数，需要发动机特性曲线
        self._gamma0 = 0.0  # 飞机初始速度方向与水平面的夹角
        self._Tdeltat = 0.0  # 发动机推力系数
        # 飞机气动参数
        self._CL0 = 0.0
        self._CLa = 0.0
        self._CLM = 0.0
        self._CLdeltae = 0.0  # 升降舵偏转对升力的影响
        self._CD0 = 0.0
        self._CDa = 0.0
        self._CDM = 0.0
        self._Cm0 = 0.0
        self._Cma = 0.0
        self._Cmq = 0.0
        self._Cmadot = 0.0
        self._CmM = 0.0
        self._Cmdeltae = 0.0
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
        # 操纵力矩
        self._Xdeltat = 0.0
        self._Zdeltae = 0.0
        self._Mdeltae = 0.0
        self._Mdeltat = 0.0
        # 稳定判据
        self._a0 = 1.0
        self._a1 = 1.0
        self._a2 = 1.0
        self._a3 = 1.0
        self._Stability = True

    def setu(self, u):
        self._u = u

    def setv(self, v):
        self._v = v

    def setzt(self, zt):
        self._zt

    def setTV(self, TV):
        self._TV = TV

    def setgamma0(self, gamma0):
        self._gamma0 = gamma0

    def setTdeltat(self, Tdeltat):
        self._Tdeltat = Tdeltat

    def setCL0(self, CL0):
        self._CL0 = CL0

    def setCLa(self, CLa):
        self._CLa = CLa

    def setCLM(self, CLM):
        self._CLM = CLM

    def setCLdeltae(self, CLdeltae):
        self._CLdeltae = CLdeltae

    def setCD0(self, CD0):
        self._CD0 = CD0

    def setCDa(self, CDa):
        self._CDa = CDa

    def setCDM(self, CDM):
        self._CDM = CDM

    def setCm0(self, Cm0):
        self._Cm0 = Cm0

    def setCma(self, Cma):
        self._Cma = Cma

    def setCmq(self, Cmq):
        self._Cmq = Cmq

    def setCmadot(self, Cmadot):
        self._Cmadot = Cmadot

    def setCmM(self, CmM):
        self._CmM = CmM

    def setCmdeltae(self, Cmdeltae):
        self._Cmdeltae = Cmdeltae

    def printInfo(self):
        print('*'*80)
        print("飞机参数")
        print('-'*80)
        print(f"质量:  {self._m:>10} kg")
        print(f"翼展:  {self._l:>10} m")
        print(f"翼面积:{self._Sw:>10} m^2")
        print(
            f"转动惯量:Ix :{self._Ix:>10} kg*m^2 Iy:{self._Iy:>10} kg*m^2 Iz:{self._Iz:>10} kg*m^2")
        print(
            f"惯性积:Ixy:{self._Ixy:>13} m^4 Iyz:{self._Iyz:>12} m^4 Izx:{self._Izx:>12} m^4")
        print(f"平均气动弦长cA:             {self._cA} m")
        print(f"翼展变形率lambda:{self._u}\n翼角变形率rho:{self._v}")
        print(f"零迎角升力系数CL0:          {self._CL0}")
        print(f"迎角升力系数CLa:            {self._CLa}")
        print(f"马赫数升力系数CLM:          {self._CLM}")
        print(f"零迎角阻力系数CD0:          {self._CD0}")
        print(f"迎角阻力系数CDa:            {self._CDa}")
        print(f"马赫数阻力系数CDM:          {self._CDM}")
        print(f"迎角力矩系数Cma:            {self._Cma}")
        print(f"地面轴角速度力矩系数Cmq:     {self._Cmq}")
        print(f"机体轴角速度力矩系数Cmadot:  {self._Cmadot}")
        print(f"马赫数力矩系数CmM:           {self._CmM}")
        print(f"发动机推力系数Tdeltat:       {self._Tdeltat}")
        print(f"推力对速度偏导数:            {self._TV}")
        print(f"gamma0:                    {self._gamma0}")
        print(f"推力矢量到机体x轴的距离zt:   {self._zt} m")
        print(f"稳定性:                    {self._Stability}")
        print('-'*80)
        print('*'*80)

    def computeDV(self, V0, rho, M0):
        self._DV = (1/V0)*((1/2)*rho*(V0**2)) * \
            self._Sw*(2*self._CD0+M0*self._CDM)
        # print('*'*20)
        # print("computeDV")
        # print('*'*20)
        # print(f"V0:{V0}")
        # print(f"rho:{rho}")
        # print(f"Sw:{self._Sw}")
        # print(f"CD0:{self._CD0}")
        # print(f"M0:{M0}")
        # print(f"CDM:{self._CDM}")
        # print(f"DV:{self._DV}")

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
        print('*'*20)
        print("computeLV")
        print('*'*20)
        print(f"V0:{V0}")
        print(f"rho:{rho}")
        print(f"Sw:{self._Sw}")
        print(f"CL0:{self._CL0}")
        print(f"M0:{M0}")
        print(f"CLM:{self._CLM}")
        print(f"LV:{self._LV}")

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
        # print('*'*20)
        # print("computeXV")
        # print('*'*20)
        # print(f"DV:{self._DV}")
        # print(f"TV:{self._TV}")
        # print(f"m:{self._m}")
        # print(f"XV:{self._XV}")

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

    def computeZV(self, V0):
        self._ZV = self._LV/(self._m*V0)
        print('*'*20)
        print("computeZV")
        print('*'*20)
        print(f"LV:{self._LV}")
        print(f"m:{self._m}")
        print(f"V0:{V0}")

    def computeZa(self, V0):
        self._Za = self._La/(self._m*V0)

    def computeMV(self, V0):
        self._MV = -(V0*(self._MaV+self._TV*self._zt)/self._Iy)

    def computeMa(self):
        self._Ma = -self._Maa/self._Iy

    def computeMadot(self):
        self._Madot = -self._Maadot/self._Iy

    def computeMq(self):
        self._Mq = -self._Maq/self._Iy
        # print('*'*20)
        # print("computeMq")
        # print('*'*20)
        # print(f"Maq:{self._Maq}")
        # print(f"Iy:{self._Iy}")

    def computeXdeltat(self, V0):
        self._Xdeltat = -(self._Tdeltat/(self._m*V0))

    def computeZdeltae(self, V0):
        self._Zdeltae = self._Ldeltae/(self._m*V0)

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

    # 一定要在计算完大导数之后再使用该函数判断稳定性
    def judgeStability(self):
        self._a3 = self._Mq+self._Za+self._Madot+self._XV
        self._a2 = self._XV*(self._Mq+self._Za+self._Madot) + \
            self._Za*self._Mq+self._Ma-self._Xa*self._ZV
        self._a1 = self._XV*(self._Za*self._Mq+self._Ma)+self._Xtheta * \
            (self._ZV*self._Madot-self._MV) - \
            self._Xa*(self._ZV*self._Mq+self._MV)
        self._a0 = self._Xtheta*(self._ZV*self._Ma-self._MV*self._Za)
        if ((self._a0 > 0) and (self._a1 > 0) and (self._a2 > 0) and (self._a3 > 0)):
            self._Stability = True
        else:
            self._Stability = False

    def printBigDerivative(self):
        print('*'*80)
        print("气动导数")
        print('-'*80)
        print(f"XV:{self._XV:<25} Xa:{self._Xa:<25} Xtheta:{self._Xtheta:<25}")
        print(f"ZV:{self._ZV:<25} Za:{self._Za:<25} MV:{self._MV:<25}")
        print(f"Ma:{self._Ma:<25} Madot:{self._Madot:<22} Mq:{self._Mq:<25}")
        print('-'*80)
        print('*'*80)

    def getXV(self):
        return self._XV

    def getXa(self):
        return self._Xa

    def getXtheta(self):
        return self._Xtheta

    def getZV(self):
        return self._ZV

    def getZa(self):
        return self._Za

    def getMV(self):
        return self._MV

    def getMa(self):
        return self._Ma

    def getMadot(self):
        return self._Madot

    def getMq(self):
        return self._Mq

    def getStability(self):
        self.judgeStability()
        return self._Stability


# # 飞机物理参数
# m  # 质量
# l  # 翼展
# Sw  # 翼面积
# Ix  # 滚转转动惯量
# Iy  # 俯仰转动惯量
# Iz  # 偏航转动惯量
# Ixy  # 惯性积
# Iyz
# Izx
# cA  # 平均气动弦长
# u = 0.0  # 翼展变形率
# v = 0.0  # 翼角变形率
# zt = 0.0  # 发动机推力向量到机体 x 轴的距离
# TV = 0.0  # 发动机推力 T 对飞行速度 V 的偏导数，需要发动机特性曲线
# gamma0 = 0.0  # 飞机初始速度方向与水平面的夹角
# Tdeltat = 0.0  # 发动机推力系数
# # 飞机气动参数
# CL0 = 0.2
# CLa = 3.0
# CLM = 0.2
# CLdeltae = 0.0  # 升降舵偏转对升力的影响
# CD0 = 0.03
# CDa = 0.2
# CDM = 0.05
# Cm0 = 0.0
# Cma = -0.5
# Cmq = -7.0
# Cmadot = -3.0
# CmM = -0.06
# Cmdeltae = 0.0
# # 中间量
# DV = 0.0
# Da = 0.0
# LV = 0.0
# La = 0.0
# Ldeltae = 0.0
# MaV = 0.0
# Maa = 0.0
# Maadot = 0.0
# Maq = 0.0
# Madeltae = 0.0
# # 大导数
# XV = 0.0
# Xa = 0.0
# Xtheta = 0.0
# ZV = 0.0
# Za = 0.0
# MV = 0.0
# Ma = 0.0
# Madot = 0.0
# Mq = 0.0
# # 操纵力矩
# Xdeltat = 0.0
# Zdeltae = 0.0
# Mdeltae = 0.0
# Mdeltat = 0.0
g = 9.8  # 重力加速度
M = 285.0  # 马赫数对应声速
M0 = 0.9  # 初始马赫数
# V0 = 266.0  # 初始速度
V0 = M0*M
h = 11000.0  # 飞行高度
rho = 0.0371  # 空气密度
if __name__ == "__main__":
    print('*'*80)
    print("Plane类测试")
    print('*'*80)
    a = 285.0  # 空速？
    a0 = 3.62  # 初始迎角
    print(f"V0:{V0}")
    # 设置飞机物理参数
    m = 9000.0
    l = 10.0
    Sw = 27.95
    Ix = 0.0
    Iy = 7447.0
    Iz = 0.0
    Ixy = 0.0  # 惯性积
    Iyz = 0.0
    Izx = 0.0
    cA = 3.097
    u = 0.0  # 翼展变形率
    v = 0.0  # 翼角变形率
    zt = 0.0
    TV = 0.0
    gamma0 = 0.0
    Tdeltat = 0.0  # 发动机推力系数

    # 气动导数，之后用拟合后的函数取代
    CL0 = 0.246
    CLa = 3.9
    CLM = 0.23
    CLdeltae = 0.0  # 升降舵偏转对升力的影响
    CD0 = 0.0306
    CDa = 0.284
    CDM = 0.055
    Cm0 = 0.0
    Cma = -0.562
    Cmq = -7.06
    Cmadot = -2.792
    CmM = -0.0654
    Cmdeltae = 0.0  # 升降舵对俯仰力矩的系数

    # 初始化飞机
    testPlane = Plane(m, l, Sw, Ix, Iy, Iz, Ixy, Iyz, Izx, cA)
    testPlane.setu(u=u)
    testPlane.setv(v=v)
    testPlane.setzt(zt=zt)
    testPlane.setTV(TV=TV)
    testPlane.setgamma0(gamma0=gamma0)
    testPlane.setTdeltat(Tdeltat=Tdeltat)
    testPlane.setCL0(CL0=CL0)
    testPlane.setCLa(CLa=CLa)
    testPlane.setCLM(CLM=CLM)
    testPlane.setCLdeltae(CLdeltae=CLdeltae)
    testPlane.setCD0(CD0=CD0)
    testPlane.setCDa(CDa=CDa)
    testPlane.setCDM(CDM=CDM)
    testPlane.setCm0(Cm0=Cm0)
    testPlane.setCma(Cma=Cma)
    testPlane.setCmq(Cmq=Cmq)
    testPlane.setCmadot(Cmadot=Cmadot)
    testPlane.setCmM(CmM=CmM)
    testPlane.setCmdeltae(Cmdeltae=Cmdeltae)
    # 求大导数
    testPlane.computeBigDerivative(V0, rho, M0, g)
    # 打印大导数
    testPlane.printBigDerivative()
    testPlane.judgeStability()
    testPlane.printInfo()
