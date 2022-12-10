class Plane(object):
    def __init__(self, m, l, Sw, Ix, Iy, Iz, Ixy, Iyz, Izx, cA):
        # 客观世界参数与初始条件
        self._g = 9.8  # 重力加速度
        self._M = 340.0  # 马赫数对应声速
        self._M0 = 0.0  # 初始马赫数
        self._V0 = 0.0
        self._h = 0.0  # 飞行高度
        self._rho = 0.0371  # 空气密度
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
        self._alpha0 = 0.0
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
        self._CDa2 = 0.0
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
        self._wnp = 0.0
        self._zetap = 0.0
        self._wnsp = 0.0
        self._zetasp = 0.0
        self._pStablity = True
        self._spStability = True
        self._Stability = True

    def setg(self, g):
        self._g = g

    def setM(self, M):
        self._M = M

    def setM0(self, M0):
        self._M0 = M0

    def setV0(self, V0):
        self._V0 = V0

    def seth(self, h):
        self._h = h

    def setrho(self, rho):
        self._rho = rho

    def setalpha0(self, alpha0):
        self._alpha0 = alpha0

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

    def setCDa2(self, CDa2):
        self._CDa2 = CDa2

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
        print(f"质量:  {self._m:>15} kg")
        print(f"翼展:  {self._l:>15} m")
        print(f"翼面积:{self._Sw:>15} m^2")
        print(f"重力加速度g:  {self._g:>10} kg*m/s^2")
        print(f"当地马赫数M:  {self._M:>10} m/s")
        print(f"初始马赫数M0: {self._M0:>10}")
        print(f"初始速度:     {self._V0:>10} m/s")
        print(f"飞行高度h:    {self._h:>10} m")
        print(f"空气密度rho:  {self._rho:>10}")
        print(f"初始迎角:     {self._alpha0:>10}")
        print(f"升降舵偏转对升力的影响CLdeltae: {self._CLdeltae:>10}")
        print(f"升降舵对俯仰力矩的系数Cmdeltae: {self._Cmdeltae:>10}")
        print(
            f"转动惯量 Ix :{self._Ix:>10} kg*m^2 Iy:{self._Iy:>10} kg*m^2 Iz:{self._Iz:>10} kg*m^2")
        print(
            f"惯性积   Ixy:{self._Ixy:>10} m^4    Iyz:{self._Iyz:>9} m^4    Izx:{self._Izx:>9} m^4")
        print(f"平均气动弦长cA:             {self._cA} m")
        print(f"翼展变形率lambda:           {self._u}")
        print(f"翼角变形率rho:              {self._v}")
        print(f"零迎角升力系数CL0:          {self._CL0}")
        print(f"迎角升力系数CLa:            {self._CLa}")
        print(f"马赫数升力系数CLM:          {self._CLM}")
        print(f"零迎角阻力系数CD0:          {self._CD0}")
        print(f"迎角阻力系数CDa:            {self._CDa}")
        print(f"马赫数阻力系数CDM:          {self._CDM}")
        print(f"迎角力矩系数Cma:            {self._Cma}")
        print(f"地面轴角速度力矩系数Cmq:    {self._Cmq}")
        print(f"机体轴角速度力矩系数Cmadot: {self._Cmadot}")
        print(f"马赫数力矩系数CmM:          {self._CmM}")
        print(f"发动机推力系数Tdeltat:      {self._Tdeltat}")
        print(f"推力对速度偏导数:           {self._TV}")
        print(f"gamma0:                     {self._gamma0}")
        print(f"推力矢量到机体x轴的距离zt:  {self._zt} m")
        print(f"稳定性:                     {self._Stability}")
        print('-'*80)
        print('*'*80)

    def computeDV(self):
        self._DV = (1/self._V0)*((1/2)*self._rho*(self._V0**2)) * \
            self._Sw*(2*self._CD0+self._M0*self._CDM)
        # print('*'*20)
        # print("computeDV")
        # print('*'*20)
        # print(f"V0:{self._V0}")
        # print(f"rho:{self._rho}")
        # print(f"Sw:{self._Sw}")
        # print(f"CD0:{self._CD0}")
        # print(f"M0:{self._M0}")
        # print(f"CDM:{self._CDM}")
        # print(f"DV:{self._DV}")

    def computeDa(self,):
        self._Da = self._CDa*(1/2)*self._rho*self._V0**2*self._Sw
        # print('*'*20)
        # print("computeDa")
        # print('*'*20)
        # print(f"CDa:{self._CDa}")
        # print(f"rho:{rho}")
        # print(f"V0:{V0}")
        # print(f"Sw:{self._Sw}")
        # print(f"Da:{self._Da}")

    def computeLV(self):
        self._LV = (1/self._V0)*((1/2)*self._rho*(self._V0**2)) * \
            self._Sw*(2*self._CL0+self._M0*self._CLM)
        # print('*'*20)
        # print("computeLV")
        # print('*'*20)
        # print(f"V0:{self._V0}")
        # print(f"rho:{self._rho}")
        # print(f"Sw:{self._Sw}")
        # print(f"CL0:{self._CL0}")
        # print(f"M0:{self._M0}")
        # print(f"CLM:{self._CLM}")
        # print(f"LV:{self._LV}")

    def computeLa(self):
        self._La = self._CLa*(1/2)*self._rho*self._V0**2*self._Sw

    def computeLdeltae(self):
        self._Ldeltae = self._CLdeltae*(1/2)*self._rho*self._V0**2*self._Sw

    def computeMaV(self):
        self._MaV = ((2*self._Cm0+self._M0*self._CmM)/self._V0) * \
            ((1/2)*self._rho*self._V0**2)*self._cA*self._Sw

    def computeMaa(self):
        self._Maa = ((1/2)*self._rho*self._V0**2)*self._cA*self._Sw*self._Cma

    def computeMaadot(self):
        self._Maadot = ((1/2)*self._rho*self._V0**2)*(self._cA **
                                                      2/(2*self._V0))*self._Sw*self._Cmadot

    def computeMaq(self):
        self._Maq = ((1/2)*self._rho*self._V0**2)*(self._cA **
                                                   2/(2*self._V0))*self._Sw*self._Cmq
        # print('*'*20)
        # print("computeMaq")
        # print('*'*20)
        # print(f"rho:{self._rho}")
        # print(f"V0:{self._V0}")
        # print(f"cA:{self._cA}")
        # print(f"Sw:{self._Sw}")
        # print(f"Cmq:{self._Cmq}")

    def computeMadeltae(self):
        self._Madeltae = ((1/2)*self._rho*self._V0**2)*(self._cA**2 /
                                                        (2*self._V0))*self._Sw*self._Cmdeltae

    def computeXV(self):
        self._XV = (self._DV-self._TV)/self._m
        # print('*'*20)
        # print("computeXV")
        # print('*'*20)
        # print(f"DV:{self._DV}")
        # print(f"TV:{self._TV}")
        # print(f"m:{self._m}")
        # print(f"XV:{self._XV}")

    def computeXa(self):
        self._Xa = self._Da/(self._m*self._V0)-self._g/self._V0
        # print('*'*20)
        # print("computeXa")
        # print('*'*20)
        # print(f"Da:{self._Da}")
        # print(f"m:{self._m}")
        # print(f"V0:{self._V0}")
        # print(f"g:{self._g}")
        # print(f"Xa:{self._Xa}")

    def computeXtheta(self):
        self._Xtheta = self._g/self._V0

    def computeZV(self):
        self._ZV = self._LV/(self._m*self._V0)
        # print('*'*20)
        # print("computeZV")
        # print('*'*20)
        # print(f"LV:{self._LV}")
        # print(f"m:{self._m}")
        # print(f"V0:{self._V0}")

    def computeZa(self):
        self._Za = self._La/(self._m*self._V0)

    def computeMV(self):
        self._MV = -(self._V0*(self._MaV+self._TV*self._zt)/self._Iy)

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

    def computeXdeltat(self):
        self._Xdeltat = -(self._Tdeltat/(self._m*self._V0))

    def computeZdeltae(self):
        self._Zdeltae = self._Ldeltae/(self._m*self._V0)

    def computeMdeltae(self):
        self._Mdeltae = -self._Madeltae/self._Iy

    def computeMdeltat(self):
        self._Mdeltat = -(self._Tdeltat*self._zt)/self._Iy

    def computeBigDerivative(self):
        self.computeDV()
        self.computeDa()
        self.computeLV()
        self.computeLa()
        self.computeLdeltae()
        self.computeMaV()
        self.computeMaa()
        self.computeMaadot()
        self.computeMaq()
        self.computeMadeltae()
        self.computeXV()
        self.computeXa()
        self.computeXtheta()
        self.computeXdeltat()
        self.computeZV()
        self.computeZa()
        self.computeZdeltae()
        self.computeMV()
        self.computeMadot()
        self.computeMa()
        self.computeMq()
        self.computeMdeltae()
        self.computeMdeltat()

    def computewnp(self):
        self._wnp = (-self._ZV*self._g/self._V0)**0.5

    def computezetap(self):
        self._zetap = -self._XV/(2*abs(self._wnp))

    def computewnsp(self):
        self._wnsp = (self._ZV*self._Mq/self._V0-self._Ma)**0.5

    def computezetasp(self):
        self._zetasp = -(self._Mq+self._Madot+self._Za /
                         self._V0)/(2*abs(self._wnsp))

    # 一定要在计算完大导数之后再使用该函数判断稳定性
    def judgeStability(self):
        self._a3 = self._Mq+self._Za+self._Madot+self._XV
        self._a2 = self._XV*(self._Mq+self._Za+self._Madot) + \
            self._Za*self._Mq+self._Ma-self._Xa*self._ZV
        self._a1 = self._XV*(self._Za*self._Mq+self._Ma)+self._Xtheta * \
            (self._ZV*self._Madot-self._MV) - \
            self._Xa*(self._ZV*self._Mq+self._MV)
        self._a0 = self._Xtheta*(self._ZV*self._Ma-self._MV*self._Za)
        self.computewnp()
        self.computezetap()
        self.computewnsp()
        self.computezetasp()

        def isStable(zeta):
            if(zeta**2 > 0.0):
                return True
            else:
                return False
        self._pStablity = isStable(self._zetap)
        self._spStablity = isStable(self._zetasp)
        if((self._pStablity == True) and (self._spStability == True)):
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

    def transform(self, u, v):
        self._u = u
        self._v = v
        x = self._u
        y = self._v
        self._CL0 = -0.041*x**2 + 0.001*x*y - 0.019*y**2 + 0.106*x + 0.002*y + 0.076
        self._CLa = -0.007*x**2 - 0.000*x*y - 0.004*y**2 + 0.010*x - 0.001*y + 0.101
        self._CD0 = +0.002*x**2 + 0.000*x*y + 0.001*y**2 + 0.002*x - 0.000*y + 0.016
        self._CDa = -0.000*x**2 + 0.000*x*y - 6.666*y**2 + 0.000*x - 8.333*y + 0.000
        self._CDa2 = +0.000*x**2 + 3.079*x*y + 2.226*y**2 - 0.000*x - 2.436*y + 0.000
        self._Cm0 = +0.123*x**2 - 0.085*x*y + 0.039*y**2 - 0.337*x - 0.091*y + 0.246
        self._Cma = +0.011*x**2 - 0.040*x*y - 6.666*y**2 - 0.023*x - 0.043*y - 0.046
        self.computeBigDerivative()
        self.judgeStability()

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


if __name__ == "__main__":
    print('*'*80)
    print("Plane类测试")
    print('*'*80)
    # 设置客观世界参数
    g = 9.8  # 重力加速度
    M = 285.0  # 马赫数对应声速
    M0 = 0.9  # 初始马赫数
    V0 = M0*M
    h = 11000.0  # 飞行高度
    rho = 0.0371  # 空气密度
    alpha0 = 3.62  # 初始迎角
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
    CLdeltae = 0.0119  # 升降舵偏转对升力的影响
    CD0 = 0.0306
    CDa = 0.284
    CDM = 0.055
    Cm0 = 0.0
    Cma = -0.562
    Cmq = -7.06
    Cmadot = -2.792
    CmM = -0.0654
    Cmdeltae = -0.0385  # 升降舵对俯仰力矩的系数

    # 初始化飞机
    testPlane = Plane(m, l, Sw, Ix, Iy, Iz, Ixy, Iyz, Izx, cA)
    testPlane.setg(g=g)
    testPlane.setM(M=M)
    testPlane.setM0(M0=M0)
    testPlane.setV0(V0=V0)
    testPlane.seth(h=h)
    testPlane.setrho(rho=rho)
    testPlane.setalpha0(alpha0=alpha0)
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
    testPlane.computeBigDerivative()
    # 打印大导数
    testPlane.printBigDerivative()
    testPlane.judgeStability()
    testPlane.printInfo()
    testPlane.transform(0.5, 1)
    testPlane.printBigDerivative()
    testPlane.printInfo()
