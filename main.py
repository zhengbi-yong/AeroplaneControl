from plane import Plane

if __name__ == "__main__":
    # 设置客观世界参数
    g = 9.8  # 重力加速度
    M = 340.0  # 马赫数对应声速
    M0 = 0.088  # 初始马赫数
    V0 = M0*M
    h = 2000.0  # 飞行高度
    rho = 0.0371  # 空气密度
    alpha0 = 0.0  # 初始迎角
    # 设置飞机物理参数
    m = 1247.0
    l = 10.18
    Sw = 17.09
    Ix = 1420.9
    Iy = 4067.5
    Iz = 4786.0
    Ixy = 0.0  # 惯性积
    Iyz = 0.0
    Izx = 0.0
    cA = 1.74
    u = 0.0  # 翼展变形率
    v = 0.0  # 翼角变形率
    zt = 0.0
    TV = 0.0
    gamma0 = 0.0
    Tdeltat = 0.0  # 发动机推力系数
    # 气动导数，之后用拟合后的函数取代
    # CL0 = 0.246
    # CLa = 3.9
    # CLM = 0.23
    CLdeltae = 0.0119  # 升降舵偏转对升力的影响
    # CD0 = 0.0306
    # CDa = 0.284
    # CDM = 0.055
    # Cm0 = 0.0
    # Cma = -0.562
    # Cmq = -7.06
    # Cmadot = -2.792
    # CmM = -0.0654
    Cmdeltae = -0.0385  # 升降舵对俯仰力矩的系数

    # 初始化飞机
    TransformablePlane = Plane(m, l, Sw, Ix, Iy, Iz, Ixy, Iyz, Izx, cA)
    TransformablePlane.setg(g=g)
    TransformablePlane.setM(M=M)
    TransformablePlane.setM0(M0=M0)
    TransformablePlane.setV0(V0=V0)
    TransformablePlane.seth(h=h)
    TransformablePlane.setrho(rho=rho)
    TransformablePlane.setalpha0(alpha0=alpha0)
    TransformablePlane.setzt(zt=zt)
    TransformablePlane.setTV(TV=TV)
    TransformablePlane.setgamma0(gamma0=gamma0)
    TransformablePlane.setTdeltat(Tdeltat=Tdeltat)
    TransformablePlane.transform(u, v)
    # TransformablePlane.setCL0(CL0=CL0)
    # TransformablePlane.setCLa(CLa=CLa)
    # TransformablePlane.setCLM(CLM=CLM)
    # TransformablePlane.setCLdeltae(CLdeltae=CLdeltae)
    # TransformablePlane.setCD0(CD0=CD0)
    # TransformablePlane.setCDa(CDa=CDa)
    # TransformablePlane.setCDM(CDM=CDM)
    # TransformablePlane.setCm0(Cm0=Cm0)
    # TransformablePlane.setCma(Cma=Cma)
    # TransformablePlane.setCmq(Cmq=Cmq)
    # TransformablePlane.setCmadot(Cmadot=Cmadot)
    # TransformablePlane.setCmM(CmM=CmM)
    # TransformablePlane.setCmdeltae(Cmdeltae=Cmdeltae)
    # 求大导数
    TransformablePlane.computeBigDerivative()
    # 打印大导数
    TransformablePlane.printBigDerivative()
    TransformablePlane.judgeStability()
    TransformablePlane.printInfo()
    # TransformablePlane.transform(0.5, 1)
    TransformablePlane.printBigDerivative()
    TransformablePlane.printInfo()
