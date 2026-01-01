import iapws.iapws97 as iapws
import matplotlib.pyplot as plt

# 给定参数
Pel = 200  # 额定功率MW
p0 = 12.74  # 新蒸汽压力MPa
t0 = 550 + 273.15   # 新蒸汽温度K
pc = 0.0042    # 凝汽器压力kPa

# 其他参数
n = 3000    # 汽轮机转速
eta_m = 0.995   # 机械效率
eta_g = 0.985   # 发电效率
eta_ri = 0.88   # 相对内效率
C1 = 0.02   # 余速损失系数%
C2 = 15     # 系数kJ/kg

# 损失估算
delta_p0 = 0.05 * p0
delta_pc = 0.04 * pc


p01 = p0 - delta_p0
pc1 = pc + delta_pc
Point_0 = iapws.IAPWS97(P=p0, T=t0)
Point_01 = iapws.IAPWS97(P=p01, h=Point_0.h)
Point_ct = iapws.IAPWS97(P=pc, s=Point_0.s)
delta_Ht = Point_0.h - Point_ct.h
delta_Hi = delta_Ht * eta_ri
Point_c2 = iapws.IAPWS97(P=pc1, h=Point_0.h-delta_Hi)
delta_hc2 = delta_Ht * C1
Point_c3 = iapws.IAPWS97(P=pc1, h=Point_c2.h-delta_hc2)
Point_41 = iapws.IAPWS97(h=(Point_0.h + Point_c3.h) / 2, s=(Point_01.s + Point_c3.s) / 2)
Point_4 = iapws.IAPWS97(P=Point_41.P, h=Point_41.h - C2)

# Point 0
Temp1 = iapws.IAPWS97(P=p0, h=Point_0.h-50)
Temp2 = iapws.IAPWS97(P=p0, h=Point_0.h+100)
x = [Temp1.s, Point_0.s, Temp2.s]
y = [Point_0.h - 50, Point_0.h, Point_0.h + 100]
plt.plot(x, y, '-')    # 等压线

x = [Point_0.s - 0.1, Point_0.s + 1]
y = [Point_0.h, Point_0.h]
plt.plot(x, y, '-')    # 等焓线

# Point 1
Temp1 = iapws.IAPWS97(P=p01, h=Point_0.h-50)
Temp2 = iapws.IAPWS97(P=p01, h=Point_0.h+100)
x = [Temp1.s, Point_01.s, Temp2.s]
y = [Point_0.h - 50, Point_0.h, Point_0.h + 100]
plt.plot(x, y, '-')    # 等压线

# Point 2t & 2
Temp1 = iapws.IAPWS97(P=pc, h=Point_ct.h-50)
Temp2 = iapws.IAPWS97(P=pc, h=Point_ct.h+50)
x = [Temp1.s, Point_0.s, Temp2.s]
y = [Point_ct.h - 50, Point_ct.h, Point_ct.h + 50]
plt.plot(x, y, '-')    # 背压线
x = [Point_0.s - 0.1, Point_0.s + 0.1]
y = [Point_ct.h, Point_ct.h]
plt.plot(x, y, '-')    # 辅助等焓线
x = [Point_0.s, Point_0.s]
y = [Point_0.h, Point_ct.h]
plt.plot(x, y, 'o-')    # 等熵焓降线
Temp1 = iapws.IAPWS97(P=pc1, h=Point_c2.h-50)
Temp2 = iapws.IAPWS97(P=pc1, h=Point_c2.h+50)
x = [Temp1.s, Point_c2.s, Temp2.s]
y = [Point_c2.h - 50, Point_c2.h, Point_c2.h + 50]
plt.plot(x, y, '-')
x = [Point_c2.s - 0.1, Point_c2.s + 0.1]
y = [Point_c2.h, Point_c2.h]
plt.plot(x, y, '-')    # 辅助等焓线
x = [Point_c2.s, Point_c2.s]
y = [Point_0.h, Point_c2.h]
plt.plot(x, y, '-')    # 实际焓降线

# Point 3
x = [Point_c2.s, Point_c3.s]
y = [Point_c2.h, Point_c3.h]
plt.plot(x, y, 'o-', linewidth=3)
x = [Point_c3.s - 0.1, Point_c3.s + 0.1]
y = [Point_c3.h, Point_c3.h]
plt.plot(x, y, '-')    # 辅助等焓线

# Point 41
x = [Point_01.s, Point_c3.s]
y = [Point_0.h, Point_c3.h]
plt.plot(x, y, 'o-')
Temp1 = iapws.IAPWS97(P=Point_41.P, h=Point_41.h-50)
Temp2 = iapws.IAPWS97(P=Point_41.P, h=Point_41.h+50)
x = [Temp1.s, Point_41.s, Temp2.s]
y = [Point_41.h - 50, Point_41.h, Point_41.h + 50]
plt.plot(x, y, '-')    # 等压线

# Point 4
x = [Point_01.s, Point_4.s, Point_c3.s]
y = [Point_0.h, Point_4.h, Point_c3.h]
plt.plot(x, y, 'o-', linewidth=3)    # 等压线


def draw():
    plt.text(6.5, 2050, "pc={:.4f}MPa".format(pc))
    plt.text(6.9, 2000, "hc={:.4f}kJ/kg".format(Point_c3.h))
    plt.text(7.3, 2110, "Dhc2={:.4f}kJ/kg".format(delta_hc2))
    plt.text(7.3, 2250, "pc1={:.6f}MPa".format(pc1))
    plt.text(7.2, 2800, "DHi={:.4f}kJ/kg".format(delta_Hi))
    plt.text(6.5, 2650, "DHt={:.4f}kJ/kg".format(delta_Ht))
    plt.text(6.4, 3300, "h0={:.4f}kJ/kg".format(Point_0.h))
    plt.text(6.4, 3520, "p0={:.2f}MPa".format(Point_0.P))
    plt.text(6.8, 3520, "p01={:.3f}MPa".format(Point_01.P))
    plt.xlabel('s[kJ/(kg.K)]')
    plt.ylabel('h[kJ/kg]')
    plt.show()


draw()
