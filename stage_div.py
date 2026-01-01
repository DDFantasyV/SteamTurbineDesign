import numpy as np
import iapws.iapws97 as iapws
from stage import GT, G0, ln, alpha_1, h2, miu_n, u, c1, Point_2, Point_1t
from main import Point_c3, n
import pandas as pd

# 选定参数
omega_m = 0.075  # 第一压力级平均反动度
e1 = 1    # 第一压力级部分进汽度
fai = 0.97     # 第一压力级喷嘴速度系数
theta = 3.01  # 径高比
xaz = 0.552    # 末级速比
alpha_z = 90
alpha = 0.065     # 估取重热系数

delta_Gsg = 0.015 * G0 / 3.6       # 高压缸平衡环漏汽，估计为1.5%
G1 = GT - delta_Gsg  # 进入第一压力级的流量
x1 = u / c1   # 第一级速比
xa1 = x1 * fai * (1 - omega_m) ** 0.5
d1 = ((60 * G1 * xa1) / (e1 * np.pi ** 2 * n * ln * 10 ** (-3) * Point_1t.rho * miu_n * (1 - omega_m) ** 0.5 * np.sin(alpha_1 * np.pi / 180))) ** 0.5      # 第一压力级直径
d1 = 0.928
# Gz = G1 - np.sin(alpha_1 * np.pi / 180) * G0 / 3.6
# cz = (2 * 0.02 * delta_Ht * 1000) ** 0.5
# dz = ((Gz * Point_c3.v * theta) / (np.pi * cz * np.sin(alpha_z * np.pi / 180))) ** 0.5 * 1000    # 末级直径
dz = 2

d2 = (dz - d1) * 0.02 + d1     # 第二级直径
d3 = (dz - d1) * 0.04 + d1     # 第三级直径
d4 = (dz - d1) * 0.06 + d1     # 第四级直径
d5 = (dz - d1) * 0.10 + d1     # 第五级直径
d6 = (dz - d1) * 0.16 + d1     # 第六级直径
d7 = (dz - d1) * 0.26 + d1     # 第七级直径
d8 = (dz - d1) * 0.42 + d1     # 第八级直径
d9 = (dz - d1) * 0.60 + d1     # 第九级直径
d10 = (dz - d1) * 0.80 + d1     # 第十级直径

xa2 = (xaz - xa1) * 0.02 + xa1     # 第二级速比
xa3 = (xaz - xa1) * 0.04 + xa1     # 第三级速比
xa4 = (xaz - xa1) * 0.06 + xa1     # 第四级速比
xa5 = (xaz - xa1) * 0.10 + xa1     # 第五级速比
xa6 = (xaz - xa1) * 0.16 + xa1     # 第六级速比
xa7 = (xaz - xa1) * 0.26 + xa1     # 第七级速比
xa8 = (xaz - xa1) * 0.42 + xa1     # 第八级速比
xa9 = (xaz - xa1) * 0.60 + xa1     # 第九级速比
xa10 = (xaz - xa1) * 0.80 + xa1     # 第十级速比

delta_h_t1 = 12.3245 * (d1 / xa1) ** 2   # 第一级理想比焓降
delta_h_t2 = 12.3245 * (d2 / xa2) ** 2   # 第二级理想比焓降
delta_h_t3 = 12.3245 * (d3 / xa3) ** 2   # 第三级理想比焓降
delta_h_t4 = 12.3245 * (d4 / xa4) ** 2   # 第四级理想比焓降
delta_h_t5 = 12.3245 * (d5 / xa5) ** 2   # 第五级理想比焓降
delta_h_t6 = 12.3245 * (d6 / xa6) ** 2   # 第六级理想比焓降
delta_h_t7 = 12.3245 * (d7 / xa7) ** 2   # 第七级理想比焓降
delta_h_t8 = 12.3245 * (d8 / xa8) ** 2   # 第八级理想比焓降
delta_h_t9 = 12.3245 * (d9 / xa9) ** 2   # 第九级理想比焓降
delta_h_t10 = 12.3245 * (d10 / xa10) ** 2     # 第十级理想比焓降
delta_h_tz = 12.3245 * (dz / xaz) ** 2   # 末级理想比焓降

delta_h = (delta_h_t1 + delta_h_t2 + delta_h_t3 + delta_h_t4 + delta_h_t5 + delta_h_t6 + delta_h_t7 + delta_h_t8 + delta_h_t9 + delta_h_t10 + delta_h_tz) / 11  # 平均理想焓降

z = (1 + alpha) * (h2 - Point_c3.h) / delta_h
z = 14
print(h2 - Point_c3.h)

alpha_ = 3.9 * 10 ** (-4) * (h2 - Point_c3.h) / 7 * z / (z + 1)
delta_alpha = abs((alpha_ - alpha) / alpha_)
if delta_alpha < 0.01:
     print("重热系数校核通过")


d2_re =  (dz - d1) * 0.01 + d1
d3_re =  (dz - d1) * 0.018 + d1
d4_re =  (dz - d1) * 0.038 + d1
d5_re =  (dz - d1) * 0.06 + d1
d6_re =  (dz - d1) * 0.092 + d1
d7_re =  (dz - d1) * 0.15 + d1
d8_re =  (dz - d1) * 0.22 + d1
d9_re =  (dz - d1) * 0.32 + d1
d10_re =  (dz - d1) * 0.42 + d1
d11_re =  (dz - d1) * 0.53 + d1
d12_re =  (dz - d1) * 0.65 + d1
d13_re =  (dz - d1) * 0.80 + d1

xa2_re = (xaz - xa1) * 0.01 + xa1
xa3_re = (xaz - xa1) * 0.018 + xa1
xa4_re = (xaz - xa1) * 0.038 + xa1
xa5_re = (xaz - xa1) * 0.06 + xa1
xa6_re = (xaz - xa1) * 0.092 + xa1
xa7_re = (xaz - xa1) * 0.15 + xa1
xa8_re = (xaz - xa1) * 0.22 + xa1
xa9_re = (xaz - xa1) * 0.32 + xa1
xa10_re = (xaz - xa1) * 0.42 + xa1
xa11_re = (xaz - xa1) * 0.53 + xa1
xa12_re = (xaz - xa1) * 0.65 + xa1
xa13_re = (xaz - xa1) * 0.80 + xa1

delta_h_t1_re = 12.3245 * (d1 / xa1) ** 2
delta_h_t2_re = 12.3245 * (d2_re / xa2_re) ** 2
delta_h_t3_re = 12.3245 * (d3_re / xa3_re) ** 2
delta_h_t4_re = 12.3245 * (d4_re / xa4_re) ** 2
delta_h_t5_re = 12.3245 * (d5_re / xa5_re) ** 2
delta_h_t6_re = 12.3245 * (d6_re / xa6_re) ** 2
delta_h_t7_re = 12.3245 * (d7_re / xa7_re) ** 2
delta_h_t8_re = 12.3245 * (d8_re / xa8_re) ** 2
delta_h_t9_re = 12.3245 * (d9_re / xa9_re) ** 2
delta_h_t10_re = 12.3245 * (d10_re / xa10_re) ** 2
delta_h_t11_re = 12.3245 * (d11_re / xa11_re) ** 2
delta_h_t12_re = 12.3245 * (d12_re / xa12_re) ** 2
delta_h_t13_re = 12.3245 * (d13_re / xa13_re) ** 2
delta_h_tz_re = 12.3245 * (dz / xaz) ** 2

delta_h_t1_re_ = 12.3245 * (d1 / xa1) ** 2 - 4.336004795
delta_h_t2_re_ = 12.3245 * (d2_re / xa2_re) ** 2 - 4.336004795
delta_h_t3_re_ = 12.3245 * (d3_re / xa3_re) ** 2 - 4.336004795
delta_h_t4_re_ = 12.3245 * (d4_re / xa4_re) ** 2 - 4.336004795
delta_h_t5_re_ = 12.3245 * (d5_re / xa5_re) ** 2 - 4.336004795
delta_h_t6_re_ = 12.3245 * (d6_re / xa6_re) ** 2 - 4.336004795
delta_h_t7_re_ = 12.3245 * (d7_re / xa7_re) ** 2 - 4.336004795
delta_h_t8_re_ = 12.3245 * (d8_re / xa8_re) ** 2 - 4.336004795
delta_h_t9_re_ = 12.3245 * (d9_re / xa9_re) ** 2 - 4.336004795
delta_h_t10_re_ = 12.3245 * (d10_re / xa10_re) ** 2 - 4.336004795
delta_h_t11_re_ = 12.3245 * (d11_re / xa11_re) ** 2 - 4.336004795
delta_h_t12_re_ = 12.3245 * (d12_re / xa12_re) ** 2 - 4.336004795
delta_h_t13_re_ = 12.3245 * (d13_re / xa13_re) ** 2 - 4.336004795
delta_h_tz_re_ = 12.3245 * (dz / xaz) ** 2 - 4.336004795

zhijing = '直径'
subi = '速比'
bihanjiang = '比焓降'
tiaozhenghoubihanjiang = '调整后比焓降'

df = pd.DataFrame(columns = ['级序号' , '1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14'])
new_row1 = pd.DataFrame({'级序号': [zhijing] , '1': [d1] , '2': [d2_re] , '3': [d3_re] ,'4': [d4_re] , '5': [d5_re] , '6': [d6_re] , '7': [d7_re] , '8': [d8_re] , '9': [d9_re] , '10': [d10_re] , '11': [d11_re] , '12': [d12_re] , '13': [d13_re] , '14': [dz]}) 
new_row2 = pd.DataFrame({'级序号': [subi] , '1': [xa1] , '2': [xa2_re] , '3': [xa3_re] ,'4': [xa4_re] , '5': [xa5_re] , '6': [xa6_re] , '7': [xa7_re] , '8': [xa8_re] , '9': [xa9_re] , '10': [xa10_re] , '11': [xa11_re] , '12': [xa12_re] , '13': [xa13_re] , '14': [xaz]})    
new_row3 = pd.DataFrame({'级序号': [bihanjiang] , '1': [delta_h_t1_re] , '2': [delta_h_t2_re] , '3': [delta_h_t3_re] ,'4': [delta_h_t4_re] , '5': [delta_h_t5_re] , '6': [delta_h_t6_re] , '7': [delta_h_t7_re] , '8': [delta_h_t8_re] , '9': [delta_h_t9_re] , '10': [delta_h_t10_re] , '11': [delta_h_t11_re] , '12': [delta_h_t12_re] , '13': [delta_h_t13_re] , '14': [delta_h_tz_re]})  
new_row4 = pd.DataFrame({'级序号': [tiaozhenghoubihanjiang] , '1': [delta_h_t1_re_] , '2': [delta_h_t2_re_] , '3': [delta_h_t3_re_] ,'4': [delta_h_t4_re_] , '5': [delta_h_t5_re_] , '6': [delta_h_t6_re_] , '7': [delta_h_t7_re_] , '8': [delta_h_t8_re_] , '9': [delta_h_t9_re_] , '10': [delta_h_t10_re_] , '11': [delta_h_t11_re_] , '12': [delta_h_t12_re_] , '13': [delta_h_t13_re_] , '14': [delta_h_tz_re_]})
df = pd.concat([df, new_row1, new_row2, new_row3, new_row4], ignore_index=True)
print(df)  

#详细计算（第一压力级）

e = 0.89    # 部分进汽度
psi = 0.928  # 动叶速度系数
beta_2 = 20.5     # 动叶出汽角 19-22
miu_1 = 0
a = 1.2     # 经验系数 单列级1.2
K1 = 1.2   # 经验系数 1.0-1.3
sn = 6      # 喷嘴组数
Be = 0.15   # 与级类型有关的系数 单列级一般0.15
ec = 0  # 装有护罩弧段长度与整个圆周长度之比
Ce = 0.012   # 与级类型有关的系数 单列级一般0.012

G = G1  # 进入压力级进气量
delta_h_n = (1-omega_m) * delta_h_t1_re_  # 喷嘴的滞止理想比焓降
c1t = (2 * delta_h_n * 1000) ** 0.5
c1 = fai * c1t      # 喷嘴出口汽流速度
h1t = h2 - delta_h_n
Point_pzt = iapws.IAPWS97(h=h1t, s=Point_2.s)
epsilon_n = Point_pzt.P / Point_2.P     # 喷嘴压比
An = (G * Point_pzt.v) / (miu_n * c1t) * 10 ** 4    # 喷嘴出口面积
ca = (2 * delta_h_t1_re_ * 1000) ** 0.5     # 级的假想速度
u1 = ca * xa1     # 级的圆周速度
ln_ = An * 100 / (e * np.pi * d1 * 1000 * np.sin(alpha_1 * np.pi / 180))    # 喷嘴高度
ln_ = 32  #为方便制造，取整数32
dhn = (1 - fai ** 2) * delta_h_n   # 喷嘴损失
h1 = h1t + dhn   # 喷嘴出口比焓值
Point_pz = iapws.IAPWS97(h=h1, P=Point_pzt.P)
w1 = (c1 ** 2 + u1 ** 2 - 2 * c1 * u1 * np.cos(alpha_1 * np.pi / 180)) ** 0.5     # 动叶进口汽流相对速度
beta_1 = np.arcsin((c1 * np.sin(alpha_1 * np.pi / 180)) / w1) * 180 / np.pi     # 进汽角 25-40
dhw1 = w1 ** 2 / 2000

w2t = (2 * omega_m * delta_h_t1_re_ * 1000 + w1 ** 2) ** 0.5
w2 = psi * w2t      # 动叶出口相对速度
h2t = h1 - omega_m * delta_h_t1_re_
Point_dy = iapws.IAPWS97(h=h2t, s=Point_pz.s)
Ab = (G * Point_dy.v * 10 ** 4) / (miu_n * w2t)     # 动叶出口面积
lb = ln_ + 1 + 1     # 动叶高度
c2 = (w2 ** 2 + u1 ** 2 - 2 * w2 * u1 * np.cos(beta_2 * np.pi / 180)) ** 0.5
alpha_2 = np.arcsin((w2 * np.sin(beta_2 * np.pi / 180)) / c2) * 180 / np.pi
dhb = (1 - psi ** 2) * w2t ** 2 / 2000      # 动叶损失
h2 = h2t + dhb      # 动叶出口比焓值
dhc2 = c2 ** 2 / 2000   # 余速损失
dhu = dhn + dhb + dhc2  # 轮周损失
delta_hu = delta_h_t1_re_ - dhu     # 轮周有效比焓降
eta_u = delta_hu / (delta_h_t1_re_ - miu_1 * dhc2)  # 轮周效率
Pu1 = u1 * (c1 * np.cos(alpha_1 * np.pi / 180) + c2 * np.cos(alpha_2 * np.pi / 180)) / 1000
eta_u_ = Pu1 / (delta_h_t1_re_ - miu_1 * dhc2)
deta_u = abs((eta_u_ - eta_u) / eta_u)
if deta_u < 0.01:
    print("轮周效率校核通过")

dhl = a / ln_ * delta_hu    # 叶高损失
dh_theta = 0.7 * (lb / d1 / 1000) ** 2 * (delta_h_t1_re_ - miu_1 * dhc2)    # 扇形损失
v = (Point_pz.v + Point_dy.v) / 2
delta_Pf = K1 * (u / 100) ** 3 * (d1 * 1000 * 10 ** (-3)) ** 2 / v
dhf = delta_Pf / G     # 叶轮摩擦损失
zeta_w = Be / e * (1 - e - ec / 2) * xa1 ** 3    # 鼓风损失
zeta_s = Ce * sn * xa1 / (e * d1 * 1000 * 10 ** (-3))   # 斥汽损失
zeta_e = zeta_w + zeta_s
dhe = zeta_e * (delta_h_t1_re_ - miu_1 * dhc2)  # 部分进汽损失
sigma_dh = dhl + dh_theta + dhf + dhe   # 级内各项损失之和
delta_hi = delta_hu - sigma_dh  # 级的有效比焓降
eta_i = delta_hi / (delta_h_t1_re_ - miu_1 * dhc2)  # 级效率
Pi = G * delta_hi  # 级的内功率
