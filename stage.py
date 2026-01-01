import numpy as np
import iapws.iapws97 as iapws
from main import Pel, delta_Ht, eta_ri, eta_m, eta_g, Point_01

# 选定参数
m = 1.2     # 抽汽回热
delta_D0 = 0.03     # 漏汽量
delta_h0_t = 74.1   # 调节级比焓降
xa = 0.408  # 调节级速比
dm = 1000   # 调节级平均直径
omega_m = 0.08  # 平均反动度
e = 0.89    # 部分进汽度
alpha_1 = 16.0     # 喷嘴出汽角 13-17
beta_2 = 20.5     # 动叶出汽角 19-22
fai = 0.97  # 喷嘴速度系数
miu_n = 0.97    # 喷嘴流量系数
psi = 0.928  # 动叶速度系数
miu_b = 0.964   # 动叶流量系数
miu_1 = 0


D0 = 3600 * Pel * m / (delta_Ht * eta_ri * eta_m * eta_g) / (1 - delta_D0)
G0 = D0
delta_Gv1 = 0.01 * D0
GT = (G0 - delta_Gv1) / 3.6
Gn = GT     # 调节级喷嘴流量
delta_h0_n = (1 - omega_m) * delta_h0_t     # 喷嘴的滞止理想比焓降
c1t = (2 * delta_h0_n * 1000) ** 0.5
c1 = fai * c1t      # 喷嘴出口汽流速度
h1t = Point_01.h - delta_h0_n
Point_1t = iapws.IAPWS97(h=h1t, s=Point_01.s)
epsilon_n = Point_1t.P / Point_01.P     # 喷嘴压比
An = (Gn * Point_1t.v) / (miu_n * c1t) * 10 ** 4    # 喷嘴出口面积
ca = (2 * delta_h0_t * 1000) ** 0.5     # 级的假想速度
u = ca * xa     # 级的圆周速度
ln = An * 100 / (e * np.pi * dm * np.sin(alpha_1 * np.pi / 180))    # 喷嘴高度
ln = 24
dhn = (1 - fai ** 2) * delta_h0_n   # 喷嘴损失
h1 = Point_1t.h + dhn   # 喷嘴出口比焓值
Point_1 = iapws.IAPWS97(h=h1, P=Point_1t.P)
w1 = (c1 ** 2 + u ** 2 - 2 * c1 * u * np.cos(alpha_1 * np.pi / 180)) ** 0.5     # 动叶进口汽流相对速度
beta_1 = np.arcsin((c1 * np.sin(alpha_1 * np.pi / 180)) / w1) * 180 / np.pi     # 进汽角 25-40
dhw1 = w1 ** 2 / 2000

w2t = (2 * omega_m * delta_h0_t * 1000 + w1 ** 2) ** 0.5
w2 = psi * w2t      # 动叶出口相对速度
h2t = h1 - omega_m * delta_h0_t
Point_2 = iapws.IAPWS97(h=h2t, s=Point_1.s)
Ab = (Gn * Point_2.v * 10 ** 4) / (miu_n * w2t)     # 动叶出口面积
lb = ln + 1 + 1     # 动叶高度
# beta_2 = np.arcsin(Ab * 100 / (dm * lb * np.pi)) * 180 / np.pi
c2 = (w2 ** 2 + u ** 2 - 2 * w2 * u * np.cos(beta_2 * np.pi / 180)) ** 0.5
alpha_2 = np.arcsin((w2 * np.sin(beta_2 * np.pi / 180)) / c2) * 180 / np.pi
dhb = (1 - psi ** 2) * w2t ** 2 / 2000      # 动叶损失
h2 = h2t + dhb      # 动叶出口比焓值
dhc2 = c2 ** 2 / 2000   # 余速损失
dhu = dhn + dhb + dhc2  # 轮周损失
delta_hu = delta_h0_t - dhu     # 轮周有效比焓降
eta_u = delta_hu / (delta_h0_t - miu_1 * dhc2)  # 轮周效率
Pu1 = u * (c1 * np.cos(alpha_1 * np.pi / 180) + c2 * np.cos(alpha_2 * np.pi / 180)) / 1000
eta_u_ = Pu1 / (delta_h0_t - miu_1 * dhc2)
deta_u = abs((eta_u_ - eta_u) / eta_u)

print(alpha_1)
print(beta_1)
print(beta_2)

if (13 <= alpha_1 <= 17
        and 25 <= beta_1 <= 40
        and 19 <= beta_2 <= 22):
    print("角度校核通过")
if deta_u < 0.01:
    print("轮周效率校核通过")


a = 1.2     # 经验系数 单列级1.2
K1 = 1.2   # 经验系数 1.0-1.3
sn = 6      # 喷嘴组数
Be = 0.15   # 与级类型有关的系数 单列级一般0.15
ec = 0  # 装有护罩弧段长度与整个圆周长度之比
Ce = 0.012   # 与级类型有关的系数 单列级一般0.012


dhl = a / ln * delta_hu    # 叶高损失
dh_theta = 0.7 * (lb / dm) ** 2 * (delta_h0_t - miu_1 * dhc2)    # 扇形损失
v = (Point_1.v + Point_2.v) / 2
delta_Pf = K1 * (u / 100) ** 3 * (dm * 10 ** (-3)) ** 2 / v
dhf = delta_Pf / Gn     # 叶轮摩擦损失
zeta_w = Be / e * (1 - e - ec / 2) * xa ** 3    # 鼓风损失
zeta_s = Ce * sn * xa / (e * dm * 10 ** (-3))   # 斥汽损失
zeta_e = zeta_w + zeta_s
dhe = zeta_e * (delta_h0_t - miu_1 * dhc2)  # 部分进汽损失
sigma_dh = dhl + dh_theta + dhf + dhe   # 级内各项损失之和
h2_ = h2 + sigma_dh
Point_2_ = iapws.IAPWS97(h=h2_, P=Point_2.P)

delta_hi = delta_hu - sigma_dh  # 级的有效比焓降
eta_i = delta_hi / (delta_h0_t - miu_1 * dhc2)  # 级效率
Pi = Gn * delta_hi  # 级的内功率