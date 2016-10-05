import numpy as np


def stringMethod():
    N_str = 11  # ストリングの「区間」の数
    N = 100
    d_old = 100000  # 収束条件に使用
    count = 0
    dt = 0.001
    s_sum = 0.0  # 地点間の距離の計算に仕様
    tol = N_str ** (-4)
    phi = np.zeros((N_str + 1, 9))  # 配位空間上の点
    phi_old = np.zeros((N_str, 9))
