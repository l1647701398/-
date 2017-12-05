# coding=utf-8

'''
作者:Li Jianghua 学号：3117001132 班级：机械7003
程序:多项式曲线拟合算法
'''
import matplotlib.pyplot as plt
import math
import numpy as np

def lagrange(x, y, _x):
    L = 0
    n = len(x)
    for i in np.arange(n):
        wx = 1
        a = 0
        b = 0
        for j in np.arange(n):
            if j != i:
                a = _x - x[j]
                b = x[i] - x[j]
                wx = wx * a / b
        # print wx
        L = L + y[i] * wx

    # for i in np.arange(n):
    #     di=1.0
    #     for j in np.arange(n):
    #         if j!=i:
    #             a=x[i]-x[j]
    #             print a
    #             di=di*a
    #     d.append(di)
    #     print d
    # a=0
    # for k in np.arange(n):
    #     a= _x-x[k]
    #     w= a*w
    #     a= y[k]/(a*d[k])
    #     L= L+a
    return L


def newton(x, yy, _x):
    n = len(yy)
    for k in np.arange(1, n):
        for i in np.arange(k, n):
            i = n + k - i - 1
            yy[i] = (yy[i] - yy[i - 1]) / (x[i] - x[i - k])
    # print y
    _yy = yy[n - 1]
    for j in np.arange(n - 1):
        j = n - 2 - j
        _yy = _yy * (_x - x[j]) + yy[j]

    return _yy


# TSS
def tss(a, b, c, d):
    n = len(b)
    u = np.zeros(n)
    y = np.zeros(n)
    u[0] = b[0]
    y[0] = d[0]
    l = np.zeros(n)
    x = np.zeros(n)
    for k in np.arange(1, n):
        l[k] = a[k] / u[k - 1]
        u[k] = b[k] - l[k] * c[k - 1]
        y[k] = d[k] - l[k] * y[k - 1]
    x[n - 1] = y[n - 1] / u[n - 1]
    print
    for k in np.arange(n - 1):
        k = n - 2 - k

        x[k] = (y[k] - c[k] * x[k + 1]) / u[k]
    return x


def findk(x, _x):
    k = 1

    n = len(x)
    for i in np.arange(1, n - 1):

        if _x <= x[i]:
            k = i
            return k
        else:
            k = i + 1

    return k


# splinem
def splinem(x, y, r0, d0, un, dn):
    n = len(x)
    M = np.zeros(n)
    h = np.zeros(n)
    a = np.zeros(n)
    b = np.zeros(n)
    c = np.zeros(n)
    for i in np.arange(n):
        M[i] = y[i]

    for k in [1, 2]:
        for j in np.arange(k, n):
            j = n + k - j - 1
            M[j] = (M[j] - M[j - 1]) / (x[j] - x[j - k])

    h[1] = x[1] - x[0]
    for i in np.arange(1, n - 1):
        h[i + 1] = x[i + 1] - x[i]
        c[i] = h[i + 1] / (h[i] + h[i + 1])
        a[i] = 1 - c[i]
        b[i] = 2
        M[i] = 6 * M[i + 1]

    M[0] = d0
    M[n - 1] = dn
    c[0] = r0
    b[0] = 2
    a[n - 1] = un
    b[n - 1] = 2
    # TSS
    M = tss(a, b, c, M)
    return M


# 3次样条插值
def evasoline(x, y, _x, r0, d0, un, dn):
    M = splinem(x, y, r0, d0, un, dn)
    _y = np.zeros(len(_x))
    for i in np.arange(len(_x)):
        k = findk(x, _x[i])
        h = x[k] - x[k - 1]
        x_ = x[k] - _x[i]
        x__ = _x[i] - x[k - 1]

        _y[i] = (M[k - 1] * (x_ ** 3) / 6 + M[k] * (x__ ** 3) / 6 + (y[k - 1] - M[k - 1] * (h ** 2) / 6) * x_ + (
                y[k] - M[k] * (h ** 2) / 6) * x__) / h
        # print _y[i]
    return _y


# 3次样条插值导数
def d_evasoline(x, y, _x, M):
    _dy = 0
    k = findk(x, _x)
    h = x[k] - x[k - 1]
    x_ = x[k] - _x
    x__ = _x - x[k - 1]
    _dy = -M[k - 1] * (x_ ** 2) / (2 * h) + M[k] * (x__ ** 2) / (2 * h) + (y[k] - y[k - 1]) / h - (
            M[k] - M[k - 1]) * h / 6
    _ddy = M[k - 1] * (x_) / h + M[k] * x__ / h
    # print _dy
    # return _dy
    return math.sqrt(1 + _dy ** 2)


# 画出拟合后的曲线
if __name__ == '__main__':
    fig = plt.figure()
    axes = fig.add_subplot(111)

    # 生成曲线上的各个点
    x = np.arange(21)
    y = [9.01, 8.96, 7.96, 7.97, 8.02, 9.05, 10.13, 11.18, 12.26, 13.28, 13.32, 12.61, 11.29, 10.22, 9.15, 7.9, 7.95,
         8.86, 9.81, 10.8, 10.93]
    yy = np.array(y)

    _x = np.arange(0, 20.009, 0.01)

    M = splinem(x, y, r0=0, d0=0, un=0, dn=0)
    # _x=20.0
    # _y=-lagrange(x, y, _x)
    # _y=-newton(x, yy, _x)
    _y = -evasoline(x, yy, _x, r0=0, d0=0, un=0, dn=0)
    # print _y
    flag1 = True
    flag2 = True
    b = 20.0
    a = 0.0
    e = 10 ** (-10)
    h1 = b - a
    T2n = 0
    Tn = h1 * (d_evasoline(x, yy, a, M) + d_evasoline(x, yy, b, M)) / 2

    while flag2:
        h2 = h1 / 2
        S = 0
        xx = a + h2
        while flag1:
            S = S + d_evasoline(x, yy, xx, M)
            xx = xx + h1
            if xx >= b:
                flag1 = False
        T2n = (Tn + h1 * S) / 2

        if abs(T2n - Tn) < e:

            I = T2n
            flag2 = False
        else:
            Tn = T2n
            h1 = h2
            flag2 = True
        flag1 = True
    print "所需光缆长度近似值为：", I

    axes.plot(_x, _y, color='r', linestyle='-', marker='*', label='data source')
    axes.plot(x, -np.array(y), color='k', linestyle='-', marker='*', label='data fitting')
    axes.plot(x, np.zeros_like(x), color='g', linestyle='--', marker='*', label='water level')
    axes.legend()
    plt.title(u"河道深度图")
    plt.show()
