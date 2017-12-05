# coding=utf-8

'''
作者:Li Jianghua(黎江华) 学号：3117001132 班级：机械7003
程序:非线性方程求解
'''
from math import *
import numpy as np
def g(x,theta):
    return cos(x*sin(theta))

def f(x):
    flag1 = True
    flag2 = True
    b = pi
    a = 0.0
    e = 10 ** (-5)
    h1 = b - a
    Tn = h1 * (g(x,a) + g(x,b)) / 2
    T2n = 0
    while flag2:
        h2 = h1 / 2
        S = 0
        xx = a + h2
        while flag1:
            S = S + g(x,xx)
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
    return I/pi


# 测试程序
if __name__ == '__main__':
    #N=input("请输入最大迭代次数：")
    N=50
    i=0
    x=0
    x0 = 2.5
    e=10**(-10)
    while i<N:
        i=i+1
        x1=x0-f(x0)
        x2=x1-f(x1)
        x=x2-(x2-x1)**2/((x2-x1)-(x1-x0))
        print x,x-x2
        if abs(x-x2)<=e:
            break
        else:
            x0=x

    if i>=N:
        print "Nmax次迭代后仍不收敛"
    else:
        print "方程的解为：",x

