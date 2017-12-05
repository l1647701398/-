# coding=utf-8

'''
作者:Li Jianghua（黎江华） 学号：3117001132 班级：机械7003
程序:线性方程组求解
'''
import matplotlib.pyplot as plt
import numpy as np
import math
from struct import *

def GAUSS(A,b):
    m,n=A.shape
    x=np.zeros(n)
    for k in range(n-1):
        if A[k][k]==0:
            print "无法求解"
            return -1
        for i in range(k+1,n):
            aik=A[i][k]/A[k][k]
            for j in range(k+1,n):
                A[i][j]=A[i][j]-aik*A[k][j]
            b[i]=b[i]-aik*b[k]


    x[n-1]=b[n-1]/A[n-1][n-1]
    for k in range(n-2,-1,-1):
        S=b[k]
        for j in range(k+1,n):
            S=S-A[k][j]*x[j]
            x[k]=S/A[k][k]

    return x

# 测试程序
if __name__ == '__main__':
    file = open('dat20172.dat', 'rb')
    s = file.read(8)
    s = s[::-1]
    id = ''.join([hex(ord(c)).replace('0x', '') for c in s[4:8]])
    ver = ''.join([hex(ord(c)).replace('0x', '0') for c in s[:4]])
    print "id:",id ,'ver:' ,ver
    s=file.read(4)
    #print ' '.join([hex(ord(c)).replace('0x', '0') for c in s[:4]])
    n,= unpack("i",s[:4])
    s = file.read(4)
    q, = unpack("i", s[:4])
    s = file.read(4)
    p, = unpack("i", s[:4])
    print "方程阶数:",n,"带状矩阵上带宽:",q,"带状矩阵下带宽:",p
    A=np.zeros([n,n])
    b = np.zeros(n)
    if int(ver)==102:
        for i in range(n):
            for j in range(n):
                s=file.read(4)
                A[i][j],=unpack('f',s)

        for k in range(n):
            s = file.read(4)
            b[k], = unpack('f', s)
    if int(ver)==202:
        AA=np.zeros([n,p+q+1])
        for i in range(n):
            for j in range(q+p+1):
                s=file.read(4)
                AA[i][j],=unpack('f',s)
        for k in range(n):
            s=file.read(4)
            b[k],=unpack('f',s)
        for i in range(n):
            if i<p:
                for j in range(i+q):
                    A[i][j]=AA[i][p+j]
            elif i>=p and i<=n-q:
                for j in range(p+q+1):
                    A[i][i-p+j]=AA[i][j]
            elif i>n-q:
                for j in range(n-i+p):
                    A[i][i-p+j]=AA[i][j]



    #print b
    X=GAUSS(A,b)
    print X


