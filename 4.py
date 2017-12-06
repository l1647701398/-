# coding=utf-8

'''
作者:Li Jianghua（黎江华） 学号：3117001132 班级：机械7003
程序:线性方程组求解
'''
from numpy import *
import math
from struct import *

#普通矩阵LU分解
def LU(A):
    m,n=A.shape
    UL=zeros([m,n])
    for j in range(n):
        UL[0][j]=A[0][j]
    for k in range(1,n):
        UL[k][0]=A[k][0]/UL[0][0]
    for i in range(1,n):
        for j in range(i,n):
            sum=0
            for t in range(i):
                sum=sum+UL[i][t]*UL[t][j]
            UL[i][j]=A[i][j]-sum

        for k in range(i+1,n):
            sum=0
            for t in range(i):
                sum=sum+UL[k][t]*UL[t][i]
            UL[k][i]=(A[k][i]-sum)/UL[i][i]

    return UL

#带状矩阵LU分解
def LLUU(A,p,q):
    n, m = A.shape
    UULL = zeros([n, m])
    for j in range(q + 1):
        UULL[0][j+p] = A[0][j+p]
    for k in range(1, p + 1):
        UULL[k][p-k] = A[k][p-k] / UULL[0][p]
    for i in range(1, n):
        for j in range(i, min(n, i + q + 1)):
            sum = 0
            for t in range(max(0, i - p , j - q ), i):
                sum = sum + UULL[i][t+p-i] * UULL[t][j+p-t]
            UULL[i][j+p-i] = A[i][j+p-i] - sum

        for k in range(i + 1, min(n, i + p + 1)):
            sum = 0
            for t in range(max(0, k - p , i - q ), i):
                sum = sum + UULL[k][t+p-k] * UULL[t][i+p-t]
            UULL[k][i+p-k] = (A[k][i+p-k] - sum) / UULL[i][i+p-i]

    return UULL

def TSSGAUSS(UULL,p,q,b):
    n,m=UULL.shape
    y=zeros(n)
    y[0]=b[0]
    for k in range(1,n):

        sum=0
        for i in range(max(0,k-p),k):
            sum=sum+UULL[k][i+p-k]*y[i]
        y[k]=b[k]-sum

    x=zeros(n)
    x[n-1]=y[n-1]/UULL[n-1][p]
    for k in range(n-2,-1,-1):
        sum=0
        for j in range(k+1,min(n,k+q+1)):
            sum=sum+UULL[k][j+p-k]*x[j]
        x[k]=(y[k]-sum)/UULL[k][k+p-k]
    return x


# 测试程序
if __name__ == '__main__':
    file = open('dat20171.dat', 'rb')
    s = file.read(8)
    s = s[::-1]
    id = ''.join([hex(ord(c)).replace('0x', '') for c in s[4:8]])
    ver = ''.join([hex(ord(c)).replace('0x', '0') for c in s[:4]])
    print "id:",id ,'ver:' ,int(ver)
    s=file.read(4)
    n,= unpack("i",s[:4])
    s = file.read(4)
    q, = unpack("i", s[:4])
    s = file.read(4)
    p, = unpack("i", s[:4])
    print "方程阶数:",n,"带状矩阵上带宽:",q,"带状矩阵下带宽:",p
    AA = zeros([n, p + q + 1])
    b = zeros(n)
    if int(ver)==102:
        A = zeros([n, n])
        for i in range(n):
            for j in range(n):
                s=file.read(4)
                A[i][j],=unpack('f',s)

        for k in range(n):
            s = file.read(4)
            b[k], = unpack('f', s)

        for i in range(n):
            for j in range(max(0,p-i),min(p+q+1,n+p-i)):
                AA[i][j]=A[i][i+j-p]

    if int(ver)==202:
        for i in range(n):
            for j in range(q+p+1):
                s=file.read(4)
                AA[i][j],=unpack('f',s)
        for k in range(n):
            s=file.read(4)
            b[k],=unpack('f',s)

    UULL=LLUU(AA, p, q)#UULL中存放L和U的对应值（去除零元素）
    X=TSSGAUSS(UULL, p, q, b)
    print X



