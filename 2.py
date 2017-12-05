# coding=utf-8

'''
作者:Li Jianghua 学号：3117001132 班级：机械7003
程序:多项式曲线拟合和指数形式最小二乘算法
'''
import matplotlib.pyplot as plt
import numpy as np
import math

#正交化方法求解最小二乘法
def LSS(G):
    m, n = G.shape
    n=n-1
    summ = 0
    w=np.zeros(m)
    for k in np.arange(n):
        summ = 0
        for j in np.arange(k, m):
            summ =summ+(G[j][k])**2
        sigema = -np.sign(G[k][k])*math.sqrt(summ)
        w[k]=G[k][k]-sigema
        for j in np.arange(k+1,m):
            w[j]=G[j][k]
        berta = sigema * w[k]
        G[k][k]=sigema
        for j in np.arange(k+1,n+1):
            summ=0
            for i in np.arange(k,m):
                summ = summ + w[i]*G[i][j]
            t=summ/berta
            for i in np.arange(k,m):
                G[i][j]=G[i][j]+t*w[i]

    x=np.zeros(n)

    x[n-1]=G[n-1][n]/G[n-1][n-1]
    for i in np.arange(n-1):
        i = n-i-2
        summ =0
        for j in np.arange(i,n):
            summ=summ+G[i][j]*x[j]
        x[i]=(G[i][n]-summ)/G[i][i]
    rho=0
    for i in np.arange(n,m):

        rho=rho +G[i][n]**2
    return x,rho

def px(a,n,_x):
    sum=0
    for i in np.arange(n):
        sum=sum+a[i]*pow(_x,i)

    return sum

# 画出拟合后的曲线
if __name__ == '__main__':
    fig = plt.figure()
    axes = fig.add_subplot(111)

    # 生成曲线上的各个点
    x = np.arange(25)
    y = [15,14,14,14,14,15,16,18,20,20,23,25,28,31,34,31,29,27,25,24,22,20,18,17,16]


    yy=np.array(y)
    N=len(x)
    order = input("请输入插值多项式阶数：")
    G=np.zeros([N,order+1])
    for i in np.arange(N):
        #G[i][order]=y[i]#多项式形式的近似函数
        G[i][order] = np.log(y[i])#自然对数e 指数形式的近似函数
    for j in range(order):
        for i in range(N):
            G[i][j]= pow(x[i],j)
    a, e = LSS(G)
    #插值点
    _x=np.arange(0,24,0.01)
    _y=px(a,order,_x)
    _y=np.exp(_y)
    axes.plot(_x, _y, color='r', linestyle='-', marker='*',label='data fitting')
    axes.plot(x, np.array(y),  color='k', linestyle='', marker='*',label='data source')
    axes.plot(x, np.zeros_like(x), color='g', linestyle='--', marker='*',label='water level')
    axes.legend()
    plt.title("the temperature")
    plt.show()
    print "这一天的平均气温为：", np.mean(y)
