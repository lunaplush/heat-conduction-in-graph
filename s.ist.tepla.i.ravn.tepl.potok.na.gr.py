from math import *

lambda_=1
c=1
ro=1
ksi=1

Tl=1790             #temperatura liquidus
Ts=1740             #temperatura solidus
#Ts=1790             #temperatura solidus
Tkr=1773

n=12
delta_x=0.1/n
temp=0              #temperatura
t=[]                #temperatura tabliza
life=5

delta_tau=5     # в секундах
# условие устойчивости delta_tau <= c_min * ro_min * delta_x**2 / (2 * lam_max)
delta_tau_max = 592 * 660 * delta_x**2 / (2 * 35)     #  = 0.5581714285714287

table_heatcapacity=[[975,660],[1031,635],[1032,610],[1079,592],[1727,696],[1757,720],[1758,743],[1773,748],[1797,824],[1873,846]]
table_ro=[[975,7700],[1079,7500],[1727,7330],[1757,7280],[1758,7265],[1773,7230],[1797,7020],[1873,6600]]  # плотность
table_lambda=[[975,33],[1079,28],[1757,33],[1797,34.8],[1873,34.8]]    # теплопроводность
table_initial_T=[[0,1823],[0.05,1820],[0.07,1812],[0.09,Tl],[0.095,Ts],[0.1,1700]]

def linealinterpol(x1,x2,y1,y2,x):
    return ((y2-y1)/(x2-x1))*(x-x1)+y1    

def table_interpol(t,table):
    i=0
    while t>table[i][0]:
        i=i+1
    x1=table[i-1][0]
    y1=table[i-1][1]
    x2=table[i][0]
    y2=table[i][1]
    return linealinterpol(x1,x2,y1,y2,t)

def c(temp):
    return table_interpol(temp,table_heatcapacity)

def lam(temp):
    return table_interpol(temp,table_lambda)

def ro(temp):
    return table_interpol(temp,table_ro)


# решаем уравнение теплопроводности
# задаем начальные условия
x = [i*delta_x for i in range(0,n+1)] # узлы сетки
#def f(x):
#    return 50000000*(x**2-0.06*x+0.0008)*(x**2-0.14*x+0.0048)+1500
def T0(x):
    return table_interpol(x,table_initial_T)
t=[[T0(x[i]) for i in range(0,n+1)]]

def next_t(t_im1,t_i,t_ip1):
    # bez istocnika
    #r=(delta_tau/(c(t_i)*ro(t_i)*(delta_x**2)))*((lam(t_i)-lam(t_im1))*(t_i-t_im1)+lam(t_i)*(t_ip1-2*t_i+t_im1))+t_i

    #s istochnikom
    u = 0               #isto4nik / stok  tepla
    r=delta_tau*((lam(t_im1)*t_im1 - lam(t_im1)*t_i - lam(t_i)*t_i + lam(t_i)*t_ip1)/(c(t_i)*ro(t_i)*(delta_x**2)) + u) + t_i
    return r


for j in range(0,life+1):
    q = 1000000*(1.12*exp(-3.12*j*delta_tau)+0.503)
    q =100000
    t.append([0 for i in range(0,n+1)])
    t[j+1][0]=t[j][0]   
    t[j+1][n]=t[j][n]
    for i in range(1,n): t[j+1][i]=next_t(t[j][i-1],t[j][i],t[j][i+1]) 

    #t[j+1][0]= 47500.0*delta_x/lam(t[j][0]) + t[j+1][1]   #gran usl II-go roda

#gran usl II-go roda, teplovoj potok =0
    #t[j+1][0]= t[j+1][1]   
#gran usl II-go roda, teplovoj potok =q
#    q = (1.12*exp(-3.13*i*delta_tau)+0.503)*(10**6)

    t[j+1][0]= q*delta_x/lam(t[j][0]) + t[j+1][1]
    t[j+1][n]= -q*delta_x/lam(t[j][n]) + t[j+1][n-1]
    
from numpy import *
from math import *
import matplotlib.pyplot as plt

def table2d_interpol(table,stroka,tmp):
    x1=floor(tmp)
    if x1==tmp:
        return table[stroka][x1]
    x2=x1+1
    y1=table[stroka][x1]
    y2=table[stroka][x2]
    return linealinterpol(x1,x2,y1,y2,tmp)

def g(tmp):
    return table2d_interpol(t,i,tmp)


for i in range(0,life+1):
    tmp = linspace(0, len(t[i])-1, 33)  # на отрезке 0 и 0.1 100 точек 
    y=[]
    for x in tmp:
        y.append(g(x))
    plt.plot(tmp*delta_x, y)

plt.show()



