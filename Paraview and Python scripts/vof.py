# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import math
import random

def find_alpha(V,M1,M2,M3):
    
    m1 = min(M1,M2,M3)
    m3 = max(M1,M2,M3)
    m2 = M1+M2+M3-m1-m3
    
    mag=math.sqrt(m1*m1+m2*m2+m3*m3)
    
    m1=m1/mag
    m2=m2/mag
    m3=m3/mag
    
    eps=1.0e-12
    
    m12=m1+m2
    m=min(m12,m3)
    
    V1=m1**2/max(6.0*m2*m3,eps)
    V2=V1+(m2-m1)/(2.0*m3)
    V31=(m3**2.0*(3.0*m12-m3)+m1**2.0*(m1-3.0*m3)+m2**2.0*(m2-3.0*m3))/(6.0*m1*m2*m3)
    V32=m12/(2.0*m3)
    
    if abs(m-m3)<eps:
        V3=V31
    else:
        V3=V32

    
    if V1 > V:
        alpha=(6.0*m1*m2*m3*V)**(1.0/3.0)
    elif V1<=V<V2:
        alpha=0.5*(m1+math.sqrt(m1**2.0+8.0*m2*m3*(V-V1)))
    elif  V2<=V<V3:
        alpha = cubic_root(-1.0, 3.0*m12, -3.0*(m1**2+m2**2), m1**3+m2**3-6*m1*m2*m3*V)
    else:
        if abs(V3-V31)<eps:
            alpha = cubic_root(-2.0, 3.0, -3*(m1**2+m2**2+m3**2), m1**3+m2**3+m3**3-6*m1*m2*m3*V)
        else:
            alpha = m3*V + m12/2
            
    return alpha

def cubic_root(a3,a2,a1,a0):
    a2=a2/a3
    a1=a1/a3
    a0=a0/a3
    
    p=a1/3.0-a2**2.0/9.0
    q=(a1*a2-3.0*a0)/6.0-a2**3.0/27.0
    
    theta = math.acos(q/math.sqrt(-p**3.0))/3.0
    
    return math.sqrt(-p)*(math.sqrt(3.0)*math.sin(theta)-math.cos(theta))-a2/3.0

def get_V(alpha,M1,M2,M3):
    
    m1 = min(M1,M2,M3)
    m3 = max(M1,M2,M3)
    m2 = M1+M2+M3-m1-m3
    
    mag=math.sqrt(m1*m1+m2*m2+m3*m3)
    
    m1=m1/mag
    m2=m2/mag
    m3=m3/mag
    
    eps=1.0e-12

    m12=m1+m2
    m=min(m12,m3)
    
    V1=m1**2/max(6.0*m2*m3,eps)

    if alpha < m1:
        V = alpha**3 / (6*m1*m2*m3)
    elif m1 <= alpha < m2:
        V = alpha*(alpha-m1)/(2*m2*m3) + V1
    elif m2 <= alpha < m:
        V = alpha**2*(3*m12-alpha)+m1**2*(m1-3*alpha)+m2**2*(m2-3*alpha)
        V = V / (6*m1*m2*m3)
    else:
        if abs(m-m3)<eps:
            V = alpha**2*(3-2*alpha)+m1**2*(m1-3*alpha)+m2**2*(m2-3*alpha)+m3**2*(m3-3*alpha)
            V = V / (6*m1*m2*m3)
        else:
            V = (2*alpha-m12)/(2*m3)
    
    return V

for i in range(1000000):
    alpha=random.random()
    if alpha > 0.5 : alpha = 1 - alpha
    m1=random.random()
    m2=random.random()
    m3=random.random()
    V=get_V(alpha,m1,m2,m3)
    try:
        alpha2=find_alpha(V,m1,m2,m3)
    except:
        print("Error at:",V,alpha)
    # print(V,alpha,alpha2)
    if abs(alpha-alpha2)>1.0e-12:
        print("Error at:",V,alpha,abs(alpha-alpha2))
        
    