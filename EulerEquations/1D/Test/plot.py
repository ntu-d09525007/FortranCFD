# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 17:32:58 2020

@author: luisl
"""
import matplotlib.pyplot as plt

files=["FineGridSolution_001.plt","CRWENO-LD_001.plt",
       "CRWENO_001.plt","OCRWENO-LD_001.plt"]

legend=["Fine Grid Solution","CRWENO-lD","CRWENO","OCRWENO-LD"]

index = 0
x, y = [], []
        
for file in files:
    with open(file,"r") as f:
        for _ in range(1):
            next(f)
        xd,yd=[],[]
        for line in f:
            values = [ float(s) for s in line.split() ]
            xd.append(values[0])
            yd.append(values[1])
        x.append(xd)
        y.append(yd)
        index+=1

for xd, yd, name in zip(x,y,legend):
    plt.plot(xd,yd,label=name)
         
plt.xlabel('x')
plt.ylabel('density')
plt.legend()
plt.show()
