#!/usr/bin/env python
# coding: utf-8

# In[37]:


import os
import pandas as pd
import matplotlib.pyplot as plt


# In[51]:


data=[]
cnt=0
while os.path.isfile(f'dataset_{cnt:06d}.csv'):
    with open(f'dataset_{cnt:06d}.csv') as f:
        df=pd.read_csv (f,usecols= ['Points:0','Points:2'])
        df=df.rename(columns={'Points:0': 'x', 'Points:2': 'z'})
        data.append(df)
        cnt=cnt+1
fig, ax = plt.subplots()
for i in range(0,len(data),5):
    ax.plot(data[i]['x'],data[i]['z'],'*')
plt.show()


# In[53]:


data=[]
cnt=0
while os.path.isfile(f'../CLSVOF/dataset_{cnt:06d}.csv'):
    with open(f'../CLSVOF/dataset_{cnt:06d}.csv') as f:
        df=pd.read_csv (f,usecols= ['Points:0','Points:2'])
        df=df.rename(columns={'Points:0': 'x', 'Points:2': 'z'})
        data.append(df)
        cnt=cnt+1
fig, ax = plt.subplots()
for i in range(0,len(data),5):
    ax.plot(data[i]['x'],data[i]['z'],'*')
plt.show()


# In[ ]:




