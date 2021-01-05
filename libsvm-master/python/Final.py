# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 10:17:37 2021

@author: Luis
"""

import csv
from svmutil import *
import numpy as np

x_raw=[]
with open('train.csv') as f:
    contents = csv.reader(f)
    n=0
    for line in contents:
        n=n+1
        if n==1 : types=line
        if n>1 : 
            dat={types[i]:line[i] for i in range(len(types))}
            x_raw.append(dat)

x_new=[]
cancell=[]
adr=[]
for x in x_raw:
    
    xx=[]
    
    if x['is_canceled'] == '1':
        cancell.append(1)
    else:
        cancell.append(-1)
    
    adr.append(float(x['adr']))
    
    # Hotel type
    if x['hotel']=='Resort Hotel':
        xx.append(1)
    else:
        xx.append(-1)
        
    # Deposit Type
    if x['deposit_type'] == 'No Deposit':
        xx.extend([1,-1,-1])
    elif x['deposit_type'] == 'Non refund':
        xx.extend([-1,1,-1])
    else:
        xx.extend([-1,-1,1])
        
    # Company
    if x['company'] != '':
        xx.append(-1)
    else:
        xx.append(1)
        
    # Agent
    if x['agent'] != '':
        xx.append(-1)
    else:
        xx.append(1)
    
    # Market Segment
    if 'Online' in x['market_segment'] :
        xx.extend([1,-1,-1,-1])
    elif 'Offline' in x['market_segment']:
        xx.extend([-1,1,-1,-1])
    elif 'Direct' in x['market_segment']:
        xx.extend([-1,-1,1,-1])
    else:
        xx.extend([-1,-1,-1,1])
        
    #Distribution Channel
    if 'T' in x['distribution_channel']:
        xx.extend([1,-1,-1])
    elif 'Direct' in x['distribution_channel']:
        xx.extend([-1,1,-1])
    else:
        xx.extend([-1,-1,1])
        
    # Customer Type
    if x['customer_type'] == 'Transient':
        xx.extend([1,-1,-1,-1])
    elif x['customer_type'] == 'Transient-party':
        xx.extend([-1,1,-1,-1])
    elif x['customer_type'] == 'Contract':
        xx.extend([-1,-1,1,-1])
    else:
        xx.extend([-1,-1,-1,1])
    
    # Booking Time
    xx.append(float(x['lead_time'])/365)
    xx.append(float(x['arrival_date_week_number']))
    xx.append(float(x['arrival_date_day_of_month']))
    
    # Members
    xx.append(float(x['adults']))
    if x['children'] != '':
        xx.append(float(x['children']))
    else:
        xx.append(0)
    if x['children'] != '':
        xx.append(float(x['babies']))
    else:
        xx.append(0)
        
    # Meals
    if x['meal']=='BB':
        xx.extend([1,-1,-1])
    elif x['meal'] =='HB':
        xx.extend([-1,1,-1])
    else:
        xx.extend([-1,-1,1])
        
    # Room Type
    if x['reserved_room_type'] == 'A':
        xx.extend([1,-1,-1,-1])
    elif x['reserved_room_type'] == 'D':
        xx.extend([-1,1,-1,-1])
    elif x['reserved_room_type'] == 'E':
        xx.extend([-1,-1,1,-1])   
    else:
        xx.extend([-1,-1,-1,1])
        
    if x['assigned_room_type'] == 'A':
        xx.extend([1,-1,-1,-1])
    elif x['assigned_room_type'] == 'D':
        xx.extend([-1,1,-1,-1])
    elif x['assigned_room_type'] == 'E':
        xx.extend([-1,-1,1,-1])   
    else:
        xx.extend([-1,-1,-1,1])
    
    # Stay length
    xx.append(float(x['stays_in_weekend_nights']))
    xx.append(float(x['stays_in_week_nights']))
    
    # Previous Booking Informations
    if x['previous_cancellations']==1:
        xx.append(1)
    else:
        xx.append(-1)
    
    if x['previous_bookings_not_canceled']==1:
        xx.append(1)
    else:
        xx.append(-1)
        
    # Special Requests
    xx.append(float(x['total_of_special_requests']))
    xx.append(float(x['required_car_parking_spaces']))
        
    x_new.append(xx)
    
m = svm_train(cancell[0:1000], x_new[0:1000], '-s 0 -c 1 -t 2 -g 0.001')
p_label, p_acc, p_val = svm_predict(cancell[0:1000], x_new[0:1000], m, '-q')    
    
    