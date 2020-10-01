# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 20:55:56 2020

@author: Luis
"""
import os

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)

dirname='pov_new'
createFolder('./'+dirname+'/')

for i in range(1000):
    
    flag1=flag2=flag3=0
    try:        
        with open(str(i)+".pov","r") as fp: content=fp.readlines()
        fp.close()
        fn = open('./'+dirname+'/'+str(i)+".pov","w")
        fn.write('#include "my setting.inc" \n\n')
    
        for line in content:
            if flag1==0 and "mesh2 {" in line: flag1=1
            if flag1==1 and "matrix" in line: flag2=1
            if flag1*flag2==1 and ">" in line: flag3=1
            if flag1==1:fn.write(line)
            if flag3==1:
                fn.write("material {milk}")
                fn.write("}")
                break
        fn.close()
    except:
        break

            
                
        
                    
                
