# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 10:15:48 2018

@author: aalco
"""

s = open('contig45_predictedgene_takahashii.txt','r')
seq = s.read()

def fastahandler(string):
    string_list = string.splitlines()
    if ('>' in string_list[0]) == True:    
        string_list.remove(string_list[0])
    stringline = ''.join(string_list)
    stringline.replace('\n','')
    return stringline
x = fastahandler(seq)
if '\n' in x :
    print('true')
