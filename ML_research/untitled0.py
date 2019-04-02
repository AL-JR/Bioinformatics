# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 13:41:18 2018

@author: aalco
"""

def solution(S, K):
    # write your code in Python 3.6
    mod_S = ''.join(S.upper().split('-'))
    mod_S_len = len(mod_S)
    new_S = ''

    if mod_S_len > K:
        
        right = K
        left = 0
        
        while right <= mod_S_len:
            
            string_chunk = mod_S[left:right]
            
            
            if right >= mod_S_len:
                new_S += mod_S[left:mod_S_len]
            else:
                new_S += string_chunk + '-'
            right  += K
            left += K
        new_S += mod_S[left:mod_S_len]
    else:
        new_S = mod_S
        
    return new_S

s = '9jhV'

print(solution(s,1))