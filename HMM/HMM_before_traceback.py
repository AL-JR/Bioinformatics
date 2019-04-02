# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 16:16:07 2018

@author: aalco
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 11:24:21 2018

@author: aalco
"""
from argparse import ArgumentParser
from numpy import zeros
from decimal import *

#if __name__ == '__main__':
#    parser = ArgumentParser(description = 'Problem Set 3: Viternii Algorithmn')
#    parser.add_argument('-s', dest = 'soluble_file', required = True , help = "Training set for soluble sequences")
#    parser.add_argument('-t', dest = 'membrane_file', required = True , help = "Training set for transmembrane sequences")
#    parser.add_argument('-st', dest = 'states_file', required = True , help = "Training set for state sequences")
#    args = parser.parse_args()
#    
def adddict(word,seqdict):
    if word in seqdict:
        seqdict[word] += 1
    elif word not in seqdict:
        seqdict[word] = 1
    return
def emissionp(seq_string):
    word_length = 1
    seq_string_length = len(seq_string)
    left_reading_frame_index = 0
    right_reading_frame_index = left_reading_frame_index + word_length
    word_dict = {}
    freq_dict = {}
    word_counter = 0.0
    
    while right_reading_frame_index <= seq_string_length:
        word = seq_string[left_reading_frame_index:right_reading_frame_index]
    
        adddict(word,word_dict)
    
        left_reading_frame_index+= 1
        right_reading_frame_index+=1
        word_counter += 1
    for key in word_dict:
        adddict(key,freq_dict)
        freq_dict[key] = word_dict[key]/word_counter
    for key in word_dict:
        print(key + " Frequencies: " + str(word_dict[key]/word_counter))
    print("Word Length: " + str(word_length))
    print(sum(freq_dict.values()))
    print('*****^^^^^^^^^^^')

    return freq_dict
def transmissionp(state_list):
    word_dict = {}
    freq_dict = {}
    S_counter = 0.0
    T_counter = 0.0
    for seq in state_list:
        
        seq_string = seq
        seq_string = seq_string.replace('\n','')
        seq_string_length = len(seq_string)
        left_reading_frame_index = 0
        right_reading_frame_index = 2
#        word_dict = {}
#        word_counter = 0.0
        
        while right_reading_frame_index <= seq_string_length:
            word = seq_string[left_reading_frame_index:right_reading_frame_index]
            if word[0] == 'S':
                S_counter += 1
            elif word[0] == 'T':
                T_counter += 1
            adddict(word,word_dict)
            left_reading_frame_index+= 1
            right_reading_frame_index+=1
#            word_counter += 1
            
    
    for key in word_dict:
        adddict(key,freq_dict)
        if key[0] == 'S':       
            freq_dict[key] = word_dict[key]/S_counter
        elif key[0] == 'T':
            freq_dict[key] = word_dict[key]/T_counter


    for key in word_dict:
        print(key + " Frequencies: " + str(freq_dict[key]))
#        print(key + " Count: " + str(word_dict[key]))
    print("Word Length: " + '2')
    
    
    assert (freq_dict['SS'] + freq_dict['ST'] == 1) and (freq_dict['TT'] + freq_dict['TS']) == 1
    return freq_dict


def startp(state_list):
    start_dict = {'S':0,'T':0}
    start_freq = {}
    for line in state_list:
        if line[0] == 'S':
            start_dict['S']+= 1
        if line[0] == 'T':
            
            start_dict['T'] += 1
    for key in start_dict:
        adddict(key,start_freq)
        start_freq[key] = start_dict[key]/len(state_list)


    
    assert start_freq['S'] + start_freq['T'] == 1 
    return start_freq
def direction(route1,route2,max_route,row_num):
    if row_num == 0:
        if route1 == max_route: return 'SS'
        elif route2 == max_route: return 'TS'
    elif row_num == 1:
        if route1 == max_route:
            return 'TT'
        elif route2 == max_route: 
            return 'ST'
    
    return
def state_step(soluble_max,membrane_max):
    if soluble_max > membrane_max:
        return 'S'
    else:
        return 'T'
#########################

soluble_file = open('soluble_sequences.txt','r')
membrane_file = open('transmembrane_sequences.txt','r')
state_file = open('state_sequences.txt', 'r') 
soluble_text = soluble_file.read()
membrane_text = membrane_file.read()  
state_text = state_file.readlines()
soluble_file.close() , membrane_file.close(), state_file.close()


soluble_text = soluble_text.replace('\n','')
membrane_text = membrane_text.replace('\n','')

soluble_freq = emissionp(soluble_text)
membrane_freq = emissionp(membrane_text)
transmission_freq = transmissionp(state_text)
start_freq = startp(state_text)
dir_list = []
path_list = []
seq = 'STFQFCAIPQVMAIATLALVF'
seq_len = len(seq)
virt_mtx = zeros([2 , seq_len], float)
seq_index = 0 
for amino in seq:
    
    print(seq_index)
    if seq_index == 0:
        for i in range(2):
            if i == 0:
                p = 1 * start_freq['S'] * soluble_freq[amino]
                virt_mtx[i,seq_index] = p
            elif i == 1:
                p = 1 * start_freq['T'] * membrane_freq[amino]
                virt_mtx[i, seq_index] = p
                
    elif seq_index != 0:
        for i in range(2):
            if i == 0:
                p_SS = virt_mtx[i, seq_index - 1] * transmission_freq['SS'] * soluble_freq[amino]
                p_TS = virt_mtx[i + 1, seq_index - 1] * transmission_freq['TS'] * soluble_freq[amino]
                soluble_max =  max(p_SS, p_TS)
                dir_list.append(direction(p_SS,p_TS,p_max,i))
                virt_mtx[i, seq_index] = soluble_max
                
            elif i == 1:
                p_TT = virt_mtx[i, seq_index - 1] * transmission_freq['TT'] * membrane_freq[amino]
                p_ST = virt_mtx[i - 1, seq_index - 1] * transmission_freq['ST'] * membrane_freq[amino]
                membrane_max = max(p_TT, p_ST)

                dir_list.append(direction(p_TT,p_ST,p_max,i) )
                virt_mtx[i, seq_index] = membrane_max
        path_list.append(state_step(soluble_max,membrane_max))
    seq_index += 1
    
    
    
hmm = ''.join(path_list)
print(hmm)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    