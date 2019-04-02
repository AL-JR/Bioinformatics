# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 11:24:21 2018

@author: aalco
"""
from argparse import ArgumentParser
from numpy import zeros, log

#Argparse
if __name__ == '__main__':
    parser = ArgumentParser(description = 'Problem Set 3: Viterbii Algorithmn')
    parser.add_argument('-s', dest = 'soluble_file', required = True , help = "Training set for soluble sequences")
    parser.add_argument('-t', dest = 'membrane_file', required = True , help = "Training set for transmembrane sequences")
    parser.add_argument('-st', dest = 'states_file', required = True , help = "Training set for state sequences")
    parser.add_argument('-x', dest = 'seq', required = True, help = "Sequence whose HMM path is to be determined")
    args = parser.parse_args()

#Methods    
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
        print(key + " Frequencies: " + str(freq_dict[key]))
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
        
        while right_reading_frame_index <= seq_string_length:
            word = seq_string[left_reading_frame_index:right_reading_frame_index]
            if word[0] == 'S':
                S_counter += 1
            elif word[0] == 'T':
                T_counter += 1
            adddict(word,word_dict)
            left_reading_frame_index+= 1
            right_reading_frame_index+=1
            
    
    for key in word_dict:
        adddict(key,freq_dict)
        if key[0] == 'S':       
            freq_dict[key] = word_dict[key]/S_counter
        elif key[0] == 'T':
            freq_dict[key] = word_dict[key]/T_counter

    
    print('Transmission Probabilities')
    for key in word_dict:
        print(key + " Frequencies: " + str(freq_dict[key]))

    
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

    
    print('Start Probabilities')
    for key in start_dict:
        print(key + " Frequencies: " + str(start_freq[key]))
    
    assert start_freq['S'] + start_freq['T'] == 1 
    return start_freq
def direction(route1,route2,max_route,row_num):
    if row_num == 0:
        if route1 == max_route: 
           
            return 'SS'
        elif route2 == max_route: 
            return 'TS'
    elif row_num == 1:
        if route1 == max_route:
            return 'TT'
        elif route2 == max_route: 
            return 'ST'
    
    return
def getnextstateandcoord(last_tup, dir_list):
    r = last_tup[0]
    c = last_tup[1]
    tup_index = (2*c + r) - 2
    transition = dir_list[tup_index]
    next_state = transition[0]
    current_state = transition[1]
    if next_state == 'S' and current_state == 'S':
        next_coord = (0,c - 1)
    elif next_state == 'S' and current_state == 'T':
        next_coord = (r-1,c-1)
    elif next_state == 'T' and current_state == 'T':
        next_coord = (r, c - 1)
    elif next_state == 'T' and current_state == 'S':
        next_coord = (r + 1, c - 1)
    return next_state,next_coord
def traceback(score_list, path_list):
    last_S = virt_mtx[0,seq_len-1]
    last_T = virt_mtx[1,seq_len-1]
    if last_S > last_T:
        path_list.append('S')
        score_list.append(last_S)
        start_coord = (0, seq_len - 1)
    else:
        path_list.append('T')
        score_list.append(last_T)
        start_coord = (1, seq_len - 1)
    next_coord = start_coord
    
    for i in range(seq_len-1,0,-1):
    
        next_state,next_coord = getnextstateandcoord(next_coord, dir_list)
    
        
        path_list.append(next_state)
        score_list.append(virt_mtx[next_coord[0],next_coord[1]])
    hmm = ''.join(path_list)
    hmm = hmm[::-1]
    return hmm, score_list
def score(score_list):
    x = 1
    for i in range(len(score_list)):
        x+= score_list[i]
    return x
def fastahandler(string):
    string_list = string.splitlines()
    if ('>' in string_list[0]) == True:    
        string_list.remove(string_list[0])
    stringline = ''.join(string_list)
    stringline.replace('\n','')
    stringline.upper()
    return stringline


#Main Code
soluble_file = open(args.soluble_file,'r')
membrane_file = open(args.membrane_file,'r')
state_file = open(args.states_file, 'r') 
seq_file = open(args.seq,'r')
seq_text = seq_file.read()
soluble_text = soluble_file.read()
membrane_text = membrane_file.read()  
state_text = state_file.readlines()
soluble_file.close() , membrane_file.close(), state_file.close(), #seq_file.close()


soluble_text = soluble_text.replace('\n','')
membrane_text = membrane_text.replace('\n','')
seq = fastahandler(seq_text)
dir_list = []
score_list = []
path_list = []
seq_len = len(seq)
virt_mtx = zeros([2 , seq_len], float)
print('Soluble Emission Probabilities')
soluble_freq = emissionp(soluble_text)
print(' ')
print('Tansmembrane Emission Probabilities')
membrane_freq = emissionp(membrane_text)
print(' ')
transmission_freq = transmissionp(state_text)
print(' ')
start_freq = startp(state_text)
print(' ')
seq_index = 0 
for amino in seq:
    
    if seq_index == 0:
        for i in range(2):
            if i == 0:
                p = log(1 * start_freq['S'] * soluble_freq[amino])
                virt_mtx[i,seq_index] = p
            elif i == 1:
                p = log(1 * start_freq['T'] * membrane_freq[amino])
                virt_mtx[i, seq_index] = p
                
    elif seq_index != 0:
        for i in range(2):
            if i == 0:
                p_SS = virt_mtx[i, seq_index - 1] + log( transmission_freq['SS'] * soluble_freq[amino])
                p_TS = virt_mtx[i + 1, seq_index - 1]  + log(transmission_freq['TS'] * soluble_freq[amino])

                soluble_max =  max(p_SS, p_TS)
                dir_list.append(direction(p_SS,p_TS,soluble_max,i))
                virt_mtx[i, seq_index] = soluble_max
                
            elif i == 1:
                p_TT = virt_mtx[i, seq_index - 1] +  log(transmission_freq['TT'] * membrane_freq[amino])

                p_ST = virt_mtx[i - 1, seq_index - 1] + log(transmission_freq['ST'] * membrane_freq[amino])


                membrane_max = max(p_TT, p_ST)

                dir_list.append(direction(p_TT,p_ST,membrane_max,i) )
                virt_mtx[i, seq_index] = membrane_max
    seq_index += 1
   
hmm_path, score_list = traceback(score_list,path_list)
path_score = score(score_list)
print('Hidden Markov Path: ')
print(hmm_path)
print(' ')
print('Path Score: ' + str(path_score))  
print(' ')
print('Scoring Matrix: ')    
print(virt_mtx)
    
    
    
    
    
    
    
    
    
    
    
    