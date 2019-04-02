# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 11:24:21 2018

@author: aalco
"""


def adddict(word,seqdict):
    if word in seqdict:
        seqdict[word] += 1
    elif word not in seqdict:
        seqdict[word] = 1
    return
def missionp(seq_string,word_length):
    seq_string_length = len(seq_string)
    left_reading_frame_index = 0
    right_reading_frame_index = left_reading_frame_index + word_length
    word_dict = {}
    word_counter = 0.0
    
    while right_reading_frame_index <= seq_string_length:
        word = seq_string[left_reading_frame_index:right_reading_frame_index]
    
        adddict(word,word_dict)
    
        left_reading_frame_index+= 1
        right_reading_frame_index+=1
        word_counter += 1
    
    for key in word_dict:
        print(key + " Frequencies: " + str(word_dict[key]/word_counter))
    print("Word Length: " + str(word_length))
    return word_dict

soluble_file = open('soluble_sequences.txt','r')
membrane_file = open('transmembrane_sequences.txt','r')
state_file = open('state_sequences.txt', 'r')
soluble_text = soluble_file.read()
membrane_text = membrane_file.read()  
state_text = state_file.read()
soluble_file.close() , membrane_file.close(), state_file.close()

soluble_text = soluble_text.replace('\n','')
membrane_text = membrane_text.replace('\n','')
state_text = state_text.replace('\n','')

soluble_dict = missionp(soluble_text, 1)
membrane_dict = missionp(membrane_text, 1)
state_dict = missionp(state_text,2)

