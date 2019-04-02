# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 15:42:12 2018

@author: aalco
"""
from argparse import ArgumentParser

#Argparse
#if __name__ == '__main__':
#    parser = ArgumentParser(description = 'Problem Set 4: Clustering')
#    parser.add_argument('-p', dest = 'protein_file', required = True , help = "file contaning protein sequences and their respective labels")
#    args = parser.parse_args()
#    
#    
#    
#    
def protdict():
    standard = 'ARNDCQEGHILKMFPSTWYV'
    protd = {}
    for acid in standard:
        protd[acid] = 0
    return protd
def adddict(word,seqdict):
    if word in seqdict:
        seqdict[word] += 1
    elif word not in seqdict:
        seqdict[word] = 1
    return
def aafreq(string,protd):
    aaseq = string.replace('\n','')
    word_length = 1
    seq_string_length = len(aaseq)
    left_reading_frame_index = 0
    right_reading_frame_index = left_reading_frame_index + word_length
    
    freq_dict = {}
    word_counter = 0.0
    
    while right_reading_frame_index <= seq_string_length:
        word = aaseq[left_reading_frame_index:right_reading_frame_index]

        protd[word]+= 1
    
        left_reading_frame_index+= 1
        right_reading_frame_index+=1
        word_counter += 1
    for key in protd:
        adddict(key,freq_dict)
        freq_dict[key] = round(protd[key]/word_counter,3)

    return freq_dict
def vecstring(label,freqs):
    vec_str = label 
    for key in aa_order:
        if key != 'V':
            aa_freq = '\t' + str(freqs[key]) 
        else: 
            aa_freq = '\t' + str(freqs[key])  + '\n'
        
        vec_str += aa_freq
    return vec_str
            
def veclist(key_list,protd_list):
    vector_list = []
    i = 0
    while i < len(key_list):
        label = key_list[i]
        freqs = protd_list[i]
        vector = vecstring(label,freqs)
        
        vector_list.append(vector)
        i+=1
    return vector_list

def morpheustext(vlist):
    morpheus_txt = 'PROTEIN\tA\tR\tN\tD\tC\tQ\tE\tG\tH\tI\tL\tK\tM\tF\tP\tS\tT\tW\tY\tV\n'
    i = 0
    while i < len(vlist):
        morpheus_txt += vlist[i]
        
        i +=1

    return morpheus_txt


#prot_file = open(args.protein_file, 'r')
prot_file = open('test.txt','r')
prot_text = prot_file.read()
prot_file.close()
aa_order = 'ARNDCQEGHILKMFPSTWYV'
prot_chunks = prot_text.split('>')
prot_chunks.remove('')

protd_list = []
key_list = []
i = 0
while i < len(prot_chunks):
    prot_lines = prot_chunks[i].splitlines()
    prot_chunks[i] = prot_lines
    plabel = prot_chunks[i][0]
    key_list.append(plabel)
    i+=1
i = 0
j = 1
temp = ''
while i < len(key_list):
    while j < len(prot_chunks[i]):
        temp+= prot_chunks[i][j]
        j+=1
        
    protd = protdict()
    prot_freq = aafreq(temp,protd)
    protd_list.append(prot_freq)
    temp =''
    i+=1
    j = 1
    
vlist = veclist(key_list,protd_list)



x = morpheustext(vlist)
morph_txt = open('Morpheus_text.txt','w')
morph_txt.write(x)
morph_txt.close()

print('Process Completed')






