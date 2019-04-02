# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 03:18:41 2018

@author: aalco
"""

import numpy as np
from argparse import ArgumentParser

if __name__ == '__main__':
 parser = ArgumentParser(description = "Smith-Waterman Local Alignment")
 parser.add_argument("-p",dest = 'pen' , default = -2 , type = int , help = "gap penalty value")
 parser.add_argument("-b", dest = 'blosum_file', required = True , help = "File containing Blosum matrix")
 parser.add_argument("-s1", dest = 'seq1', required =True, help = 'sequence on y axis')
 parser.add_argument("-s2", dest = 'seq2', required =True, help = 'sequence on x axis')
 args = parser.parse_args()
 
def fastahandler(string):
    string_list = string.splitlines()
    if ('>' in string_list[0]) == True:    
        string_list.remove(string_list[0])
    stringline = ''.join(string_list)
    stringline.replace('\n','')
    return stringline

def blosumhandler(blosum_lines):
    mod_lines = []
    for line in blosum_lines:
        x = line.replace('\n','')
        x = x.split(' ')
        if line == blosum_lines[1]:
            blos_aa = ''.join(mod_lines[0])

        for thing in x:
            if thing == '':
                x.remove(thing)
        mod_lines.append(x)

    blos_dim = len(blos_aa)
    
    blosum_mtx = np.zeros((blos_dim,blos_dim),int)
    
    i_list = 1
    j_list = 1
    while i_list <= blos_dim:
        while j_list <= blos_dim:
            list_score = mod_lines[i_list][j_list]
            
            blosum_mtx[i_list-1][j_list-1] = list_score
            j_list += 1
        
        i_list += 1
        j_list = 1
    return blosum_mtx,blos_aa



def getscore(i,j,mtx):
    return mtx[i][j]

def matrix_init(seq1,seq2):
    dim1 = len(seq1) +1
    dim2 = len(seq2) +1
    return np.zeros((dim1,dim2))
def getblosumscore(firstA,secondA , blos_aa_order):
    i_blos = blos_aa_order.index(firstA)
    j_blos = blos_aa_order.index(secondA)
    
    return blosum_mtx[i_blos][j_blos]


def swscore(i,j,aligment_mtx,penalty,dir_list):

    up_score = getscore(i-1,j,alignment_mtx) + penalty
    left_score = getscore(i,j-1,alignment_mtx) +penalty
    diag_score = getscore(i-1,j-1,alignment_mtx)  + getblosumscore(seq1[i-1],seq2[j-1], blos_aa) 
    
    score_list = [up_score,left_score,diag_score]
    max_index = score_list.index(max(score_list))
    score = score_list[max_index]
    
    if max_index == 0:
        dir_list.append('U')
    elif max_index == 1:
        dir_list.append('L')
    elif max_index == 2:
        dir_list.append('D')
    
    if score < 0:
        return 0
    else:
        return score
def getdirection(tup):
    new_i = tup[0] -1
    new_j = tup[1] -1
    dir_index =  new_i * len(seq2) + new_j
    return dir_list[dir_index]
def nextcoord(direction, tup):
    i = tup[0]
    j = tup[1]
    if direction == 'U':
        return (i-1,j)
    elif direction == 'L':
        return (i,j-1)
    elif direction == 'D':
        return (i-1,j-1)
    return
def getAA(seq1,seq2,tup_coord,direction):
    s1_ind = tup_coord[0] -1 
    s2_ind = tup_coord[1] -1
    if direction == 'D':
        seq1_aa = seq1[s1_ind]
        seq2_aa = seq2[s2_ind]
    elif direction == 'U':
        seq2_aa = '-'
        seq1_aa = seq1[s1_ind]
    elif direction == 'L':
        seq1_aa = '-'
        seq2_aa = seq2[s2_ind]
    
    return (seq1_aa,seq2_aa)
def alignmentadjuster(seq1,seq2,coord):
    r = coord[0]
    c = coord[1]
    a = ''
    b = ''
    up_string = a + seq1
    down_string = b + seq2
    if r == 0:
        for k in range(c):
            a+= ' '
    elif c == 0:
        for k in range(r):
            b += ' '
    return up_string,down_string

def alignmentformatter(alignment_stringlist,coord):
    up_string = ''
    down_string = ''
    for tup in alignment_stringlist:
        up_string += tup[0]
        down_string += tup[1]
       
    
    r = coord[0]
    c = coord[1]
    a = ''
    b = ''
    if r == 0:
        for k in range(c):
            a+= ' '
    elif c == 0:
        for k in range(r):
            b += ' '
    up_string = a + up_string
    down_string = b + down_string
    alignment_string = 'Seq1: ' + up_string[::-1] + '\n' + 'Seq2: ' + down_string[::-1]
    print(alignment_string)
    return 


blosum_file = open(args.blosum_file, 'r')
blosum_lines = blosum_file.readlines()
blosum_file.close()

blosum_mtx,blos_aa = blosumhandler(blosum_lines)
seq1_file = open(args.seq1,'r')
seq1_text = seq1_file.read()
seq1 = fastahandler(seq1_text)
seq2_file = open(args.seq2,'r')
seq2_text = seq2_file.read()
seq2 = fastahandler(seq2_text)
seq1_file.close()
seq2_file.close()
penalty = args.pen
seq1.upper()
seq2.upper()           
    

alignment_mtx = matrix_init(seq1,seq2)
dir_list = []
i_mtx = 1
j_mtx = 1
max_list = [(0,0),0]

while i_mtx <= len(seq1):
    while j_mtx <= len(seq2):
        cell_score = swscore(i_mtx,j_mtx,alignment_mtx,penalty,dir_list)
        
        alignment_mtx[i_mtx][j_mtx] = cell_score
        if max_list[1] < cell_score:
            max_list = [(i_mtx,j_mtx),cell_score]
        j_mtx += 1
    i_mtx +=1   
    j_mtx = 1


max_dim = max(max_list[0][0],max_list[0][1])
coord = max_list[0]

alignment_stringlist = []
for k in range( max_dim, 0 ,-1):
    direction = getdirection(coord)
    AA = getAA(seq1,seq2,coord,direction)
    alignment_stringlist.append(AA)
    coord = nextcoord(direction,coord)
    if coord[0] == 0 or coord[1] == 0 :
        break
print('Score of alignment: ' + str(max_list[1]))
print('Alignment String:')
alignmentformatter(alignment_stringlist,coord)

    