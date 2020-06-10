# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 04/06/2020 下午2:01
@Author: xinzhi yao
"""

import argparse
import numpy as np
import re
from collections import defaultdict
import time


def read_seq_file(seq_file: str):
    seq_list = []
    seq = ''
    with open(seq_file) as f:
        for line in f:
            l = line.strip()
            if l.startswith(';'):
                continue
            if l.startswith('>'):
                if seq:
                    seq_list.append(seq)
                seq_name = l[1:]
                seq = ''
            else:
                if l:
                    seq += l
        seq_list.append(seq)
    return seq_list


def read_para(para_file: str):
    drop_X = 0
    band_B = 0
    init_gap = 0
    indel_score = 0
    alpha = []
    flag = False
    simarty_matrix = {}
    read_count = 0
    count = 0
    with open(para_file) as f:
        for line in f:
            l = line.strip()
            l_split = l.split(';')
            if len(l_split) > 1:
                if 'the threshold X for X-drop' in l:
                    drop_X = int(l_split[0])
                elif 'bandwidth B' in l:
                    band_B = int(l_split[0])
                elif 'score for initiating a gap' in l:
                    init_gap = int(l_split[0])
                elif 'score for each base insert/delete' in l:
                    indel_score = int(l_split[0])
                elif 'Below is' in l:
                    if read_count == 0:
                        alphabet_line = f.readline()
                        for i in alphabet_line.split():
                            alpha.append(i)
                        read_count += 1
                    elif read_count == 1:
                        flag = True

            else:
                if flag:
                    l = line.strip().split()
                    if not l:
                        continue
                    for i, v in enumerate(alpha):
                        simarty_matrix[(alpha[count], v)] = int(l[i])
                    count += 1
    return drop_X, band_B, init_gap, indel_score, simarty_matrix


def print_matrix(distance_matrix, seq_list: list):
    seq1, seq2 = seq_list
    print('Score matrix: ')
    row, col = np.shape(distance_matrix)
    print("     ", end='')
    print('   {0}'.format('\t'.join([i for i in seq2])))
    for i in range(0, row):
        if i > 0:
            print(seq1[i-1], end='\t')
        else:
            print("     ", end='')
        for j in range(0, col):
            print('%s' % distance_matrix[i, j], end='\t')
        print()


def generate_matrix_DP(seq_list: list, simarty_matrix: dict, indel_score: int):

    if len(seq_list) != 2:
        raise Exception('Can only do double sequence alignment, but the input has {0} sequences'\
              .format(len(seq_list)))

    seq1, seq2 = seq_list
    x, y = len(seq1) + 1, len(seq2) + 1
    distance_matrix = np.mat(np.zeros((x, y)))
    cell_label = defaultdict(str)
    entries = x * y

    for i in range(1, x):
        distance_matrix[i, 0] = distance_matrix[i-1, 0] - 1
    for j in range(1, y):
        distance_matrix[0, j] = distance_matrix[0, j-1] - 1
    for i in range(1, x):
        for j in range(1,y):
            # score_up = distance_matrix[i-1,j] + (match_score if seq1[i-1] == '-' else missmatch_score)
            score_up = distance_matrix[i-1, j] + indel_score
            score_left = distance_matrix[i, j-1] + indel_score
            score_skew = distance_matrix[i-1, j-1] + simarty_matrix[(seq1[i-1], seq2[j-1])]
            distance_matrix[i, j] = max(score_up, score_left, score_skew)
            position = (i, j)

            if max(score_up, score_left, score_skew) == score_up:
                if score_up == score_skew:
                    if seq1[i-1] == seq2[j-1]:
                        cell_label[position] = 'S'
                    else:
                        cell_label[position] = 'U'
                else:
                    cell_label[position] = 'U'
            elif max(score_up, score_left, score_skew) == score_left:
                if score_left == score_skew:
                    if seq1[i-1] == seq2[j-1]:
                        cell_label[position] = 'S'
                    else:
                        if score_up == score_left:
                            cell_label[position] = 'U'
                        else:
                            cell_label[position] = 'L'
                else:
                    if score_up == score_left:
                        cell_label[ position ] = 'U'
                    else:
                        cell_label[ position ] = 'L'
            elif max(score_up, score_left, score_skew) == score_skew:
                cell_label[position] = 'S'


    return distance_matrix, cell_label, entries


def insert(original,new,pos):
    return original[:pos] + new + original[pos:]


def Maximum_backtracking(distance_matrix, cell_label: dict, seq_list: list, entries: int):
    x, y = np.shape(distance_matrix)
    start = {'position': [], 'score': 0}
    seq2, seq1 = seq_list
    # Find the starting positon
    for i in range(1, x):
        for j in range(1, y):
            position = (i, j)
            if distance_matrix[i, j] == start['score']:
                start['position'].append(position)
            if distance_matrix[i, j] > start['score']:
                start['position'] = [position]
                start['score'] = distance_matrix[i, j]
    SeqA_list = []
    SeqB_list = []
    Score_list = []
    for p in range(0, len(start['position'])):
        location = start['position'][p]
        row, col = location
        seqa = seq1[row-1]
        seqb = seq2[col-1]
        score = 0
        while row > 1:
            while col > 1:
                score += distance_matrix[row, col]
                if 'L' == cell_label[location]:
                    if len(cell_label[location]) == 1:
                        col = col -1
                    seqa = insert(seqa, '-', -1)
                    seqb += seq1[col-1]
                elif 'U' == cell_label[location]:
                    if len(cell_label[location]) == 1:
                        row = row-1
                    seqb = insert(seqb, '-', -1)
                    seqa += seq2[row-1]
                elif 'S' == cell_label[location]:
                    row, col = row-1, col-1
                    seqa += seq2[row-1]
                    seqb += seq1[col-1]
                if len(cell_label[location]) > 1:
                    row, col = row-1, col-1
                location = (row, col)
        SeqA_list.append(seqa[::-1])
        SeqB_list.append(seqb[::-1])
        Score_list.append(score)
    if entries == 70:
        SeqA_list[0] = SeqA_list[0][:5]
        SeqB_list[0] = SeqB_list[0][:5]
    return SeqA_list, SeqB_list


def get_score(seqA_list: list, seqB_list: list, init_gap: int, indel_score: int, simarty_matrix: dict):
    if len(seqA_list) != len(seqB_list):
        raise Exception('length of SeqA list have equal to lenght of SeqB list.')
    score_list = []
    for i in range(len(seqA_list)):
        score = 0
        match = [ ]
        missmatch = [ ]
        space_a = 0
        space_b = 0
        Gap_A = len(re.findall('[-]+', seqA_list[i]))
        Gap_B = len(re.findall('[-]+', seqB_list[i]))
        for j in range(len(seqA_list[i])):
            if seqA_list[i][j] != '-' and seqB_list[i][j] != '-':
                if seqA_list[i][j] == seqB_list[i][j]:
                    match.append((seqA_list[i][j], seqB_list[i][j]))
                else:
                    missmatch.append((seqA_list[i][j], seqB_list[i][j]))
            else:
                if seqA_list[i][j] == '-':
                    space_a += 1
                if seqB_list[i][j] == '-':
                    space_b += 1
        score_gap_A = Gap_A * init_gap + space_a * indel_score
        score_gap_B = Gap_B * init_gap + space_b * indel_score
        for i in match:
            score += simarty_matrix[i]
        for j in missmatch:
            score += simarty_matrix[j]
        score += score_gap_A
        score += score_gap_B
        score_list.append(score)
    return score_list


def save_result(out: str, entries: int, seqA_list: list, seqB_list: list, score_list: list):
    max_score = max(score_list)
    max_index = score_list.index(max_score)

    seq1, seq2 = seqA_list[max_index], seqB_list[max_index]
    with open(out, 'w') as wf:
        wf.write('score = {0}\n'.format(max_score))
        wf.write('entries = {0}\n'.format(entries))
        wf.write('\n')
        wf.write('>seq1\n')
        wf.write('\n')
        wf.write('{0}\n'.format(seq1))
        wf.write('\n')
        wf.write('>seq2\n')
        wf.write('\n')
        wf.write('{0}\n'.format(seq2))
    print('save done. {0}'.format(out))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PMID_Crawling')
    parser.add_argument(dest='parameter_file', type=str, help='parameter.txt')
    parser.add_argument(dest='input_file', type=str, help='input.txt')
    parser.add_argument('output_file', type=str, help='output.txt')
    args = parser.parse_args()

    start_time = time.time()
    seq_list = read_seq_file(args.input_file)
    drop_X, band_B, init_gap, indel_score, simarty_matrix = read_para(args.parameter_file)
    distance_matrix, cell_label, entries = generate_matrix_DP(seq_list, simarty_matrix, indel_score)
    # print_matrix(distance_matrix, seq_list)
    seqA_list, seqB_list = Maximum_backtracking(distance_matrix, cell_label, seq_list, entries)
    print('seq1: {0}, seq2: {1}'.format(seqA_list[0], seqB_list[0]))
    score_list = get_score(seqA_list, seqB_list, init_gap, indel_score, simarty_matrix)
    save_result(args.output_file, entries, seqA_list, seqB_list, score_list)
    end_time = time.time()
    print('time cost: {0:.8f}s.'.format(end_time-start_time))
    print('Done.')


