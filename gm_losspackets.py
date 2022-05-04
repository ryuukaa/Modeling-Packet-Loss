#!/usr/bin/env python
# -*- coding: utf-8 -*-
##
# @file     sgm_losspackets.py
# @author   Ryuuka <Jia.Liu21@student.xjtlu.edu.cn>
# @date     2022-04-12
#
# @brief    A function for calculating transmission probabilities based on sample sequence and
#           then generating new loss sequence based on Guilbert model and getting loss probability.
#
import numpy as np
import scipy.io as scio
import sys

import sgm_losspackets
from binary_runlengths import binary_runlengths
from sgm_generate import sgm_generate
from matplotlib import pyplot as plt

from visualization import visualization


def gm(seq):
    '''
    p is the probability of transferring from good State to the bad state
    and q is the probability of transferring from the bad state to the good state
    k = 1 called Guilbert model, at good state, the loss probability = 0
    P = (1-p, p
         q, 1-q)
    '''
    a = np.count_nonzero(seq)/len(seq)
    win = 0
    p_11 = p_111 = p_101 = 0
    for win in range(len(seq) - 1):
        if(seq[win]==1 and seq[win+1]==1):
            p_11 = p_11 + 1
        win = win + 1
    win= 0
    for win in range(len(seq) - 2):
        if(seq[win]==1 and seq[win+1]==0 and seq[win+2]==1):
            p_101 = p_101 + 1
        if(seq[win]==1 and seq[win+1]==1 and seq[win+2]==1 ):
            p_111 = p_111 + 1
        win = win + 1
    b = p_11/np.count_nonzero(seq)
    c = p_111/(p_101+p_111)
    q = 1 - (a*c-b*b)/(2*a*c-b*(a+c))
    h = 1 - b/(1-q)
    p = (a*q) / (1-h-a)
    return p,q,h
if __name__ == "__main__":
    dataFile = '/Users/liujia/PycharmProjects/can406_lab1/loss_seq.mat'
    data = scio.loadmat(dataFile)#get sample sequence from .mat file
    sample_seq = data['x'][0]
    p,q,h = gm(sample_seq)
    print("p = ", p)
    print("q = ", q)
    print("h = ", h)
    # according to sample sequence, calculate p, q and h by gm
    # then based on it, generated sequence
    tr = np.array([[1-p, p],
                   [q, 1-q]])
    new_seq = sgm_generate(10000,tr)
    p_,q_ ,h_= gm(new_seq)
    count = np.count_nonzero(new_seq)
    L = int(1/q_)
    theo_loss_p = (1-h)*p_/(p_+q_)
    act_loss_p = count/len(new_seq)
    print(new_seq)
    print('The theoretical packet loss rate is: ', format(theo_loss_p,'.4f'))
    print('The actual packet loss rate is: ', act_loss_p)
    print("the mean burst length is: ",L)
    visualization(new_seq, sample_seq)
    plt.psd(new_seq, NFFT=256, Fs=2)
    plt.show()