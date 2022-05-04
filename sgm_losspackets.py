#!/usr/bin/env python
# -*- coding: utf-8 -*-
##
# @file     sgm_losspackets.py
# @author   Ryuuka <Jia.Liu21@student.xjtlu.edu.cn>
# @date     2022-04-12
#
# @brief    A function for calculating transmission probabilities based on sample sequence and
#           then generating new loss sequence based on simple Guilbert model and getting loss probability.
#
import numpy as np
import scipy.io as scio
import sys
from binary_runlengths import binary_runlengths
from sgm_generate import sgm_generate
from matplotlib import pyplot as plt

from visualization import visualization


def sgm(seq):
    '''
    p is the probability of transferring from good State to the bad state
    and q is the probability of transferring from the bad state to the good state
    k = 1, h = 0, at good state, the loss probability = 0, at bad state, loss probability = 1
    P = (1-p, p
         q, 1-q)
    '''
    win = 0
    n_01 = n_10 = n_0 = n_1 = 0
    for win in range(len(seq) - 1):
        if(seq[win]==0 and seq[win+1]==1):
            n_01 = n_01 + 1
        if(seq[win]==1 and seq[win+1]==0):
            n_10 = n_10 + 1
        win = win + 1
    n_1 = np.count_nonzero(seq)
    n_0 = len(seq) - n_1
    p = n_01/n_0
    q = n_10/n_1
    return p,q
if __name__ == "__main__":
    dataFile = '/Users/liujia/PycharmProjects/can406_lab1/loss_seq.mat'
    data = scio.loadmat(dataFile)#get sample sequence from .mat file
    sample_seq = data['x'][0]
    p,q = sgm(sample_seq)
    print("p = ", p)
    print("q = ", q)
    #according to sample sequence, calculate p and q by sgm
    #then based on it, generated sequence
    tr = np.array([[1-p, p],
                   [q, 1-q]])
    new_seq = sgm_generate(10000,tr)
    p_,q_ = sgm(new_seq)
    print(p/(p+q))
    count = np.count_nonzero(new_seq)
    theo_loss_p = p_/(p_+q_)
    act_loss_p = count/len(new_seq)
    L = int(1/q)
    print(new_seq)
    print('The theoretical packet loss rate is: ', format(theo_loss_p,'.4f'))
    print('The actual packet loss rate is: ', act_loss_p)
    print("the mean burst length is: ",L)
    #visualization
    visualization(new_seq, sample_seq)
    plt.psd(new_seq, NFFT=256, Fs=2, sides = 'twosided')
    plt.show()