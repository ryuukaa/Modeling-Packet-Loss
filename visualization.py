#!/usr/bin/env python
# -*- coding: utf-8 -*-
##
# @file     visualization.py
# @author   Ryuuka <Jia.Liu21@student.xjtlu.edu.cn>
# @date     2022-04-12
#
# @brief    A function for visualization and drawing histograms to compare sample and generated sequence
#
from binary_runlengths import binary_runlengths
from matplotlib import pyplot as plt

def visualization(seq1, seq2):
    len_zero = binary_runlengths(seq1)[0]
    len_one = binary_runlengths(seq1)[1]
    slen_zero = binary_runlengths(seq2)[0]
    slen_one = binary_runlengths(seq2)[1]
    plt.hist(len_zero,10,histtype='step', stacked=True, label='sgm generated sequence')
    plt.hist(slen_zero,10,histtype='step', stacked=True, label='sample sequence')
    plt.title('run lengths for zero of binary sequences')
    plt.legend()
    plt.show()
    plt.hist(len_one, 10, histtype='step', stacked=True, label='sgm generated sequence')
    plt.hist(slen_one, 10, histtype='step', stacked=True, label='sample sequence')
    plt.title('run lengths for one of binary sequences')
    plt.legend()
    plt.show()