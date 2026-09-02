"""
Version 1.0
Authors: ecsLab
Data: 2025.01.05

计算序列的有序率eta
"""

import random
from .. import ladderpath as lp



def get_omega_max(seqs, display=False):
    """
    计算omega_max，全A的序列，并且按原来的各个targets的长度切割
    """
    allAs = []
    for i in seqs:
        allAs.append('A'*len(i))
    lpjson = lp.get_ladderpath(allAs, save_file_name=None, show_version=False, estimate_eta=False)
    if display:
        print('omega_max:', lpjson['ladderpath-index'], lpjson['order-index'], lpjson['size-index'])
    return lpjson['order-index']


def get_omega_min(seqs, n, display=False):
    """
    计算omega_min，shuffle n 次，并且按原来的各个targets的长度切割
    """
    merged = list(''.join(seqs))
    omega = 0
    for _ in range(n):
        random.shuffle(merged)
        merged_shuffled = ''.join(merged)
        newseqs = []
        i = 0
        for s in seqs:
            newseqs.append(merged_shuffled[i : i+len(s)])
            i += len(s)
        lpjson = lp.get_ladderpath(newseqs, save_file_name=None, show_version=False, estimate_eta=False)
        if display:
            print(lpjson['ladderpath-index'], lpjson['order-index'], lpjson['size-index'])
        omega += lpjson['order-index']
    omega = omega / n
    return omega


def get_eta(seqs, n_min=10, lpjson=None, display=False):
    if lpjson is None:
        lpjson = lp.get_ladderpath(seqs, save_file_name=None, show_version=False, estimate_eta=False)
    omega = lpjson['order-index']

    omega_max = get_omega_max(seqs, display)
    omega_min = get_omega_min(seqs, n_min, display)
    if display:
        print(f"omega: {omega}, omega_max: {omega_max}, omega_min: {omega_min}")

    eta = (omega - omega_min) / (omega_max - omega_min)
    return eta


# 关于eta的计算：max，全相同序列的omega（不再按各target长度切割）
def _cal_omega_max_AllIdentical(aStringOrLen):
    # This function only works for a single sequence
    if isinstance(aStringOrLen, str):
        temp = bin(len(aStringOrLen))[2:] # Calculate the omega_max of a sequence of length S
    elif isinstance(aStringOrLen, int):
        temp = bin(aStringOrLen)[2:]
    else: # 报错
        print('!!! Wrong: input type wrong.')
        return
    lpIdx_AllIdentical = len(temp) + temp.count('1') - 1
    if isinstance(aStringOrLen, str):
        return len(aStringOrLen) - lpIdx_AllIdentical
    else:
        return aStringOrLen - lpIdx_AllIdentical
