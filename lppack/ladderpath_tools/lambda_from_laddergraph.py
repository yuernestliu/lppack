"""
Version 2.0.1
Authors: ecsLab (Jingwen Zhang et al.)
Data: 2025.06.05

利用梯图来计算其中两个targets(leafs)之间的梯径距离
"""

import ladderpath as lp


def get_target_ID(two_leaf_strs0, lpjson):
    # 找到subset中几个target的序号
    strs_subset = two_leaf_strs0.copy()
    ids=[]
    for ID, info in lpjson["targets"].items():
        if info[2] in strs_subset:
            ids.append(ID)
            strs_subset.remove(info[2])
    assert len(ids) == len(two_leaf_strs0)
    return ids

def separate_basic_blocks(COMP):
    # 把基本单位和梯元分开，返回梯元的id（保留重复）和基本单位的使用个数
    new_COMP=[]
    block_count=0
    for i in COMP:
        if type(i)==str:
            block_count+=len(i)
        else:
            new_COMP.append(i)
    return new_COMP,block_count

def if_id_in_COMP(COMP):
    # 查看是不是所有的梯元都划完了
    for i in COMP:
        if type(i)==int:
            return True
    return False


def get(leaves1or2, lpjson): # leaves1or2 could only be 1 or 2 leaf
    ids = None
    if len(leaves1or2) == 1 or len(leaves1or2) == 2:
        lp.fill_lpjson_STR(lpjson)
        if isinstance(leaves1or2[0], int):
            if leaves1or2[0] in lpjson['targets']:
                if len(leaves1or2) == 1:
                    ids = leaves1or2
                else: # len(leaves1or2) == 2
                    if leaves1or2[1] in lpjson['targets']:
                        ids = leaves1or2
        else:
            ids = get_target_ID(leaves1or2, lpjson)
    
    if ids is None:
        print('!!!Wrong: input wrong -> leaves1or2')
        return None

    basic_block_num=0
    comp=[]
    ladderon_dic={}
    for i in ids:
        comp+=lpjson["targets"][i][0]
    pool=set()

    while if_id_in_COMP(comp):
        ids,basic_blocks=separate_basic_blocks(comp)
        comp=[]
        basic_block_num+=basic_blocks
        id_set=set()
        for i in ids:
            if i in ladderon_dic:
                ladderon_dic[i]+=1
            else:
                ladderon_dic[i]=1
            if i not in id_set and i not in pool:
                comp+=lpjson["ladderons"][i][0]
                pool.add(i)
                id_set.add(i)

    for i in comp:
        basic_block_num+=len(i)
    count=basic_block_num

    for ladderon,times in ladderon_dic.items():
        count+=(times-1)
    return count
