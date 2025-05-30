"""
Version 2.0.2
Authors: ecsLab (Yu Liu, Lanxin Ma, et al.)
Data: 2025.04.04

Based on ladderpath standard JSON format (V1.0.2.20250404_Alpha)
"""

import ladderpath as lp

# 利用梯径进行压缩
def compress(lpjson, display=False, SEP='@'):
    try:
        if lpjson['input_type'] != 'list':
            print('Warning: Abort. The input_type for targets must be list, otherwise we cannot compress it. Cannot be dict.')
            return
    except:
        print('Warning: ladderpath standard JSON format has been updated to V1.0.2.20250404_Alpha, or above.')
        print('         What you used is an old version; It will be abandoned in future versions.')
        print('')

    targets, lds = lpjson['targets'], lpjson['ladderons']

    ld_incl = [] # 记录所有牵扯到的梯元
    for infoList in targets.values():
        for i in infoList[0]:
            if isinstance(i, int) and (i not in ld_incl):
                ld_incl.append(i) # 首先将targets中的梯元ID全部按顺序放入（初始化）
    # 读jd2.json，有：ld_incl = [4, 0, 2, 3]

    i = 0
    while i < len(ld_incl):
        infoList0 = lds[ld_incl[i]][0] # 对应的下一层的那个梯元
        for j in infoList0:
            if isinstance(j, int) and (j not in ld_incl):
                ld_incl.append(j)
        i += 1
    ld_incl.sort(reverse=True) # 遍历之后，得到 所有牵扯到的梯元
    # ld_incl = [4, 3, 2, 1, 0]

    convert_dict = {val:i for i, val in enumerate(ld_incl)} #梯元的ID要转换
    # convert_dict = {4: 0, 3: 1, 2: 2, 1: 3, 0: 4}

    compressed_list = [] #从小梯元排起
    for i in ld_incl:
        infoList0 = lds[i][0]
        compressed_list.append(SEP)
        for j in infoList0:
            if isinstance(j, int):
                compressed_list.append(convert_dict[j])
            else:
                compressed_list.append(j) # 得到所有牵涉到的梯元的压缩后的序列
    # compressed_list = ['@', 'BC', '@', 'EF', '@', 'DCD', '@', 'D', 0, '@', 3, 3]

    all_compressed_list = [] #targets的梯元组成序列   
    targetslist = lp.reconstruct_targets_list(lpjson)
    for targetID in targetslist:
        all_compressed_list.append(SEP)
        all_compressed_list.extend(convert_dict[i] if isinstance(i, int) else i for i in targets[targetID][0])
    # all_compressed_list = ['@', 'A', 0, 4, 2, 1, 1, '@', 4, 'A', '@', 0, 2, 'A', 1]

    all_compressed_list.extend(compressed_list)
    all_compressed_list.insert(0, len(targetslist)) #将targets的具体数量记录在最前面
    # all_compressed_list = [3, '@', 'A', 0, 4, 2, 1, 1, '@', 4, 'A', '@', 0, 2, 'A', 1, '@', 'BC', '@', 'EF', '@', 'DCD', '@', 'D', 0, '@', 3, 3]

    if display:
        print(all_compressed_list)
        
    return all_compressed_list



# 解压函数
def decompressed(compressed_list, display=False):
    nTargets = compressed_list[0]  # compressed_list的第一个元素是targets的数量
    SEP = compressed_list[1]
    current_idx = 1  # 从compresses_list的第二个元素开始
    ladderons = {}  # 存储梯元的字典
    basic_building_blocks = set()  # 基本字母组成
    final_strs = []  # 最终解压缩的序列

    # 得到basic_building_blocks，知道字母组成
    for item in compressed_list:
        if isinstance(item, str) and item != SEP:
            basic_building_blocks.update(item)
    # print("basic_building_blocks:", list(basic_building_blocks))

    # 遍历列表找到第 nTargets+1 个 SEP 的位置 前nTargets个@分别是各个target压缩表示，第nTargets+1个@开始才是具体的梯元
    sep_count = 0 #储存便历过程中@出现次数
    for i in range(1, len(compressed_list)):
        if compressed_list[i] == SEP:
            sep_count += 1
        # 当找到第 nTargets+1 个 SEP 时，设置 current_idx 到这个 SEP 之后的位置
            if sep_count == nTargets + 1:
                current_idx = i + 1
                break

    # 找到梯元开始的位置后开始得到具体梯元
    while current_idx < len(compressed_list):
        if compressed_list[current_idx] == SEP:
            current_idx += 1
        ladderon_content = []
        while current_idx < len(compressed_list) and compressed_list[current_idx] != SEP:
            ladderon_content.append(compressed_list[current_idx])
            current_idx += 1
        ladderon_index = len(ladderons)
        ladderons[ladderon_index] = ladderon_content
        current_idx += 1  # 跳过下一个 SEP 或结束

    # 输出结果
    # print("ladderons:", ladderons)

    # 重置current_idx
    current_idx = 1
    
    # 得到目标字符串
    while current_idx < len(compressed_list) and len(final_strs) < nTargets:
        if compressed_list[current_idx] == SEP:
            current_idx += 1
        target_str = []
        while current_idx < len(compressed_list) and compressed_list[current_idx] != SEP:
            item = compressed_list[current_idx]
            if isinstance(item, int):  # 如果是梯元编号，需要解压梯元
                target_str.append(decompress_ladderon(item, ladderons))
            else:  # 如果是字母，直接作为基本构建块
                target_str.append(item)
            current_idx += 1
        final_strs.append(''.join(target_str))
        current_idx += 1  # 跳过下一个SEP

    if display:
        print(final_strs)
    return final_strs


#获取梯元具体指代字母，主要是递归部分
def decompress_ladderon(ladderon_idx, ladderons):
    ladderon = ladderons[ladderon_idx]  # 获取梯元内容
    decompressed_str = []
    for item in ladderon:
        if isinstance(item, int):
            decompressed_str.append(decompress_ladderon(item, ladderons))  # 递归解压梯元
        else:
            decompressed_str.append(item)  # 基本块直接添加
    return ''.join(decompressed_str)

