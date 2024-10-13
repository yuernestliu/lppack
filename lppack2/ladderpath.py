"""
Version 2.0.3
Authors: Dinger & ecsLab (Yu Liu, Jingwen Zhang, et al.)
Data: 2024.09.28

梯径的标准化格式 输出, V1.0.1.20240928_Alpha
(& 清空/填满缓存 & 读取JSON文档 & 偏序多重集表示 & 梯图表示)
"""


import re
import json
import heapq
from time import time
import random
from collections import defaultdict, Counter
from dataclasses import dataclass
import graphviz
from typing import List, Dict, Tuple


@dataclass
class Pattern:
    s: str
    update: int                   ## Pattern更新的次数
    idxs: List[Tuple[int, int]]   ## Pattern.s出现的上家梯元的id以及起始位置
    def key(self) -> int:         ## Pattern.s在所有上家中出现的次数
        return len(self.idxs)
    def __lt__(self, other):      ## 出现次数多的Pattern会在堆中排前面
        # for big heap
        return self.key() > other.key()

ID_t = int  ## 定义类型别名
class Ladderon:
    def __init__(self, id, s):
        ## 梯元的基本信息
        self.ID: ID_t = id                                   ## 序号
        self.STR: str = s                                    ## 内容
        self.POS: Dict[ID_t, List[int]] = defaultdict(list)  ## 在上级梯元中的位置信息{ID:[起始位置们]}
        self.COMP: List[Tuple[int, ID_t]] = []               ## 构成表示[基础组件ID，梯元ID]
    def make_ref(self) -> 'LadderonRef':
        ## 生成平凡子串：最大子串
        return LadderonRef(self, 0, len(self.STR), self.STR)
    def asdict(self) -> dict:
        ## 把梯元属性生成一个字典
        return {"ID": self.ID, "STR": self.STR, "POS": self.POS, "COMP": self.COMP}
    def __repr__(self):
        ## 把梯元信息的字典形式写成JSON文件
        return json.dumps(self.asdict())

@dataclass
class LadderonRef:
    ## 一个Ladderon的子串
    ladderon: Ladderon  ## 谁的子串
    start: int          ## 在self.ladderon中的起始位置
    end: int            ## 在self.ladderon中的结束位置(前闭后开)
    s: str              ## 子串的内容（字符序列）
    def slice(self, start, end) -> 'LadderonRef':
        ## 通过首尾位置的挪动，生成self.ladderon的一个新子串
        return LadderonRef(self.ladderon, self.start + start, self.start + end, self.s[start:end])


def find_components(ss: List[LadderonRef]) -> List[Dict[str, List[Tuple[int, int]]]]:
    last_cs = defaultdict(list)
    for idx_s, s in enumerate(ss):
        if s is None:
            continue
        for idx_c, c in enumerate(s.s):
            ## 初始化：last_cs的键值对为所有的{基础单元：[(上家梯元id,起始位置),...]}
            last_cs[c].append((idx_s, idx_c))
    result = []
    while True:
        result.append(last_cs)
        new_cs = defaultdict(list)
        for c, idxs in last_cs.items():
            l = len(c)
            cs = defaultdict(list)
            for idx_s, idx_c in idxs:
                if idx_c + l >= len(ss[idx_s].s):
                    ## 如果超过梯元的长度，跳过
                    continue
                ## 将每个子串往后扩展1，继承扩展前的梯元上家id和起始位置
                cs[ss[idx_s].s[idx_c:idx_c+l+1]].append((idx_s, idx_c))
            for k,v in cs.items():
                if len(v) <= 1:
                    ## 判断：如果扩展子串k不是一个重复结构，跳过
                    continue
                if len(set(i for i,_ in v)) == 1:
                    ## 如果两个位置上的扩展子串k都在同一个梯元里面，并且有重叠，也跳过
                    idx_cs = [j for _,j in v]
                    if max(idx_cs) - min(idx_cs) < l + 1:
                        continue
                ## 保留合法的扩展子串k及其信息v
                new_cs[k].extend(v)
        if new_cs:
            ## 得到了扩展后更长的一组子串
            last_cs = new_cs
        else:
            ## 如果没有更长的重复结构，结束循环
            break
    ## 最后result是所有基础单元和重复结构的字典，值是对应的[(上家梯元id,起始位置),...]
    return result

def find_components_with_c(ss: List[LadderonRef], luts: List[Dict[int, List[int]]], p: Pattern) -> List[Tuple[int, int]]:
    if p.update >= len(luts):
        ## 如果p.update大于或等于luts中梯元的个数，说明不需要进一步扩展或更新。
        assert p.update == len(luts)
        ## 确保p.update与luts的长度相等，如果相等，直接返回当前模式p的匹配位置p.idxs。
        return p.idxs
    ## 如果p.update小于luts的长度
    cp = re.compile(p.s)
    new_idxs = []
    grouped_idxs = defaultdict(list)
    for idx_s, idx_c in p.idxs:
        ## 将同一个id中的不同出现位置整理到一起
        grouped_idxs[idx_s].append(idx_c)
    for idx_s, idx_cs in grouped_idxs.items():
        need_update = False
        ## 如果这个上家是luts中id从p.update往后的某个梯元的上家，则需要更新
        for l in range(p.update, len(luts)):
            if idx_s in luts[l]:
                need_update = True
                break
        if need_update:
            idx_ss = {idx_s}
            for l in range(p.update, len(luts)):
                next_idx_ss = set()
                for i in idx_ss:
                    if i in luts[l]:
                        next_idx_ss.update(luts[l][i])
                    else:
                        next_idx_ss.add(i)
                idx_ss = next_idx_ss
            for i in idx_ss:
                ## 在target里面找
                for it in cp.finditer(ss[i].s):
                    new_idxs.append((i, it.start()))
        else:
            ## 如果不需要更新的话，直接把整理好的id变回原本[（ID，开始位置）,...]的形式
            for idx_c in idx_cs:
                new_idxs.append((idx_s, idx_c))
    if len(new_idxs) <= 1:
        ## 不是重复结构，返回空值
        return None
    if len(set(i for i,_ in new_idxs)) == 1:
        ## 如果所有的重复结构都在同一个梯元中出现，且有重叠，也返回空值
        idx_cs = [j for _,j in new_idxs]
        if max(idx_cs) - min(idx_cs) < len(p.s):
            return None
    # print(f"update {p.s} in {ss} with {luts} from {p.idxs} to {new_idxs}")
    return new_idxs

def find_ladderpath_of_level0(ss: List[LadderonRef], components: List[Ladderon]) -> int:
    base_components: Dict[str, Ladderon] = {}
    for s in ss:
        if s is None: continue
        for i, c in enumerate(s.s):
            if c not in base_components:
                base_components[c] = Ladderon(len(components), c)
                components.append(base_components[c])
            base_components[c].POS[s.ladderon.ID].append(s.start + i)
            s.ladderon.COMP.append((s.start + i, base_components[c].ID))
    return sum(0 if s is None else len(s.s) for s in ss)

def find_ladderpath_with_cs(ss: List[LadderonRef], p: Pattern, components: List[Ladderon]) -> Tuple[List[LadderonRef], Dict[int, List[int]], int]:
    # print(f"replace {p.s} in {ss} with {p.idxs}")
    m_idx_cs = defaultdict(list)
    for idx_s, idx_c in p.idxs:
        m_idx_cs[idx_s].append(idx_c)
    replaced_cnt = 0
    lut = {}
    new_ladderon = Ladderon(len(components), p.s)
    for idx_s, idx_cs in m_idx_cs.items():
        s: LadderonRef = ss[idx_s]
        new_s: List[LadderonRef] = []
        lut[idx_s] = []
        last_i = 0
        for i in sorted(idx_cs):
            assert i + len(p.s) <= len(s.s)
            assert s.s[i:i+len(p.s)] == p.s
            if i < last_i:
                continue
            if i > last_i:
                new_s.append(s.slice(last_i, i))
            replaced_cnt += 1
            new_ladderon.POS[s.ladderon.ID].append(s.start + i)
            s.ladderon.COMP.append((s.start + i, new_ladderon.ID))
            last_i = i + len(p.s)
        if last_i < len(s.s):
            new_s.append(s.slice(last_i, len(s.s)))
        if not new_s:
            ss[idx_s] = None
        else:
            ss[idx_s] = new_s[0]
            lut[idx_s].append(idx_s)
        ss.extend(new_s[1:])
        lut[idx_s].extend(list(range(len(ss) - len(new_s) + 1, len(ss))))
    ss.append(new_ladderon.make_ref())
    c_pos = len(ss) - 1
    for idx_s in m_idx_cs:
        lut[idx_s].append(c_pos)
    # print(f"got {ss} with {lut}")
    assert replaced_cnt > 1
    components.append(new_ladderon)
    return ss, lut, replaced_cnt - 1

def find_ladderpath(ss: List[str]) -> Tuple[int, int, int, List[Ladderon]]:
    str_lens = [len(s) for s in ss]   ## 提取输入targets的长度列表
    total_str_len = sum(str_lens)     ## 计算规模度
    ## 第一次输出：target个数、最长字符串长度、规模度
    # print(f'{len(ss)} strings, max_str_len: {max(str_lens)}, total_str_len: {total_str_len}')
    total_start_time = time()
    ## 先把target标上序号变成Ladderon类
    components = [Ladderon(i, s) for i, s in enumerate(ss)]
    ## 把ss改成LadderonRef类的列表
    ss = [LadderonRef(l, 0, len(l.STR), l.STR) for l in components]

    lp = 0
    luts = []
    component_start_time = time()
    ## 找到所有的重复结构及其关系
    css = find_components(ss)
    component_end_time = time()
    need_recompute = False
    level = len(css) - 1
    new_level = True
    # print(f'Start from {level} level, components time: {component_end_time - component_start_time}s')
    while level > 0:
        if need_recompute:
            luts = []
            component_start_time = time()
            last_duration = component_start_time - component_end_time
            css = find_components(ss)
            component_end_time = time()
            # print(f'after {last_duration}s, Recompute components using {component_end_time - component_start_time}s')
            need_recompute = False
            assert len(css) - 1 <= level
            if level != len(css) - 1:
                new_level = True
                # if level_component_cnt > 0:
                #     print(f'processed {level} level, {len(ss)} strings, {level_component_cnt} components')
            level = len(css) - 1
            if level <= 0:
                break
        if new_level:
            level_component_cnt = 0
            # process_start_time = time()
        new_level = False
        cs = [Pattern(c, 0, idxs) for c, idxs in css[level].items()]
        heapq.heapify(cs)
        while True:
            if not cs:
                new_level = True
                break
            max_cs = None
            updated_cs: List[Pattern] = []
            while max_cs is None or max_cs.key() < cs[0].key():
                next_max = heapq.heappop(cs)
                idxs = find_components_with_c(ss, luts, next_max)
                if idxs is not None:
                    new_c = Pattern(next_max.s, len(luts), idxs)
                    updated_cs.append(new_c)
                    if max_cs is None or len(idxs) > len(max_cs.idxs):
                        max_cs = new_c
                if len(cs) == 0:
                    break
            if time() - component_end_time > component_end_time - component_start_time:
                need_recompute = True
            if max_cs is None:
                new_level = True
                break
            for c in updated_cs:
                if c is not max_cs:
                    heapq.heappush(cs, c)
            level_component_cnt += 1
            # if level_component_cnt % 1000 == 0:
            #     now_time = time()
            #     print(f'processing {level} level, {len(ss)} strings, {level_component_cnt} components using {now_time - process_start_time}s')
            #     process_start_time = now_time
            ss, lut, new_lp = find_ladderpath_with_cs(ss, max_cs, components)
            lp += new_lp
            luts.append(lut)
            if need_recompute:
                break
        if new_level:
            # if level_component_cnt > 0:
            #     print(f'processed {level} level, {len(ss)} strings, {level_component_cnt} components')
            level -= 1
    lp += find_ladderpath_of_level0(ss, components)
    duration = time() - total_start_time
    # print(f'LP: {lp}, order: {total_str_len - lp}, size: {total_str_len}, time: {duration}s')
    for l in components:
        l.COMP = [i for _, i in sorted(l.COMP)]
    
    return lp, total_str_len, duration, components



# 输出 JSON 格式
def get_ladderpath(targets0: List[str], 
    info='V1.0.1.20240928_Alpha', 
    estimate_eta = False, # 是否需要估算eta
    estimate_eta_para = [10, 'global'], # 计算eta时，需要几次以计算omega_min，采用哪种方法
    save_file_name=None, 
    show_version=True) -> dict: #输出是梯径的JSON标准格式

    if not valid_input(targets0): # 判断targets0是否在合法的
        print('Wrong: Input is not valid!!!')
        return

    if show_version:
        print('This version of the ladderpath JSON format is:', info)

    ss, duplications_info = uniquenize(targets0) # 去重

    lp, total_str_len, duration, components = find_ladderpath(ss) #确保ss中没有重复
    length = len(ss)       # 字符串个数
    # 标准格式初始化
    data = {
        "info": info, 
        "ladderpath-index": lp, 
        "order-index": total_str_len - lp, 
        "size-index": total_str_len, 
        "eta" : None,

        "ladderons": {}, 
        "basic_building_blocks": [], 
        "targets": {},
        "duplications_info": {},
        "eta_info": {
            "omega_max_AllIdentical": None,
            "omega_max_Sorted": None,
            "omega_min_Shuffle_list": [],
            "omega_min_LocalDist_list": [],
            "omega_min_EvenDist_list": []
            }
        }

    
    # 筛选出所有的基础单元
    num_of_basic_blocks = 0
    for item in components:
        if len(item.STR) == 1:
            data["basic_building_blocks"].append(item.STR)
            num_of_basic_blocks += 1

    # 处理所有的target和梯元
    for i in range(len(components) - num_of_basic_blocks):
        ladderon_info = components[i]
        id = ladderon_info.ID
        # 将原本的ID转换为标准ID
        if id < length:
            new_id = -id - 1
        else:
            new_id = id - length
        string = ladderon_info.STR
        pos = ladderon_info.POS
        comp = ladderon_info.COMP
        new_comp = []
        
        for j in comp:
            if len(components[j].STR) == 1:
                # 将组成列表中的基础元件由ID替换成对应的字符
                basic_block = components[j].STR
                if len(new_comp) >= 1 and type(new_comp[-1]) == str:
                    # 相邻字符连接成串
                    new_comp[-1] += basic_block
                else:
                    new_comp.append(basic_block)
            else:
                # 将原本的ID转换为标准ID，同样分梯元和target处理
                if j < length:
                    new_comp.append(-j - 1)
                else:
                    new_comp.append(j - length)
        
        if not pos:
            # 判断是target的条件
            data["targets"][new_id] = [new_comp, len(string), "", 1] #暂时REP设为1

        else:
            # 如果是梯元，则继续处理pos
            new_pos = {}
            for index, situ in pos.items():
                # 对pos中的ID进行转换
                if index < length:
                    index = -int(index) - 1
                else:
                    index = int(index) - length
                new_pos[index] = situ
            data["ladderons"][new_id] = [new_comp, len(string), "", new_pos]

    if duplications_info: # 处理targets中有重复的情况
        temp_dup_info = {}
        for s, val in duplications_info.items():
            thisTargetID = val[0]  # val[0]是target的ID
            extra_duplications = len(val)-2
            thisInfoList = data["targets"][thisTargetID] #该target的infoList

            thisInfoList[3] += extra_duplications
            data["ladderpath-index"] += extra_duplications
            data["size-index"] += (extra_duplications * len(s))

            if len(thisInfoList[0]) > 1: #表明 该s还不是梯元，则要将其放入梯元(否则 该s已经既是target也是梯元）
                newLadderonID = len(data['ladderons'])
                data['ladderons'][newLadderonID] = thisInfoList[:3] + [{ thisTargetID:[0] }]
                thisInfoList[0] = [newLadderonID]

                for key, infoList in data['ladderons'].items(): #将连接至targets的ID改成ladderon的ID
                    if key != newLadderonID and thisTargetID in infoList[3]:
                        infoList[3][newLadderonID] = infoList[3].pop(thisTargetID)
            temp_dup_info[val[0]] = val[1:]
        data["order-index"] = data["size-index"] - data["ladderpath-index"]
        data['duplications_info'] = temp_dup_info
    
    if estimate_eta: # 计算eta，这里默认使用global方法: omega_max_AllIdentical & omega_min_Shuffle_list
        if estimate_eta_para[1] == 'global':
            combined_str = ''.join(targets0)
            omega_max = cal_omega_max_AllIdentical(combined_str)
            data['eta_info']['omega_max_AllIdentical'] = omega_max

            temp_list = cal_omega_min_Shuffle_list(combined_str, estimate_eta_para[0])
            data['eta_info']['omega_min_Shuffle_list'] = temp_list
            omega_min = min(temp_list)

            data['eta'] = (data['order-index'] - omega_min) / (omega_max - omega_min)
        else:
            print('Warning: estimate_eta_para[1] is not global, abort function.')


    # 将处理好的格式导入JSON文件
    if save_file_name:
        save_file_name += '.json' if not save_file_name.endswith('.json') else ''
        save_ladderpath_json(data, save_file_name)

    return data


# 将targets序列去重
def uniquenize(strs):
    counter = Counter(strs)
    if len(strs) != len(counter):
        strs_unique = []  # 最后要输出的没有重复序列的targets
        duplications_info = {}  # 最后要输出的重复序列的信息 Dict[str:List(int)]。字典的value(List)第一entry是这个str的ID，后面的entries是其str在原始序列中出现的所有位置
        for i, s in enumerate(strs):
            if counter[s] == 1:
                strs_unique.append(s)
            else:
                if s in duplications_info:
                    duplications_info[s].append(i)
                else:
                    duplications_info[s] = [-len(strs_unique)-1, i]
                    strs_unique.append(s)
        return strs_unique, duplications_info
    else:
        return strs, {}



# =============== 对梯径JSON格式的 标准后续处理 ===============
# 储存JSON文件
def save_ladderpath_json(data, save_file_name):
    with open(save_file_name, "w", encoding="utf-8", newline='') as f:
        json.dump(data, f, ensure_ascii=False, indent=4)


# 读取JSON文件，梯径标准格式（会将dict的key转换成数字）
def load_ladderpath_json(filename):
    with open(filename, 'r', encoding='utf-8') as file:
        lpjson = json.load(file, object_pairs_hook=convert_keys)
    return lpjson

def convert_keys(pairs):
    return {int(k) if is_integer(k) else k: v for k, v in pairs}

def is_integer(s):
    if s.startswith('-'):  # 检查是否以负号开头
        return s[1:].isdigit()  # 如果是负数，检查负号后的部分是否为数字
    return s.isdigit()  # 如果不是负数，直接检查是否为数字




# 清空JSON中每个梯元STR的显式表示（清空缓存）
def clear_lpjson_STR(lpjson):
    for infoList in lpjson['ladderons'].values():
        infoList[2] = ''
    for infoList in lpjson['targets'].values():
        infoList[2] = ''


# 填满JSON中每个梯元STR的显式表示（填满缓存）
def fill_lpjson_STR(lpjson):
    for ldID in reversed(lpjson['ladderons'].keys()): #先从小的梯元填满
        lpjson['ladderons'][ldID][2] = fill(ldID, lpjson)
    for infoList in lpjson['targets'].values():
        if infoList[2]:
            assert len(infoList[2]) == infoList[1] # 判断必要条件
            pass
        else:
            thisstr = ''
            for i_COMP in infoList[0]:
                if isinstance(i_COMP, int):
                    thisstr += lpjson['ladderons'][i_COMP][2]
                else:
                    thisstr += i_COMP
            infoList[2] = thisstr

def fill(ldID, lpjson):
    if lpjson['ladderons'][ldID][2]:
        assert len(lpjson['ladderons'][ldID][2]) == lpjson['ladderons'][ldID][1] # 判断必要条件
        return lpjson['ladderons'][ldID][2]
    thisstr = ''
    for i_COMP in lpjson['ladderons'][ldID][0]:
        if isinstance(i_COMP, int):
            thisstr += fill(i_COMP, lpjson)
        else:
            thisstr += i_COMP
    lpjson['ladderons'][ldID][2] = thisstr
    return thisstr



# 梯径的偏序多重集表示形式 partially_ordered_multiset
def POM_from_JSON(lpjson, display_str=False):
    temp = bbb_multiplicity(lpjson) # basic_building_blocks的计数
    pom_final_pool = {0: temp} #最终的偏序多重集 dict表示

    fill_lpjson_STR(lpjson)
    lds = lpjson['ladderons']

    level_known = {} #每个梯元位于哪一层
    for ldID in reversed(lds.keys()): #先从小的梯元开始
        find_level(ldID, level_known, lds)

    each_level_comp = defaultdict(list) #每层中包含了哪些梯元
    for ldID, level in level_known.items():
        each_level_comp[level].append(ldID)
    if len(each_level_comp) > 0: #如果不是完全由basic building blocks 构成
        assert len(each_level_comp) == max(each_level_comp.keys()) #确保没有空的层

    for i in range(1, len(each_level_comp)+1):
        levelinfo = {}
        for ldID in each_level_comp[i]:
            levelinfo[lds[ldID][2]] = sum(len(val) for val in lds[ldID][3].values())-1 #获得梯元的重数
        original_order = list(levelinfo.keys())
        temp = sorted(levelinfo.items(), key=lambda item: (-item[1], -original_order.index(item[0]))) #排序 重数从大到小，字符串按照ID倒序
        pom_final_pool[i] = dict(temp)

    for ID, pos in lpjson['duplications_info'].items(): # 处理targets中有重复的情况
        for val_strs in pom_final_pool.values():
            if lpjson['targets'][ID][2] in val_strs: #lpjson['targets'][ID][2]是ladderon STR本身
                val_strs[ lpjson['targets'][ID][2] ] += (len(pos) - 1) # 更改的是pom_final_pool

    pom_str = ''
    if display_str:
        pom_str = str_partially_ordered_multiset(pom_final_pool)
        print(pom_str)

    return pom_final_pool, pom_str
    # final_pool是偏序多重集的字典表示
    # {0: {'A': 3, 'D': 3, 'C': 2, 'B': 1, 'E': 1, 'F': 1},
    #  1: {'EF': 2, 'BC': 2, 'DCD': 1},
    #  2: {'DBC': 1},
    #  3: {'DBCDBC': 1}}

def find_level(x, level_known, lds): #找到每个梯元相应的层级（递归）
    if x in level_known:
        return level_known[x]
    else:
        infoList0 = lds[x][0]
        if len(infoList0) == 1 and isinstance(infoList0[0], str):
            level_known[x] = 1
            return 1
        else:
            comp_levels = [] # 所有组分位于的层级
            for i in infoList0:
                if isinstance(i, int):
                    comp_levels.append( 1 + find_level(i, level_known, lds) )
            level_known[x] = max(comp_levels)
            return level_known[x]

def str_partially_ordered_multiset(pom): #转化成字符串形式
    pom_str = '{ '
    # 使用列表生成式简化拼接
    pom_str += ' // '.join(
        ', '.join(f'{ld}({multi})' if multi != 1 else f'{ld}' for ld, multi in infoList.items()) 
        for infoList in pom.values()
    )
    pom_str += ' }'
    return pom_str



# 展示3个指标
def disp3index(lpjson):
    index3 = (lpjson['ladderpath-index'], lpjson['order-index'], lpjson['size-index'])
    print(f'( Ladderpath-index:{index3[0]},  Order-index:{index3[1]},  Size-index:{index3[2]} )')
    return index3



# 如果eta存在，则取出；如果不存在，则计算
def get_eta(lpjson, estimate_eta_para=[10, 'global'], use_as_update=False):
    # use_as_update=False 这个参数不允许外部设为 True
    if (lpjson['eta'] is not None) and (not use_as_update):
        return lpjson['eta']
    else:
        fill_lpjson_STR(lpjson)
        targets_list = reconstruct_targets_list(lpjson)
        targets_len, combined_str_list = [], []
        for ID in targets_list:
            temp = lpjson['targets'][ID]
            targets_len.append( temp[1] ) # 长度的list
            combined_str_list.append( temp[2] )
        combined_str = ''.join(combined_str_list)

        if use_as_update:
            omega_max = lpjson['eta_info']['omega_max_AllIdentical']
        else:
            omega_max = cal_omega_max_AllIdentical(combined_str)
            lpjson['eta_info']['omega_max_AllIdentical'] = omega_max

        # 计算omega_min: 计算N_omega_min次，取最小的
        omega_min_list = cal_omega_min_Shuffle_list(combined_str, estimate_eta_para[0])
        if use_as_update:
            lpjson['eta_info']['omega_min_Shuffle_list'].extend( omega_min_list )
        else:
            lpjson['eta_info']['omega_min_Shuffle_list'] = omega_min_list
        omega_min = min(lpjson['eta_info']['omega_min_Shuffle_list'])
        
        eta = (lpjson['order-index'] - omega_min) / (omega_max - omega_min)
        lpjson['eta'] = eta
        return eta

def update_eta(lpjson, estimate_eta_para=[10, 'global']): #更新eta值，即多增加N_omega_min次计算omega_min
    if lpjson['eta'] is None:
        print('Wrong & Abort: There is no eta value. get_eta() first and then update_eta().')
        return
    else:
        get_eta(lpjson, estimate_eta_para, use_as_update=True)



# 画梯图
def ellipse_len(seq): # 画梯图的associated函数
    # length = np.sqrt(len(seq))/2
    length = (len(seq)**(1/3) )/2
    return length

def draw_laddergraph(lpjson, show_longer_than = 0, style = "ellipse", 
    warning_n_ladderons_to_show = 500,
    rankdir = "BT", color = "grey",
    save_fig_name = None, figformat = "pdf", cleanGVfile=True):
    # Draw the laddergraph.
    # "show_longer_than": When the length of the ladderon > show_longer_than, this ladderon will be displayed.
    #     Note that "show_longer_than" should always be >= 1, and the basic building blocks are also omitted.
    # "style" dictates how the laddergraph is displayed. It can either be 
    #     "ellipse" (the sequence won't be displayed, but the size of the ellipse is positively related to the length of the sequence),
    #     or "box" (the sequence will be displayed),
    #     or "box-OnlyShowTargetID" (the sequence will be displayed, but only show the target's ID).
    # "warning_n_ladderons_to_show": The largest number of ladderons allowed to draw.
    # "rankdir" is the order of the nodes, should be BT (from bottom to top), TB, LR, RL.
    # "color" can be "grey", "red", "#808080", etc.
    # "save_fig_name" is the file name of the figure.
    # "figformat" is the format of the exported figure, which could be "pdf", "png", etc.
    # "cleanGVfile" 是否删除render图片的.gv格式文件
    
    n_ladderons_to_show = 0
    for ldID, infoList in lpjson['ladderons'].items():
        if len(infoList[2]) > show_longer_than:
            n_ladderons_to_show += 1
    if n_ladderons_to_show > warning_n_ladderons_to_show:
        print(f'Warning: too many ladderons (>{n_ladderons_to_show}) to draw!!!')


    if style not in ['ellipse', 'box', 'box-OnlyShowTargetID']:
        print('! Wrong: the parameter \'style\' can only be either \'ellipse\', \'box\', \'box-OnlyShowTargetID\'')
        return

    onlyshowtargetid = False
    if style == "box-OnlyShowTargetID":
        style = "box"
        onlyshowtargetid = True
    
    ladderons = lpjson['ladderons']    
    fill_lpjson_STR(lpjson) # 先将STR填满缓存

    targets = lpjson['targets']
    targets_2ID = {val[2]: ID for ID, val in targets.items()} # targets_strs按输入来，不一定是-1，-2，-3...
    ladderons_2ID = {val[2]: ID for ID, val in ladderons.items()}
    bbbs_2ID = { val:f'b{i}' for i, val in enumerate(lpjson['basic_building_blocks'])} #按bbb List中的顺序来，赋予其ID b1, b2,...

    # initialize the graph
    g = graphviz.Digraph()
    g.attr(rankdir=rankdir)
    if style == "box": # a detailed version, showing ladderons length
        g.attr('node', shape='box')
        for target, ID in targets_2ID.items(): # draw all target sequences
            if target in ladderons_2ID: # 说明这个target也是梯元
                if onlyshowtargetid:
                    g.node(str(ID), label = str(ID) + ' (&)')
                else:
                    g.node(str(ID), label = target + ' (&)')
            else:
                if onlyshowtargetid:
                    g.node(str(ID), label = str(ID))
                else:
                    g.node(str(ID), label = target)
        if show_longer_than == 0: # display basic building blocks or not
            for bbb, IDstr in bbbs_2ID.items():
                g.node(IDstr, label = bbb, shape = 'hexagon', color=color)

    else: # style == "ellipse":
        g.attr('node', shape='ellipse')
        for target, ID in targets_2ID.items(): # draw all target sequences
            length = ellipse_len(target)
            if target in ladderons_2ID:
                IDanother = ladderons_2ID[target]
                templabel = f'{ID}({IDanother})'
                g.node(str(ID), label = templabel, **{'width':str(length),'height':str(length/2)})
            else:
                g.node(str(ID), label = str(ID), **{'width':str(length),'height':str(length/2)})
        if show_longer_than == 0: # display basic building blocks or not
            length = ellipse_len( next(iter(bbbs_2ID)) )
            for bbb, IDstr in bbbs_2ID.items():
                g.node(IDstr, label = bbb, **{'width':str(length),'height':str(length/2)}, \
                        shape = 'hexagon', color = color)

    for ldID, infoList in ladderons.items():
        dict_temp = dict()
        for linkedTo, positions in infoList[3].items(): # 遍历POS
            if len(infoList[2]) > show_longer_than:
                if style == "box":
                    g.node(str(ldID), label = infoList[2], style = 'filled', color=color)
                else: # style == "ellipse":
                    lenk = ellipse_len(infoList[2])
                    g.node(str(ldID), label = str(ldID), **{'width':str(lenk),'height':str(lenk/2)}, \
                            style = 'filled', color = color)
                dict_temp[str(ldID), str(linkedTo)] = len(positions)

        for id_linkedto, times in dict_temp.items():
            for _ in range(times):
                g.edge(id_linkedto[0], id_linkedto[1], color = color)

    for ID, pos in lpjson['duplications_info'].items(): # 处理targets中有重复的情况
        for _ in range(len(pos)-1):
            g.edge(str( lpjson['targets'][ID][0][0] ), str(ID), color = color)

    if show_longer_than == 0: # display the links from basic building blocks
        for ldID, infoList in ladderons.items():
            for comp0 in infoList[0]: # COMP
                if isinstance(comp0, str):
                    for bbb in comp0:
                        g.edge(bbbs_2ID[bbb], str(ldID), color = color)
        for ldID, infoList in targets.items():
            if infoList[2] not in ladderons: #如果既是target又是ladderon，这里就不画了，因为上面画过了
                for comp0 in infoList[0]: # COMP
                    if isinstance(comp0, str):
                        for bbb in comp0:
                            g.edge(bbbs_2ID[bbb], str(ldID), color = color)

    if save_fig_name:
        g.render(filename=save_fig_name, format=figformat, cleanup=cleanGVfile)
    return g


# ================ 辅助函数 ================
# 判断最初的input strs是否合法
def valid_input(targets0):
    for s in targets0:
        if len(s) <= 1: #不能是最基本单元
            print('Wrong: There are strings with length 1 or 0!!!')
            return False
    return True


# basic_building_blocks的计数
def bbb_multiplicity(lpjson):
    bbb_multi = defaultdict(int)
    for item_name in ['ladderons', 'targets']:
        for infoList in lpjson[item_name].values():
            for i in infoList[0]:
                if isinstance(i, str):
                    # 将字母计数结果添加到 bbb_multi 中
                    for letter, count in Counter(i).items():
                        bbb_multi[letter] += count
    temp = dict(sorted(bbb_multi.items(), key=lambda item: (-item[1], lpjson['basic_building_blocks'].index(item[0])) )) # 排序，重数从大到小
    return temp


# 重构targets的编号序列
def reconstruct_targets_list(lpjson):
    targetslist = []
    for dupID, pos in lpjson['duplications_info'].items():
        if len(targetslist) <= max(pos):
            targetslist.extend([None]* (max(pos)-len(targetslist)+1) )
        for j in pos:
            targetslist[j] = dupID
    i, IDtoFill = 0, 0
    while IDtoFill > -len(lpjson['targets']):
        if i >= len(targetslist):
            targetslist.append(None)
        if targetslist[i] is None:
            IDtoFill -= 1
            targetslist[i] = IDtoFill
        else:
            IDtoFill = min(IDtoFill, targetslist[i])
        i += 1
    return targetslist


# 关于eta的计算：max，全相同序列的omega（不再按各target长度切割）
def cal_omega_max_AllIdentical(combined_str):
    # This function only works for a single sequence
    temp = bin(len(combined_str))[2:] # Calculate the omega_max of a sequence of length S
    lpIdx_AllIdentical = len(temp) + temp.count('1') - 1
    return len(combined_str) - lpIdx_AllIdentical

# 关于eta的计算：max，排序后按照一整个序列combined_str计算，（不再按各target长度切割）
def cal_omega_max_Sorted(combined_str):
    sorted_combined_str = ''.join(sorted(combined_str)) # 字符串排序
    lpindex, sizeindex, _, _ = find_ladderpath([sorted_combined_str])
    return sizeindex - lpindex


# 关于eta的计算：min，随机shuffle原序列（不再按各target长度切割）
def cal_omega_min_Shuffle_list(combined_str, N_omega_min):
    omega_min_list = []
    for _ in range(N_omega_min):
        temp = list(combined_str)
        random.shuffle(temp)
        shuffled_combined_str = ''.join(temp) # 将combined_str打乱

        # strs_new, index = [], 0 # 这是需要切割的方式
        # for length in targets_len:
        #     strs_new.append(shuffled_combined_str[index:index+length])
        #     index += length
        # lpindex, sizeindex, _, _ = find_ladderpath(strs_new)
        # omega_min_list.append( sizeindex - lpindex )

        lpindex, sizeindex, _, _ = find_ladderpath([shuffled_combined_str])
        omega_min_list.append( sizeindex - lpindex )
    return omega_min_list

# 关于eta的计算：min，按照原序列的概率分布
def cal_omega_min_LocalDist_list():
    pass

# 关于eta的计算：min，按照均匀概率分布
def cal_omega_min_EvenDist_list():
    pass
