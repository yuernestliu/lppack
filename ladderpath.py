"""
Version 3.0.0
Authors: ecsLab
Data: 2026.06.12

CUDA and macm backends use Unicode character coordinates for non-ASCII input.
"""

VERSION_ladderpathpy = '3.0.0'
VERSION_lpjson = '3.0.0'


import json
import heapq
from time import time
import random
from collections import defaultdict, Counter
from dataclasses import dataclass
from typing import List, Dict, Tuple


@dataclass
class Pattern:
    s: str
    update: int
    idxs: List[Tuple[int, int]]
    def key(self) -> int:
        return len(self.idxs)
    def __lt__(self, other):
        # for big heap
        return self.key() > other.key()

ID_t = int
class Ladderon:
    def __init__(self, id, s):

        self.ID: ID_t = id
        self.STR: str = s
        self.POS: Dict[ID_t, List[int]] = defaultdict(list)
        self.COMP: List[Tuple[int, ID_t]] = []
    def make_ref(self) -> 'LadderonRef':

        return LadderonRef(self, 0, len(self.STR), self.STR)
    def asdict(self) -> dict:

        return {"ID": self.ID, "STR": self.STR, "POS": self.POS, "COMP": self.COMP}
    def __repr__(self):

        return json.dumps(self.asdict())

@dataclass
class LadderonRef:

    ladderon: Ladderon
    start: int
    end: int
    s: str
    def slice(self, start, end) -> 'LadderonRef':

        return LadderonRef(self.ladderon, self.start + start, self.start + end, self.s[start:end])


def find_components(ss: List[LadderonRef]) -> List[Dict[str, List[Tuple[int, int]]]]:
    last_cs = defaultdict(list)
    for idx_s, s in enumerate(ss):
        if s is None:
            continue
        for idx_c, c in enumerate(s.s):

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

                    continue

                cs[ss[idx_s].s[idx_c:idx_c+l+1]].append((idx_s, idx_c))
            for k,v in cs.items():
                if len(v) <= 1:

                    continue
                if len(set(i for i,_ in v)) == 1:

                    idx_cs = [j for _,j in v]
                    if max(idx_cs) - min(idx_cs) < l + 1:
                        continue

                new_cs[k].extend(v)
        if new_cs:

            last_cs = new_cs
        else:

            break

    return result

def find_all(a: str, b: str):
    index = 0
    while True:
        index = a.find(b, index)
        if index == -1:
            return
        yield index
        index += 1

def find_components_with_c(ss: List[LadderonRef], luts: List[Dict[int, List[int]]], p: Pattern) -> List[Tuple[int, int]]:
    if p.update >= len(luts):

        assert p.update == len(luts)

        return p.idxs

    new_idxs = []
    grouped_idxs = defaultdict(list)
    for idx_s, idx_c in p.idxs:

        grouped_idxs[idx_s].append(idx_c)
    for idx_s, idx_cs in grouped_idxs.items():
        need_update = False

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

                for it in find_all(ss[i].s, p.s):
                    new_idxs.append((i, it))
        else:

            for idx_c in idx_cs:
                new_idxs.append((idx_s, idx_c))
    if len(new_idxs) <= 1:

        return None
    if len(set(i for i,_ in new_idxs)) == 1:

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
    str_lens = [len(s) for s in ss]
    total_str_len = sum(str_lens)

    # print(f'{len(ss)} strings, max_str_len: {max(str_lens)}, total_str_len: {total_str_len}')
    total_start_time = time()

    components = [Ladderon(i, s) for i, s in enumerate(ss)]

    ss = [LadderonRef(l, 0, len(l.STR), l.STR) for l in components]

    lp = 0
    luts = []
    component_start_time = time()

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




def get_ladderpath(targets0,
    fill_ladderons_STR = False,
    save_file_name=None, 
    show_version=True,
    mode='cpu') -> dict:

    
    if not isinstance(targets0, (list, dict)):
        raise TypeError('Input ERROR: Input must be a list or dict!')
    input_type = 'list' if isinstance(targets0, list) else 'dict'
    if not isinstance(mode, str):
        raise TypeError("Parameter ERROR: mode must be one of 'cpu', 'gpu-cuda', or 'macm'.")
    mode = mode.lower()
    if mode not in ('cpu', 'gpu-cuda', 'macm'):
        raise ValueError("Parameter ERROR: mode must be one of 'cpu', 'gpu-cuda', or 'macm'.")


    if input_type == 'list':
        if len(set(targets0)) != len(targets0):
            raise ValueError('Input ERROR: There are duplications in the list. Use dict instead.')
        ss = targets0
        duplications_info = {}
    else: # input_type == 'dict'
        ss = list(targets0)
        duplications_info = {}
        for idx, (_, multi) in enumerate(targets0.items(), start=1):
            if multi > 1:
                duplications_info[-idx] = multi


    for s0 in ss:
        if len(s0) < 2:
            raise ValueError('Input ERROR: Single letters cannot be in the targets.')

    use_unicode_positions = any(not s.isascii() for s in ss)


    if mode == 'gpu-cuda':
        from . import ladderpath_cuda as lp_cuda

        data = lp_cuda.get_ladderpath(ss, ch=use_unicode_positions, show_version=show_version)
        data['input_type'] = input_type

    elif mode == 'macm':
        from . import ladderpath_macm as lp_macm

        data = lp_macm.get_ladderpath(ss, ch=use_unicode_positions, show_version=show_version)
        data['input_type'] = input_type

    else:
        if show_version:
            print('lppack version:', VERSION_ladderpathpy)
            print('lpjson version:', VERSION_lpjson)

        lp, total_str_len, duration, components = find_ladderpath(ss)
        length = len(ss)

        data = {
            "info": VERSION_lpjson, 
            "ladderpath-index": lp, 
            "order-index": total_str_len - lp, 
            "size-index": total_str_len, 

            "basic_building_blocks": [], 
            "ladderons": {}, 
            "targets": {},
            "duplications_info": {},
            "input_type": input_type
            }

        

        num_of_basic_blocks = 0
        for item in components:
            if len(item.STR) == 1:
                data["basic_building_blocks"].append(item.STR)
                num_of_basic_blocks += 1


        for i in range(len(components) - num_of_basic_blocks):
            ladderon_info = components[i]
            id = ladderon_info.ID

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

                    basic_block = components[j].STR
                    if len(new_comp) >= 1 and type(new_comp[-1]) == str:

                        new_comp[-1] += basic_block
                    else:
                        new_comp.append(basic_block)
                else:

                    if j < length:
                        new_comp.append(-j - 1)
                    else:
                        new_comp.append(j - length)
            
            if not pos:

                data["targets"][new_id] = [new_comp, len(string), "", 1]

            else:

                new_pos = {}
                for index, situ in pos.items():

                    if index < length:
                        index = -int(index) - 1
                    else:
                        index = int(index) - length
                    new_pos[index] = situ
                data["ladderons"][new_id] = [new_comp, len(string), "", new_pos]


    if duplications_info:
        if isinstance(next(iter(duplications_info.values())), int):
            input_is_dict = True
            temp_dup_info = duplications_info
        else:
            input_is_dict = False
            temp_dup_info = {}

        for s, val in duplications_info.items():
            if input_is_dict:
                thisTargetID = s
                extra_duplications = val - 1
            else: 
                thisTargetID = val[0]
                extra_duplications = len(val)-2
                temp_dup_info[val[0]] = val[1:]

            thisInfoList = data["targets"][thisTargetID]

            thisInfoList[3] += extra_duplications
            data["ladderpath-index"] += extra_duplications
            data["size-index"] += (extra_duplications * thisInfoList[1])

            if _for_duplications_info_ladderonExist(thisInfoList[0]):
                newLadderonID = len(data['ladderons'])
                data['ladderons'][newLadderonID] = thisInfoList[:3] + [{ thisTargetID:[0] }]

                for key, infoList in data['ladderons'].items():
                    if key != newLadderonID and thisTargetID in infoList[3]:
                        infoList[3][newLadderonID] = infoList[3].pop(thisTargetID)

                data["targets"][thisTargetID][0] = [newLadderonID]

        data["order-index"] = data["size-index"] - data["ladderpath-index"]
        data['duplications_info'] = temp_dup_info


    if fill_ladderons_STR:
        fill_lpjson_STR(data)


    if save_file_name:
        save_file_name += '.json' if not save_file_name.endswith('.json') else ''
        save_ladderpath_json(data, save_file_name)

    return data


def _for_duplications_info_ladderonExist(thisInfoList_0): # thisInfoList_0 is thisInfoList[0]

    if len(thisInfoList_0) > 1:
        return True
    if isinstance(thisInfoList_0[0], int):
        return False
    else:
        if len(thisInfoList_0[0]) > 1:
            return True
    return False




def save_ladderpath_json(data, save_file_name):
    with open(save_file_name, "w", encoding="utf-8", newline='') as f:
        json.dump(data, f, ensure_ascii=False, indent=4)



def load_ladderpath_json(filename):
    with open(filename, 'r', encoding='utf-8') as file:
        lpjson = json.load(file, object_pairs_hook=convert_keys)
    return lpjson

def convert_keys(pairs):
    return {int(k) if is_integer(k) else k: v for k, v in pairs}

def is_integer(s):
    if s.startswith('-'):
        return s[1:].isdigit()
    return s.isdigit()





def clear_lpjson_STR(lpjson):
    for infoList in lpjson['ladderons'].values():
        infoList[2] = ''
    for infoList in lpjson['targets'].values():
        infoList[2] = ''



def fill_lpjson_STR(lpjson):
    for ldID in reversed(lpjson['ladderons'].keys()):
        lpjson['ladderons'][ldID][2] = fill(ldID, lpjson)
    for infoList in lpjson['targets'].values():
        if infoList[2]:
            assert len(infoList[2]) == infoList[1]
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
        assert len(lpjson['ladderons'][ldID][2]) == lpjson['ladderons'][ldID][1]
        return lpjson['ladderons'][ldID][2]
    thisstr = ''
    for i_COMP in lpjson['ladderons'][ldID][0]:
        if isinstance(i_COMP, int):
            thisstr += fill(i_COMP, lpjson)
        else:
            thisstr += i_COMP
    lpjson['ladderons'][ldID][2] = thisstr
    return thisstr




def POM_from_JSON(lpjson, display_str=False):
    temp = bbb_multiplicity(lpjson)
    pom_final_pool = {0: temp}

    fill_lpjson_STR(lpjson)
    lds = lpjson['ladderons']

    level_known = {}
    for ldID in reversed(lds.keys()):
        find_level(ldID, level_known, lds)

    each_level_comp = defaultdict(list)
    for ldID, level in level_known.items():
        each_level_comp[level].append(ldID)
    if len(each_level_comp) > 0:
        assert len(each_level_comp) == max(each_level_comp.keys())

    for i in range(1, len(each_level_comp)+1):
        levelinfo = {}
        for ldID in each_level_comp[i]:
            levelinfo[lds[ldID][2]] = sum(len(val) for val in lds[ldID][3].values())-1
        original_order = list(levelinfo.keys())
        temp = sorted(levelinfo.items(), key=lambda item: (-item[1], -original_order.index(item[0])))
        pom_final_pool[i] = dict(temp)

    for ID, pos in lpjson['duplications_info'].items():
        if isinstance(pos, int):
            thismulti = pos - 1
        else:
            thismulti = len(pos) - 1

        for val_strs in pom_final_pool.values():
            if lpjson['targets'][ID][2] in val_strs:
                val_strs[ lpjson['targets'][ID][2] ] += thismulti

    pom_str = ''
    if display_str:
        pom_str = str_partially_ordered_multiset(pom_final_pool)
        print(pom_str)

    return pom_final_pool, pom_str

    # {0: {'A': 3, 'D': 3, 'C': 2, 'B': 1, 'E': 1, 'F': 1},
    #  1: {'EF': 2, 'BC': 2, 'DCD': 1},
    #  2: {'DBC': 1},
    #  3: {'DBCDBC': 1}}

def find_level(x, level_known, lds):
    if x in level_known:
        return level_known[x]
    else:
        infoList0 = lds[x][0]
        if len(infoList0) == 1 and isinstance(infoList0[0], str):
            level_known[x] = 1
            return 1
        else:
            comp_levels = []
            for i in infoList0:
                if isinstance(i, int):
                    comp_levels.append( 1 + find_level(i, level_known, lds) )
            level_known[x] = max(comp_levels)
            return level_known[x]

def str_partially_ordered_multiset(pom):
    pom_str = '{ '

    pom_str += ' // '.join(
        ', '.join(f'{ld}({multi})' if multi != 1 else f'{ld}' for ld, multi in infoList.items()) 
        for infoList in pom.values()
    )
    pom_str += ' }'
    return pom_str




def disp3index(lpjson):
    index3 = (lpjson['ladderpath-index'], lpjson['order-index'], lpjson['size-index'])
    print(f'( Ladderpath-index:{index3[0]},  Order-index:{index3[1]},  Size-index:{index3[2]} )')
    return index3



def ellipse_len(seq):
    # length = np.sqrt(len(seq))/2
    length = (len(seq)**(1/3) )/2
    return length

def draw_laddergraph(lpjson, show_longer_than = 0, 
    style = "ellipse", 
    TargetNames_UserDefined = None,
    warning_n_ladderons_to_show = 500,
    rankdir = "BT", 
    color = "grey", 
    figsize = None,
    save_fig_name = None, 
    figformat = "pdf", 
    cleanGVfile=True,
    render_fig=True,
    show_node_size=True):



    fill_lpjson_STR(lpjson)

    from . import helpers_laddergraph

    g = helpers_laddergraph.draw_laddergraph(lpjson, show_longer_than=show_longer_than,
        style = style, 
        TargetNames_UserDefined = TargetNames_UserDefined,
        warning_n_ladderons_to_show = warning_n_ladderons_to_show,
        rankdir = rankdir, color = color, figsize = figsize,
        save_fig_name = save_fig_name, figformat = figformat, cleanGVfile=cleanGVfile,
        render_fig=render_fig,
        show_node_size=show_node_size)

    return g




def bbb_multiplicity(lpjson):
    bbb_multi = defaultdict(int)
    for item_name in ['ladderons', 'targets']:
        for infoList in lpjson[item_name].values():
            for i in infoList[0]:
                if isinstance(i, str):

                    for letter, count in Counter(i).items():
                        bbb_multi[letter] += count
    temp = dict(sorted(bbb_multi.items(), key=lambda item: (-item[1], lpjson['basic_building_blocks'].index(item[0])) ))
    return temp



def reconstruct_targets_list(lpjson):
    targetslist = []
    for target_id, info_list in lpjson['targets'].items():
        repetition = info_list[3] if len(info_list) > 3 else 1
        targetslist.extend([target_id] * int(repetition))
    return targetslist
