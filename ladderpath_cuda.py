"""
Version 3.0.0
Authors: Hu Zhu
Data: 2026.06.12

CUDA GPU backend for Ladderpath.
This module receives unique targets; the public entry point handles multiplicity.
"""

import json
import math
import time
from collections import defaultdict, Counter
from typing import List
from numba import cuda
import numpy as np
import os
from numba import njit
import pydivsufsort as pyd
os.environ['NUMBA_CUDA_LOW_OCCUPANCY_WARNINGS'] = '0'


VERSION_ladderpathpy_cuda = '3.0.0'
VERSION_lpjson = '3.0.0'


@cuda.jit
def cal_delect_gpu(d_suffix_array,d_lcp,now_pos,now_pos_del,now_len):
    idx =cuda.grid(1)
    if idx < len(d_suffix_array):

        over_len = d_lcp[idx]
        this_start = d_suffix_array[idx]
        this_end = this_start + over_len - 1

        for pos in now_pos_del:
            if over_len < 2:
                d_lcp[idx]=0
                return
            start = pos
            end = pos + now_len-1

            if end < this_start or start > this_end:
                continue
            elif start <= this_start:
                d_lcp[idx] = 0
                return
            else:
                over = start - this_start

            over_len = min(over,over_len)

        for pos in now_pos:
            if over_len < 2:
                d_lcp[idx]=0
                return
            start = pos
            end = pos + now_len-1

            if end < this_start or start > this_end:
                continue
            elif start > this_start:
                over = start - this_start
            else:
                over = end - this_start + 1

            over_len = min(over, over_len)
        d_lcp[idx] = over_len


@cuda.jit
def init_lcp_for_s_list(d_lcp,d_sa,end_list):
    idx,idy = cuda.grid(2)
    if idx < len(d_lcp):
        if idy < len(end_list)-1:
            start = d_sa[idx]
            break_point=end_list[idy+1]
            if start <break_point:
                cuda.atomic.min(d_lcp,idx,break_point-start)

# @cuda.jit
# def init_lcp_for_s_list_one(d_lcp,d_sa,end_list):
#     idx,idy = cuda.grid(1)
#     if idx < len(d_lcp):
#         for i in range(1,len(end_list)):
#             start = d_sa[idx]
#             valid_len=d_lcp[idx]
#             break_start=end_list[i-1]
#             break_end=end_list[i]
#             if break_start <= start < break_end:
#                 d_lcp[idx]=min(valid_len,break_end-start)


@cuda.jit
def cal_del_gpu_same_len(d_suffix_array,d_last_select,max_len,d_valid):
    idx =cuda.grid(1)
    if idx < len(d_suffix_array):
        if d_valid[idx]:
            this_start = d_suffix_array[idx]
            for select_start in d_last_select:
                if select_start+max_len <= this_start or this_start+max_len <= select_start:
                    continue
                d_valid[idx]=False
                return


@cuda.jit
def cal_del_gpu_same_len_v1(d_suffix_array, d_last_select, max_len, d_valid):
    idx = cuda.grid(1)
    if idx >= len(d_suffix_array):
        return

    if not d_valid[idx]:
        return

    this_start = d_suffix_array[idx]
    n_select = len(d_last_select)
    if n_select == 0:
        return


    low = 0
    high = n_select
    while low < high:
        mid = (low + high) // 2
        if d_last_select[mid] < this_start:
            low = mid + 1
        else:
            high = mid

    if low > 0:
        select_start_prev = d_last_select[low - 1]
        if (this_start - select_start_prev) < max_len:
            d_valid[idx] = False
            return

    if low < n_select:
        select_start_next = d_last_select[low]
        if (select_start_next - this_start) < max_len:
            d_valid[idx] = False
            return


@cuda.jit
def count_non_overlap_gpu(d_suffix_array,d_valid_flag_sa,d_str_pos_in_sa,max_len,d_non_overlap_len):
    idx =cuda.grid(1)
    if idx < len(d_str_pos_in_sa):
        if d_non_overlap_len[idx] < 2:
            return
        start_end=d_str_pos_in_sa[idx]
        start = start_end[0]
        end = start_end[1]

        valid_len=0
        last_pos=-1
        for index in range(start,end):
            if d_valid_flag_sa[index]:
                now_pos=d_suffix_array[index]
                if last_pos == -1 or now_pos >= last_pos + max_len:
                    last_pos=now_pos
                    valid_len+=1
        d_non_overlap_len[idx]=valid_len

@cuda.jit
def get_max_non_overlap_pos(d_sort_sa,
                            d_valid_flag,
                            d_str_pos_in_sa,
                            max_index,
                            max_len,
                            d_non_overlap_pos_out):

    idx =cuda.grid(1)
    if idx == 0 :
        start_end=d_str_pos_in_sa[max_index]
        start = start_end[0]
        end = start_end[1]

        valid_len=0
        last_pos=-1
        for index in range(start,end):
            if d_valid_flag[index]:
                now_pos=d_sort_sa[index]
                if last_pos == -1 or now_pos >= last_pos + max_len:
                    last_pos=now_pos
                    d_non_overlap_pos_out[valid_len]=now_pos
                    valid_len+=1

@cuda.jit
def find_breakpoints_gpu(d_text, d_suffixes, d_breakpoints, max_len):
    idx = cuda.grid(1)
    n_suffixes = d_suffixes.shape[0]

    if idx < n_suffixes - 1:
        d_breakpoints[idx] = False
        for i in range(max_len):
            c1 = d_text[d_suffixes[idx] + i]
            c2 = d_text[d_suffixes[idx + 1] + i]
            if c1 != c2:
                d_breakpoints[idx] = True
                return

def find_most_frequent_longest_substring(suffix_array, lcp, join_s_int, end_list,max_str,first_pos=False,all_pos=False):
    len_s=len(join_s_int)
    threads_per_block = 256
    d_text = cuda.to_device(join_s_int)

    lcp[1:] = np.maximum(lcp[1:], lcp[:-1])

    sum_strs=len(end_list)
    if sum_strs > 1:
        # blocks_per_grid = math.ceil(len(suffix_array) / threads_per_block)
        threads_block = (32, 16)
        blocks_grid_x = math.ceil(len_s / threads_block[0])
        blocks_grid_y = math.ceil(sum_strs / threads_block[1])
        blocks_grid = (blocks_grid_x, blocks_grid_y)

        d_lcp = cuda.to_device(lcp)
        d_sa = cuda.to_device(suffix_array)
        d_end_list = cuda.to_device(end_list)

        init_lcp_for_s_list[blocks_grid,threads_block](d_lcp,d_sa, d_end_list)

        lcp=d_lcp.copy_to_host()

    lp=0
    remain_len=len_s
    max_sub=min(math.floor(len_s/2), max_str)

    select_ld_info_muti=[]

    now_id=0
    ld_info_list=[]
    while True:
        mask = lcp > 1
        lcp = lcp[mask]
        lcp = np.where(lcp > max_sub, max_sub, lcp)
        if lcp.size == 0:
            break
        suffix_array = suffix_array[mask]

        while True:
            max_len = lcp.max()
            if max_len < 2:
                break


            temp_mask= lcp == max_len
            temp_sa=suffix_array[temp_mask]
            lcp = np.where(lcp==max_len, lcp - 1, lcp)
            if temp_sa.size < 2:

                selected = []
                continue

            n_suffixes = len(temp_sa)
            d_suffixes = cuda.to_device(temp_sa)
            d_breakpoints = cuda.device_array(n_suffixes - 1, dtype=np.bool_)

            blocks_per_grid = math.ceil((n_suffixes - 1) / threads_per_block)
            find_breakpoints_gpu[blocks_per_grid, threads_per_block](
                d_text, d_suffixes, d_breakpoints, max_len
            )

            breakpoints = d_breakpoints.copy_to_host()  # shape: (n_suffixes-1,)

            break_indices = np.where(breakpoints)[0]
            if break_indices.size == n_suffixes-1:

                selected = []
                continue

            group_starts = np.concatenate([[0], break_indices + 1])
            group_ends = np.concatenate([break_indices + 1, [n_suffixes]])

            new_sa=np.array([], dtype=np.int64)


            now_start=0
            start_end=[]
            count_sub_str=0

            for start, end in zip(group_starts, group_ends):
                group_size = end - start
                if group_size > 1:
                    sort_pos = np.sort(temp_sa[start:end])
                    new_sa = np.concatenate((new_sa, sort_pos))
                    start_end.append([now_start, now_start + group_size])
                    now_start += group_size
                    count_sub_str += 1

            d_sort_sa=cuda.to_device(new_sa)
            d_valid_flag=cuda.to_device(np.ones(now_start, dtype=np.bool_))
            d_str_pos_in_sa=cuda.to_device(np.array(start_end,dtype=np.int64))
            d_non_overlap_len=cuda.to_device(np.full(count_sub_str,2, dtype=np.int64))

            blocks_per_grid_sa = math.ceil(now_start / threads_per_block)
            blocks_per_grid_substr = math.ceil(count_sub_str / threads_per_block)

            selected=[]
            selected_del=[]
            while True:

                count_non_overlap_gpu[blocks_per_grid_substr, threads_per_block](
                    d_sort_sa, d_valid_flag, d_str_pos_in_sa, max_len, d_non_overlap_len
                )


                non_overlap_len=d_non_overlap_len.copy_to_host()
                max_index=np.argmax(non_overlap_len)
                max_count=int(non_overlap_len[max_index])

                if max_count < 2:

                    break


                d_non_overlap_pos_out = cuda.device_array(max_count, dtype=np.int64)


                get_max_non_overlap_pos[1, 128](
                    d_sort_sa, d_valid_flag, d_str_pos_in_sa,
                    max_index, max_len,
                    d_non_overlap_pos_out
                )


                non_overlap_pos = d_non_overlap_pos_out.copy_to_host()


                cal_del_gpu_same_len_v1[blocks_per_grid_sa, threads_per_block](
                    d_sort_sa, d_non_overlap_pos_out, max_len, d_valid_flag
                )

                select_pos=non_overlap_pos[0]
                delect_pos=non_overlap_pos[1:]

                selected.append(select_pos)
                selected_del.extend(delect_pos)

                if first_pos:
                    select_ld_info_muti.append([select_pos,select_pos+max_len,max_count])
                elif all_pos:
                    ld_info_list.append([non_overlap_pos,max_len,now_id])
                    now_id += 1

                lp += max_count - 1
                remain_len = remain_len - max_len * (max_count - 1)
                max_sub =math.floor(min(max(min(remain_len - max_len, max_len), max_len / 2, (remain_len - max_len) / 2),max_str))

            if selected:
                break


        if max_len<2:
            break

        d_lcp = cuda.to_device(lcp)
        d_suffix_array = cuda.to_device(suffix_array)
        blocks_per_grid = math.ceil(len(suffix_array) / threads_per_block)

        d_select=cuda.to_device(selected)
        d_select_del=cuda.to_device(selected_del)
        cal_delect_gpu[blocks_per_grid, threads_per_block](d_suffix_array,d_lcp,d_select,d_select_del,max_len)

        lcp=d_lcp.copy_to_host()

    lp=lp+remain_len

    if first_pos:
        return lp,select_ld_info_muti
    elif all_pos:
        return lp,ld_info_list
    else:
        return lp

# @cuda.jit
# def find_father_id(all_start,all_end,ladder_start,father_end_id,father_id):
#     idx = cuda.grid(1)

#     if idx < ladder_start.shape[0]:
#         start = ladder_start[idx]
#         limit = father_end_id[idx]

#         res = -1
#         for k in range(limit, -1, -1):

#             if all_start[k] <= start < all_end[k]:
#                 res = k


#         father_id[idx] = res

# def find_father_id_chr(father_id_all_chr,ladder_start,ladder_len,ladder_id):

#     fa_id_chr=np.array(father_id_all_chr)
#     for i in range(len(ladder_start)):
#         start = ladder_start[i]
#         end = start + ladder_len[i]
#         ld_id=ladder_id[i]
#         fa_id_chr[start:end] = ld_id
#     return fa_id_chr


def find_ladderpath(s_list, model=0, ch=False):
    ch = ch or any(not s.isascii() for s in s_list)


    end=0
    join_s=''
    end_list=[0]
    max_str=0
    for s in s_list:
        len_s=len(s)
        if len_s>max_str:
            max_str=len_s
        end += len_s
        join_s += s
        end_list.append(end)

    if ch:
        join_s_int = np.array([ord(c) for c in join_s], dtype=np.int32)
    else:
        join_s_int=join_s.encode('utf-8')

    sa=pyd.divsufsort(join_s_int)
    lcps=pyd.kasai(join_s_int, sa)

    if not ch:
        join_s_int=np.frombuffer(join_s_int, dtype=np.uint8)

    if model==0:
        return find_most_frequent_longest_substring(sa, lcps, join_s_int, end_list, max_str)

    elif model==1:
        lp,first_ld_pos=find_most_frequent_longest_substring(sa, lcps, join_s_int, end_list, max_str,first_pos=True)
        ladder_dic={}
        for ladder_info in first_ld_pos:
            start_pos=ladder_info[0]
            end_pos=ladder_info[1]
            multi=ladder_info[2]
            ladder_dic[join_s[start_pos:end_pos]]=multi
        return ladder_dic

    else:
        join_id_list = []
        now_id = -1
        ta_ld_dic = {}
        now_start = 0
        for s in s_list:
            join_id_list.extend([now_id] * len(s))
            ta_ld_dic[now_id] = [now_start, {}, {},len(s)]
            now_id -= 1
            now_start += len(s)

        join_id_list = np.array(join_id_list)
        valid_flag = np.ones(len(join_id_list), dtype=bool)

        lp, ld_info_list = find_most_frequent_longest_substring(sa, lcps, join_s_int, end_list, max_str, all_pos=True)

        for ld_info in ld_info_list:
            start_poses = ld_info[0]
            ld_len = ld_info[1]
            ld_id = ld_info[2]

            component_dic = defaultdict(list)

            select_start_pos = start_poses[0]
            for start_pos in start_poses:
                father_id = join_id_list[start_pos]
                father_info = ta_ld_dic[father_id]
                start_pos_in_f = int(start_pos - father_info[0])
                component_dic[int(father_id)].append(start_pos_in_f)
                father_info[1][start_pos_in_f] = ld_id
                if start_pos == select_start_pos:
                    join_id_list[start_pos:start_pos + ld_len] = ld_id
                else:
                    valid_flag[start_pos:start_pos + ld_len] = False
            ta_ld_dic[ld_id] = [start_poses[0], {}, component_dic , ld_len]

        basic_blocks = set()
        for i, father_id in enumerate(join_id_list):
            if valid_flag[i]:
                c = join_s[i]
                basic_blocks.update(c)
                father_info = ta_ld_dic[father_id]
                start_pos_in_f = i - father_info[0]
                father_info[1][start_pos_in_f] = c

        return lp,end, ta_ld_dic,basic_blocks



# def get_lp(s_list,ch=False):

#     return find_ladderpath(s_list,0,ch)

# def get_lp_with_w(s_list,ch=False):

#     return find_ladderpath(s_list,1,ch)


def get_ladderpath(ss, ch=False, show_version=True):


    if show_version:
        print('lppack CUDA backend version:', VERSION_ladderpathpy_cuda)
        print('lpjson version:', VERSION_lpjson)

    lp,end, ta_ld_dic,basic_blocks = find_ladderpath(ss, model=2, ch=ch)


    data = {
        "info": VERSION_lpjson,
        "ladderpath-index": int(lp),
        "order-index": int(end - lp),
        "size-index": int(end),

        "basic_building_blocks": sorted(basic_blocks),
        "ladderons": {},
        "targets": {},
        "duplications_info": {}
    }

    target_dic={}
    ladder_dic={}
    for i,info in ta_ld_dic.items():
        comp=[]
        sub_dict=info[1]
        sorted_key=sorted(sub_dict.keys())
        last_comp=sub_dict[sorted_key[0]]

        for key in sorted_key[1:]:
            this_comp=sub_dict[key]
            if isinstance(last_comp, str) and isinstance(this_comp, str):
                last_comp =last_comp+this_comp
            else:
                comp.append(last_comp)
                last_comp=this_comp
        comp.append(last_comp)
        if i < 0:
            target_dic[i] = [comp,info[3],'',1]
        else:
            ladder_dic[i] = [comp,int(info[3]),'',info[2]]

    data['targets'] = target_dic
    data['ladderons'] = ladder_dic

    return data
