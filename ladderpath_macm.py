"""
Version 3.0.0
Authors: ecsLab
Data: 2026.06.12

Apple Silicon Mac GPU backend for ladderpath.

This file mirrors the CUDA backend structure in ladderpath_cuda.py, but ports the
CUDA kernels to MLX custom Metal kernels. It intentionally does not fall back to
CPU when `mode="macm"` is requested.
"""

import math
import platform
from collections import defaultdict

import numpy as np
import pydivsufsort as pyd


VERSION_ladderpathpy_macm = "3.0.0"
VERSION_lpjson = "3.0.0"
THREADS_PER_GROUP = 256


def _import_mlx():
    try:
        import mlx.core as mx
    except ImportError as exc:
        raise RuntimeError(
            "macm backend requires MLX. Install it in an arm64 Apple Silicon "
            "Python environment, for example: pip install mlx"
        ) from exc
    return mx


def is_macm_available():
    """Return True when the current platform is macOS."""
    return platform.system() == "Darwin"


def is_native_arm64_python():
    """Return True when the current Python process is native arm64."""
    return platform.machine() == "arm64"


def is_apple_gpu_backend_available():
    """Return True when MLX can see an Apple Metal GPU backend."""
    if not is_macm_available():
        return False
    try:
        mx = _import_mlx()
    except RuntimeError:
        return False
    return bool(mx.metal.is_available())


def _require_macm():
    if not is_macm_available():
        raise RuntimeError("macm backend requires macOS on Apple Silicon Mac.")
    mx = _import_mlx()
    if not mx.metal.is_available():
        raise RuntimeError("macm backend requires an available MLX Metal GPU backend.")
    return mx


def _to_numpy(mx, arr):
    mx.eval(arr)
    return np.array(arr)


def _mx_scalars_to_ints(mx, *arrs):
    mx.eval(*arrs)
    return tuple(int(np.array(arr)) for arr in arrs)


def _mx_int64(mx, arr):
    return mx.array(np.asarray(arr, dtype=np.int64), dtype=mx.int64)


def _mx_params(mx, *values):
    return mx.array(np.asarray(values, dtype=np.int64), dtype=mx.int64)


def _threadgroup(n):
    return (min(THREADS_PER_GROUP, max(1, int(n))), 1, 1)


def _mx_text(mx, arr):
    arr = np.asarray(arr)
    if arr.dtype == np.uint8:
        return mx.array(arr, dtype=mx.uint8)
    return mx.array(arr.astype(np.int32, copy=False), dtype=mx.int32)


def _kernel_init_lcp_for_s_list(mx):
    source = """
        long n = params[0];
        long end_n = params[1];
        uint idx = thread_position_in_grid.x;
        if (idx >= n) {
            return;
        }

        long valid_len = lcp[idx];
        long start = sa[idx];
        for (long j = 1; j < end_n; ++j) {
            long break_point = end_list[j];
            if (start < break_point) {
                long candidate = break_point - start;
                if (candidate < valid_len) {
                    valid_len = candidate;
                }
            }
        }
        lcp_out[idx] = valid_len;
    """
    return mx.fast.metal_kernel(
        name="lp_macm_init_lcp_for_s_list",
        input_names=["lcp", "sa", "end_list", "params"],
        output_names=["lcp_out"],
        source=source,
    )


def _kernel_find_breakpoints(mx):
    source = """
        long n_suffixes = params[0];
        long max_len = params[1];
        uint idx = thread_position_in_grid.x;
        if (idx >= n_suffixes - 1) {
            return;
        }

        breakpoints[idx] = false;
        for (long i = 0; i < max_len; ++i) {
            long p1 = suffixes[idx] + i;
            long p2 = suffixes[idx + 1] + i;
            if (text[p1] != text[p2]) {
                breakpoints[idx] = true;
                return;
            }
        }
    """
    return mx.fast.metal_kernel(
        name="lp_macm_find_breakpoints",
        input_names=["text", "suffixes", "params"],
        output_names=["breakpoints"],
        source=source,
    )


def _kernel_count_non_overlap(mx):
    source = """
        long n_groups = params[0];
        long max_len = params[1];
        uint idx = thread_position_in_grid.x;
        if (idx >= n_groups) {
            return;
        }

        long start = str_pos_in_sa[idx * 2];
        long end = str_pos_in_sa[idx * 2 + 1];
        long valid_len = 0;
        long last_pos = -1;

        for (long index = start; index < end; ++index) {
            if (valid_flag[index]) {
                long now_pos = sort_sa[index];
                if (last_pos == -1 || now_pos >= last_pos + max_len) {
                    last_pos = now_pos;
                    valid_len += 1;
                }
            }
        }
        non_overlap_len[idx] = valid_len;
    """
    return mx.fast.metal_kernel(
        name="lp_macm_count_non_overlap",
        input_names=["sort_sa", "valid_flag", "str_pos_in_sa", "params"],
        output_names=["non_overlap_len"],
        source=source,
    )


def _kernel_get_max_non_overlap_pos(mx):
    source = """
        long max_index = params[0];
        long max_len = params[1];
        long max_count = params[2];
        uint idx = thread_position_in_grid.x;
        if (idx != 0) {
            return;
        }

        long start = str_pos_in_sa[max_index * 2];
        long end = str_pos_in_sa[max_index * 2 + 1];
        long valid_len = 0;
        long last_pos = -1;

        for (long index = start; index < end; ++index) {
            if (valid_flag[index]) {
                long now_pos = sort_sa[index];
                if (last_pos == -1 || now_pos >= last_pos + max_len) {
                    last_pos = now_pos;
                    if (valid_len < max_count) {
                        non_overlap_pos_out[valid_len] = now_pos;
                    }
                    valid_len += 1;
                }
            }
        }
    """
    return mx.fast.metal_kernel(
        name="lp_macm_get_max_non_overlap_pos",
        input_names=["sort_sa", "valid_flag", "str_pos_in_sa", "params"],
        output_names=["non_overlap_pos_out"],
        source=source,
    )


def _kernel_cal_del_same_len(mx):
    source = """
        long n = params[0];
        long n_select = params[1];
        long max_len = params[2];
        uint idx = thread_position_in_grid.x;
        if (idx >= n) {
            return;
        }

        if (!valid[idx]) {
            valid_out[idx] = false;
            return;
        }

        long this_start = suffix_array[idx];
        if (n_select == 0) {
            valid_out[idx] = true;
            return;
        }

        long low = 0;
        long high = n_select;
        while (low < high) {
            long mid = (low + high) / 2;
            if (last_select[mid] < this_start) {
                low = mid + 1;
            } else {
                high = mid;
            }
        }

        if (low > 0) {
            long select_start_prev = last_select[low - 1];
            if ((this_start - select_start_prev) < max_len) {
                valid_out[idx] = false;
                return;
            }
        }

        if (low < n_select) {
            long select_start_next = last_select[low];
            if ((select_start_next - this_start) < max_len) {
                valid_out[idx] = false;
                return;
            }
        }

        valid_out[idx] = true;
    """
    return mx.fast.metal_kernel(
        name="lp_macm_cal_del_same_len",
        input_names=["suffix_array", "last_select", "valid", "params"],
        output_names=["valid_out"],
        source=source,
    )


def _kernel_update_lcp_after_selection(mx):
    source = """
        long n = params[0];
        long n_selected = params[1];
        long n_selected_del = params[2];
        long now_len = params[3];
        uint idx = thread_position_in_grid.x;
        if (idx >= n) {
            return;
        }

        long over_len = lcp[idx];
        long this_start = suffix_array[idx];
        long this_end = this_start + over_len - 1;

        for (long i = 0; i < n_selected_del; ++i) {
            if (over_len < 2) {
                lcp_out[idx] = 0;
                return;
            }
            long start = selected_del[i];
            long end = start + now_len - 1;
            if (end < this_start || start > this_end) {
                continue;
            } else if (start <= this_start) {
                lcp_out[idx] = 0;
                return;
            } else {
                long over = start - this_start;
                if (over < over_len) {
                    over_len = over;
                }
            }
        }

        for (long i = 0; i < n_selected; ++i) {
            if (over_len < 2) {
                lcp_out[idx] = 0;
                return;
            }
            long start = selected[i];
            long end = start + now_len - 1;
            if (end < this_start || start > this_end) {
                continue;
            } else if (start > this_start) {
                long over = start - this_start;
                if (over < over_len) {
                    over_len = over;
                }
            } else {
                long over = end - this_start + 1;
                if (over < over_len) {
                    over_len = over;
                }
            }
        }
        lcp_out[idx] = over_len;
    """
    return mx.fast.metal_kernel(
        name="lp_macm_update_lcp_after_selection",
        input_names=["suffix_array", "lcp", "selected", "selected_del", "params"],
        output_names=["lcp_out"],
        source=source,
    )


_KERNELS = {}


def _get_kernel(mx, name, factory):
    kernel = _KERNELS.get(name)
    if kernel is None:
        kernel = factory(mx)
        _KERNELS[name] = kernel
    return kernel


def _mlx_init_lcp_for_s_list(mx, lcp, suffix_array, end_list):
    if len(lcp) == 0:
        return lcp
    kernel = _get_kernel(mx, "init_lcp_for_s_list", _kernel_init_lcp_for_s_list)
    outputs = kernel(
        inputs=[
            _mx_int64(mx, lcp),
            _mx_int64(mx, suffix_array),
            _mx_int64(mx, end_list),
            _mx_params(mx, len(lcp), len(end_list)),
        ],
        grid=(len(lcp), 1, 1),
        threadgroup=_threadgroup(len(lcp)),
        output_shapes=[(len(lcp),)],
        output_dtypes=[mx.int64],
    )
    return _to_numpy(mx, outputs[0])


def _mlx_find_breakpoints_device(mx, text_mx, suffixes, max_len):
    n_suffixes = len(suffixes)
    if n_suffixes <= 1:
        return np.empty(0, dtype=np.bool_)
    kernel = _get_kernel(mx, "find_breakpoints", _kernel_find_breakpoints)
    outputs = kernel(
        inputs=[text_mx, _mx_int64(mx, suffixes), _mx_params(mx, n_suffixes, max_len)],
        grid=(n_suffixes - 1, 1, 1),
        threadgroup=_threadgroup(n_suffixes - 1),
        output_shapes=[(n_suffixes - 1,)],
        output_dtypes=[mx.bool_],
    )
    return _to_numpy(mx, outputs[0]).astype(np.bool_, copy=False)


def _mlx_count_non_overlap_device(mx, sort_sa_mx, valid_flag_mx, str_pos_in_sa_mx, n_groups, max_len):
    kernel = _get_kernel(mx, "count_non_overlap", _kernel_count_non_overlap)
    outputs = kernel(
        inputs=[sort_sa_mx, valid_flag_mx, str_pos_in_sa_mx, _mx_params(mx, n_groups, max_len)],
        grid=(n_groups, 1, 1),
        threadgroup=_threadgroup(n_groups),
        output_shapes=[(n_groups,)],
        output_dtypes=[mx.int64],
    )
    return outputs[0]


def _mlx_get_max_non_overlap_pos_device(mx, sort_sa_mx, valid_flag_mx, str_pos_in_sa_mx, max_index, max_len, max_count):
    kernel = _get_kernel(mx, "get_max_non_overlap_pos", _kernel_get_max_non_overlap_pos)
    outputs = kernel(
        inputs=[
            sort_sa_mx,
            valid_flag_mx,
            str_pos_in_sa_mx,
            _mx_params(mx, max_index, max_len, max_count),
        ],
        grid=(1, 1, 1),
        threadgroup=(1, 1, 1),
        output_shapes=[(max_count,)],
        output_dtypes=[mx.int64],
    )
    return outputs[0]


def _mlx_cal_del_same_len_device(mx, suffix_array_mx, last_select_mx, max_len, valid_mx, n_suffixes, n_select):
    kernel = _get_kernel(mx, "cal_del_same_len", _kernel_cal_del_same_len)
    outputs = kernel(
        inputs=[suffix_array_mx, last_select_mx, valid_mx, _mx_params(mx, n_suffixes, n_select, max_len)],
        grid=(n_suffixes, 1, 1),
        threadgroup=_threadgroup(n_suffixes),
        output_shapes=[(n_suffixes,)],
        output_dtypes=[mx.bool_],
    )
    return outputs[0]


def _mlx_update_lcp_after_selection(mx, suffix_array, lcp, selected, selected_del, now_len):
    if len(suffix_array) == 0:
        return lcp
    selected = np.asarray(selected, dtype=np.int64)
    selected_del = np.asarray(selected_del, dtype=np.int64)
    kernel = _get_kernel(mx, "update_lcp_after_selection", _kernel_update_lcp_after_selection)
    outputs = kernel(
        inputs=[
            _mx_int64(mx, suffix_array),
            _mx_int64(mx, lcp),
            _mx_int64(mx, selected),
            _mx_int64(mx, selected_del),
            _mx_params(mx, len(suffix_array), len(selected), len(selected_del), now_len),
        ],
        grid=(len(suffix_array), 1, 1),
        threadgroup=_threadgroup(len(suffix_array)),
        output_shapes=[(len(suffix_array),)],
        output_dtypes=[mx.int64],
    )
    return _to_numpy(mx, outputs[0])


def find_most_frequent_longest_substring(suffix_array, lcp, join_s_int, end_list, max_str, first_pos=False, all_pos=False):
    mx = _require_macm()

    suffix_array = np.asarray(suffix_array, dtype=np.int64)
    lcp = np.asarray(lcp, dtype=np.int64)
    join_s_int = np.asarray(join_s_int)
    end_list = np.asarray(end_list, dtype=np.int64)

    len_s = len(join_s_int)
    join_s_mx = _mx_text(mx, join_s_int)
    lcp[1:] = np.maximum(lcp[1:], lcp[:-1])

    if len(end_list) > 1:
        lcp = _mlx_init_lcp_for_s_list(mx, lcp, suffix_array, end_list)

    lp = 0
    remain_len = len_s
    max_sub = min(math.floor(len_s / 2), max_str)

    select_ld_info_muti = []
    now_id = 0
    ld_info_list = []

    while True:
        mask = lcp > 1
        lcp = lcp[mask]
        np.minimum(lcp, max_sub, out=lcp)
        if lcp.size == 0:
            break
        suffix_array = suffix_array[mask]

        while True:
            max_len = int(lcp.max())
            if max_len < 2:
                break

            temp_mask = lcp == max_len
            temp_sa = suffix_array[temp_mask]
            lcp[temp_mask] -= 1
            if temp_sa.size < 2:
                selected = []
                continue

            breakpoints = _mlx_find_breakpoints_device(mx, join_s_mx, temp_sa, max_len)
            break_indices = np.where(breakpoints)[0]
            if break_indices.size == len(temp_sa) - 1:
                selected = []
                continue

            group_starts = np.concatenate([[0], break_indices + 1])
            group_ends = np.concatenate([break_indices + 1, [len(temp_sa)]])

            sorted_suffix_groups = []
            group_ranges = []
            n_group_positions = 0

            for start, end in zip(group_starts, group_ends):
                group_size = end - start
                if group_size > 1:
                    sorted_positions = np.sort(temp_sa[start:end])
                    sorted_suffix_groups.append(sorted_positions)
                    group_ranges.append([n_group_positions, n_group_positions + group_size])
                    n_group_positions += group_size

            grouped_suffixes = np.concatenate(sorted_suffix_groups).astype(np.int64, copy=False)
            group_ranges = np.asarray(group_ranges, dtype=np.int64)
            n_groups = len(group_ranges)
            grouped_suffixes_mx = _mx_int64(mx, grouped_suffixes)
            group_ranges_mx = _mx_int64(mx, group_ranges)
            valid_flag_mx = mx.ones((n_group_positions,), dtype=mx.bool_)

            selected = []
            selected_del = []
            while True:
                non_overlap_len_mx = _mlx_count_non_overlap_device(
                    mx, grouped_suffixes_mx, valid_flag_mx, group_ranges_mx, n_groups, max_len
                )
                max_index_mx = mx.argmax(non_overlap_len_mx)
                max_count_mx = mx.max(non_overlap_len_mx)
                max_index, max_count = _mx_scalars_to_ints(mx, max_index_mx, max_count_mx)

                if max_count < 2:
                    break

                non_overlap_pos_mx = _mlx_get_max_non_overlap_pos_device(
                    mx, grouped_suffixes_mx, valid_flag_mx, group_ranges_mx, max_index, max_len, max_count
                )
                valid_flag_mx = _mlx_cal_del_same_len_device(
                    mx, grouped_suffixes_mx, non_overlap_pos_mx, max_len, valid_flag_mx, n_group_positions, max_count
                )
                non_overlap_pos = _to_numpy(mx, non_overlap_pos_mx)

                select_pos = int(non_overlap_pos[0])
                delect_pos = non_overlap_pos[1:]

                selected.append(select_pos)
                selected_del.extend(int(pos) for pos in delect_pos)

                if first_pos:
                    select_ld_info_muti.append([select_pos, select_pos + max_len, max_count])
                elif all_pos:
                    ld_info_list.append([non_overlap_pos, max_len, now_id])
                    now_id += 1

                lp += max_count - 1
                remain_len = remain_len - max_len * (max_count - 1)
                max_sub = math.floor(
                    min(
                        max(
                            min(remain_len - max_len, max_len),
                            max_len / 2,
                            (remain_len - max_len) / 2,
                        ),
                        max_str,
                    )
                )

            if selected:
                break

        if max_len < 2:
            break

        lcp = _mlx_update_lcp_after_selection(mx, suffix_array, lcp, selected, selected_del, max_len)

    lp = lp + remain_len

    if first_pos:
        return lp, select_ld_info_muti
    if all_pos:
        return lp, ld_info_list
    return lp


def find_ladderpath(s_list, model=0, ch=False):
    _require_macm()
    ch = ch or any(not s.isascii() for s in s_list)

    end = 0
    join_s = ""
    end_list = [0]
    max_str = 0
    for s in s_list:
        len_s = len(s)
        if len_s > max_str:
            max_str = len_s
        end += len_s
        join_s += s
        end_list.append(end)

    if ch:
        join_s_int = np.array([ord(c) for c in join_s], dtype=np.int32)
    else:
        join_s_int = join_s.encode("utf-8")

    sa = pyd.divsufsort(join_s_int)
    lcps = pyd.kasai(join_s_int, sa)

    if not ch:
        join_s_int = np.frombuffer(join_s_int, dtype=np.uint8)

    if model == 0:
        return find_most_frequent_longest_substring(sa, lcps, join_s_int, end_list, max_str)

    if model == 1:
        lp, first_ld_pos = find_most_frequent_longest_substring(
            sa, lcps, join_s_int, end_list, max_str, first_pos=True
        )
        ladder_dic = {}
        for ladder_info in first_ld_pos:
            start_pos = ladder_info[0]
            end_pos = ladder_info[1]
            multi = ladder_info[2]
            ladder_dic[join_s[start_pos:end_pos]] = multi
        return ladder_dic

    join_id_list = []
    now_id = -1
    ta_ld_dic = {}
    now_start = 0
    for s in s_list:
        join_id_list.extend([now_id] * len(s))
        ta_ld_dic[now_id] = [now_start, {}, {}, len(s)]
        now_id -= 1
        now_start += len(s)

    join_id_list = np.array(join_id_list)
    valid_flag = np.ones(len(join_id_list), dtype=bool)

    lp, ld_info_list = find_most_frequent_longest_substring(
        sa, lcps, join_s_int, end_list, max_str, all_pos=True
    )

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
        ta_ld_dic[ld_id] = [start_poses[0], {}, component_dic, ld_len]

    basic_blocks = set()
    for i, father_id in enumerate(join_id_list):
        if valid_flag[i]:
            c = join_s[i]
            basic_blocks.update(c)
            father_info = ta_ld_dic[father_id]
            start_pos_in_f = i - father_info[0]
            father_info[1][start_pos_in_f] = c

    return lp, end, ta_ld_dic, basic_blocks


def get_ladderpath(ss, ch=False, show_version=True):

    _require_macm()

    if show_version:
        print("lppack Apple Silicon backend version:", VERSION_ladderpathpy_macm)
        print("lpjson version:", VERSION_lpjson)
        if not is_native_arm64_python():
            print("Warning: current Python is not native arm64; performance may be limited.")

    lp, end, ta_ld_dic, basic_blocks = find_ladderpath(ss, model=2, ch=ch)

    data = {
        "info": VERSION_lpjson,
        "ladderpath-index": int(lp),
        "order-index": int(end - lp),
        "size-index": int(end),
        "basic_building_blocks": sorted(basic_blocks),
        "ladderons": {},
        "targets": {},
        "duplications_info": {},
    }

    target_dic = {}
    ladder_dic = {}
    for i, info in ta_ld_dic.items():
        comp = []
        sub_dict = info[1]
        sorted_key = sorted(sub_dict.keys())
        last_comp = sub_dict[sorted_key[0]]

        for key in sorted_key[1:]:
            this_comp = sub_dict[key]
            if isinstance(last_comp, str) and isinstance(this_comp, str):
                last_comp = last_comp + this_comp
            else:
                comp.append(last_comp)
                last_comp = this_comp
        comp.append(last_comp)
        if i < 0:
            target_dic[i] = [comp, info[3], "", 1]
        else:
            ladder_dic[i] = [comp, int(info[3]), "", info[2]]

    data["targets"] = target_dic
    data["ladderons"] = ladder_dic

    return data
