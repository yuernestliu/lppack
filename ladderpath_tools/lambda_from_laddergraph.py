"""
Version 2.1.0
Authors: ecsLab (Jingwen Zhang et al.)
    revised by Xiaojun Hu
    updated by codex
Date: 2026.03.31

利用梯图来计算其中若干 targets (leaves) 的 lambda。

当前版本相较于 2.0.2 增加了几项显式约束：
1. `lpjson` 必须来自 `lp.get_ladderpath()` 或 `lp.load_ladderpath_json()`。
   不允许直接用 `json.load()` 读取梯径 JSON。
2. 当前版本不支持带 duplicates 的 `lpjson`，即 `duplications_info` 必须为空。
3. `leaves` 必须是非空列表，且元素类型必须全为 int 或全为 str。
"""

from .. import ladderpath as lp


def get_target_ID(two_leaf_strs0, lpjson):
    # 找到subset中几个target的序号
    strs_subset = two_leaf_strs0.copy()
    ids = []
    for ID, info in lpjson["targets"].items():
        if info[2] in strs_subset:
            ids.append(ID)
            strs_subset.remove(info[2])
    if len(ids) != len(two_leaf_strs0):
        raise ValueError(
            "输入的 target 字符串未能全部在 lpjson['targets'] 中匹配到。"
        )
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
        if type(i) == int:
            return True
    return False


def _validate_lpjson_ready(lpjson):
    if not isinstance(lpjson, dict):
        raise ValueError("lpjson 必须是 dict。")

    if "targets" not in lpjson or "ladderons" not in lpjson:
        raise ValueError("lpjson 缺少 'targets' 或 'ladderons' 字段。")

    if not lpjson["targets"]:
        raise ValueError("lpjson['targets'] 不能为空。")

    target_key = next(iter(lpjson["targets"].keys()))
    ladderon_key = next(iter(lpjson["ladderons"].keys())) if lpjson["ladderons"] else None

    if not isinstance(target_key, int) or (
        ladderon_key is not None and not isinstance(ladderon_key, int)
    ):
        raise ValueError(
            "lpjson 的 key 类型不正确。"
            "请使用 lp.get_ladderpath(...) 或 lp.load_ladderpath_json(...) 获取 lpjson，"
            "不要直接使用 json.load()。"
        )

    if lpjson.get("duplications_info"):
        raise ValueError(
            "当前 lambda_from_laddergraph.get() 不支持带 duplicates 的 lpjson。"
            "请先去重，或使用不含 duplications_info 的梯径对象。"
        )


def _validate_leaves(leaves):
    if not isinstance(leaves, list) or len(leaves) == 0:
        raise ValueError("leaves 必须是非空 list。")

    first = leaves[0]
    if isinstance(first, int):
        if not all(isinstance(x, int) for x in leaves):
            raise ValueError("leaves 中的元素类型必须一致，不能混用 int 和 str。")
    elif isinstance(first, str):
        if not all(isinstance(x, str) for x in leaves):
            raise ValueError("leaves 中的元素类型必须一致，不能混用 int 和 str。")
    else:
        raise ValueError("leaves 只允许是 int 列表或 str 列表。")

    if len(set(leaves)) != len(leaves):
        raise ValueError(
            "当前 lambda_from_laddergraph.get() 不支持重复 leaves。"
            "请先去重或明确 duplicates 的语义后再计算。"
        )


# def get(leaves1or2, lpjson): # leaves1or2 could only be 1 or 2 leaf
#     ids = None
#     if len(leaves1or2) == 1 or len(leaves1or2) == 2:
#         lp.fill_lpjson_STR(lpjson)
#         if isinstance(leaves1or2[0], int):
#             if leaves1or2[0] in lpjson['targets']:
#                 if len(leaves1or2) == 1:
#                     ids = leaves1or2
#                 else: # len(leaves1or2) == 2
#                     if leaves1or2[1] in lpjson['targets']:
#                         ids = leaves1or2
#         else:
#             ids = get_target_ID(leaves1or2, lpjson)
    
#     if ids is None:
#         print('!!!Wrong: input wrong -> leaves1or2')
#         return None

#     basic_block_num=0
#     comp=[]
#     ladderon_dic={}
#     for i in ids:
#         comp+=lpjson["targets"][i][0]
#     pool=set()

#     while if_id_in_COMP(comp):
#         ids,basic_blocks=separate_basic_blocks(comp)
#         comp=[]
#         basic_block_num+=basic_blocks
#         id_set=set()
#         for i in ids:
#             if i in ladderon_dic:
#                 ladderon_dic[i]+=1
#             else:
#                 ladderon_dic[i]=1
#             if i not in id_set and i not in pool:
#                 comp+=lpjson["ladderons"][i][0]
#                 pool.add(i)
#                 id_set.add(i)

#     for i in comp:
#         basic_block_num+=len(i)
#     count=basic_block_num

#     for ladderon,times in ladderon_dic.items():
#         count+=(times-1)
#     return count

    
def get(leaves, lpjson):  # 支持任意数量的 leaf
    _validate_leaves(leaves)
    _validate_lpjson_ready(lpjson)

    lp.fill_lpjson_STR(lpjson)

    # 判断是字符串形式还是id形式
    if isinstance(leaves[0], int):
        for leaf in leaves:
            if leaf not in lpjson['targets']:
                raise ValueError(f"输入的 target ID {leaf} 不在 lpjson['targets'] 中。")
        ids = leaves
    else:
        ids = get_target_ID(leaves.copy(), lpjson)

    basic_block_num = 0
    comp = []
    ladderon_dic = {}
    for i in ids:
        comp += lpjson["targets"][i][0]
    pool = set()

    while if_id_in_COMP(comp):
        ids_temp, basic_blocks = separate_basic_blocks(comp)
        comp = []
        basic_block_num += basic_blocks
        id_set = set()
        for i in ids_temp:
            if i in ladderon_dic:
                ladderon_dic[i] += 1
            else:
                ladderon_dic[i] = 1
            if i not in id_set and i not in pool:
                comp += lpjson["ladderons"][i][0]
                pool.add(i)
                id_set.add(i)

    for i in comp:
        basic_block_num += len(i)
    count = basic_block_num

    for ladderon, times in ladderon_dic.items():
        count += (times - 1)
    return count
